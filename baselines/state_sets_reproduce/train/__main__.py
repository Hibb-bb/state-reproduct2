import json
import os
from pathlib import Path
import pickle
import shutil
import re
from os.path import join, exists
from typing import Dict, List, Optional

import hydra
import torch
import sys

import lightning.pytorch as pl
from lightning.pytorch.loggers import CSVLogger, WandbLogger
from lightning.pytorch.callbacks import ModelCheckpoint
from omegaconf import DictConfig, OmegaConf
from lightning.pytorch.plugins.precision import MixedPrecision

from cell_load.utils.modules import get_datamodule

sys.path.append("/home/mohsen/projects/state-sets-reproduce/src/")

from state_sets_reproduce.models import (
    # scGPTForPerturbationModel,
    CPAPerturbationModel,
    # SCVIPerturbationModel,
    # LowRankLinearModel,
    # GEARSPerturbationModel,
)
from state_sets_reproduce.callbacks import BatchSpeedMonitorCallback

import logging


logger = logging.getLogger(__name__)
torch.set_float32_matmul_precision("medium")


def _get_train_data_dir_from_toml(toml_config_path: Optional[str]) -> Optional[str]:
    """Resolve training data directory from TOML config. Returns None if not applicable."""
    if not toml_config_path or not os.path.isfile(toml_config_path):
        return None
    try:
        import toml
        toml_cfg = toml.load(toml_config_path)
    except Exception:
        return None
    datasets = toml_cfg.get("datasets") or {}
    training = toml_cfg.get("training") or {}
    for name, split in training.items():
        if split == "train" and name in datasets:
            data_dir = datasets[name]
            if isinstance(data_dir, str) and os.path.isdir(data_dir):
                return data_dir
    return None


def _get_var_dims_from_metadata_and_h5ad(
    toml_config_path: Optional[str],
    dm,
) -> Optional[Dict]:
    """
    Build var_dims from metadata.json (n_genes, optionally gene_names) and first h5ad
    in the training data dir. Used when embed_key is null so we avoid cell_load's
    get_dim_for_obsm (which looks in obsm/ and fails for raw .X data).
    """
    data_dir = _get_train_data_dir_from_toml(toml_config_path)
    if not data_dir:
        return None
    meta_path = os.path.join(data_dir, "metadata.json")
    if not os.path.isfile(meta_path):
        return None
    try:
        with open(meta_path, "r") as f:
            meta = json.load(f)
    except Exception:
        return None
    n_genes = meta.get("n_genes")
    gene_names = meta.get("gene_names")  # optional list
    # Find first h5ad in data dir (for n_genes / gene_names if missing from metadata)
    h5ad_path = None
    for f in sorted(os.listdir(data_dir)):
        if f.endswith(".h5ad"):
            h5ad_path = os.path.join(data_dir, f)
            break
    if h5ad_path is not None and (n_genes is None or gene_names is None):
        try:
            import anndata
            adata = anndata.read_h5ad(h5ad_path, backed="r")
            if n_genes is None:
                n_genes = adata.n_vars
            if gene_names is None:
                gene_names = list(adata.var.index)
            try:
                adata.close()
            except AttributeError:
                pass
        except Exception as e:
            logger.warning("Could not read h5ad for var_dims: %s", e)
            if n_genes is None or gene_names is None:
                return None
    if n_genes is None or gene_names is None:
        return None
    n_genes = int(n_genes)
    if isinstance(gene_names, list) and len(gene_names) != n_genes:
        logger.warning("metadata n_genes=%d but gene_names length=%d; using len(gene_names)", n_genes, len(gene_names))
        n_genes = len(gene_names)
    pert_dim = len(dm.pert_onehot_map)
    batch_dim = len(dm.batch_onehot_map)
    pert_names = list(dm.pert_onehot_map.keys())
    return {
        "input_dim": n_genes,
        "output_dim": n_genes,
        "gene_dim": n_genes,
        "hvg_dim": n_genes,
        "pert_dim": pert_dim,
        "batch_dim": batch_dim,
        "gene_names": gene_names,
        "pert_names": pert_names,
    }


def get_lightning_module(
    model_type: str,
    data_config: dict,
    model_config: dict,
    training_config: dict,
    var_dims: dict,
):
    """Create model instance based on config."""
    # combine the model config and training config
    module_config = {**model_config, **training_config}
    module_config["embed_key"] = data_config["embed_key"]
    module_config["output_space"] = data_config["output_space"]
    module_config["gene_names"] = var_dims["gene_names"]

    print(f"Found {len(module_config['gene_names'])} genes in the data")
    module_config["batch_size"] = training_config["batch_size"]

    if data_config["output_space"] == "gene":
        gene_dim = var_dims["hvg_dim"]
    else:
        gene_dim = var_dims["gene_dim"]

    if model_type.lower() == "lowranklinear":
        return LowRankLinearModel.from_pretrained_embeddings(
            input_dim=var_dims["input_dim"],
            output_dim=var_dims["output_dim"],
            pert_dim=var_dims["pert_dim"],
            gene_dim=gene_dim,
            pert_names=var_dims["pert_names"],
            **module_config,
        )
    elif model_type.lower() == "cpa":
        return CPAPerturbationModel(
            input_dim=var_dims["input_dim"],
            output_dim=var_dims["output_dim"],
            pert_dim=var_dims["pert_dim"],
            gene_dim=gene_dim,
            **module_config,
        )
    elif model_type.lower() == "scvi":
        return SCVIPerturbationModel(
            input_dim=var_dims["input_dim"],
            gene_dim=gene_dim,
            hvg_dim=var_dims["hvg_dim"],
            output_dim=var_dims["output_dim"],
            pert_dim=var_dims["pert_dim"],
            batch_dim=var_dims["batch_dim"],
            **module_config,
        )
    elif (
        model_type.lower() == "scgpt-chemical" or model_type.lower() == "scgpt-genetic"
    ):
        pretrained_path = module_config["pretrained_path"]
        assert pretrained_path is not None, "pretrained_path must be provided for scGPT"

        model_dir = Path(pretrained_path)
        model_config_file = model_dir / "args.json"
        model_file = model_dir / "best_model.pt"

        model = scGPTForPerturbationModel(
            ntoken=module_config["ntoken"],
            n_drug_tokens=module_config[
                "n_perts"
            ],  # only used for chemical perturbations
            vocab=module_config["vocab"],
            gene_names=var_dims["gene_names"],
            d_model=module_config["d_model"],
            nhead=module_config["nhead"],
            d_hid=module_config["d_hid"],
            nlayers=module_config["nlayers"],
            nlayers_cls=module_config["n_layers_cls"],
            n_cls=1,
            dropout=module_config["dropout"],
            pad_token_id=module_config["pad_token_id"],
            pad_value=module_config["pad_value"],
            pert_pad_id=module_config["pert_pad_id"],
            do_mvc=module_config["do_MVC"],
            cell_emb_style=module_config["cell_emb_style"],
            mvc_decoder_style=module_config["mvc_decoder_style"],
            use_fast_transformer=module_config["use_fast_transformer"],
            lr=module_config["lr"],
            step_size_lr=module_config["step_size_lr"],
            include_zero_gene=module_config["include_zero_gene"],
            embed_key=module_config["embed_key"],
            perturbation_type=module_config["perturbation_type"],
            pert_names=var_dims["pert_names"],
            control_pert=module_config["control_pert"],
        )

        load_param_prefixes = module_config["load_param_prefixes"]

        if load_param_prefixes is not None:
            model_dict = model.model.state_dict()
            pretrained_dict = torch.load(model_file)
            pretrained_dict = {
                k: v
                for k, v in pretrained_dict.items()
                if any(
                    [
                        k.startswith(prefix)
                        for prefix in module_config["load_param_prefixes"]
                    ]
                )
            }
            for k, v in pretrained_dict.items():
                print(f"Loading params {k} with shape {v.shape}")

            model_dict.update(pretrained_dict)
            model.model.load_state_dict(model_dict)
        else:
            try:
                model.model.load_state_dict(torch.load(model_file))
                print(f"Loading all model params from {model_file}")
            except:
                # only load params that are in the model and match the size
                model_dict = model.model.state_dict()
                pretrained_dict = torch.load(model_file)
                pretrained_dict = {
                    k: v
                    for k, v in pretrained_dict.items()
                    if k in model_dict and v.shape == model_dict[k].shape
                }
                for k, v in pretrained_dict.items():
                    print(f"Loading params {k} with shape {v.shape}")

                model_dict.update(pretrained_dict)
                model.model.load_state_dict(model_dict)

        return model
    elif model_type.lower() == "gears":
        return GEARSPerturbationModel(
            input_dim=var_dims["input_dim"],
            output_dim=var_dims["output_dim"],
            pert_dim=var_dims["pert_dim"],
            gene_dim=gene_dim,
            pert_names=var_dims["pert_names"],
            **module_config,
        )
    else:
        raise ValueError(f"Unknown model type: {model_type}")


def get_latest_step_checkpoint(directory):
    # Get all checkpoint files
    files = os.listdir(directory)

    # Extract step numbers using regex, excluding files with 'val_loss'
    step_numbers = []
    for f in files:
        if f.startswith("step=") and "val_loss" not in f:
            # Extract the number between 'step=' and '.ckpt'
            match = re.search(r"step=(\d+)(?:-v\d+)?\.ckpt", f)
            if match:
                step_numbers.append(int(match.group(1)))

    if not step_numbers:
        raise ValueError("No checkpoint files found")

    # Get the maximum step number
    max_step = max(step_numbers)

    # Construct the checkpoint path
    checkpoint_path = join(directory, f"step={max_step}.ckpt")

    return checkpoint_path


def get_loggers(
    output_dir: str,
    name: str,
    wandb_project: str,
    wandb_entity: str,
    local_wandb_dir: str,
    use_wandb: bool = False,
    cfg: dict = None,
) -> List:
    """Set up logging to local CSV and optionally WandB."""
    # Always use CSV logger
    csv_logger = CSVLogger(save_dir=output_dir, name=name, version=0)
    loggers = [csv_logger]

    # Add WandB if requested
    if use_wandb:
        wandb_logger = WandbLogger(
            name=name,
            project=wandb_project,
            entity=wandb_entity,
            dir=local_wandb_dir,
            tags=cfg["wandb"].get("tags", []) if cfg else [],
        )
        if cfg is not None:
            wandb_logger.experiment.config.update(cfg)
        loggers.append(wandb_logger)

    return loggers


def get_checkpoint_callbacks(
    output_dir: str, name: str, val_freq: int, ckpt_every_n_steps: int
) -> List[ModelCheckpoint]:
    """Create checkpoint callbacks based on validation frequency."""
    checkpoint_dir = join(output_dir, name, "checkpoints")
    callbacks = []

    # Save best checkpoint based on validation loss
    best_ckpt = ModelCheckpoint(
        dirpath=checkpoint_dir,
        filename="step={step}-val_loss={val_loss:.4f}",
        save_last="link",  # Will create last.ckpt symlink to best checkpoint
        monitor="val_loss",
        mode="min",
        save_top_k=1,  # Only keep the best checkpoint
        every_n_train_steps=val_freq,
    )
    callbacks.append(best_ckpt)

    # Also save periodic checkpoints (without affecting the "last" symlink)
    periodic_ckpt = ModelCheckpoint(
        dirpath=checkpoint_dir,
        filename="{step}",
        save_last=False,  # Don't create/update symlink
        every_n_train_steps=ckpt_every_n_steps,
        save_top_k=-1,  # Keep all periodic checkpoints
    )
    callbacks.append(periodic_ckpt)

    return callbacks


@hydra.main(version_base=None, config_path="../configs", config_name="config")
def train(cfg: DictConfig) -> None:
    """Main training function."""
    # Convert config to YAML for logging
    cfg_yaml = OmegaConf.to_yaml(cfg, resolve=True)
    cfg = OmegaConf.to_container(cfg, resolve=True)
    print(cfg_yaml)

    # Setup output directory
    run_output_dir = join(cfg["output_dir"], cfg["name"])
    if os.path.exists(run_output_dir) and cfg["overwrite"]:
        print(f"Output dir {run_output_dir} already exists, overwriting")
        shutil.rmtree(run_output_dir)
    os.makedirs(run_output_dir, exist_ok=True)

    # Set up wandb directory if needed
    if cfg["use_wandb"]:
        os.makedirs(cfg["wandb"]["local_wandb_dir"], exist_ok=True)

    # Set random seeds
    if "train_seed" in cfg["training"]:
        pl.seed_everything(cfg["training"]["train_seed"])

    # if the provided pert_col is drugname_drugconc, hard code the value of control pert
    # this is because it's surprisingly hard to specify a list of tuples in the config as a string
    if cfg["data"]["kwargs"]["pert_col"] == "drugname_drugconc":
        cfg["data"]["kwargs"]["control_pert"] = "[('DMSO_TF', 0.0, 'uM')]"

    cfg["model"]["kwargs"]["control_pert"] = cfg["data"]["kwargs"]["control_pert"]

    # Initialize data module. this is backwards compatible with previous configs
    try:
        sentence_len = cfg["model"]["cell_set_len"]
    except KeyError:
        if cfg["model"]["name"].lower() in [
            "cpa",
            "gears",
            "scvi",
            "lowranklinear",
        ] or cfg["model"]["name"].lower().startswith("scgpt"):
            if (
                "cell_sentence_len" in cfg["model"]["kwargs"]
                and cfg["model"]["kwargs"]["cell_sentence_len"] > 1
            ):
                sentence_len = cfg["model"]["kwargs"]["cell_sentence_len"]
                cfg["training"]["batch_size"] = 1
            else:
                sentence_len = 1
        else:
            sentence_len = cfg["model"]["kwargs"]["transformer_backbone_kwargs"][
                "n_positions"
            ]

    if (
        cfg["model"]["name"].lower().startswith("scgpt")
    ):  # scGPT uses log-normalized expression
        # cfg["data"]["kwargs"]["hvg_names_uns_key"] = "hvg_names" if cfg["data"]["kwargs"]["train_task"] != "replogle" else None # TODO: better to not hardcode this

        model_dir = Path(cfg["model"]["kwargs"]["pretrained_path"])

        vocab_file = model_dir / "vocab.json"

        vocab = json.load(open(vocab_file, "r"))
        cfg["model"]["kwargs"]["pad_token_id"] = vocab["<pad>"]
        for s in cfg["model"]["kwargs"]["special_tokens"]:
            if s not in vocab:
                vocab[s] = len(vocab)

        cfg["model"]["kwargs"]["vocab"] = vocab
        cfg["model"]["kwargs"]["ntoken"] = len(vocab)
        cfg["model"]["kwargs"]["d_model"] = cfg["model"]["kwargs"]["embsize"]

        logger.info(f"Added vocab and hvg_names_uns_key to data kwargs for scGPT")

    dm = get_datamodule(
        name=cfg["data"]["name"],
        kwargs=cfg["data"]["kwargs"],
        batch_size=cfg["training"]["batch_size"],
        cell_sentence_len=sentence_len,
    )

    dm.setup()

    if cfg["model"]["name"].lower() in ["cpa", "scvi", "lowranklinear", "gears"] or cfg[
        "model"
    ]["name"].lower().startswith("scgpt"):
        cfg["model"]["kwargs"]["n_cell_types"] = len(dm.cell_type_onehot_map)
        cfg["model"]["kwargs"]["n_perts"] = len(dm.pert_onehot_map)
        cfg["model"]["kwargs"]["n_batches"] = len(dm.batch_onehot_map)

    if (
        cfg["model"]["name"].lower() == "gears"
        and cfg["data"]["kwargs"]["pert_col"].lower() == "cytokine"
    ):
        logger.info("Using parse coexpression graph for GEARS")
        cfg["model"]["kwargs"][
            "coexpression_graph_path"
        ] = "/large_storage/ctc/userspace/mohsen/state_revisions/gears_prep/parse_coexpression.csv"
    elif (
        cfg["model"]["name"].lower() == "gears"
        and cfg["data"]["kwargs"]["pert_col"].lower() == "drugname_drugconc"
    ):
        logger.info("Using tahoe coexpression graph for GEARS")
        cfg["model"]["kwargs"][
            "coexpression_graph_path"
        ] = "/large_storage/ctc/userspace/mohsen/state_revisions/gears_prep/tahoe_coexpression.csv"

    if cfg["model"]["name"].lower() == "lowranklinear":
        if (
            cfg["model"]["kwargs"]["pert_emb"] == "identity"
        ):  # Use the identity matrix as the perturbation embeddings (one-hot encoding)
            cfg["model"]["kwargs"]["pert_emb_path"] = "identity"
        elif (
            cfg["model"]["kwargs"]["pert_emb"] == "scgpt"
        ):  # scGPT: Genetic perturbation data
            cfg["model"]["kwargs"][
                "pert_emb_path"
            ] = f"/large_storage/goodarzilab/userspace/mohsen/VCI-models/scGPT/scGPT_human/gene_embeddings.h5"
        elif (
            cfg["model"]["kwargs"]["pert_emb"] == "tahoe_rdkit"
        ):  # Tahoe: Chemical perturbation data
            cfg["model"]["kwargs"][
                "pert_emb_path"
            ] = "/large_storage/goodarzilab/userspace/mohsen/VCI/tahoe/tahoe_rdkit_embs.h5"
        elif (
            cfg["model"]["kwargs"]["pert_emb"] == "gears_norman"
        ):  # Extract GEARS perturbation embeddings from the trained GEARS on Norman2019 dataset
            cfg["model"]["kwargs"][
                "pert_emb_path"
            ] = "/large_storage/goodarzilab/userspace/mohsen/VCI-models/GEARS/gears_norman.h5"
        else:
            raise ValueError(
                f"Unknown perturbation embedding: {cfg['model']['kwargs']['pert_emb']}"
            )

        if (
            cfg["model"]["kwargs"]["gene_emb"] == "training_data"
        ):  # Use the training data as the gene embeddings
            # 1. Perform PCA on the training data
            raise NotImplementedError("PCA on training data is not implemented yet")
        elif (
            cfg["model"]["kwargs"]["gene_emb"] == "gears_norman"
        ):  # Extract GEARS gene embeddings from the trained GEARS on Norman2019 dataset
            cfg["model"]["kwargs"][
                "gene_emb_path"
            ] = "/large_storage/goodarzilab/userspace/mohsen/VCI-models/GEARS/gears_norman.h5"
        elif (
            cfg["model"]["kwargs"]["gene_emb"] == "scgpt"
        ):  # Extract scGPT's vocabulary embeddings
            cfg["model"]["kwargs"][
                "gene_emb_path"
            ] = f"/large_storage/goodarzilab/userspace/mohsen/VCI-models/scGPT/scGPT_human/gene_embeddings.h5"
        else:
            raise ValueError(f"Unknown gene embedding: {cfg['model']['gene_emb']}")

    with open(join(run_output_dir, "config.yaml"), "w") as f:
        # f.write()
        new_cfg_yaml = OmegaConf.to_yaml(cfg, resolve=True)
        f.write(new_cfg_yaml)

    # When embed_key is null (raw .X), cell_load get_var_dims() uses get_dim_for_obsm(embed_key) and fails.
    # Build var_dims from metadata.json (+ first h5ad for gene_names / n_genes if missing) instead.
    dm_var_dims = None
    if cfg["data"]["kwargs"].get("embed_key") is None:
        dm_var_dims = _get_var_dims_from_metadata_and_h5ad(
            cfg["data"]["kwargs"].get("toml_config_path"),
            dm,
        )
        if dm_var_dims is not None:
            logger.info(
                "Using var_dims from metadata.json + h5ad (embed_key=null): n_genes=%d",
                dm_var_dims["input_dim"],
            )
    if dm_var_dims is None:
        # dm.embed_key = None
        dm_var_dims = dm.get_var_dims()

    print(dm.embed_key)
    dm.embed_key = "X"

    # if cfg["model"]["name"].lower() == "gears":
    # if cfg["data"]["kwargs"]["pert_col"].lower() == "drugname_drugconc":
    #     logger.info("Using tahoe gene names for the baselines ")
    #     path = "/large_storage/ctc/userspace/mohsen/state_revisions/gears_prep/tahoe_gene_names.txt"
    #     with open(path, "r") as f:
    #         gene_names = f.readlines()
    #     gene_names = [gene.strip() for gene in gene_names]
    #     dm_var_dims["gene_names"] = gene_names
    #     dm_var_dims["gene_dim"] = len(gene_names)
    # elif cfg["data"]["kwargs"]["pert_col"].lower() == "cytokine":
    #     logger.info("Using parse gene names for the baselines ")
    #     path = "/large_storage/ctc/userspace/mohsen/state_revisions/gears_prep/parse_hvg_names.txt"
    #     with open(path, "r") as f:
    #         gene_names = f.readlines()
    #     gene_names = [gene.strip() for gene in gene_names]
    #     dm_var_dims["gene_names"] = gene_names
    #     dm_var_dims["gene_dim"] = len(gene_names)

    # Create model
    model = get_lightning_module(
        cfg["model"]["name"],
        cfg["data"]["kwargs"],
        cfg["model"]["kwargs"],
        cfg["training"],
        dm_var_dims,
    )

    # Set up logging
    loggers = get_loggers(
        output_dir=cfg["output_dir"],
        name=cfg["name"],
        wandb_project="a",
        wandb_entity="a",
        local_wandb_dir="a",
        use_wandb=False,
        cfg=cfg,
    )

    # If using wandb, store the run path in a text file for eval
    # that matches the old train_lightning.py logic
    for lg in loggers:
        if isinstance(lg, WandbLogger):
            wandb_info_path = os.path.join(run_output_dir, "wandb_path.txt")
            with open(wandb_info_path, "w") as f:
                f.write(lg.experiment.path)
            break

    # Set up callbacks
    ckpt_callbacks = get_checkpoint_callbacks(
        cfg["output_dir"],
        cfg["name"],
        8000,
        cfg["training"].get("ckpt_every_n_steps", 4000),
    )
    # Add BatchSpeedMonitorCallback to log batches per second to wandb
    batch_speed_monitor = BatchSpeedMonitorCallback()
    callbacks = ckpt_callbacks + [batch_speed_monitor]

    logger.info("Loggers and callbacks set up.")

    if cfg["model"]["name"].lower().startswith("scgpt"):
        plugins = [
            MixedPrecision(
                precision="bf16-mixed",
                device="cuda",
            )
        ]
    else:
        plugins = []

    if torch.cuda.is_available():
        accelerator = "gpu"
    else:
        accelerator = "cpu"

    # Decide on trainer params
    trainer_kwargs = dict(
        accelerator=accelerator,
        devices=1,
        max_steps=cfg["training"].get("max_steps", -1),  # for normal models
        max_epochs=cfg["training"].get("max_epochs", -1),
        check_val_every_n_epoch=None,
        val_check_interval=9000,
        logger=loggers,
        plugins=plugins,
        callbacks=callbacks,
        gradient_clip_val=cfg["training"].get("gradient_clip_val", None),
    )

    if cfg["model"]["name"].lower() == "cpa":
        trainer_kwargs["gradient_clip_val"] = 0

    # If it's SimpleSum, override to do exactly 1 epoch, ignoring `max_steps`.
    if (
        cfg["model"]["name"].lower() == "celltypemean"
        or cfg["model"]["name"].lower() == "globalsimplesum"
    ):
        trainer_kwargs["max_epochs"] = 1  # do exactly one epoch
        # delete max_steps to avoid conflicts
        del trainer_kwargs["max_steps"]

    # Build trainer
    trainer = pl.Trainer(**trainer_kwargs)

    # Load checkpoint if exists (skip if overwrite is True)
    checkpoint_path = None
    if not cfg.get("overwrite", False):
        checkpoint_path = join(ckpt_callbacks[0].dirpath, "last.ckpt")
        if not exists(checkpoint_path):
            checkpoint_path = None
        else:
            logging.info(f"!! Resuming training from {checkpoint_path} !!")
    else:
        logger.info("Training from scratch (overwrite=True)")

    logger.info("Starting trainer fit.")

    # Train
    trainer.fit(
        model,
        datamodule=dm,
        ckpt_path=checkpoint_path,
    )

    # at this point if checkpoint_path does not exist, manually create one
    checkpoint_path = join(ckpt_callbacks[0].dirpath, "final.ckpt")
    if not exists(checkpoint_path):
        trainer.save_checkpoint(checkpoint_path)

    # save the data_module
    with open(join(run_output_dir, "data_module.pkl"), "wb") as f:
        pickle.dump(dm, f)


if __name__ == "__main__":
    train()
