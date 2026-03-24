#!/usr/bin/env python3
"""
Sanity-check outputs from Marson *control* generation (predict.sh marson_controls).

Expects get_predictions to have written under --generation_dir:
  adata_pred.h5ad, adata_real.h5ad (and optional *_gene.h5ad).

Example:
  python scripts/check_control_generation_legit.py \\
    /mnt/experiments/cpa/cpa_marson/marson/generation_controls
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np


def _load(path: Path) -> ad.AnnData:
    if not path.is_file():
        raise FileNotFoundError(path)
    return ad.read_h5ad(path)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "generation_dir",
        type=Path,
        help="Directory containing adata_pred.h5ad / adata_real.h5ad",
    )
    p.add_argument(
        "--control-token",
        default="NTC",
        help="Expected perturbation label for controls in obs['pert_name'] (if present).",
    )
    args = p.parse_args()
    gdir: Path = args.generation_dir

    pred_p = gdir / "adata_pred.h5ad"
    real_p = gdir / "adata_real.h5ad"

    try:
        adata_p = _load(pred_p)
        adata_r = _load(real_p)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    if adata_p.n_obs != adata_r.n_obs:
        print(
            f"ERROR: n_obs mismatch pred={adata_p.n_obs} real={adata_r.n_obs}",
            file=sys.stderr,
        )
        return 1

    if adata_p.obs_names.tolist() != adata_r.obs_names.tolist():
        print("ERROR: obs_names differ between pred and real", file=sys.stderr)
        return 1

    if adata_p.shape != adata_r.shape:
        print(
            f"ERROR: shape mismatch pred={adata_p.shape} real={adata_r.shape}",
            file=sys.stderr,
        )
        return 1

    if np.isnan(adata_p.X.data).any() if hasattr(adata_p.X, "data") else np.isnan(adata_p.X).any():
        print("ERROR: NaNs in adata_pred.X", file=sys.stderr)
        return 1
    if np.isnan(adata_r.X.data).any() if hasattr(adata_r.X, "data") else np.isnan(adata_r.X).any():
        print("ERROR: NaNs in adata_real.X", file=sys.stderr)
        return 1

    if "pert_name" in adata_p.obs.columns:
        uniques = adata_p.obs["pert_name"].astype(str).unique()
        if len(uniques) > 1:
            print(
                f"WARN: multiple pert_name values in pred obs (expected single control): {sorted(uniques)[:20]}...",
                file=sys.stderr,
            )
        if args.control_token not in set(uniques):
            print(
                f"WARN: control token {args.control_token!r} not in pert_name uniques: {sorted(uniques)[:20]}",
                file=sys.stderr,
            )

    for optional in ("adata_pred_gene.h5ad", "adata_real_gene.h5ad"):
        op = gdir / optional
        if op.is_file():
            aux = _load(op)
            if aux.n_obs != adata_p.n_obs:
                print(
                    f"ERROR: {optional} n_obs={aux.n_obs} != pred n_obs={adata_p.n_obs}",
                    file=sys.stderr,
                )
                return 1

    print(
        f"OK: {gdir} — n_obs={adata_p.n_obs}, n_vars={adata_p.n_vars}, "
        f"obs columns: {list(adata_p.obs.columns)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
