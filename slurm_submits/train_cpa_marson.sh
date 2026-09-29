#!/bin/bash
# CPA on Marson using the shared NTC basal pool from marson.toml, matched by donor
# and timepoint. Streams the train_hvg shards through cell-load.
#
# Equivalent of the old:  sh /work/baselines/scripts/setup_and_run.sh
# minus the install step, since the venv is already built at baselines/.venv.
#
#SBATCH --job-name=cpa_marson
#SBATCH --partition=h100-reserved
#SBATCH --qos=mid
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=384G
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/mnt/main0/home/sudarshan.babu/cell_bench_runs/cpa_state_ntcpool/logs/train_%j.out
#SBATCH --error=/mnt/main0/home/sudarshan.babu/cell_bench_runs/cpa_state_ntcpool/logs/train_%j.out

set -euo pipefail

REPO=/mnt/main0/home/sudarshan.babu/state_reproduce_marson/baselines
OUT_ROOT=/mnt/main0/home/sudarshan.babu/cell_bench_runs/cpa_state_ntcpool

mkdir -p "$OUT_ROOT/logs"

# train.sh reads these; without CELL_LOAD_REPO it silently falls back to the PyPI
# cell-load, which lacks the `barcode` kwarg the Marson invocation passes.
export CELL_LOAD_REPO=/mnt/main0/home/sudarshan.babu/cell-load
export OUTPUT_DIR_BASE="$OUT_ROOT"

# overwrite=false so a requeued job resumes from the newest checkpoint rather than
# restarting; find_resume_checkpoint falls back to step=*.ckpt if last.ckpt is missing.
export OVERWRITE=false
export NUM_WORKERS=24

# Leave TMPDIR at the default /tmp, which is a 28 TB local xfs disk. Pointing it at
# the NFS output volume OOM-killed job 3610063: multiprocessing puts its shared-memory
# heap arenas under TMPDIR and unlinks them while still mapped, and dirty pages of an
# unlinked NFS file cannot be reclaimed, so RSS climbed to 402 GB over 171 arenas.

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export HYDRA_FULL_ERROR=1
export TOKENIZERS_PARALLELISM=false

# W&B is off by default here; the harness takes use_wandb from its config, and
# train.sh does not thread a flag through, so leave the key unset to stay offline.
export WANDB_MODE="${WANDB_MODE:-offline}"
export WANDB_DIR="$OUT_ROOT/wandb"
export WANDB_CACHE_DIR="$OUT_ROOT/wandb_cache"
mkdir -p "$WANDB_DIR" "$WANDB_CACHE_DIR"

cd "$REPO"

echo "job=$SLURM_JOB_ID node=$SLURMD_NODENAME start=$(date -Is)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

srun --unbuffered bash scripts/train.sh cpa marson marson
