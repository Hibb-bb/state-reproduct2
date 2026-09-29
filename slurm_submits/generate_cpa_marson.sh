#!/bin/bash
# CPA generation on Marson: perturbed test cells, then the NTC pool.
# Both passes use the same basal pool as training, matched by donor and timepoint.
#
# The job may be submitted before training reaches $CKPT. It waits for a stable
# checkpoint before starting inference.
#
# Both passes must come from the same checkpoint, since the delta metric subtracts
# the generated control mean from the generated perturbed mean.
#
#SBATCH --job-name=cpa_marson_gen
#SBATCH --partition=h100-reserved
#SBATCH --qos=mid
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --time=06:00:00
#SBATCH --open-mode=append
#SBATCH --output=/mnt/main0/home/sudarshan.babu/cell_bench_runs/cpa_state_ntcpool/logs/gen_%j.out
#SBATCH --error=/mnt/main0/home/sudarshan.babu/cell_bench_runs/cpa_state_ntcpool/logs/gen_%j.out

set -euo pipefail

REPO=/mnt/main0/home/sudarshan.babu/state_reproduce_marson/baselines
OUT_ROOT=/mnt/main0/home/sudarshan.babu/cell_bench_runs/cpa_state_ntcpool

CKPT="${CKPT:-step=200000.ckpt}"
CKPT_PATH="$OUT_ROOT/cpa_marson/marson/checkpoints/$CKPT"

mkdir -p "$OUT_ROOT/logs"

export CELL_LOAD_REPO=/mnt/main0/home/sudarshan.babu/cell-load
export OUTPUT_DIR_BASE="$OUT_ROOT"

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export HYDRA_FULL_ERROR=1
export TOKENIZERS_PARALLELISM=false

export WANDB_MODE="${WANDB_MODE:-offline}"
export WANDB_DIR="$OUT_ROOT/wandb"
export WANDB_CACHE_DIR="$OUT_ROOT/wandb_cache"
mkdir -p "$WANDB_DIR" "$WANDB_CACHE_DIR"

cd "$REPO"

echo "job=$SLURM_JOB_ID node=$SLURMD_NODENAME ckpt=$CKPT start=$(date -Is)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# Wait for training to finish writing the selected checkpoint.
until [ -f "$CKPT_PATH" ]; do sleep 60; done
prev=-1
while size=$(stat -c %s "$CKPT_PATH") && [ "$size" != "$prev" ]; do
    prev=$size
    sleep 60
done
echo "checkpoint ready: $CKPT_PATH ($prev bytes) at $(date -Is)"

GEN_ROOT="$OUT_ROOT/cpa_marson/marson"

# Each pass is skipped if it already produced a prediction file, so a rerun after one
# pass fails does not redo several hours of work on 2.4M test cells.
run_pass() {
    local label=$1 dataset=$2 dir=$3
    if [ -f "$GEN_ROOT/$dir/adata_pred.h5ad" ]; then
        echo "=== $label: already present at $dir/adata_pred.h5ad, skipping ==="
        return 0
    fi
    echo "=== $label ==="
    srun --unbuffered bash scripts/predict.sh cpa "$dataset" "$dataset" "$CKPT"
}

# Perturbed test cells (162 test_hvg shards) -> generation/
run_pass "pass 1/2: perturbed (marson)" marson generation

# Shared NTC pool -> generation_controls/
run_pass "pass 2/2: controls (marson_controls)" marson_controls generation_controls

echo "done=$(date -Is)"
