#!/bin/bash
# Install environment and run Marson *control* generation (marson_generation_controls.toml → ctrl_for_gen_hvg).
# Usage: ./scripts/install_and_gen_control.sh [CHECKPOINT]
#   CHECKPOINT: optional, e.g. step=8000.ckpt or absolute path (default in predict.sh: step=196000.ckpt)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASELINES_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$BASELINES_DIR"

if [ $# -ge 1 ]; then
    export CHECKPOINT="$1"
fi

echo "=========================================="
echo "Installing environment and running Marson control generation"
echo "Data TOML: ${BASELINES_DIR}/marson_generation_controls.toml"
echo "Output dir: /mnt/experiments/cpa/cpa_marson/marson/generation_controls/"
if [ -n "${CHECKPOINT:-}" ]; then
    echo "Checkpoint: ${CHECKPOINT}"
fi
echo "=========================================="

if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    uv venv .venv --python 3.12
fi

echo "Activating virtual environment..."
. .venv/bin/activate

echo "Installing requirements..."
uv pip install -r requirements.txt \
    --extra-index-url https://download.pytorch.org/whl/cu128 \
    --index-strategy unsafe-best-match

echo "Installing torch-scatter from source..."
uv pip install git+https://github.com/rusty1s/pytorch_scatter.git

echo "Installing baselines package..."
uv pip install -e .

echo "=========================================="
echo "Installation complete!"
echo "=========================================="

echo "=========================================="
echo "Running control generation (predict.sh marson_controls)..."
echo "=========================================="

bash scripts/predict.sh cpa marson_controls marson_controls ${CHECKPOINT:+$CHECKPOINT}

echo "=========================================="
echo "Control generation complete!"
echo "=========================================="
