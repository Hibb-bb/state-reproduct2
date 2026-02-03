#!/bin/bash
# Script to install packages and run training
# Usage: ./scripts/install_and_train.sh [MODEL_NAME] [DATASET_NAME] [FOLD_ID]

set -e

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASELINES_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$BASELINES_DIR"

# Parse arguments
MODEL_NAME=${1:-cpa}
DATASET_NAME=${2:-replogle}
FOLD_ID=${3:-1}

echo "=========================================="
echo "Installing packages and running training"
echo "Model: $MODEL_NAME"
echo "Dataset: $DATASET_NAME"
echo "Fold: $FOLD_ID"
echo "=========================================="

# Create virtual environment if it doesn't exist
if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    uv venv .venv --python 3.12
fi

# Activate virtual environment
echo "Activating virtual environment..."
. .venv/bin/activate

# Install requirements using uv
echo "Installing requirements..."
uv pip install -r requirements.txt

# Install torch-scatter (needed but commented out in requirements.txt)
echo "Installing torch-scatter..."
uv pip install torch-scatter --no-build-isolation

# Install the package itself in editable mode
echo "Installing baselines package..."
uv pip install -e .

echo "=========================================="
echo "Installation complete!"
echo "=========================================="

# Generate .hepg2.toml (always regenerate to ensure it's up to date)
echo "Generating .hepg2.toml file..."
.venv/bin/python rep_toml.py

echo "=========================================="
echo "Starting training..."
echo "=========================================="

# Run the training script
HYDRA_FULL_ERROR=1  ./scripts/train.sh cpa replogle 1 

echo "=========================================="
echo "Training complete!"
echo "=========================================="

