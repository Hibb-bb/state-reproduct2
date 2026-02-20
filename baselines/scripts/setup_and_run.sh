#!/bin/bash
# Script to install packages and run training for Marson dataset
# Usage: ./setup_and_run.sh [MODEL_NAME] [DATASET_NAME]
# Defaults: MODEL_NAME=cpa, DATASET_NAME=marson
# Example: ./setup_and_run.sh cpa marson

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASELINES_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Parse arguments with defaults for Marson
MODEL_NAME=${1:-cpa}
DATASET_NAME=${2:-marson}
FOLD_ID=$3

echo "=========================================="
echo "Setting up environment and installing packages for Marson"
echo "Model: $MODEL_NAME"
echo "Dataset: $DATASET_NAME"
echo "=========================================="

# Create virtual environment if it doesn't exist
VENV_DIR="$BASELINES_DIR/.venv"
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment at $VENV_DIR..."
    python3 -m venv "$VENV_DIR"
fi

# Activate virtual environment
echo "Activating virtual environment..."
source "$VENV_DIR/bin/activate"

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip setuptools wheel

# Install requirements
echo "Installing requirements from requirements.txt..."
pip install -r "$BASELINES_DIR/requirements.txt"

# Install baselines package in editable mode
echo "Installing baselines package in editable mode..."
cd "$BASELINES_DIR"
pip install -e .

echo ""
echo "=========================================="
echo "Package installation complete!"
echo "=========================================="
echo ""
echo "=========================================="
echo "Starting training..."
echo "=========================================="

# Run the training script
"$BASELINES_DIR/scripts/train.sh" "$MODEL_NAME" "$DATASET_NAME" "$FOLD_ID"

