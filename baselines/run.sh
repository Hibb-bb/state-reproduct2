#!/bin/bash
#SBATCH --account=p32234
#SBATCH --job-name=cpa
#SBATCH --nodes=1
#SBATCH --partition=gengpu
#SBATCH --gres=gpu:a100:1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=10:00:00
#SBATCH --output=slurm-cpa.out
#SBATCH --error=slurm-cpa.err

set -e

# Get the directory where this script is located
# Use $0 for sh/bash compatibility
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Create virtual environment if it doesn't exist
if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    uv venv .venv --python 3.12
fi

# Activate virtual environment (use . instead of source for sh compatibility)
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

# Download replogle data if it doesn't exist
if [ ! -f "/mnt/experiments/cpa/replogle/replogle.h5ad" ]; then
    echo "Downloading replogle dataset..."
    .venv/bin/python download_rep.py
fi

# Generate .hepg2.toml if it doesn't exist (needed for replogle fold 1)
if [ ! -f ".hepg2.toml" ]; then
    echo "Generating .hepg2.toml file..."
    .venv/bin/python rep_toml.py
fi

# Run the training script
HYDRA_FULL_ERROR=1 sh scripts/train.sh cpa replogle 1

