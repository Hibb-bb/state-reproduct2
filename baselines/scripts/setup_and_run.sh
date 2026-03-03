#!/usr/bin/env bash
# Robusted setup_and_run.sh: installs requirements, makes sure pip exists,
# installs baselines editable (may upgrade torch), then installs torch_scatter
# matching the *final* torch in the venv (wheel -> fallback build-from-source),
# then runs training.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASELINES_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Optional: remove existing venv for a clean reinstall
if [ "${1:-}" = "--clean" ]; then
    shift
    VENV_DIR="$BASELINES_DIR/.venv"
    if [ -d "$VENV_DIR" ]; then
        echo "Removing existing virtual environment at $VENV_DIR..."
        rm -rf "$VENV_DIR"
    fi
fi

MODEL_NAME=${1:-cpa}
DATASET_NAME=${2:-marson}
FOLD_ID=${3:-}

echo "=========================================="
echo "Setting up environment and installing packages for Marson"
echo "Model: $MODEL_NAME"
echo "Dataset: $DATASET_NAME"
echo "=========================================="

# find uv
UV=""
if command -v uv >/dev/null 2>&1; then
    UV=uv
elif [ -n "${UV_HOME:-}" ] && [ -x "$UV_HOME/uv" ]; then
    UV="$UV_HOME/uv"
elif [ -x "$HOME/.local/bin/uv" ]; then
    UV="$HOME/.local/bin/uv"
fi
if [ -z "$UV" ]; then
    echo "Error: uv not found. Install it with: curl -LsSf https://astral.sh/uv/install.sh | sh" >&2
    exit 1
fi

VENV_DIR="$BASELINES_DIR/.venv"
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment at $VENV_DIR with uv (Python 3.12)..."
    "$UV" venv "$VENV_DIR" --python 3.12
fi

echo "Activating virtual environment..."
# shellcheck disable=SC1091
. "$VENV_DIR/bin/activate"

# Bootstrap pip if missing
if ! python -m pip --version >/dev/null 2>&1; then
    echo "pip not found in venv — attempting to bootstrap pip with ensurepip..."
    if python -m ensurepip --upgrade >/dev/null 2>&1; then
        echo "Bootstrapped pip via ensurepip."
    else
        echo "ensurepip unavailable; attempting to bootstrap pip via get-pip.py..."
        TMP_GET_PIP="/tmp/get-pip-$$.py"
        if command -v curl >/dev/null 2>&1; then
            curl -sS https://bootstrap.pypa.io/get-pip.py -o "$TMP_GET_PIP"
        elif command -v wget >/dev/null 2>&1; then
            wget -q -O "$TMP_GET_PIP" https://bootstrap.pypa.io/get-pip.py
        else
            echo "Neither curl nor wget available to fetch get-pip.py. Please install pip manually." >&2
            exit 1
        fi
        python "$TMP_GET_PIP"
        rm -f "$TMP_GET_PIP"
    fi
fi

# Ensure pip/setuptools/wheel up-to-date
python -m pip install --upgrade pip setuptools wheel

# pip wrapper
pip_cmd() { python -m pip "$@"; }

# Install requirements EXCLUDING torch_scatter (we'll install torch_scatter later to match final torch)
# If your requirements.txt explicitly pins torch_scatter, consider removing it there; here we filter it out.
REQ_FILE="$BASELINES_DIR/requirements.txt"
TMP_REQ="/tmp/req-$$.txt"
# Filter out lines that reference torch-scatter (case-insensitive)
grep -i -vE '^\s*#|torch[-_.]?scatter' "$REQ_FILE" > "$TMP_REQ" || true

echo "Installing requirements (excluding torch_scatter) from requirements.txt..."
pip_cmd install -r "$TMP_REQ"

# Install baselines in editable mode (this may upgrade torch to baselines' required version)
echo "Installing baselines package in editable mode (may upgrade torch)..."
cd "$BASELINES_DIR"
pip_cmd install -e .

# At this point torch (and CUDA wheels) may have been upgraded.
# Install torch_scatter now — it must match the *installed* torch.

echo "Installing torch_scatter matching the installed torch..."

# remove any existing torch_scatter to avoid mismatched binary left behind
pip_cmd uninstall -y torch-scatter || true
set +e
pip_cmd cache purge || true
set -e

# Detect installed torch version and cuda
PYTORCH_BASE_VER="$(python - <<'PY'
import torch, re
v = torch.__version__
base = re.split(r'\+', v)[0]
print(base)
PY
)"

PYTORCH_CUDA_RAW="$(python - <<'PY'
import torch
c = getattr(torch.version, "cuda", None)
print("" if c is None else c)
PY
)"

if [ -n "$PYTORCH_CUDA_RAW" ]; then
    CUDA_NODOT="$(echo "$PYTORCH_CUDA_RAW" | tr -d '.')"
    PYG_CUDA_STR="cu${CUDA_NODOT}"
else
    PYG_CUDA_STR="cpu"
fi

echo "Detected torch version: ${PYTORCH_BASE_VER}"
if [ -n "$PYTORCH_CUDA_RAW" ]; then
    echo "Detected torch CUDA version: ${PYTORCH_CUDA_RAW} -> wheel tag: ${PYG_CUDA_STR}"
else
    echo "No CUDA detected in torch (assuming CPU build)."
fi

WHEEL_INDEX_URL="https://data.pyg.org/whl/torch-${PYTORCH_BASE_VER}+${PYG_CUDA_STR}.html"
echo "Attempting to install prebuilt wheel from: ${WHEEL_INDEX_URL}"

set +e
pip_cmd install --no-cache-dir --only-binary torch-scatter -f "${WHEEL_INDEX_URL}" torch-scatter
INSTALL_EXIT=$?
set -e

if [ $INSTALL_EXIT -eq 0 ]; then
    echo "Successfully installed torch-scatter from wheel index."
else
    echo "No prebuilt wheel available (exit $INSTALL_EXIT). Building torch_scatter from source against installed torch..."

    # help cmake find torch
    CMAKE_PREFIX_PATH="$(python - <<'PY'
import torch
print(torch.utils.cmake_prefix_path)
PY
)"
    export CMAKE_PREFIX_PATH
    echo "Set CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}"

    if [ -z "${CUDA_HOME:-}" ] && [ -d "/usr/local/cuda" ]; then
        export CUDA_HOME="/usr/local/cuda"
        echo "Set CUDA_HOME=${CUDA_HOME}"
    fi

    pip_cmd install --no-build-isolation --verbose torch_scatter
fi

# cleanup TMP_REQ
rm -f "$TMP_REQ" || true

echo "Installing any remaining local editable requirements (if needed)..."
# (optional) re-run any post-editable installs if your workflow needs them

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
_NAME" "$FOLD_ID"