#!/bin/bash

# Detect and use Python from virtual environment if available
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASELINES_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
VENV_PYTHON="$BASELINES_DIR/.venv/bin/python"

if [ -f "$VENV_PYTHON" ]; then
    PYTHON_CMD="$VENV_PYTHON"
elif command -v python3 >/dev/null 2>&1; then
    PYTHON_CMD="python3"
elif command -v python >/dev/null 2>&1; then
    PYTHON_CMD="python"
else
    echo "Error: Python not found. Please ensure Python is installed or activate the virtual environment."
    exit 1
fi

# Add local cell-load repo to PYTHONPATH to use local version instead of installed package
CELL_LOAD_REPO="/mnt/sudarshan/cell-load"
if [ -d "$CELL_LOAD_REPO/src" ]; then
    export PYTHONPATH="$CELL_LOAD_REPO/src:$PYTHONPATH"
    echo "Using local cell-load repo from: $CELL_LOAD_REPO/src"
else
    echo "Warning: Local cell-load repo not found at $CELL_LOAD_REPO/src, using installed package"
fi

MODEL_NAME=$1
DATASET_NAME=$2
FOLD_ID=$3
if [ $# -eq 4 ]; then
    CKPT=$4
else
    CKPT=""
fi

# Define output directory (matching train.sh)
OUTPUT_DIR_BASE="/mnt/experiments/cpa"

# Define test tasks for each fold (matching train.sh structure)
if [ "$DATASET_NAME" = "replogle" ]; then
    OUTPUT_DIR="${OUTPUT_DIR_BASE}/${MODEL_NAME}_replogle_v2/fold${FOLD_ID}/"
    if [ -z "$CKPT" ]; then
        CKPT="final.ckpt"
    fi
elif [ "$DATASET_NAME" = "tahoe" ]; then
    OUTPUT_DIR="${OUTPUT_DIR_BASE}/${MODEL_NAME}_tahoe/tahoe_generalization/"
    if [ -z "$CKPT" ]; then
        CKPT="last.ckpt"
    fi
    if [ "$MODEL_NAME" = "lrlm" ]; then
        CKPT="final.ckpt"
    fi
elif [ "$DATASET_NAME" = "parse" ]; then
    OUTPUT_DIR="${OUTPUT_DIR_BASE}/${MODEL_NAME}_parse/${FOLD_ID}/"
    if [ -z "$CKPT" ]; then
        CKPT="last.ckpt"
    fi
elif [ "$DATASET_NAME" = "xaira" ]; then
    OUTPUT_DIR="${OUTPUT_DIR_BASE}/${MODEL_NAME}_xaira/${FOLD_ID}/"
    if [ -z "$CKPT" ]; then
        CKPT="final.ckpt"
    fi
elif [ "$DATASET_NAME" = "marson" ]; then
    OUTPUT_DIR="${OUTPUT_DIR_BASE}/${MODEL_NAME}_marson/marson/"
    GENERATION_DIR="${OUTPUT_DIR}generation"
    DATA_TOML="${BASELINES_DIR}/marson_generation.toml"
    if [ -z "$CKPT" ]; then
        CKPT="last.ckpt"
    fi
fi

echo "Generating Predictions for $MODEL_NAME on $DATASET_NAME (fold: $FOLD_ID)"
echo "Output directory: $OUTPUT_DIR"

PREDICT_ARGS="--output_dir ${OUTPUT_DIR} --checkpoint ${CKPT}"
[ -n "${GENERATION_DIR:-}" ] && PREDICT_ARGS="${PREDICT_ARGS} --generation_dir ${GENERATION_DIR}"
[ -n "${DATA_TOML:-}" ] && PREDICT_ARGS="${PREDICT_ARGS} --data_toml ${DATA_TOML}"
$PYTHON_CMD -m state_sets_reproduce.train.get_predictions $PREDICT_ARGS 
