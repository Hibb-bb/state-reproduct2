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

MODEL_NAME=$1
DATASET_NAME=$2
FOLD_ID=$3
if [ $# -eq 4 ]; then
    CKPT=$4
else
    CKPT=""
fi

# Define output directory (matching train.sh)
OUTPUT_DIR_BASE="outputs"

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
fi

echo "Generating Predictions for $MODEL_NAME on $DATASET_NAME (fold: $FOLD_ID)"
echo "Output directory: $OUTPUT_DIR"

$PYTHON_CMD -m state_sets_reproduce.train.get_predictions \
    --output_dir ${OUTPUT_DIR} \
    --checkpoint ${CKPT} 
