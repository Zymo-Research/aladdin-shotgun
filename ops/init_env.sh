#!/bin/bash

SCRIPT_DIR="./ops"
CONDA_PREFIX="./.conda"

source "$SCRIPT_DIR/conda_helpers.sh"


if ! command -v conda &> /dev/null
then
    echo
    echo "ERROR: 'conda' command could not be found."
    echo
    exit
fi

echo

echo "INITIALISING THE CONDA ENVIRONMENT..."
echo
conda_env_init "$CONDA_PREFIX" "$SCRIPT_DIR" "$1"

echo "INSTALLING PRE-COMMIT HOOKS FOR GIT..."
conda run --live-stream --prefix "$CONDA_PREFIX" pre-commit install -t pre-commit -t commit-msg
echo
