#!/bin/bash
# Activate boltzgen environment and run FCGR4 design workflow

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Activate boltzgen environment
# Adjust this path based on your conda/virtualenv setup
if [ -d "$CONDA_PREFIX" ] && [ -n "$CONDA_PREFIX" ]; then
    echo "Using existing conda environment: $CONDA_DEFAULT_ENV"
else
    # Try to activate boltzgen conda environment
    if command -v conda &> /dev/null; then
        echo "Activating boltzgen conda environment..."
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate boltzgen || {
            echo "Warning: Could not activate boltzgen conda environment"
            echo "Please activate it manually: conda activate boltzgen"
        }
    else
        echo "Warning: conda not found. Please activate boltzgen environment manually."
    fi
fi

# Run the workflow script
echo "Running FCGR4 design workflow..."
python3 fcgr4_design_workflow.py



