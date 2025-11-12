#!/bin/bash
#
# Fine-Grained Binder Design Pipeline
# ====================================
#
# Simple wrapper for the binder optimization pipeline
#
# Usage:
#   ./optimize_binder.sh config/cxcr4_optimization.yaml
#

set -e

# Check if config file is provided
if [ $# -eq 0 ]; then
    echo "ERROR: No configuration file provided"
    echo ""
    echo "Usage: $0 <config.yaml>"
    echo ""
    echo "Examples:"
    echo "  $0 configs/cxcr4_optimization_example.yaml"
    echo ""
    exit 1
fi

CONFIG=$1

# Check if config exists
if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Configuration file not found: $CONFIG"
    exit 1
fi

# Run the pipeline
python3 scripts/optimize_binder.py "$CONFIG"
