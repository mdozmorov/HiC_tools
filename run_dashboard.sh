#!/usr/bin/env bash
# Run the Hi-C Tools Dashboard locally

set -e

if command -v uv &>/dev/null; then
    echo "Using uv..."
    [ -d ".venv" ] || uv venv
    source .venv/bin/activate
    uv pip install -r requirements.txt
else
    echo "Using pip..."
    pip install -r requirements.txt
fi

streamlit run dashboard.py
