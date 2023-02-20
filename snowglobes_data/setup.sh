#! /bin/bash

# Set up directory structure for a Python package
mkdir -p src/snowglobes_data
cp -r ../LICENSE ../detector_configurations.dat ../channels ../effic ../smear ../xscns src/snowglobes_data/

# Build package using configurations in pyproject.toml
pip install --upgrade build
python -m build
