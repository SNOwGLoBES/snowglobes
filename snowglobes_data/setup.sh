#! /bin/bash

# Set up directory structure for a Python package
mkdir -p src/snowglobes_data
cp -r ../LICENSE ../detector_configurations.dat ../channels ../effic ../smear ../xscns src/snowglobes_data/

# Get version number from latest git tag
VERSION=`git describe --abbrev=0 | sed 's/v//'`
echo "__version__ = '$VERSION'" > src/snowglobes_data/__init__.py

# Create __init__ files for submodules. Without these, code like `files(snowglobes_data.channels)` raises an error in importlib.resources under Python 3.9
# See https://setuptools.pypa.io/en/latest/userguide/datafiles.html#accessing-data-files-at-runtime
touch src/snowglobes_data/__init__.py src/snowglobes_data/channels/__init__.py src/snowglobes_data/effic/__init__.py src/snowglobes_data/smear/__init__.py src/snowglobes_data/xscns/__init__.py

# Build package using configurations in pyproject.toml
pip install --upgrade build
python -m build
