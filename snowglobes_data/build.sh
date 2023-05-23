#! /bin/bash

# Set up directory structure for a Python package
mkdir -p src/snowglobes_data
cp -r ../LICENSE ../detector_configurations.dat ../channels ../effic ../smear ../xscns src/snowglobes_data/

# Get version number from environment variable (in GitHub Action); fall back to latest git tag for local installs
if [ ! $GIT_VERSION ]; then
  GIT_VERSION=`git describe --dirty | sed 's/v//'`;
fi
echo "__version__ = '$GIT_VERSION'" > src/snowglobes_data/__init__.py

# Create __init__ files for submodules. Without these, code like `files(snowglobes_data.channels)` raises an error in importlib.resources under Python 3.9
# See https://setuptools.pypa.io/en/latest/userguide/datafiles.html#accessing-data-files-at-runtime
touch src/snowglobes_data/__init__.py src/snowglobes_data/channels/__init__.py src/snowglobes_data/effic/__init__.py src/snowglobes_data/smear/__init__.py src/snowglobes_data/xscns/__init__.py

# Build package using configurations in pyproject.toml
pip install --upgrade build
python -m build
