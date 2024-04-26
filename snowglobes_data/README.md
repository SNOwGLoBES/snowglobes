# snowglobes_data

Module to make [SNOwGLoBES](https://github.com/SNOwGLoBES/snowglobes) data files available as a PyPI package.

This lets any Python package install SNOwGLoBES data files as a dependency via `pip` and access them via the [`importlib.resources`](https://docs.python.org/3/library/importlib.resources.html) module.

## Usage

Install this package manually using
```bash
pip install snowglobes_data
```
or automatically by listing it as a dependency e.g. in the `pyproject.toml` or `requirements.txt` file of another Python package.

You can then access the data files without cloning the SNOwGLoBES repo and manually hard-coding its path in the script:
```Python
from importlib.resources import files  # on Python 3.9 and later
#from importlib_resources import files  # on Python 3.8 or earlier
import snowglobes_data

# Use `files(snowglobes_data)` instead of the hard-coded path to the SNOwGLoBES repo:
with open(files(snowglobes_data).joinpath("detector_configurations.dat")) as detectors:
    for line in detectors:
        print(line.strip())
```

Note: `files` is available since Python 3.9. If you are using an older Python version, please install and use the backport [`importlib_resources`](https://pypi.org/project/importlib-resources/) instead.