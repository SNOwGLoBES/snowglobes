# Cross section scripts
A collection of various Python scripts to handle cross section files.

## Installation
All scripts are tested under Python 3.x.

Some scripts require the following additional libraries:
* [matplotlib](https://matplotlib.org/stable/users/installing.html)
* [NumPy](https://numpy.org/install/)
* [PyROOT](https://root.cern/manual/python/)
* [SciPy](https://scipy.org/install.html)

## Usage information

### compare_xs.py
Plot multiple cross section files.
Run `python compare_xs.py` for usage information.

### compute_nu_e_cross_sections.py
Generate file with cross section for neutrino-electron scattering, based on [DOI:10.1088/0954-3899/29/11/013](https://doi.org/10.1088/0954-3899/29/11/013).

### universal_xs.py
Produce a cross section file for an arbitrary A, Z, and Q. This gives a rough approximation and may only be accurate to within a factor of a few.
Run `python universal_xs.py` for usage information.

The folder `parameterization/` contains scripts used to generate this universal XS fit.