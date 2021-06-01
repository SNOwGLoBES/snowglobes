"""
from_sntools.py

Script to generate a GLoBES cross section file from a cross section implemented
in sntools (https://github.com/JostMigenda/sntools).
To use this script:
  1) Install sntools. (See link above for instructions.)
  2) If necessary, adjust the lines marked with a "TODO" comment below
  3) Run this script. (Python 3.7 or later recommended.)
"""
from math import log10
from numpy import linspace
from scipy import integrate
from sntools.interaction_channels import o16eb as c  # TODO: select interaction channel


def xs(eNu):
    eE_min, eE_max = c.bounds_eE(eNu)
    if eNu >= c.bounds_eNu[0]:
        xs_in_natural_units = integrate.quad(lambda eE: c.dSigma_dE(eNu, eE), eE_min, eE_max, points=c._opts(eNu)["points"])[0]
    else:
        xs_in_natural_units = 0
    cm2mev = 5.067731E10  # conversion factor derived from hbar*c = 197 MeV*fm = 1
    xs_in_cm2 = xs_in_natural_units / cm2mev**2
    return xs_in_cm2


def print_xs(eNu):
    e_in_GeV = eNu / 1000
    return xs(eNu) * 1e38 / e_in_GeV


with open('xscns/xs_nuebar_O16_Suzuki2018.dat', 'w') as outfile:  # TODO: update file name
    outfile.write("# Electron antineutrino-oxygen charged-current cross section 5MeV-100MeV (10^-38 cm^2/GeV)\n")  # TODO: update description
    outfile.write("# Based on arXiv:1807.02367 (calculations) and arXiv:1809.08398 (fit).\n")  # TODO: update reference
    outfile.write("# log(energy in GeV)   nu_e       nu_mu      nu_tau       nu_e_bar    nu_mu_bar  nu_tau_bar\n\n")

    for log10eNu in linspace(log10(0.005), log10(0.100), num=1001):
        eNu = 10**log10eNu * 1000  # convert log10(GeV) to MeV
        try:
            val = print_xs(eNu)
        except:
            val = 0  # if eNu is below energy threshold of the interaction
        outfile.write(f"{log10eNu:^+20.6f} {0:^10} {0:^10} {0:^10} {val:^15.5e} {0:^10} {0:^10}\n")  # TODO: swap columns if necessary
