"""
Makes an effic file a detector/nu interaction channel.
 The effic file is 1D array with the length being dependent on the detection bin number.
 Values range from 0 to 1.
"""

import numpy as np
import sys

file_path = '../channels'
chan_name = sys.argv[1]

if chan_name == None:
    print("Need a channel name !!")
    exit()

full_path = f"{file_path}/channels_{chan_name}.dat"


###############################################################################

def IntInput(text, errortext="Invalid input", minval=-1e50, maxval=1e50):
    successful_input = False
    while not successful_input:
        try:
            response = int(input(text + "-->"))
            if response >= minval and response <= maxval:
                successful_input = True
            else:
                print("Value out of bounds, please try again")
        except ValueError:
            print(errortext)
    return response


###############################################################################


def get_detector_name():
    """Asks user to select a detector from CLI.

    Returns
    -------
    detector_name: 'str'
        detector name

    """
    detectors = {
        1: 'wc100kt30prct',
        2: 'wc100kt15prct',
        3: 'ar40kt',
        4: 'scint20kt',
        5: 'halo1',
        6: 'halo2',
        7: 'novaND',
        8: 'novaFD',
        9: 'd2O',
        10: 'wc100kt30prct_he',
        11: 'ar40kt_he',
        12: 'icecube',
        13: 'ds20',
        14: 'argo',
        15: 'lz',
        16: 'xent',
        17: 'pandax'
    }

    print('Avilable detectors (for file name):')

    for key in detectors:
        print(f"{key}) {detectors[key]}")

    dect_input = IntInput(
        "Please input a number from 1 to 17", "Invalid Number!!", 1, 17)

    return detectors[dect_input]


###############################################################################


def add_curlies(file_path):
    """As the name states this methods add { at the beginning and } at the
    end of the effic file. Need {} to initilize post_smearing

    Parameter
    ---------
    file_path : 'str'
        path to channel file
    """
    with open(file_path, 'r+') as f:
        content = f.read()
        f.truncate(0)
        f.seek(0, 0)
        f.write(r"{" + content + r"}")


###############################################################################

def get_detector_effic(detector_name, E_obs):
    """Sets up the energy effic functions for all detectors present in SNOwGLoBES

    Parameters
    ----------
    detector_name: 'str'
        detector name
    E_obs: 'float'
        Observed energy (GeV)

    Returns
    -------
    detector_effic : 'lambda'
        Lambda function for detector effic

    """

    detector_effic = {
        "ar40kt": lambda x: 1.0 if (x > 5.0E-3) else 0,
        # "ds20": not_included(),
        # "argo": not_included(),
        # "wc100kt30prct": not_included(),
        # "wc100kt100prct": not_included(),
        # "scint20kt": not_included(),
        "halo1": lambda x: 0.56,
        "halo2": lambda x: 0.36,
        # "novaND":  not_included(),
        # "novaFD":  not_included(),
        # "d20":  not_included(),
        # "wc100kt30prct_he":not_included(),
        # "icecube":not_included(),
        # "ar40kt_he":not_included(),
        "lz": lambda x: 1.0 if (x > 5.0E-7) else 0,
        "ds20": lambda x: 1.0 if (x > 5.0E-5) else 0,
        "xent": lambda x: 1.0 if (x > 5.0E-7) else 0,
        "pandax": lambda x: 1.0 if (x > 5.0E-7) else 0
    }
    return detector_effic[detector_name](E_obs)


###############################################################################

# for every channel this function will make a 1D array


def make_effic_array(full_path, detector_name):
    """Makes the effic array

    Parameter
    ---------
    full_path : 'str'
        path to channel
    detector_name: 'str'
        detector name
    """
    f = open(full_path)  # open file
    lines = f.readlines()  # read the line

    for line in lines:  # loop through all lines
        line_spl = line.split(' ')  # split line by white space , makes an arr

        if line_spl[0] == 'SN_nu':  # if general ID is used skip it
            line_spl = line_spl[1:]

        chan_name = line_spl[0]  # set chan name

        paras = line_spl[5:]  # make binning parameters array

        out_bins = int(paras[3])
        lo_enrg_out = float(paras[4])  # GeV
        hi_enrg_out = float(paras[5])  # GeV

        # make 1D e_obs array
        out_enrgs = np.linspace(lo_enrg_out, hi_enrg_out, out_bins)

        # apply detector energy res to each e_out bin

        # make empty 1D arr for effics
        effic_arr = np.zeros(out_bins)

        # fill the arr with detector effics
        ind = 0
        for elem in out_enrgs:
            effic_arr[ind] = get_detector_effic(
                detector_name, elem)
            ind = ind + 1
        # Save each effic arr as a dat for each channel in file
        effic_name = f"{detector_name}/effic_{chan_name}_{detector_name}.dat"
        print(f'\nMaking {effic_name}')
        print(f'Bins: {out_bins}')
        print(f'Low Energy Limit: {lo_enrg_out}')
        print(f'High Energy Limit: {hi_enrg_out}')
        effic_arr.tofile(effic_name, sep=",", format="%s")
        # adds curlies at the beginning and end of effic file
        add_curlies(effic_name)


###############################################################################


# Running the methods
detector_name = get_detector_name()
make_effic_array(full_path, detector_name)
