import numpy as np
import sys

file_path = './channels'
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
        13: 'kamland'
    }

    print('Available detectors (for file name):')

    for key in detectors:
        print(f"{key}) {detectors[key]}")

    dect_input = IntInput(
        "Please input a number from 1 to 13", "Invalid Number!!", 1, 13)

    return detectors[dect_input]
###############################################################################


def add_curlies(file_path):
    #As the name states this methods add { at the beginning and } at the
    # end of the effic file. Need {} to initialize post_smearing
    with open(file_path, 'r+') as f:
        content = f.read()
        f.truncate(0)
        f.seek(0, 0)
        f.write(r"{"+content+r"}")


###############################################################################
# This method stores each detectors energy efficiency.
def get_detector_effic(detector_name, E_obs):
    # E_obs in GeV
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
        # These 3 are for SNOwGLoBES 1.3
 #        "lz": lambda x: 1.0 if (x > 5.0E-7) else 0,
	# "ds20": lambda x: 1.0 if (x > 5.0E-5) else 0,
	# "xent": lambda x: 1.0 if (x > 5.0E-7) else 0,
	# "pandax": lambda x: 1.0 if (x > 5.0E-7) else 0
        "kamland_ibd": lambda x: 0.6 if (1.806e-3 < x < 4e-3) else (0.79 if (x >= 4e-3) else 0),
        "kamland": lambda x: 0.79 if (x > 5.0e-7) else 0
    }
    return detector_effic[detector_name](E_obs)
###############################################################################

# for every channel this function will make a 1D array


def make_effic_array(full_path, detector_name):

    f = open(full_path) #open file
    lines = f.readlines() #read the line

    out_bins = 200
    lo_enrg_out = 0.0005
    hi_enrg_out = 0.100

    for line in lines: # loop through all lines
        line_spl = line.split(' ') #split line by white space , makes an arr

        if line_spl[0] == 'SN_nu':  #if general ID is used skip it
            line_spl = line_spl[1:]

        chan_name = line_spl[0] #set chan name


        # paras = line_spl[5:]  #make binning parameters array

        # out_bins = int(paras[3])
        # lo_enrg_out = float(paras[4])  # GeV
        # hi_enrg_out = float(paras[5])  # GeV

#     make 1D e_obs array
        out_enrgs = np.linspace(lo_enrg_out, hi_enrg_out, out_bins)

# apply detector energy res to each e_out bin
        # make empty 1D arr for effics
        effic_arr = np.zeros(out_bins)
        elem = 0
        # fill the arr with detector effics
        if detector_name == 'kamland' and chan_name == 'ibd':
            for x in out_enrgs:
                effic_arr[elem] = get_detector_effic(
                    f'{detector_name}_ibd', out_enrgs[elem])
                elem += 1
        else:
            for x in out_enrgs:
                effic_arr[elem] = get_detector_effic(
                    detector_name, out_enrgs[elem])
                elem = elem + 1
        # Save each effic arr as a dat for each channel in file
        effic_name = f"effic/effic_{chan_name}_{detector_name}.dat"
        print(f'\nMaking {effic_name}')
        print(f'Bins: {out_bins}')
        print(f'Low Energy Limit: {lo_enrg_out}')
        print(f'High Energy Limit: {hi_enrg_out}')
        effic_arr.tofile(effic_name, sep=",", format="%s")
        add_curlies(effic_name)
###############################################################################


# Running the methods
detector_name = get_detector_name()
make_effic_array(full_path, detector_name)
