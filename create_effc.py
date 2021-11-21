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
        13: 'kamland',
        14: 'kamland_le'
    }

    print('Available detectors (for file name):')

    for key in detectors:
        print(f"{key}) {detectors[key]}")

    dect_input = IntInput(
        "Please input a number from 1 to 14", "Invalid Number!!", 1, 14)

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
# This data comes from https://arxiv.org/abs/1506.01175
kl_cc_energy = np.array([1.75557861, 1.77133033, 1.78532536, 1.7993204 , 1.81331504,
       1.82879871, 1.84382111, 1.85870413, 1.8757922 , 1.89333648,
       1.91106761, 1.92754695, 1.94630165, 1.97805824, 2.01519875,
       2.05065292, 2.0927714 , 2.13407279, 2.16572011, 2.18553525,
       2.21304774, 2.24138679, 2.26925419, 2.29501132, 2.32288327,
       2.35044618, 2.38998127, 2.42876346, 2.46754353, 2.50632359,
       2.54510366, 2.58387991, 2.62265573, 2.66224469, 2.69972287,
       2.73819022, 2.77769281, 2.81645293, 2.85522876, 2.89405254,
       2.93287929, 2.97172938, 3.01059051, 3.04944528, 3.08826311,
       3.12707459, 3.16588818, 3.20470348, 3.24683902, 3.27914046,
       3.29993184, 3.33222986, 3.35161604, 3.38757575, 3.42371778,
       3.45925374, 3.49802108, 3.53754607, 3.56571596, 3.58472374,
       3.60108523, 3.61880788, 3.63743983, 3.65414153, 3.67214665,
       3.6892084 , 3.70640902, 3.72370773, 3.74029177, 3.76212924,
       3.7901366 , 3.82055471, 3.84856443, 3.88808215, 3.9268601 ,
       3.9656495 , 4.00445291, 4.04325547, 4.08204912, 4.12084107,
       4.15963599, 4.19843686, 4.23723772, 4.27605259, 4.31487   ,
       4.35369124, 4.3925095 , 4.43132861, 4.47014985, 4.50896726,
       4.53543646])
kl_cc_effic = np.array([0.5772348 , 0.58716656, 0.59648842, 0.60581028, 0.61516172,
       0.6248322 , 0.63418659, 0.6423435 , 0.65044004, 0.65825462,
       0.6656233 , 0.67355811, 0.68198578, 0.68835244, 0.69555895,
       0.70127899, 0.70121579, 0.69859874, 0.69425461, 0.6871762 ,
       0.68020618, 0.67231839, 0.66451029, 0.65614246, 0.6479883 ,
       0.64160873, 0.64501014, 0.64815316, 0.65145751, 0.65476187,
       0.65806622, 0.66166098, 0.665288  , 0.66901707, 0.67360025,
       0.67895761, 0.68482975, 0.68965065, 0.69327768, 0.69325853,
       0.69301352, 0.69099383, 0.68813519, 0.68576056, 0.68619315,
       0.68710975, 0.68786501, 0.68849121, 0.68713087, 0.68210167,
       0.6756015 , 0.67083258, 0.66520927, 0.66178515, 0.65609807,
       0.65560027, 0.65987263, 0.66404153, 0.67114419, 0.68056211,
       0.68836458, 0.69588668, 0.70406332, 0.71204832, 0.72076307,
       0.72932539, 0.73771112, 0.74630648, 0.75557006, 0.76281564,
       0.76999844, 0.77767122, 0.78467357, 0.78939562, 0.79286131,
       0.79545579, 0.79698546, 0.79857967, 0.80085148, 0.80325236,
       0.80542737, 0.80715064, 0.80887391, 0.80953238, 0.80999724,
       0.8101717 , 0.81057202, 0.81090782, 0.81108228, 0.81154714,
       0.81165289])

def kamland_cc_effic(x):
    if x < 0.00018:
        return 0
    if x > 0.00043:
        return 0.81
    return np.interp(x, kl_cc_energy, kl_cc_effic)

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
        "kamland_cc": kamland_cc_effic,
        "kamland": lambda x: 1 if (x > 0.00035) else 0,
        "kamland_le_cc": kamland_cc_effic,
        "kamland_le": lambda x: 1 if (x > 0.00035) else 0
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

        if line_spl[0] == '%':  # custom binning
            out_bins = int(line_spl[1])
            lo_enrg_out = float(line_spl[2])
            hi_enrg_out = float(line_spl[3])
            continue

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
        if detector_name in ['kamland', 'kamland_le'] and\
            chan_name in ['ibd', 'nue_C12', 'nue_C13']:
            # Use different efficiency for charged current
            for x in out_enrgs:
                effic_arr[elem] = get_detector_effic(
                    f'{detector_name}_cc', out_enrgs[elem])
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
