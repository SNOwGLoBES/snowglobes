# Script to make a smearing matrix for expected detector response
# This script incorporates both neutrino reactions as well as detector resolution
# Run with no arguments to be run interactively
# The interactive usage will make it simple and hard to mess up
# Run with a single argument, a json file to read in, to bypass interactive running

import numpy as np
import sys
import json
import collections
import code
import math
import time

# This little bit of code below is meant to ensure that input is correctly implemented
# whether this is run in python 2.x or 3.x
try:
   input = raw_input
except NameError:
   pass

###############################################################################
# Here we define the array of parameters!  All default values *must* be defined here
# I include some small comments here
params = {}
# interaction_type: 1 = nu-e elastic, 2 = nu-nuc CC, 3 = nu-nuc NC, 4 = custom (not implemented)
params["interaction_type"] = 0        # neutrino interaction type...
params["q_threshold"] = 0.0           # MeV, q-value threshold of the interaction, g.s. to g.s.
params["nc_gamma_energy"] = 0.0       # MeV, total gamma energy from excited nucleus
params["enu_min"] = 0.5               # MeV, min neutrino energy
params["enu_max"] = 100.0             # MeV, max neutrino energy
params["n_points_enu"] = 200          # how many points in neutrino energy?
params["edet_min"] = 0.5              # MeV, min detected energy
params["edet_max"] = 100.0            # MeV, max detected energy
params["n_points_edet"] = 200         # how many bins in detected energy?
params["resolution_a"] = 0.05         # either 1-sigma resolution, or A term for function
params["resolution_b"] = 0            # B term for resolution function
params["resolution_type"] = 0         # 1 = constant, 2 = function
params["neutrino_flavor"] = 0         # nue, numu, nutau, etc.
params["targetname"] = "default_name" # output filename component
params["detname"] = "default_name"    # output filename component
params["outpath"] = "./"              # path for where to write output files
params["prec"] = 8   # factor of how many more steps we use during computations, not user-settable


# Now we have a large set of functions that we'll be using a lot while the code runs.


###############################################################################
# The next few functions are set up to collect user input and validate it (a tiny bit)
def FloatInput(text, errortext = "Invalid input", minval = -1e50, maxval = 1e50):
   successful_input = False
   while not successful_input:
      try:
         response = float(input(text + "-->"))
         if response >= minval and response <= maxval:
            successful_input = True
         else: 
            print ("Value out of bounds, please try again")
      except ValueError:
         print(errortext)
   return response

###############################################################################
def IntInput(text, errortext = "Invalid input", minval = -1e50, maxval = 1e50):
   successful_input = False
   while not successful_input:
      try:
         response = int(input(text + "-->"))
         if response >= minval and response <= maxval:
            successful_input = True
         else: 
            print ("Value out of bounds, please try again")
      except ValueError:
         print(errortext)
   return response

###############################################################################
# Take in a yes/no from the user, convert to a boolean
def BoolInput(text, errortext = "I couldn't tell what you were trying to input, try just y or n"):
   successful_input = False
   bool_result = False
   while not successful_input:
      try:
         response = input("(y/n)-->")
         if response[0] == "y" or response[0] == "Y":
            successful_input = True
            bool_result = True
         if response[0] == "n" or response[0] == "N":
            successful_input = True
            bool_result = False
         if not successful_input:
            print("I couldn't tell what you were trying to input, try just y or n")
      except:
         print("I couldn't tell what you were trying to input, try just y or n")
         continue
   return bool_result

###############################################################################
# Get inputs from the user (if needed), and define all the parameters for running
def GatherUserInputs(params):
   # interaction_type: 1 = nu-e elastic, 2 = nu-nuc CC, 3 = nu-nuc NC, 4 = custom (not implemented)
   
   print("Select a neutrino interaction type:")
   print("1) neutrino-electron elastic scattering")
   print("2) neutrino-nucleus charged-current interaction")
   print("3) neutrino-nucleus neutral-current (specific inelastic gamma) interaction")
   print("4) custom interaction (code it in in the designated spot)")
   
   params["interaction_type"] = IntInput("Interaction type", 
                               "Try again, please use an integer between 1 and 4, then press <enter>",
                               1, 4)
   
   print("The default energy range and binning for the smearing matrix is:")
   print("  enu_min = 0.5 # MeV")
   print("  enu_max = 100.0 # MeV")
   print("  n_points_enu = 200")
   print("  edet_min = 0.5 # MeV")
   print("  edet_max = 100.0 # MeV")
   print("  n_points_edet = 200")
   print("Would you like to change from these defaults?")


   params["custom_binning"] = BoolInput("Change away from default binning? ")
   
   if params["custom_binning"]: # Set up a custom binning
      print("Custom binning can be useful for special cases, such as very high energy neutrinos")
      print(" or coherent scattering.  Be certain the custom values here match those set elsewhere")
      print(" for the detector, efficiency, and backgrounds files you will be using.")

      params["enu_min"] = FloatInput("Min neutrino energy (MeV)", "Invalid input", 0, 10000)
      params["enu_max"] = FloatInput("Max neutrino energy (MeV)", "Invalid input", params["enu_min"], 10000)
      params["n_points_enu"] = IntInput("Number of points in neutrino energy", "Invalid input", 10, 100000)
      params["edet_min"] = FloatInput("Min detected energy (MeV)", "Invalid input", 0, 10000)
      params["edet_max"] = FloatInput("Max detected energy (MeV)", "Invalid input", params["edet_min"], 10000)
      params["n_points_edet"] = IntInput("Number of bins in detected energy", "Invalid input", 10, 100000)

   if params["interaction_type"] == 1:
      # Input specifically for neutrino-electron elastic scattering
      # Oh, hey, don't need any!
      print("Will be using neutrino-electron elastic scattering")
      params["detname"] = input("Detector name (for filename)-->")
   elif params["interaction_type"] == 2:
      # CC nu-nuc interaction
      print("What is the Q-value (threshold) for the reaction (in MeV)?")
      print("Only count the ground state-ground state difference.")
      print("The excited state of the final nucleus will de-excite to ground state.")
      print("Don't forget to include the mass of the produced lepton!")
      params["q_threshold"] = FloatInput("g.s.-g.s. Q-value (MeV)")
      params["detname"] = input("Detector name (for filename)-->")
      params["targetname"] = input("Target (e.g. C12) name (for filename)-->")
   elif params["interaction_type"] == 3:
      # NC nu-nuc interaction
      print("What is the total gamma-ray energy produced in this interaction (in MeV)?")
      nc_gamma_energy = FloatInput("Toal gamma energy (MeV)")
      params["detname"] = input("Detector name (for filename)-->")
      params["targetname"] = input("Target (e.g. C12) name (for filename)-->")
   elif params["interaction_type"] == 4:
      # If you want to have another interaction, you should update the code here!
      print("We don't have custom interactions implemented at this time")
   else:
      print("No valid interaction was given, please contact SNOWGLOBES authors")

   # Ask about where to write the output
   print("The current location for output files is " + params["outpath"])
   print("Would you like to change it?")
   changepath = BoolInput("Change output file path? ")
   if changepath:
      outpath_is_valid = False
      while not outpath_is_valid:
         params["outpath"] = input("Enter a new output path -->")
         try:
            testfile = open(params["outpath"] + "/temptest.test", "w")
            testfile.close()
            outpath_is_valid = True
         except OSError:
            print("Invalid filepath or no permissions")
         except IOError:
            print("Invalid filepath or no permissions")
      # Ok, outpath looks good!
   # Otherwise, leave the path as is

   # ask for the neutrino flavor!
   print("Please indicate the neutrino flavor")
   print(" 0) Perform for all 6 flavors (not possible for CC interactions)")
   print(" 1) Electron neutrino")
   print(" 2) Muon neutrino")
   print(" 3) Tau neutrino")
   print(" 4) Electron antineutrino")
   print(" 5) Muon antineutrino")
   print(" 6) Tau antineutrino")
   
   # Need a single flavor if CC interaction
   if params["interaction_type"] == 2:
      min_flav = 1
   else:
      min_flav = 0

   params["neutrino_flavor"] = IntInput("Neutrino flavor", "Invalid flavor!", min_flav, 6)

   # Query about the detector resolution parameterization
   print("How will we define Gaussian detector resolution?")
   print("For perfect resolution, choose constant resolution of zero")
   print(" 1) Constant resolution (i.e. 0.1 if sigma = 10%)")
   print(" 2) sqrt[(A/Sqrt[E])^2+ B^2] (E is detected energy in MeV)")
   params["resolution_type"] = IntInput("Resolution definition", "Invalid input", 1,2)
   if params["resolution_type"] == 1:
      params["resolution_a"] = FloatInput("Energy resolution?", "Invalid input", 0, 1)
   elif params["resolution_type"] == 2:
      params["resolution_a"] = FloatInput("Parameter A?", "Invalid input", 0, 10000)
      params["resolution_b"] = FloatInput("Parameter B?", "Invalid input", 0, 10000)
   else:
      print("no valid resolution given, try again or contact authors")
      raise

   # All the parameters should be set now, move on!

   return params

###############################################################################
# A function to convert the unicode that json gives us to strings
def convert(data):
    if isinstance(data, basestring):
        return str(data)
    elif isinstance(data, collections.Mapping):
        return dict(map(convert, data.iteritems()))
    elif isinstance(data, collections.Iterable):
        return type(data)(map(convert, data))
    else:
        return data


###############################################################################
# Now we need to create arrays which we'll operate on.
# We use numpy to take advantage of fast computations
def GetInitialArrays(params):
   # We'll return two empty arrays, one of enu, and one of edet
   # First index is detected energy bin
   # Second index is true neutrino energy point

   # We're going to increase the number of edet bins by a factor of prec
   # temporarily, to improve calculation precision
   # We don't do the same for enu, because those are points, not bins

   n_points_enu = params["n_points_enu"]
   n_points_edet = params["n_points_edet"]*params["prec"]
   unsmeared = np.zeros((n_points_edet, n_points_enu), dtype=float)
   smeared = np.zeros((n_points_edet, n_points_enu), dtype=float)
   
   # Neutrino energies are evaluated at points, not in bins
   enu_min = params["enu_min"]
   enu_max = params["enu_max"]
   enu_step = (enu_max - enu_min)/float(n_points_enu - 1)
   # Use a lambda function to evaluate neutrino energies
   enu = np.fromfunction(lambda i,j: enu_min + enu_step * j, 
                         (n_points_edet, n_points_enu), dtype=float)
   # Detected energies are evaluated in bins
   edet_min = params["edet_min"]
   edet_max = params["edet_max"]
   edet_step = (edet_max - edet_min)/float(n_points_edet)
   # Use a lambda function to evaluate detected energies (bin centers)
   edet = np.fromfunction(lambda i,j: edet_min + edet_step * (i+0.5),
                          (n_points_edet, n_points_enu), dtype=float)
   return {"unsmeared": unsmeared, "smeared": smeared, "enu": enu, "edet": edet}

###############################################################################
# Function to populate the unsmeared matrix for NC gamma interactions
def GetUnsmearedNC(params, matrixset):
   edet_min = params["edet_min"]
   edet_max = params["edet_max"]
   n_points_edet = params["n_points_edet"]*params["prec"]
   edet_step = (edet_max - edet_min)/float(n_points_edet)
   edet = matrixset["edet"]
   ncgammamat = np.empty(edet.shape)
   ncgammamat.fill(params["nc_gamma_energy"])
   mid1 = np.greater_equal(ncgammamat, edet - edet_step/2.0)
   mid2 = np.less(ncgammamat, edet + edet_step/2.0)
   # matrix returned should be 1 at NC gamma energy, 0 elsewhere
   return np.multiply(mid1.astype(int), mid2.astype(int))

###############################################################################
# Function to populate the unsmeared matrix for CC nu-nuc interactions
def GetUnsmearedCC(params, matrixset):
   edet_min = params["edet_min"]
   edet_max = params["edet_max"]
   n_points_edet = params["n_points_edet"]*params["prec"]
   edet_step = (edet_max - edet_min)/float(n_points_edet)
   edet = matrixset["edet"]
   enu = matrixset["enu"]
   q_threshold = params["q_threshold"]
   # matrix is true if enu - edet >= threshold - half_bin_width
   mid1 = np.greater_equal(enu - edet, q_threshold - edet_step/2.0)
   # matrix is true if enu - edet < q_threshold + half_bin_width
   mid2 = np.less(enu - edet, q_threshold + edet_step/2.0)
   # put these two together to get the threshold reaction matrix
   return np.multiply(mid1.astype(int), mid2.astype(int))

###############################################################################
# Function to populate the unsmeared matrix for nu-e elastic interactions
# This function needs the flavor separately specified, in case we're doing more than one
def GetUnsmearedNuE(params, matrixset, flavor):
   me = 0.511 # electron mass, MeV
   s2tw = 0.231 # sin^2(theta_W)
   if flavor == 1: # nue
      gL = 2*s2tw+1
      gR = 2*s2tw
   elif flavor == 2 or flavor == 3: # numu or nutau
      gL = 2*s2tw-1
      gR = 2*s2tw
   elif flavor == 4: # nuebar
      gL = 2*s2tw
      gR = 2*s2tw+1
   elif flavor == 5 or flavor == 6: # numubar or nutaubar
      gL = 2*s2tw
      gR = 2*s2tw-1

   enu = matrixset["enu"]
   edet = matrixset["edet"]

   # if you weren't familiar with the python term "broadcasting", look it up
   # We're going to use an expression with a mixture of matrices and scalars,
   # and the whole thing will make a nice matrix... and efficiently, too!
   tmax = 2*np.divide(np.power(enu, 2),(me + 2*enu))
   newmat = np.less(edet, tmax) # mask for whether detected energy is below
                                # the maximum allowed kinematically

   # common formula for differential cross section, with arbitrary scaling
   # find it, for example, https://arxiv.org/pdf/hep-ph/0603036.pdf
   newmat2 = (gL*gL + gR*gR*np.power(1-np.divide(edet,enu),2) 
             - gL*gR*me*np.divide(edet, np.power(enu,2)))

   return np.multiply(newmat,newmat2)


###############################################################################
# Function to populate the unsmeared matrix for a custom user interaction
def GetUnsmearedCustom(params, matrixset):
   # User must provide the function here!
   print("No user function has been set yet!  Cannot run!")
   return False


###############################################################################
# Function to downsample the matrix back to the original requested dimensions
# Basically took the recipe from https://scipython.com/blog/binning-a-2d-array-in-numpy/
def Downsample(params, matrix):
   prec = params["prec"]
   # quick double check
   if matrix.shape[0]%prec != 0:
      print("illegal downsampling!")
      raise
   newdim0 = matrix.shape[0]/prec
   newdim1 = matrix.shape[1]
   shape = (newdim0, matrix.shape[0]/newdim0,
            newdim1, matrix.shape[1]/newdim1)
   return matrix.reshape(shape).mean(-1).mean(1)


###############################################################################
# We need a function to write out the smearing matrix after it is computed!
def WriteMatrix(filename, channelname, matrix):
   outfile = open(filename, "w")
   outfile.write("energy(#" + channelname + ")<\n")
   lastelem = str(len(matrix[1,:]) - 1)
   # To start the first line...
   outfile.write("@energy = ")
   for irow in range(0, matrix.shape[0]-1):
      our_line = "{0," + lastelem + "," + ",".join([x for x in matrix[irow,:].astype(str)]) + "}:\n"
      outfile.write(our_line)
   # the last line needs a semicolon instead of a colon
   our_line = "{0," + lastelem + "," + ",".join([x for x in matrix[irow,:].astype(str)]) + "};\n"
   outfile.write(our_line)
   outfile.write(">")
   outfile.close()
   return

###############################################################################
# Here is a function to normalize the edet columns to sum to one
def NormalizeEdet(matrix):
   # recall that edet is the first index!
   # Need a bit of extra work to ensure no division by zero!
   matsum = matrix.sum(axis=0)
   ncols = len(matsum)
   for icol in range(ncols):
      if matsum[icol] > 0:
         matrix[:,icol] = matrix[:,icol]/matsum[icol]
   return matrix

###############################################################################
# We need to apply a gaussian smearing to the matrix!
# Start with a simple function to get gaussian values
def GaussianProb(sigma, dev, step):
   normterm =  step/(sigma * np.sqrt(2*np.pi))
   mainterm = np.exp(-pow(dev,2)/(2*pow(sigma,2)))
   return normterm*mainterm

###############################################################################
# Here's the full function to apply the smear
def GaussianSmear(params, matrixset, matrix):
   # First get the smearing type and params
   resolution_type = params["resolution_type"]
   resolution_a = params["resolution_a"]
   resolution_b = params["resolution_b"]
   if resolution_type == 1 and resolution_a*resolution_a < 1e-12:
      # in this case, we assume that the intent is for no smearing
      return matrix
   if resolution_type < 1 or resolution_type > 2:
      print("Invalid resolution type!")
      return matrix
   # Smearing should occur ONLY in the detected energy dimension.
   # We'll be able to get some efficiency by applying the gaussian
   # smear to a full row of constant edet at the same time.
   edet = matrixset["edet"]
   if matrix.shape != edet.shape:
      print("Matrix shape incompatibility in smearing!")
      print("edet shape: ", edet.shape)
      print("matrix shape: ", matrix.shape)
      raise

   smearedmat = np.zeros(matrix.shape, dtype=float)
   # Let's get the step between matrix elements so we know
   # how many elements over we need to do the gaussian
   edet_min = params["edet_min"]
   edet_max = params["edet_max"] 
   n_points_edet = edet.shape[0]
   edet_step = (edet_max - edet_min)/float(n_points_edet)
   # Remember that the first index is detected energy bin
   # We're going to loop over detected energies
   print("Applying gaussian smearing, this may take a bit...")
   for iedet in range(n_points_edet):
      our_edet = edet[iedet,0] # edet for this one...
      # 
      if resolution_type == 1:
         gaus_res = resolution_a * our_edet
      else:
         gaus_res = math.sqrt(pow(resolution_a/math.sqrt(our_edet),2)+pow(resolution_b,2)) * our_edet
      # Assume we want to apply the gaussian for up to 5 sigma
      gaussian_extent_bins = int(math.floor(gaus_res*5/edet_step))
      edet_row_min = max(0, iedet - gaussian_extent_bins)
      edet_row_max = min(n_points_edet-1, iedet + gaussian_extent_bins)
      # Now we loop over each row, compute the smeared contribution, and add it
      # into the final matrix
      # Actually, let's not loop over each row, that's too slow.  Matrix math!
      # Yes, we can broadcast to calculate the gaussian prob on all the entries at once!
      # Doing it this way is remarkably faster than doing it in a loop or even in a 
      # list comprehension
      row_gaus_array = GaussianProb(gaus_res, edet[edet_row_min:edet_row_max+1,0]-our_edet, edet_step)
      col_gaus_array = row_gaus_array.reshape(1, len(row_gaus_array))
      row_matrix_iedet = matrix[iedet,:].reshape(len(matrix[iedet,:]),1)
      # Maybe withs sparse matrices we could make this step faster?
      matdot = (np.dot(row_matrix_iedet, col_gaus_array)).T
      smearedmat[edet_row_min:edet_row_max+1,:] += matdot
   # Finish the loop over rows
   return smearedmat



###############################################################################
###############################################################################
# Here's the actual running of the code!
if len(sys.argv) == 1: # interactive behavior
   print("Beginning interactive smearing matrix generation")
   params = GatherUserInputs(params)
elif len(sys.argv) == 2: # read from json file
   try:
      json_file = open(sys.argv[1])
   except:
      print("Couldn't open the json file at all!")
      raise
   try:
      json_str = json_file.read()
      json_data = convert(json.loads(json_str))
   except:
      print("Some problem with parsing the json file")
      raise
   for key in json_data:
      if "comment" in key:
         continue
      else:
         params[key] = json_data[key]
print params


matrixset = GetInitialArrays(params)

# params and matrixset will be dictionaries we carry around

# Let's go through the logic if we're making matrices for nu-e elastic scattering
if params["interaction_type"] == 1:
   flavorset = []
   flavornames = ["invalid", "nue", "numu", "nutau", "nuebar", "numubar", "nutaubar"]
   # decide on flavor
   if params["neutrino_flavor"] == 0:
      flavorset = [1,2,3,4,5,6]
   else:
      flavorset = [params["neutrino_flavor"]]
   # Loop over flavors
   for iflavor in flavorset:
      our_mat1 = GetUnsmearedNuE(params, matrixset, iflavor) # generate the matrix
      our_mat2 = GaussianSmear(params, matrixset, our_mat1) # apply gaussian smearing
      our_mat3 = Downsample(params, our_mat2) # downsample back to requested dims
      our_mat4 = NormalizeEdet(our_mat3) # normalize rows of detected energy
      # Assemble the names using inputs and pre-defined formulae
      path = params["outpath"] + "/"
      channame = flavornames[iflavor] + "_e_smear"
      filename = "smear_" + flavornames[iflavor] + "_e_" + params["detname"] + ".dat"
      WriteMatrix(path + filename, channame, our_mat4)

# Next, go through the logic if we're making matrices for nu-nuc CC interactions
elif params["interaction_type"] == 2:
   iflavor = params["neutrino_flavor"]
   flavornames = ["invalid", "nue", "numu", "nutau", "nuebar", "numubar", "nutaubar"]
   # decide on flavor
   our_mat1 = GetUnsmearedCC(params, matrixset) # generate the matrix
   our_mat2 = GaussianSmear(params, matrixset, our_mat1) # apply gaussian smearing
   our_mat3 = Downsample(params, our_mat2) # downsample back to requested dims
   our_mat4 = NormalizeEdet(our_mat3) # normalize rows of detected energy
   # Assemble the names using inputs and pre-defined formulae
   path = params["outpath"] + "/"
   channame = flavornames[iflavor] + "_" + params["targetname"] + "_smear"
   filename = "smear_" + flavornames[iflavor] + "_" + params["targetname"] + "_" + params["detname"] + ".dat"
   WriteMatrix(path + filename, channame, our_mat4)

# Next, go through the logic if we're making matrices for nu-nuc NC interactions
if params["interaction_type"] == 3:
   flavorset = []
   flavornames = ["invalid", "nue", "numu", "nutau", "nuebar", "numubar", "nutaubar"]
   # decide on flavor
   if params["neutrino_flavor"] == 0:
      flavorset = [1,2,3,4,5,6]
   else:
      flavorset = [params["neutrino_flavor"]]
   # Loop over flavors
   for iflavor in flavorset:
      our_mat1 = GetUnsmearedNC(params, matrixset) # generate the matrix
      our_mat2 = GaussianSmear(params, matrixset, our_mat1) # apply gaussian smearing
      our_mat3 = Downsample(params, our_mat2) # downsample back to requested dims
      our_mat4 = NormalizeEdet(our_mat3) # normalize rows of detected energy
      # Assemble the names using inputs and pre-defined formulae
      path = params["outpath"] + "/"
      channame = "nc_" + flavornames[iflavor] + "_" + params["targetname"] + "_smear"
      filename = "smear_nc_" + flavornames[iflavor] + "_" + params["targetname"] + "_" + params["detname"] + ".dat"
      WriteMatrix(path + filename, channame, our_mat4)

# Finally, go through the logic for a custom interaction
if params["interaction_type"] == 4:
   print(" This would run for a custom interaction.  Still needs to be implemented.")

# Below is some debugging stuff...
#print(matrixset)
#testmat = GetUnsmearedNC(params, matrixset)
#testmat2 = GetUnsmearedCC(params, matrixset)
#testmat3 = GetUnsmearedNuE(params, matrixset, 1)
#print testmat
#print testmat.shape
#testmat4 = Downsample(params, testmat3)
#print testmat4
#print testmat4.shape

#iden = np.identity(800)
#print("here comes the slow part!")
#s = time.time()
#idensmear = GaussianSmear(params, matrixset, iden)
#e = time.time()
#print ("That took", e-s, "seconds")
#print idensmear[0:10,3]
#print idensmear[0:5,0:5]
#code.interact("going interactive", local=globals())
