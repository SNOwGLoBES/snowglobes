#!/usr/bin/env python
# Plot the smearing matrix using pyplot
# By jba 07/17/18
# Usage: ./smear_plotter.py <path-to-smearing-file.dat> [optional log bool] [optional non-standard binning bool]
# Usage note: after the plot is made, close the plot window to end the script

import sys
import code
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

logopt = "False" # by default, don't plot logarithmically
nsbopt = "False" # by default, don't use non-standard binning

# look for input, give helpful advice if it is problematic
infile = 0
if len(sys.argv) <2 or len(sys.argv) > 4:
    print("Usage: ./smear_plotter.py <path-to-smearing-file.dat> [opt. log bool] [opt. non-standard binning bool]")
else:
    try:
        infile = open(sys.argv[1], "r")
    except:
        print("problem encountered trying to open " + sys.argv[1])
        raise
    if len(sys.argv) == 3 or len(sys.argv) == 4:
        logopt = sys.argv[2]
    if len(sys.argv) == 4:
        nsbopt = sys.argv[3]
    

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

# This little bit of code below is meant to ensure that input is correctly implemented
# whether this is run in python 2.x or 3.x
try:
   input = raw_input
except NameError:
   pass

###############################################################################
# The next couple functions are set up to collect user input and validate it (a tiny bit)
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


# Now we have the smearing matrix file open, time to parse it.
smearrows = []
flavorname = 0 # we might as well get the name from the top line
isfirst = False
islast = False
# The first row should be ignorable, as should the second, up to an @energy = {
for iline in infile.readlines():
    if "#" in iline:
        flavorname = find_between(iline, "#", ")")
        continue # this is a header line, move on
    else:
        rowvals = find_between(iline, "{", "}")
        rowvals_str = rowvals.split(",")[2:]
        smearrows.append([float(x) for x in rowvals_str]) # first two entries don't really matter for this
    if ";" in iline:
        # this is the last one
        break

print("Read in " + str(len(smearrows)) + " lines")

# Here I've hardcoded in the energy window, if you want to use a different energy window,
# make sure you check snowglobes/glb/preamble.glb to make ensure agreement
nbins_enu = 200
enumin = 0.0005 # GeV
enumax = 0.100  # GeV
nbins_edet = 200
edetmin = 0.0005 # GeV
edetmax = 0.100  # GeV

# Here's the possibility that we'll change the binning...
do_nonstandard_binning = False
# Do we plot logarithmic color scale?
if nsbopt[0] == "T" or nsbopt[0] == "t" or nsbopt[0] == "1" or nsbopt[0] == "Y" or nsbopt[0] == "y":
   do_nonstandard_binning = True
elif nsbopt[0] == "F" or nsbopt[0] == "f" or nsbopt[0] == "0" or nsbopt[0] == "N" or nsbopt[0] == "n":
   do_nonstandard_binning = False
else:
    print "Are you sure you properly set the nonstandard binning bool?"
if do_nonstandard_binning:
   print("Please enter non-standard binning information for the smearing file")
   enumin = 0.001*FloatInput("Min neutrino energy (MeV)", "Invalid input", 0, 100000)
   enumax = 0.001*FloatInput("Max neutrino energy (MeV)", "Invalid input", enumin, 100000)
   nbins_enu = IntInput("Number of points in neutrino energy", "Invalid input", 10, 100000)
   edetmin = 0.001*FloatInput("Min detected energy (MeV)", "Invalid input", 0, 100000)
   edetmax = 0.001*FloatInput("Max detected energy (MeV)", "Invalid input", edetmin, 100000)
   nbins_edet = IntInput("Number of bins in detected energy", "Invalid input", 10, 100000)


z_mesh = np.array(smearrows)

if z_mesh.shape != (nbins_enu, nbins_edet):
    print "Problem with dimensions!"

# Do we plot logarithmic color scale?
if logopt[0] == "T" or logopt[0] == "t" or logopt[0] == "1" or logopt[0] == "Y" or logopt[0] == "y":
    normpar = colors.LogNorm()
elif logopt[0] == "F" or logopt[0] == "f" or logopt[0] == "0" or logopt[0] == "N" or logopt[0] == "n":
    normpar = None
else:
    print "Are you sure you properly set a logarithm bool?"
    normpar = None

smearplot = plt.imshow(z_mesh, # plot the smearing values as an image, because python
                       cmap="CMRmap", # nice colormap 
                       interpolation="nearest", # avoid smoothing
                       origin="lower", # (0,0) is in the bottom left corner
                       extent=(enumin,enumax,edetmin,edetmax), # set axis ranges
                       norm=normpar) # log scale
plt.title("Smearing matrix")
plt.xlabel("Neutrino energy (GeV)")
plt.ylabel("Smeared reconstructed energy (GeV)")
cbar = plt.colorbar(smearplot)

plt.show()

