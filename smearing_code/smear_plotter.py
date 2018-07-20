# Plot the smearing matrix using pyplot
# By jba 07/17/18
# Usage: ./smear_plotter.py <path-to-smearing-file.dat> [optional log bool]
# Usage note: after the plot is made, close the plot window to end the script

import sys
import code
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

logopt = "False" # by default, don't plot logarithmically

# look for input, give helpful advice if it is problematic
infile = 0
if len(sys.argv) <2 or len(sys.argv) > 3:
    print("Usage: ./smear_plotter.py <path-to-smearing-file.dat> [optional log bool]")
else:
    try:
        infile = open(sys.argv[1], "r")
    except:
        print("problem encountered trying to open " + sys.argv[1])
        raise
    if len(sys.argv) == 3:
        logopt = sys.argv[2]

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


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
nbins = 200
emin = 0.0005 # GeV
emax = 0.100  # GeV


z_mesh = np.array(smearrows)

if z_mesh.shape != (nbins, nbins):
    print "Problem with dimensions!"

# Do we plot logarithmic color scale?
if logopt[0] == "T" or logopt[0] == "t" or logopt[0] == "1":
    normpar = colors.LogNorm()
elif logopt[0] == "F" or logopt[0] == "f" or logopt[0] == "0":
    normpar = None
else:
    print "Are you sure you properly set a logarithm bool?"
    normpar = None

smearplot = plt.imshow(z_mesh, # plot the smearing values as an image, because python
                       cmap="CMRmap", # nice colormap 
                       interpolation="nearest", # avoid smoothing
                       origin="lower", # (0,0) is in the bottom left corner
                       extent=(emin,emax,emin,emax), # set axis ranges
                       norm=normpar) # log scale
plt.title("Smearing matrix")
plt.xlabel("Neutrino energy (GeV)")
plt.ylabel("Smeared reconstructed energy (GeV)")
cbar = plt.colorbar(smearplot)

plt.show()

