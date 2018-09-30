#!/usr/bin/env python
# Apply a gaussian blur to a smearing matrix
# By jba 09/18/18
# Usage: ./apply_blur.py <path-to-smearing-file.dat> <path-to-output-smearing-file.dat> [opt. blur size]
# After running, the code will plot the mean and std dev.  You can use this to confirm that you haven't
# blurred too much.  If the standard deviation is getting too big, move to a smaller number of bins for blurring

import sys
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt


# look for input, give helpful advice if it is problematic
infile = 0
outfilename = 0
blur_size = 3
if len(sys.argv) <3 or len(sys.argv) > 4:
    print("Usage: ./apply_blur.py <path-to-smearing-file.dat> <path-to-output-smearing-file.dat> [opt. blur size (bins)]")
    print("Check the plotted statistics after running to decide if the blurring is reasonable")
    sys.exit()
else:
    try:
        infile = open(sys.argv[1], "r")
    except:
        print("problem encountered trying to open " + sys.argv[1])
        raise
    try:
        outfile = open(sys.argv[2], "w")
        outfile.close()
        outfilename = sys.argv[2]
    except:
        print("problem encountered trying to open " + sys.argv[2])
        raise
    if len(sys.argv) == 4:
        blur_size = int(sys.argv[3])
    

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

z_mesh = np.array(smearrows)

blurred_z_mesh = ndimage.gaussian_filter(z_mesh, blur_size)

# Now we can set up some masking to bring things like thresholds back.
z_mesh_column_filter = np.sum(z_mesh, axis=0, keepdims=True)
blurred_z_mesh = z_mesh_column_filter*blurred_z_mesh


# Now let's re-normalize the matrix, as the blurring might mess up normalization
blurred_normalized_z_mesh = NormalizeEdet(blurred_z_mesh)

# Now we write the output as requested
WriteMatrix(outfilename, flavorname, blurred_normalized_z_mesh)

# and we're done, unless we want to keep going to analyze things a bit...

###############################################################################
def WeightedAvgAndSTD(matrix):
   """
   Return the weighted average and standard deviation.
   We are assuming 0.5 to 100 MeV
   """
   edet_min = 0.5
   edet_max = 100
   edet_step = (edet_max - edet_min)/float(matrix.shape[0])
   # Use a lambda function to evaluate detected energies (bin centers)
   edet = np.fromfunction(lambda i,j: edet_min + edet_step * (i+0.5),
                          matrix.shape, dtype=float)
   average = np.ma.average(edet, weights=matrix, axis=0)
   # Fast and numerically precise:
   variance = np.ma.average((edet-average)**2, weights=matrix, axis=0)
   return (average, np.sqrt(variance))

###############################################################################
# Let's analyze the mean and sigma before and after blurring!
def PlotMeanAndSigma(matrix1, matrix2):
   mean1, std1 = WeightedAvgAndSTD(matrix1)
   mean2, std2 = WeightedAvgAndSTD(matrix2)
   plt.figure(1)
   plt.subplot(121)
   plt.plot(mean1, label="before mean")
   plt.plot(mean2, label="after mean")
   legmean = plt.legend(loc="best")
   plt.subplot(122)
   plt.plot(std1, label="before stdev")
   plt.plot(std2, label="after stdev")
   legstd = plt.legend(loc="best")
   print("Close the plot window to finish the program")
   plt.show()
   return

PlotMeanAndSigma(z_mesh, blurred_normalized_z_mesh)
