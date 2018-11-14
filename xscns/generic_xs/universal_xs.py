# Produce a cross section file for an arbitrary A, Z, and Q.
# Uses fitted shapes and magnitudes
# This is NOT guaranteed to be accurate.  In fact, it is certain to be fairly
# inaccurate.  But it has a good chance to be accurate to within a factor of a few.
# See output from running with no arguments for instructions.
# code initially written and tested with python 2.7.10, also works with 3.7
# by jba 8/9/2018

import sys
import math

# look for input, give helpful advice if it is problematic
Q = 4.3 # MeV!
A = 40 #Argon
Z = 18 # Argon
is_nuebar = False
ofname = "./default_argon.dat"
emin = 0.0002 # GeV
emax = 0.1 # GeV
nbins = 1001 # this one cannot change!
if len(sys.argv) not in [4, 5, 7, 8]:
   print("Usage: python ./universal_xs.py Z A Q [+/- for nue/nuebar] [<filename>] [emin emax]")
   print("Q must be in MeV, the optional emin and emax must be in GeV")
   print("These are very approximate cross sections, use with bucket of salt")
if len(sys.argv) in [4,5,7,8]:
   Z = int(sys.argv[1])
   A = int(sys.argv[2])
   Q = float(sys.argv[3])
if len(sys.argv) in [5,6,8]:
   if sys.argv[4] == "-":
      is_nuebar = True
if len(sys.argv) in [6,8]:
   ofname = sys.argv[5]
if len(sys.argv) < 6:
   if is_nuebar:
      ofname = "./generic_CC_xs_nuebar_A_Z_" + str(A) + "_" + str(Z) + ".dat"
   else:
      ofname = "./generic_CC_xs_nue_A_Z_" + str(A) + "_" + str(Z) + ".dat"
if len(sys.argv) == 8:
   emin = float(sys.argv[6])
   emax = float(sys.argv[7])
if len(sys.argv) == 1:
   print ("We'll run an example run for argon, as if you entered:")
   print ("python ./universal_xs.py 18 40 4.3 +")
if len(sys.argv) not in [1,4,5,6,8]:
   sys.exit()

# ok, at this point, Q, A, Z, ofname, emin, emax, is_nuebar, and nbins are all defined


# We will use the magnitude parameterization determined in
# optimize_xs_param.py
# These are average cross sections over a Michel spectrum
# We will need to do an integral over the bins later to normalize
# NOTE THAT Q is in MeV!
def GetAverageXS(A, Z, Q):
   p0 = 11.8273
   p1 = 0.54375
   p2 = 1.3478
   p3 = 1.59275
   p4 = -0.09143
   p5 = 0.75409
   return (p0*pow(A, p1) + p2*pow(A-Z, p3))*(1+p4*pow(Q, p5))

# Michel spectrum as a function of x (neutrino energy in GeV)
# No normalization applied yet
def MichelFunction(x):
   e0 = 0.0528 # (GeV)
   return max(12/pow(e0,4)*x*x*(e0-x),0)
   

# We will use the functional shape determined in parameterize_shape.py
# This shape is not yet normalized at all!
# NOTE THAT Q is in MeV! (unlike in parameterize_shape.py)
# x, the neutrino energy, is in GeV!
def XSFunction(x, Q, A):
   q0 = 26.5460
   q1 = 158.585
   q2 = -1.9699
   q3 = 1000.0
   Qp = Q * 0.001
   xp = x - Qp
   return max((1000*xp*xp*xp + (q0 + q1/A + Qp*1000*q2)*xp*xp + q3/A*xp*xp*xp*xp),0)*(xp>0)

# We must determine the average cross section over the Michel spectrum
# We'll do a numerical integral in linear space for this purpose
# 0.01 MeV steps
x0 = 0
xstep = 0.00001
integration_points = [x0 + xstep*istep for istep in range(5290)]

michel_int = sum([MichelFunction(x) for x in integration_points])

func_int = sum([MichelFunction(x)*XSFunction(x,Q,A) for x in integration_points])

uncorrected_average_xs = func_int/michel_int
correction_factor = GetAverageXS(A, Z, Q)/uncorrected_average_xs
print("Avg xs should be: ", GetAverageXS(A,Z,Q))

# Now we make the cross section file
# This file must have 1001 points, evenly distributed logarithmically
# between emin and emax

logmin = math.log10(emin)
logmax = math.log10(emax)
logstep = (logmax - logmin)/1000.
en_log_list = [logmin + x * logstep for x in range(1000)]
en_log_list.append(logmax)
en_lin_list = [pow(10,x) for x in en_log_list]

# Open the output file
try:
   outfile = open(ofname, "w")
except:
   print("Failed to properly open " + ofname + " for writing")
   raise

outfile.write("# Generic CC cross section for A=" + str(A) + " Z=" + str(Z)
               + " Q=" + str(Q) + " (10^-38 cm^2/GeV) #\n")
outfile.write("# Be aware, these are approximate, accuracy is not guaranteed! #\n")
outfile.write("# log(energy in GeV)  nu_e      nu_mu      nu_tau      nu_e_bar      nu_mu_bar    nu_tau_bar #\n")

print("The correction factor was " + str(correction_factor))

for ielin, ielog in zip(en_lin_list, en_log_list):
   xs_plain = XSFunction(ielin, Q, A)
   xs_scaled = xs_plain * correction_factor * 0.0001
   xs_divided = xs_scaled/ielin
   # We've now gone from 10^-42 cm^2 to 10^-38 cm^2/GeV
   if is_nuebar:
      xs_divided *= 0.5*(0.35 + 5.0/float(A)) # extra factor for antineutrinos

   enstr = "{:6f}".format(ielog)
   xsstr = "{:6f}".format(xs_divided)
   zerxsstr = "{:5f}".format(0.)
   if not is_nuebar:
      outfile.write("{0}  {1}  {2}  {2}  {2}  {2}  {2}\n".format(enstr, xsstr, zerxsstr))
   else:
      outfile.write("{0}  {2}  {2}  {2}  {1}  {2}  {2}\n".format(enstr, xsstr, zerxsstr))

outfile.close()

