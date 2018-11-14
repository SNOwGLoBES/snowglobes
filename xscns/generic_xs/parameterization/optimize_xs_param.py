# Try to find a function to optimize a parameterization for neutrino cross sections
# This optimizes cross section magnitudes for Michel spectrum neutrinos.
# code initially written and tested with python 2.7.10
# by jba 7/2018

import scipy.optimize as optimize
import math
import code

# Function to read in the parameters we're optimizing for
def GetData(filename = "athar_xs.dat"):
   # Note that these cross sections are in 10^-42 cm^2
   xs_data = []
   infile = open(filename, "r")
   for line in infile.readlines():
      if line[0] == "#":
         continue # skip over header lines
      if len(line) < 5:
         continue # skip over empty lines
      else:
         split_line = line.split()
         nuc_name = split_line[0]
         a = int(split_line[1])
         z = int(split_line[2])
         q = float(split_line[3])
         xs = float(split_line[4])
         xs_data.append({"name": nuc_name, "a": a, "z": z, "q": q, "xs": xs})
   return xs_data


# Here's the function of a, z, q which we'll use to estimate the xs
# This is the function which should be played around with to make good

def calc_xs(nuc_param_dict, fit_param_list):
   a = nuc_param_dict["a"]
   z = nuc_param_dict["z"]
   q = nuc_param_dict["q"]
   p0 = fit_param_list[0]
   p1 = fit_param_list[1]
   p2 = fit_param_list[2]
   p3 = fit_param_list[3]
   p4 = fit_param_list[4]
   p5 = fit_param_list[5]

   return (p0 * pow(a, p1) + p2 * pow(a-z, p3)) * ((1 + p4 * pow(q, p5)))


# Function to compute a chi2 for how well we fit a single cross section
# I'm not sure if it is a real chi2, but it should do for our purposes
def chi2(xs_datum, fit_param_list):
   true_xs = xs_datum["xs"]
   fit_xs = calc_xs(xs_datum, fit_param_list)
   return math.sqrt(pow((true_xs - fit_xs),2))/true_xs

# Sum the chi2 over all our data points
def totalchi2(fit_param_list):
   xs_data = GetData("athar_xs.dat")
   ndat = len(xs_data)
   sum_chi2 = 0
   for xs_datum in xs_data:
      sum_chi2 += chi2(xs_datum, fit_param_list)

   print(fit_param_list, sum_chi2)
   return sum_chi2

# Now we do the minimization!

initial_guess = [1, 1, 1, 1, 0.2, 1]
coeff_range = (-40, 40)
pow_range = (-5, 5)
bounds_list = [coeff_range, pow_range, coeff_range, pow_range, coeff_range, pow_range]
result = optimize.minimize(totalchi2, initial_guess, method = "SLSQP", bounds = bounds_list,
                           options = {"maxiter":10000})

print(result)


# Now let's print out how the minimization went

xs_data = GetData("athar_xs.dat")
print("nuclide  predicted xs  fit_func xs  frac_error   q_value")
for xs_datum in xs_data:
   true_xs = xs_datum["xs"]
   fit_xs = calc_xs(xs_datum, result.x)
   namestr = xs_datum["name"].ljust(10)
   truexsstr = str(true_xs).ljust(12)
   fitxsstr = "{:6f}".format(fit_xs).ljust(15)
   errorstr = "{:4f}".format(1 - fit_xs/true_xs).ljust(15)
   qstr = "{:4f}".format(xs_datum["q"])
   print(namestr + truexsstr + fitxsstr + errorstr + qstr)


