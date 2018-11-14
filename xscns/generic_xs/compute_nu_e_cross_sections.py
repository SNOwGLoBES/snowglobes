# Compute cross sections for neutrino/electron scattering
# using calculations from doi.org/10.1088/0954-3899/29/11/013
# This script has been tested with python 2.6 and 3.7

import numpy as np
import matplotlib.pyplot as plt
import math
import code

# Goodness so many constants (from the paper)

Gmu = 1.16637e-5 # GeV^-2
me  = 0.510999e-3 # GeV  the electron mass
pi = 3.14159
s2tw = 0.23
epml = 0.5 - s2tw   # epsilon_- for mu and tau
eppl = -s2tw        # epsilon_+ for mu and tau
epme = -0.5 - s2tw  # epsilon_- for e
eppe = -s2tw        # epsilon_+ for e
hbarccorr = 3.894e-28 # (hbar*c)^2/cm^2

# from eqns 16 and 17 in the paper
# this is a differential cross section that we'll integrate over
def InstantaneousXSnu(Enu, y, epm, epp):
   return 2*Gmu*Gmu*me*Enu/pi*(epm*epm + epp*epp*(1-y)*(1-y)
                               - epm*epp*me*y/Enu)

# integrate over the differential cross section to get total
def IntegratedXSnu(Enu, epm, epp):
   ymax = 1./(1 + me/(2*Enu))
   ystep = ymax/1000.
   ylist = [ystep * i for i in range(0,1001)]
   yarray = np.array(ylist)
   xsvals = InstantaneousXSnu(Enu, yarray, epm, epp) # note the efficient python broadcasting!
   return np.sum(xsvals)*ystep*hbarccorr

# from eqns 18-21 in the paper
# calculates cross section in the Enu >> me limit
def SimpleXSnu(Enu, lepnum):
   if lepnum == 1:
      return Gmu*Gmu*me*Enu/(2*pi)*(1 + 4 * s2tw + 16/3.*s2tw*s2tw)*hbarccorr
   elif lepnum == 2 or lepnum == 3:
      return Gmu*Gmu*me*Enu/(2*pi)*(1 - 4 * s2tw + 16/3.*s2tw*s2tw)*hbarccorr
   elif lepnum == 4:
      return Gmu*Gmu*me*Enu/(2*pi)*(1/3. + 4/3. * s2tw + 16/3.*s2tw*s2tw)*hbarccorr
   else:
      return Gmu*Gmu*me*Enu/(2*pi)*(1/3. - 4/3. * s2tw + 16/3.*s2tw*s2tw)*hbarccorr

# apply electroweak corrections
# This is a whole mess of stuff, but mostly equations 43 and 44
# We also compute the simple uncorrected values from equations 18-21,
# to make a correction factor
def EWCorr(Enu, lepnum):
   nocorval = 1 + 4*s2tw + 16./3. * s2tw*s2tw
   rho = 1.013
   kappatau = 1.0064
   kappamu = 0.9970
   kappae = kappamu - 0.0179
   alpha = 1/137.
   fminus = -2/3. * math.log(2*Enu/me) - 1/6. * (pi*pi -19/4.)
   fplus = -2/3. * math.log(2*Enu/me) - 1/6. * (pi*pi -43/4.)
   s2 = 0.231
   Avale = 1-2/rho - 2*kappae*s2
   Bvale = -2*kappae*s2
   Avalmu = 1-2*kappamu*s2
   Bvalmu = -2*kappamu*s2
   Avaltau = 1-2*kappatau*s2
   Bvaltau = -2*kappatau*s2

   if lepnum == 1: # e
      Aval = Avale
      Bval = Bvale
      nocorval = 1 + 4*s2tw + 16./3. * s2tw*s2tw
   elif lepnum == 2: # mu
      Aval = Avalmu
      Bval = Bvalmu
      nocorval =  1 - 4*s2tw + 16./3. * s2tw*s2tw
   elif lepnum == 3: # tau
      Aval = Avaltau
      Bval = Bvaltau
      nocorval = 1 - 4*s2tw + 16./3. * s2tw*s2tw
   elif lepnum == 4: # ebar
      Aval = Bvale
      Bval = Avale
      nocorval = 1/3.  + 4/3.*s2tw + 16./3. * s2tw*s2tw
   elif lepnum == 5: # mubar
      Aval = Bvalmu
      Bval = Avalmu
      nocorval = 1/3.  - 4/3.*s2tw + 16./3. * s2tw*s2tw
   elif lepnum == 6: # taubar
      Aval = Bvaltau
      Bval = Avaltau
      nocorval = 1/3.  - 4/3.*s2tw + 16./3. * s2tw*s2tw
   else:
      print ("Go away!")
      return -1
   #print("fminus, fplus = ", str(fminus), str(fplus))
   corrval = rho*rho* (pow(Aval,2) * (1+alpha/pi*fminus)
                       + 1/3. * pow(Bval,2) * (1+alpha/pi*fplus))
   #print (corrval, nocorval)
   return corrval/nocorval

# That's it with the equipment to calculate cross sections, let's write!

n_final = 1001 # need 1001 points in the final xs file

# change the energy range here if you must!
en_min = 0.5 # min E_nu in MeV
en_max = 200 # max E_nu in MeV

en_log_min = math.log10(en_min*0.001)
en_log_max = math.log10(en_max*0.001)

en_log_step = (en_log_max - en_log_min)/(n_final - 1.0)
en_log_list = [en_log_min + x * en_log_step for x in range(n_final)]

en_lin_list = [math.pow(10, x) * 1000.0 for x in en_log_list]


# now we'll loop over all the energy bins and all the flavors to compute the
# cross sections, along with the EW correctino and low E correction
# for validation purposes
xs_list_e  = [0]*1001
xs_list_mu  = [0]*1001
xs_list_tau  = [0]*1001
xs_list_ebar  = [0]*1001
xs_list_mubar  = [0]*1001
xs_list_taubar  = [0]*1001
ewcorr_list_e  = [0]*1001
ewcorr_list_mu  = [0]*1001
ewcorr_list_tau  = [0]*1001
ewcorr_list_ebar  = [0]*1001
ewcorr_list_mubar  = [0]*1001
ewcorr_list_taubar  = [0]*1001
lowecorr_list_e  = [0]*1001
lowecorr_list_mu  = [0]*1001
lowecorr_list_tau  = [0]*1001
lowecorr_list_ebar  = [0]*1001
lowecorr_list_mubar  = [0]*1001
lowecorr_list_taubar  = [0]*1001
for ibin in range(1001):
    ien = en_lin_list[ibin]
    if ibin%100 == 0:
       print(ibin)
    xs_list_e[ibin] = IntegratedXSnu(ien*0.001, epme, eppe)/1e-38 * EWCorr(ien*0.001, 1)
    xs_list_mu[ibin] = IntegratedXSnu(ien*0.001, epml, eppl)/1e-38 * EWCorr(ien*0.001, 2)
    xs_list_tau[ibin] = IntegratedXSnu(ien*0.001, epml, eppl)/1e-38 * EWCorr(ien*0.001, 3)
    xs_list_ebar[ibin] = IntegratedXSnu(ien*0.001, eppe, epme)/1e-38 * EWCorr(ien*0.001, 4)
    xs_list_mubar[ibin] = IntegratedXSnu(ien*0.001, eppl, epml)/1e-38 * EWCorr(ien*0.001, 5)
    xs_list_taubar[ibin] = IntegratedXSnu(ien*0.001, eppl, epml)/1e-38 * EWCorr(ien*0.001, 6)
    ewcorr_list_e[ibin] =  EWCorr(ien*0.001, 1)
    ewcorr_list_mu[ibin] =  EWCorr(ien*0.001, 2)
    ewcorr_list_tau[ibin] =  EWCorr(ien*0.001, 3)
    ewcorr_list_ebar[ibin] =  EWCorr(ien*0.001, 4)
    ewcorr_list_mubar[ibin] =  EWCorr(ien*0.001, 5)
    ewcorr_list_taubar[ibin] =  EWCorr(ien*0.001, 6)
    lowecorr_list_e[ibin] = IntegratedXSnu(ien*0.001, epme, eppe)/SimpleXSnu(ien*0.001,1)
    lowecorr_list_mu[ibin] = IntegratedXSnu(ien*0.001, epml, eppl)/SimpleXSnu(ien*0.001,2)
    lowecorr_list_tau[ibin] = IntegratedXSnu(ien*0.001, epml, eppl)/SimpleXSnu(ien*0.001,3)
    lowecorr_list_ebar[ibin] = IntegratedXSnu(ien*0.001, eppe, epme)/SimpleXSnu(ien*0.001,4)
    lowecorr_list_mubar[ibin] = IntegratedXSnu(ien*0.001, eppl, epml)/SimpleXSnu(ien*0.001,5)
    lowecorr_list_taubar[ibin] = IntegratedXSnu(ien*0.001, eppl, epml)/SimpleXSnu(ien*0.001,6)

# Just some simple printout to make sure it is working
for i in range(len(xs_list_e)):
    print (en_log_list[i], en_lin_list[i], xs_list_e[i], xs_list_mu[i])

# Now we plot the correction factors so we can be confident we're doing this right.
plt.plot(np.array(en_lin_list), np.array(ewcorr_list_e))
plt.plot(np.array(en_lin_list), np.array(ewcorr_list_mu))
plt.plot(np.array(en_lin_list), np.array(ewcorr_list_tau))
plt.plot(np.array(en_lin_list), np.array(ewcorr_list_ebar))
plt.plot(np.array(en_lin_list), np.array(ewcorr_list_mubar))
plt.plot(np.array(en_lin_list), np.array(ewcorr_list_taubar))
plt.legend(["e", "mu", "tau", "ebar", "mubar", "taubar"], loc='best')
plt.xlabel("Enu (MeV)")
plt.ylabel("EW Correction Factor")
plt.show()

plt.plot(np.array(en_lin_list), np.array(lowecorr_list_e))
plt.plot(np.array(en_lin_list), np.array(lowecorr_list_mu))
plt.plot(np.array(en_lin_list), np.array(lowecorr_list_tau))
plt.plot(np.array(en_lin_list), np.array(lowecorr_list_ebar))
plt.plot(np.array(en_lin_list), np.array(lowecorr_list_mubar))
plt.plot(np.array(en_lin_list), np.array(lowecorr_list_taubar))
plt.legend(["e", "mu", "tau", "ebar", "mubar", "taubar"], loc='best')
plt.xlabel("Enu (MeV)")
plt.ylabel("Low E Correction Factor")
#plt.semilogx()
plt.show()


# Ok, now we need to format it for output as the new file

# The energy must be represented as log10(E in GeV)
# The cross section is given as xs/E in 10^-38 cm^2/GeV
# Make sure we convert properly!

en_formatted_list = en_log_list

xs_formatted_list_e = []
xs_formatted_list_mu = []
xs_formatted_list_tau = []
xs_formatted_list_ebar = []
xs_formatted_list_mubar = []
xs_formatted_list_taubar = []
for ixs, ien in zip(xs_list_e, en_lin_list):
    xs_formatted_list_e.append(ixs/(ien/1000.))
for ixs, ien in zip(xs_list_mu, en_lin_list):
    xs_formatted_list_mu.append(ixs/(ien/1000.))
for ixs, ien in zip(xs_list_tau, en_lin_list):
    xs_formatted_list_tau.append(ixs/(ien/1000.))
for ixs, ien in zip(xs_list_ebar, en_lin_list):
    xs_formatted_list_ebar.append(ixs/(ien/1000.))
for ixs, ien in zip(xs_list_mubar, en_lin_list):
    xs_formatted_list_mubar.append(ixs/(ien/1000.))
for ixs, ien in zip(xs_list_taubar, en_lin_list):
    xs_formatted_list_taubar.append(ixs/(ien/1000.))

# We'll make a function for writing these files since we'll need to make
# several because we need separate ones for all neutrino types
# Well... at least we could.

def MakeXSFile(en_list, xs_lists, filename):
    outfile = open(filename, "w")
    outfile.write("# # nu-e elastic cross section (10^-38 cm^2/GeV) #   From Marciano et al doi.org/10.1088/0954-3899/29/11/013 #\n")
    outfile.write("# log(energy in GeV)  nu_e      nu_mu      nu_tau      nu_e_bar      nu_mu_bar    nu_tau_bar #\n\n")
    xs_formatted_list_e = xs_lists[0]
    xs_formatted_list_mu = xs_lists[1]
    xs_formatted_list_tau = xs_lists[2]
    xs_formatted_list_ebar = xs_lists[3]
    xs_formatted_list_mubar = xs_lists[4]
    xs_formatted_list_taubar = xs_lists[5]
    for i, ien in enumerate(en_list):
        enstr = "{:+6f}".format(ien)
        xsstr_e = "{:5e}".format(xs_formatted_list_e[i])
        xsstr_mu = "{:5e}".format(xs_formatted_list_mu[i])
        xsstr_tau = "{:5e}".format(xs_formatted_list_tau[i])
        xsstr_ebar = "{:5e}".format(xs_formatted_list_ebar[i])
        xsstr_mubar = "{:5e}".format(xs_formatted_list_mubar[i])
        xsstr_taubar = "{:5e}".format(xs_formatted_list_taubar[i])
        outfile.write("{0}  {1}  {2}  {3}  {4}  {5}  {6}\n".format(enstr, xsstr_e, xsstr_mu, xsstr_tau,
                                                                   xsstr_ebar, xsstr_mubar, xsstr_taubar))


xs_lists = (xs_formatted_list_e, xs_formatted_list_mu, xs_formatted_list_tau,
            xs_formatted_list_ebar, xs_formatted_list_mubar, xs_formatted_list_taubar)
MakeXSFile(en_formatted_list, xs_lists, "xs_nue_e_he.dat")


# if  you'd like to keep the code running to check parameter values, uncomment the below line!
#code.interact("hi", local=globals())
