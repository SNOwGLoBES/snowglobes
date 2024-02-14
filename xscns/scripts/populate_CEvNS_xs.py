##############################################################################################################
#	This script, hosted in /xscns/scripts/ is intended to populate the xs_coh_helm_{element}_{flavor}
#   files which give the GLoBES-required cross-sections for CEvNS interactions. The total cross-section 
#   of interaction (as a function of incident neutrino energy) is calculated for each isotope, and the
#   weighted average is calculated using the molar fraction as the weighting factor. The weighted
#   averate is devided by the incident neutrino energy in GeV. The result is printed to the required
#   outfile for use in a SNOwGLoBES simulation.

#	Note the particular format of the cross-section - because of GLoBES, it is REQUIRED that the
#	cross-section be reported as (total cross section)/(incident neutrino energy in GeV) - the header
#	information in each file reports these units, but no documentation is available elsewhere.

#	If you need to change the energy range in the cross-sections file, you may do so by altering the
#   variable incidentenergies in Generate(). This should not be necessary in the course of normal operations.
##############################################################################################################


# Imports
import numpy as np
import warnings
import sys
warnings.filterwarnings('ignore')

##############################################################################################################
# Get the form factor from the command line
ffname = "helm" # Default form factor name

# If not enough inputs:
if len(sys.argv) != 2:
	print("Usage: python3.6 populate_CEvNS_xs.py [form factor name]")

# If enough inputs:
if len(sys.argv) == 2:
	ffname = str(sys.argv[1])
	ffname = ffname.lower()

##############################################################################################################
# Calculate the differential cross-section of interaction for one (incident energy, detected energy) pair
# for a specified element and isotope. Return the differential cross section of interaction for use in later
# functions.
def Differential(FormFactor,Er,Ei,N,M,Z):
    A = Z+N # Nucleon numbers
    
    # Constants -----------------------------------------------
    hbar = 4.136e-18/(2*np.pi) #keV*s
    c = 2.988e10 #cm/s
    G0 = 1.166e-17 #keV^-2
    Gf2 = G0**2 * (hbar*c)**6 #Fermi coupling constant squared
    # ----------------------------------------------------------

    kinetic = 2 - 2*Er/Ei + (Er/Ei)**2 - M*Er/(Ei*Ei)
    
    Qw = N - (1-4*0.231)*Z # weak neutral factor
    
    q = np.sqrt(2 * M * Er) # momentum transfer
    
    # Helm Form Factor ----------------------------------------------
    if FormFactor == "helm":

        rn0 = 5.79e-6 #1/keV
        s = 4.57e-6 #1/keV
        rn = rn0 * A**(1/3)
    
        bessel = np.sin(q*rn)/(q*rn)**2 - np.cos(q*rn)/(q*rn) # first order spherical Bessel Function in q*rn
        F = 3 * bessel* np.exp(-0.5 * (q*s)**2)/(q*rn) # form factor

        # Differential cross-section, full expression
        dxsecdEr = Gf2/(8*np.pi) * Qw**2 * M * kinetic * F**2 * (hbar*c)**(-4)
        return(dxsecdEr)
    # ----------------------------------------------------------

    # Klein-Nystrand Form Factor -------------------------------
    elif FormFactor == "klein-nystrand":
        rn0 = 5.79e-6 #1/kev
        s = 4.57e-6 #1/kev
        RA = rn0 * A**(1/3)
        ak = 2*s

        bessel = np.sin(q*RA)/(q*RA)**2 - np.cos(q*RA)/(q*RA)
        F = 3 * (bessel)/(q*RA) * (1)/(1 + q**2 * ak**2)

        dxsecdEr = Gf2/(8*np.pi) * Qw**2 * M * kinetic * F**2 * (hbar*c)**(-4)
        return(dxsecdEr)
    # ----------------------------------------------------------

    # Horowitz/Numerical Form Factor ---------------------------
    elif FormFactor == "horowitz":
        print("Please define the numerical recipe for the Horowitz form factor on line 87 before proceeding.")
        sys.exit()
    # ----------------------------------------------------------

    # Hoferichter Form Factor ----------------------------------
    elif FormFactor == "hoferichter":
        print("Please define the Hoferichter form factor on line 95 before proceeding.")
        sys.exit()
    # ----------------------------------------------------------

    else:
        print("Form Factor is not defined. Please define it.")
        sys.exit()


##############################################################################################################
# Calculate the total cross section for a given incident neutrino energy.
# The range of recoil energies for CEvNS interactions is given in keV units.
def TotalCross(FormFactor, Ei, N, M, Z):
    recoilenergies = np.linspace(0,10,1000)
    
    dcrossarray = Differential(FormFactor,recoilenergies,Ei,N,M,Z)
    mask = dcrossarray > 0
    totalcross = np.sum(dcrossarray[mask]) * (recoilenergies[1] - recoilenergies[0])
    return(totalcross)


##############################################################################################################
# Use the TotalCross() function to return the total cross sections and write them to the coherent outfile
# If you wish to set a different energy range, do so in the "incidentenergies" variable.
# The numpy.logspace(start, stop, step, num, base) function returns evenly-separated values between 
# base**start and base**stop inclusive of stop.

def Generate(FormFactor,elementname,N,Z,M,weight):
    # Incident neutrino energies
    incidentenergies = np.logspace(-3.301030, -1, 1002, base=10) #GeV logspace energies
    fmtEnergies = ["{0:.6f}".format(np.log10(i)) for i in incidentenergies] # log_10 format required by GLoBES
    kevIncident = [i*1e6 for i in incidentenergies] #keV energies

    # Calculate and format total cross-sections
    TotalCrossSections = []

    for i in range(len(kevIncident)): # For each incident energy:
        Tcross_total = 0
        for j in range(len(N)): # for each isotope:
            Tcross_unweight = TotalCross(FormFactor,kevIncident[i], N[j], M[j], Z[j]) #calculate the total xsec
            Tcross_weight = Tcross_unweight * weight[j] #weight the total xsec with the molar fraction for that isotope
            Tcross_total += Tcross_weight #add the weighted cross section to the final result
        TotalCrossSections.append(Tcross_total) # Append the total cross-section to the list of total cross-sections
    
    TotalCross_divEnergy = TotalCrossSections / incidentenergies    
    fmtTotalCross = ["{:.6e}".format(i*1e38) for i in TotalCross_divEnergy] # 1e-38/[energy in GeV] required by GLoBES
    
    # Print to the outfiles
    neutrinoname = ["nue","numu","nutau","nuebar","numubar","nutaubar"]

    for flavor in neutrinoname:
        with open(f"./../xs_coh_{FormFactor}_{flavor}_{elementname}.dat","w") as writer:
            # Header information
            writer.write(f"# Neutrino-{elementname} CEvNS cross section (10^-38 cm^2/GeV, natural abundance) # \n")
            writer.write("# log(energy in GeV)  nu_e      nu_mu      nu_tau      nu_e_bar      nu_mu_bar    nu_tau_bar # \n")
            writer.write("\n")
            # Data
            for i in range(1,1002):
                writer.write(f"{fmtEnergies[i]}   {fmtTotalCross[i]}   {fmtTotalCross[i]}   {fmtTotalCross[i]}   {fmtTotalCross[i]}   {fmtTotalCross[i]}   {fmtTotalCross[i]}\n")

##############################################################################################################
# Function calls - script runs from here. You can any or all of these calls, just comment out the undesirables.

# Isotope parameter information
Ar_Ns = [18,20,22]
Ar_Zs = [18,18,18]
Ar_masses = [3.350e7, 3.536e7, 3.722e7]
Ar_weights = [0.003336,0.000629,0.996035]

Ge_Ns = [38,40,41,42,44]
Ge_Zs = [32,32,32,32,32]
Ge_masses = [6.513e7,6.699e7,6.793e7,6.886e7,7.072e7]
Ge_weights = [0.2057,0.2745,0.0775,0.3650,0.0773]

Xe_Ns = [70,72,74,75,76,77,78,80,82]
Xe_Zs = [54,54,54,54,54,54,54,54,54]
Xe_masses = [1.154e8,1.173e8,1.191e8,1.201e8,1.210e8,1.219e8,1.229e8,1.247e8,1.266e8]
Xe_weights = [0.000952,0.000890,0.019102,0.264006,0.040710,0.21232,0.269086,0.104257,0.088573]

# Function calls
Generate(ffname,"Ar",Ar_Ns,Ar_Zs,Ar_masses,Ar_weights)
Generate(ffname,"Ge",Ge_Ns,Ge_Zs,Ge_masses,Ge_weights)
Generate(ffname,"Xe",Xe_Ns,Xe_Zs,Xe_masses,Xe_weights)
