// cevns_recoil_response.cc was developed to integrate the functionality of DukeCEvNS with SNOwGLoBES. In its
// v. 1.2 form, SNOwGLoBES did not support CEvNS interactions - the DukeCEvNS code, as the name
// suggests, does.

// Although cevns_recoil_response.cc is not intended to run autonomously of SNOwGLoBES, it can for debugging
// purposes. cevns_recoil_response.cc accepts as inputs the detector isotope(s), form factor name, energy 
// range of incident neutrinos, and energy range in which the detector is sensitive to CEvNS recoils. Its 
// output is a matrix of differential cross sections as a function of recoil energy, which is passed into
// create_smearing_matrix.py in the parent SNOwGLoBES folder, normalized, and returned in a
// SNOwGLoBES-compatible format with detector smearing applied.

// v. 1.0, finalized April 8th, 2021, supports the HELM, KLEIN-NYSTRAND, and HOROWITZ (numerical) form 
// factors, which were available in DukeCEvNS at the time of development. Implementing the chiral EFT-based 
// HOFERICHTER form factor will require first integrating the form factor with the rest of the DukeCEvNS code
// and setting up the appropriate pointers on line 312 of this file.


#include <iostream>
#include <fstream>

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <math.h>
#include <map>

#include "xscns.h"
#include "DetectorResponse.h"
#include "FormFactor.h"
#include "NuFlux.h"

#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

// Begin main
int main(int argc, char * argv[]) {

// Set up the outfile ----------------------------------------------------------------

	// Make the outfile with its horribly long but searchable name, and open it
	std::ofstream outfile;
  	std::string outfilename;
	outfilename = "./recoil_response_matrices/unsmeared_CEvNS_"+std::string(argv[1])+"_"+std::string(argv[2])+"_inLOW"+std::string(argv[3])+"_inHIGH"+std::string(argv[4])+"_inBINS"+std::string(argv[5])+"_detLOW"+std::string(argv[6])+"_detHIGH"+std::string(argv[7])+"_detBINS"+std::string(argv[8])+".matrix";
	outfile.open(outfilename);

// Collect arguments, set up mixture information -------------------------------------

	// Isotope mixture information
	#include "isomaps.h"
	#include "mixtures.h"
	std::string material = argv[1];

	// Form factor
	std::string ffname = argv[2];
	
	// Set up binning parameters
	float incidentLowThresh = strtof(argv[3], NULL); //enu_min
	float incidentHighThresh = strtof(argv[4], NULL); //enu_max
	int incidentBins = strtof(argv[5], NULL); //n_points_enu
	float detectedLowThresh = strtof(argv[6], NULL); //edet_min
	float detectedHighThresh = strtof(argv[7], NULL); //edet_max
	int detectedBins = strtof(argv[8], NULL); //n_points_edet

// Validate inputs -------------------------------------------------------------------
	
	// Should have correct number of arguments
	if (argc<8) {
		std::cout << "Incompatible inputs in dukecevns/cevns_recoil_response.cc. Usage: ./smearing [target material] [enu_min] [enu_max] [n_points_enu] [edet_min] [edet_max] [n_points_edet]" << std::endl;
		exit(0);
	}

	// Energy range must have correct limit order
	if (incidentLowThresh > incidentHighThresh) {
		std::cout << "Incorrect order of incident neutrino energy limits - enu_max must be greater than or equal to enu_min." << std::endl;
		exit(0);
	}

	if (detectedLowThresh > detectedHighThresh) {
		std::cout << "Incorrect order of detected recoil energy limits - edet_max must be greater than or equal to edet_min." << std::endl;
		exit(0);
	}

	// All inputs must be positive
	if (incidentLowThresh < 0) {
		std::cout << "All inputs to cevns_recoil_response.cc must be positive." << std::endl;
		exit(0);
	}

	if (incidentHighThresh < 0) {
		std::cout << "All inputs to cevns_recoil_response.cc must be positive." << std::endl;
		exit(0);
	}

	if (incidentBins < 0) {
		std::cout << "All inputs to cevns_recoil_response.cc must be positive." << std::endl;
		exit(0);
	}

	if (detectedLowThresh < 0) {
		std::cout << "All inputs to cevns_recoil_response.cc must be positive." << std::endl;
		exit(0);
	}

	if (detectedHighThresh < 0) {
		std::cout << "All inputs to cevns_recoil_response.cc must be positive." << std::endl;
		exit(0);
	}

	if (detectedBins < 0) {
		std::cout << "All inputs to cevns_recoil_response.cc must be positive." << std::endl;
		exit(0);
	}


// Make matrices ---------------------------------------------------------------------

	// Enu_step
	float incidentStep = (incidentHighThresh - incidentLowThresh)/(incidentBins - 1);

	// Edet_step
	float detectedStep = (detectedHighThresh - detectedLowThresh)/(detectedBins);

	// Allowed Edet values
	float Edet[detectedBins];
	for (int i = 0; i < detectedBins; i++) { //Fill array with allowed Edet
		Edet[i] = detectedLowThresh + (i*detectedStep);
	}

	// Allowed Enu values
	float Enu[incidentBins];
	for (int j = 0; j < incidentBins; j++) { //Fill array with allowed Enu
		Enu[j] = incidentLowThresh + (j*incidentStep);
	}

// Set up material -------------------------------------------------------------------

	std::vector<double> fraction = molar_fraction[material];
	std::vector<std::string> isotope_component= isotopes[material];

// Set up form factor ----------------------------------------------------------------

	int noff = 0.; // Controls ffpv, ffpa, ffnv, ffna values in upcoming loop
	if (ffname == "Unity") { // No form factor
		noff = 1.; 
	}

	// Initialize constants
	double M;
	double Delta;
	int Nn,Z,A;
	//int Zdiff, Ndiff;

	// Array of pointers to form factors
	// Protons (ffp*) and neutrons (ffn*) separated
	// Axial (*a) and vector (*v) separated
	FormFactor** ffpv;
	ffpv = new FormFactor*[max_components];
	FormFactor** ffpa;
	ffpa = new FormFactor*[max_components];

	FormFactor** ffnv;
	ffnv = new FormFactor*[max_components];
	FormFactor** ffna;
	ffna = new FormFactor*[max_components];


	// Iterate over isotopes in material, set up the form factor for each and 
	// store them in memory for later use.
	int is=0;
	std::vector<std::string>::iterator v = isotope_component.begin(); //define iterator
	std::string isotope;

	double Mtot = 0; // total mass
	v = isotope_component.begin(); // isotope iterator

	double minM = 1.e10; // minimum mass
	while( v != isotope_component.end()) { // For each isotope in material,

		isotope = *v;
		std::string isoname = std::string(isotope);

		// Collect the constants
		Z = Zs[std::string(isotope)];
		Nn = Ns[std::string(isotope)];
		Delta = Deltas[std::string(isotope)];
		M = (Z+Nn)*amu - Z*me + Delta; //MeV
		A = Nn + Z;
		
		if (M<minM) {minM=M;} // Set the minimum mass to be the mass of the isotope
		Mtot += M*fraction[is]; // Add isotope mass to weighted total


		// Default constants for all form factors
		double nvrfact = 1.0;
		double narfact = 1.0;
		double pvrfact = 1.0;
		double parfact = 1.0;


    		if (ffname == "Helm") { // set up the Helm form factor
          
			// Constants
			double nvsfact = 0.9;
			double nasfact = 0.9;
			double pvsfact = 0.9;
			double pasfact = 0.9;

			// Pointers
			Helm* helmffnv= new Helm();
			ffnv[is] = helmffnv;
			helmffnv->Setsval(nvsfact);
			helmffnv->SetRfac(nvrfact);

			Helm* helmffna= new Helm();
			ffna[is] = helmffna;
			helmffna->Setsval(nasfact);
			helmffna->SetRfac(narfact);

			Helm* helmffpv= new Helm();
			ffpv[is] = helmffpv;
			helmffpv->Setsval(pvsfact);
			helmffpv->SetRfac(pvrfact);


			Helm* helmffpa= new Helm();
			ffpa[is] = helmffpa;
			helmffpa->Setsval(pasfact);
			helmffpa->SetRfac(parfact);
		} // Close Helm if-statement

		else if (ffname == "Klein-Nystrand") { // set up the Klein-Nystrand form factor

			//Constants
			double nvak = 0.7;
			double naak = 0.7;
			double pvak = 0.7;
			double paak = 0.7;
			double nvskin = 0.0;
			double naskin = 0.0;
			double pvskin = 0.0;
			double paskin = 0.0;

			//Pointers
			Klein* kleinffnv = new Klein();
			ffnv[is] = kleinffnv;
			kleinffnv->Setakval(nvak);
			kleinffnv->SetRfac(nvrfact);
			kleinffnv->Setskinfac(nvskin);

			Klein* kleinffna = new Klein();
			ffna[is] = kleinffna;
			kleinffna->Setakval(naak);
			kleinffna->SetRfac(narfact);
			kleinffna->Setskinfac(naskin);

			Klein* kleinffpv = new Klein();
			ffpv[is] = kleinffpv;
			kleinffpv->Setakval(pvak);
			kleinffpv->SetRfac(pvrfact);
			kleinffpv->Setskinfac(pvskin);

			Klein* kleinffpa = new Klein();
			ffpa[is] = kleinffpa;
			kleinffpa->Setakval(paak);
			kleinffpa->SetRfac(parfact);
			kleinffpv->Setskinfac(paskin);
		} // Close Klein if-statement
 
		else if  (ffname =="Horowitz"){ // set up the Horowitz/Numerical form factor

			//Pointers
			Horowitz* horowitzffnv = new Horowitz();
			ffnv[is] = horowitzffnv;

			Horowitz* horowitzffna = new Horowitz();
			ffna[is] = horowitzffna;

			Horowitz* horowitzffpv = new Horowitz();
			ffpv[is] = horowitzffpv;

			Horowitz* horowitzffpa = new Horowitz();
			ffpa[is] = horowitzffpa;

			// Query Horowitz files for data and set values
			std::transform(isoname.begin(), isoname.end(),isoname.begin(), ::toupper);
			std::string horowitz_filename = isoname+".FF";

			//Neutrons
			horowitzffnv->SetFFfilename(horowitz_filename.c_str());
			horowitzffnv->ReadFFfile();
			horowitzffnv->SetRfac(nvrfact);

			horowitzffna->SetFFfilename(horowitz_filename.c_str());
			horowitzffna->ReadFFfile();
			horowitzffna->SetRfac(narfact);

			//Protons (not the most accurate)
			horowitzffpv->SetFFfilename(horowitz_filename.c_str());
			horowitzffpv->ReadFFfile();
			horowitzffpv->SetRfac(pvrfact);

			horowitzffpa->SetFFfilename(horowitz_filename.c_str());
			horowitzffpa->ReadFFfile();
			horowitzffpa->SetRfac(parfact);
		}// Close Horowitz if-statement

		// else if (ffname =="Hoferichter"){ // set up the Hoferichter/Chiral EFT ff
			// Define pointers here
		//}

		if (noff == 0) { //Set some constants for non-unity form factors
	 		ffnv[is]->SetA(A);
			ffna[is]->SetA(A);
			ffpv[is]->SetA(A);
			ffpa[is]->SetA(A);
      
			ffnv[is]->SetZ(Z);
			ffna[is]->SetZ(Z);
			ffpv[is]->SetZ(Z);
			ffpa[is]->SetZ(Z);
		}


		v++; is++; // Next isotope
	} // Close loop over isotopes

// With the form factor set up and stored in memory for each isotope, we can move on
// to calculating the differential cross-section of interaction.
// -----------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------

// Loop over components in the given isotope mixture and, for each isotope, at the
// specified (incident energy, detected energy) pair, calculate the differential cross
// section. Weight by natural abundance and sum. This is the matrix entry.

// -----------------------------------------------------------------------------------
	
	// Define constants
	double hbar = 4.136e-21/(2 * M_PI); //Reduced Planck Constant, MeV*s
	double c = 2.998e10; //speed of light, cm/s
	double G0 = 1.166e-11; //MeV^-2
	double Gf2 = pow(G0,2) * pow((hbar*c),6); //Fermi coupling constant squared
	double sin2thetaW = 0.231;
	double hbarc_MeVfm = 197.327; // MeV-fm, conversion factor for MeV momentum transfer in ff

	//Array of n_points_edet rows and n_points_enu columns
	float Unsmeared[detectedBins][incidentBins];

	//Populate Unsmeared with differential cross-section values
	for (int i = 0; i < detectedBins; i++) {
		for (int j = 0; j < incidentBins; j++) { // At a particular (detected energy = i, incident energy = j) pair:
			
			// Loop over isotope components and get parameters for dCross
			// Calculate the weighted differential cross section and set it to be the (i,j)th entry
			double dCross_Total = 0; //Initial differential cross section entry at matrix location (i,j)

			// From mixtures.h
			std::vector<double> fraction = molar_fraction[material]; // Weight factor
			std::vector<std::string> isotope_component= isotopes[material]; // Component isotopes of mixture, will loop over this

			int is=0; // Isotope loop iteration index
			std::vector<std::string>::iterator v = isotope_component.begin(); // Define iteration function
			std::string isotope;

			v = isotope_component.begin(); // Call iteration function
				
			// Begin loop over isotopes
			while( v != isotope_component.end()) {
				isotope = *v;

				// Constants from mixtures.h/isomaps.h
				Z = Zs[std::string(isotope)];
    				Nn = Ns[std::string(isotope)];
    				Delta = Deltas[std::string(isotope)];
    				M = (Z+Nn)*(amu) - Z*(me) + Delta; //MeV/c^2
				mass_fraction[is] = M/Mtot*fraction[is]; //How much of the total mass is made up of this isotope

				double Qw = Nn - (1 - 4*sin2thetaW)*Z; //Weak neutral factor

				double kinetic = 2 - (2*Edet[i])/Enu[j] + pow((Edet[i]/Enu[j]),2) - (M*Edet[i])/pow(Enu[j],2); //kinetic energy transfer term

				double Q = sqrt(2 * M * Edet[i]); //Momentum transfer, MeV
				double qq = Q/hbarc_MeVfm; //Momentum transfer, fm^(-1)

				// Form Factor
				double FormFactor;
				if (noff == 0) {
					FormFactor = ffnv[is]->FFval(qq);
				}
				else {
					FormFactor = 1.;
				}

				//Calculate the unweighted differential cross section, weight it, and add it to the total differential cross section
				double dCross_Unweight = Gf2/(8*M_PI) * Qw * Qw * M * kinetic * FormFactor * FormFactor * pow((hbar*c),-1*4);
				double dCross_Weight = dCross_Unweight * fraction[is]; // Weight the differential cross section

				dCross_Total += dCross_Weight; //Sum the weighted cross sections into the total cross section


			v++; is++;
			} //Close loop over isotopes

			// Set the (i,j)-th entry of Unsmeared
			if (dCross_Total >= 0) //A physical differential cross section
				Unsmeared[i][j] = dCross_Total;	
			else //Negative differential cross sections are unphysical and should return zero
				Unsmeared[i][j] = 0;

			// Send Unsmeared to the outfile, space-delimited
			outfile << Unsmeared[i][j] << " ";

		} //Close loop over incident neutrino energies

	outfile << "\n"; //Newline in output matrix at end of row

	} //Close loop over detected neutrino energies
// -----------------------------------------------------------------------------------

// Exit main
	return 0;
}
