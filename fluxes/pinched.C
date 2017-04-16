// Program that calculates pinched spectrum flux input files for snowglobes from input file containing parameters for the Garching quasi-thermal parameterization
// K. Scholberg; adapted from code written by Nicolas Kaiser, summer 2011
// It reads the file pinched_info.dat, and outputs one flux file per line
//      number  alpha_nu_e  alpha_nubar_e  alpha_nu_x   Eavg_nu_e   Eavg_nubar_e   Eavg_nu_x  lum_nu_e lum_nubar_e lum_nu_x
// Energies are in MeV
// Flux for nux should be for one flavor; assumes all nux fluxes are equal
// Luminosity is in ergs/sec... actually it's really ergs for the usual case 
// when the file represents fluence rather than "flux" and is integrated over
// a time bin
// 
// Flux file output goes to OUTFLUXDIR if set, otherwise to working directory



#include <iostream>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string> 
#include <sstream>

using namespace std;

ofstream outfile;
ifstream infile;

// arrays for fluxes, average energies, pinching and luminosities of neutrinos : nue, nuebar, nux
// Luminosity is for a single nux flavor

double F[3];
double E0[3];
double L[3];
double alpha[3];

double E_nu;
// [E]=GeV!!

const double dist=3.08568025e22; // [dist]=cm
const double gevpererg = 624.15;
double estep=0.0002;

void write(double a, double B[]);
double phi(double E_nu, double E0, double alpha);

//---------------------------------------------------------------------------------------------------------------------------------------------------------

int main(){

	int i;
	infile.open("pinched_info.dat");
	if (!infile.good()){
	  cout << "Can't open pinched_info.dat"<<endl;
	  exit(EXIT_FAILURE);
	}

	// Loop over all lines of input files
	while ( infile.good() ){

		infile >> i; 
		infile >> alpha[0];  
		infile >> alpha[1];  
		infile >> alpha[2];  
		infile >> E0[0];
		infile >> E0[1];
		infile >> E0[2];
		infile >> L[0];
		infile >> L[1];
		infile >> L[2];
		if (!infile.good()) break;

		// E must be in GeV
		for(int k=0; k<3; k++){
			E0[k]/=1000.;
		}

		cout << "Flux " <<i << ": alpha: "<<alpha[0]<<" "<<alpha[1]<<" "<<alpha[2]<<endl;

		cout << "Flux "<< i << ": E0: "<<E0[0]<<" "<<E0[1]<<" "<<E0[2]<<endl;

		cout << "Flux "<< i << ": Luminosity: "<<L[0]<<" "<<L[1]<<" "<<L[2]<<endl;

		// Convert luminosity to GeV/s

		L[0] *= gevpererg;
		L[1] *= gevpererg;
		L[2] *= gevpererg;

		// create filename for output file and open respective file 
		//		string filename="pinched_";
		string filename;
		std::stringstream ss;

		if (getenv("OUTFLUXDIR")==NULL){
		  ss << "pinched_" << i << ".dat";
		} else	{
		  ss << getenv("OUTFLUXDIR") << "/pinched_"<<i << ".dat";
		}
		    
		filename+=ss.str();
		cout << "Output file: "<<filename<<endl;
		cout << "--------------"<<endl;
		outfile.open(filename.c_str());

		if (!outfile.good()) {
		  cout << "Outfile "<<filename<<" not opened..."<<endl;
		  cout << "Check that directory pointed to by OUTFLUXDIR environment variable exists"<<endl;
		    exit(0);
		}

		// File should have flux in the 0.0002 GeV bin, so multiply by bin size
		E_nu=0;
		for(int i=0; i<=500; i++) {
		// Create fluxes F^0
		  for(int j=0; j<3; j++){

		    if (E0[j]>0.) {
		      F[j]=1/(4*M_PI*dist*dist)*L[j]/E0[j]*phi(E_nu,E0[j],alpha[j])*estep;  
		    } else {
		      F[j]=0.;
		    }
		  }
			// Write data to output file 
		  write(E_nu,F); // energies in output file need to be in GeV!
		  E_nu+=estep;
		}	
		outfile.close();
	}
	infile.close();
	return 0; 
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------



// Function to calculate the energy spectrum
double phi(double E_nu, double E_nu0, double alpha) {
	double N=pow((alpha+1.),(alpha+1.))/(E_nu0*tgamma(alpha+1.));
	double R=N*pow((E_nu/E_nu0),alpha)*exp((-1.)*(alpha+1.)*E_nu/E_nu0); 

	return R;
}

// Function to write energy and fluxes into fluxfile 
void write(double a, double B[]){

  // First neutrinos then antinus
  outfile << setw(8) << a << "\t " ;
  outfile << setw(8) << B[0] << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << setw(8) << B[1] << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << endl;

}




