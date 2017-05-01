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

// If the argument is non-zero, the output will create output files using 
// the simplified MSW assumption,  in $OUTFLUXDIR/nh and $OUTFLUXDIR/ih




#include <iostream>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string> 
#include <sstream>

using namespace std;

ifstream infile;
ofstream outfile;
ofstream outfile_nh;
ofstream outfile_ih;

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
void write_nh(double a, double B[], double th12);
void write_ih(double a, double B[], double th12);
double phi(double E_nu, double E0, double alpha);

//---------------------------------------------------------------------------------------------------------------------------------------------------------

void usage() {
  // print usage information
  printf("\n pinched  by K. Scholberg, N. Kaiser \n");
  printf("Usage:\n");
  printf("./pinched th12[rad]\n");
  printf("If th12 absent, no MSW assumed; if present creates NH/IH directories\n");
  
}


int main(int argc, char **argv){


  if (argc > 2) {
  usage();
  exit(-1);
}

double th12=0.;
if (argc==2) {
  th12 = atof(argv[1]);
  cout << "Assuming MSW with th12= "<<th12<<" radians"<<endl;

 } else {
  cout << "No oscillations assumed"<<endl;
 }


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
		  cout << "I think you want to set OUTFLUXDIR" <<endl;
		  exit(-1);
//		  ss << "pinched_" << i << ".dat";
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


		if (fabs(th12)>0) {
		  // Make the MSW files

		    string filename_nh;
		    std::stringstream ss_nh;
		    
		    ss_nh << getenv("OUTFLUXDIR") << "/nh/pinched_"<<i << ".dat";

		    filename_nh+=ss_nh.str();
		    cout << "Output file: "<<filename_nh<<endl;
		    cout << "--------------"<<endl;
		    outfile_nh.open(filename_nh.c_str());

		    if (!outfile_nh.good()) {
		      cout << "Outfile "<<filename_nh<<" not opened..."<<endl;
		      cout << "Check that directory pointed to by OUTFLUXDIR environment variable exists"<<endl;
		      exit(0);
		    }

		    string filename_ih;
		    std::stringstream ss_ih;
		    
		    ss_ih << getenv("OUTFLUXDIR") << "/ih/pinched_"<<i << ".dat";

		    filename_ih+=ss_ih.str();
		    cout << "Output file: "<<filename_ih<<endl;
		    cout << "--------------"<<endl;
		    outfile_ih.open(filename_ih.c_str());

		    if (!outfile_ih.good()) {
		      cout << "Outfile "<<filename_ih<<" not opened..."<<endl;
		      cout << "Check that directory pointed to by OUTFLUXDIR environment variable exists"<<endl;
		      exit(0);
		    }




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
		  if (fabs(th12)>0) {
		    write_nh(E_nu,F,th12);
		    write_ih(E_nu,F,th12);
		  }

		  E_nu+=estep;
		}	
		outfile.close();

		if (fabs(th12)>0) {
		  outfile_nh.close();
		  outfile_ih.close();
		}
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

  // No flavor transitions assumed
  // Input:  B[0]: nue, B[1]: nuebar, B[2]: nux (any one of them)
  // Output flux file: first neutrinos, e, mu, tau, then antinus, ebar, mubar, taubar

  outfile << setw(8) << a << "\t " ;
  outfile << setw(8) << B[0] << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << setw(8) << B[1] << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;

  outfile << endl;

}


// Function to write energy and fluxes into fluxfile 
void write_nh(double a, double B[], double th12){

  // Input:  B[0]: nue, B[1]: nuebar, B[2]: nux (any one of them, all equal)
  // Output flux file: first neutrinos, e, mu, tau, then antinus, ebar, mubar, taubar

  // Normal ordering:
  // nue=nux0
  // nuebar = cos^2th12 nuebar0 + sin^2th12 nux0
  // Pee = p = 0
  // Peebar = pbar = cos^2th12 

  // Neutrinos:
  // numu+nutau = (1-p)nue0 + (1+p)nux0 
  //  so numu=nutau = (nue0+nux0)/2
   
  // Antineutrinos
  // numubar + nutaubar = ((1-pbar)nuebar0 + (1+pbar)nuxbar0)
  //  so numubar = nutaubar = (sin^2th12 nuebar0 + (1+cos^2th12) nuxbar0)/2
  
  //For  th12 = 0.588336
  //double s2th12 = 0.308;
  //double c2th12 = 0.692;

  double s2th12;
  double c2th12;
  s2th12 = pow(sin(th12),2);
  c2th12 = 1-s2th12;

  outfile_nh << setw(8) << a << "\t " ;
  outfile_nh << setw(8) << B[2] << "\t " ;
  outfile_nh << setw(8) << (B[0]+B[2])/2. << "\t " ;
  outfile_nh << setw(8) << (B[0]+B[2])/2. << "\t " ;
  outfile_nh << setw(8) << c2th12*B[1]+s2th12*B[2] << "\t " ;
  outfile_nh << setw(8) << ((1.-c2th12)*B[1]+(1+c2th12)*B[2])/2. << "\t " ;
  outfile_nh << setw(8) << ((1.-c2th12)*B[1]+(1+c2th12)*B[2])/2. << "\t " ;
  outfile_nh << endl;

}


// Function to write energy and fluxes into fluxfile 
void write_ih(double a, double B[], double th12){


  // Input:  B[0]: nue, B[1]: nuebar, B[2]: nux (any one of them, all equal)
  // Output flux file: first neutrinos, e, mu, tau, then antinus, ebar, mubar, taubar

  // Inverted ordering:
  // nuebar=nuxbar0
  // nue = sin^2th12 nuebar0 + cos^2th12 nux0
  // Pee = p = sin^2th12
  // Peebar = pbar = 0

  // Neutrinos
  // numu + nutau = ((1-p)nue0 + (1+p)nux0)
  //  so numu = nutau = (cos^2th12 nue0 + (1+sin^2th12) nux0)/2

  // Antineutrinos:
  // numubar+nutaubar = (1-pbar)nuebar0 + (1+pbar)nuxbar0 
  //  so numubar=nutaubar = (nuebar0+nuxbar0)/2
   
  
  // th12 = 0.588336

  //  double s2th12 = 0.308;
  //double c2th12 = 0.692;

  double s2th12;
  double c2th12;
  s2th12 = pow(sin(th12),2);
  c2th12 = 1-s2th12;

  // First neutrinos then antinus
  outfile_ih << setw(8) << a << "\t " ;
  outfile_ih << setw(8) << s2th12*B[0]+c2th12*B[2] << "\t " ;
  outfile_ih << setw(8) << ((1-s2th12)*B[0]+(1+s2th12)*B[2])/2. << "\t " ;
  outfile_ih << setw(8) << ((1-s2th12)*B[0]+(1+s2th12)*B[2])/2. << "\t " ;
  outfile_ih << setw(8) << B[2] << "\t " ;
  outfile_ih << setw(8) << (B[1]+B[2])/2. << "\t " ;
  outfile_ih << setw(8) << (B[1]+B[2])/2. << "\t " ;
  outfile_ih << endl;


}






