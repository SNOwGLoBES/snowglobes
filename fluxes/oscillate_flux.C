// Program that applies MSW flavor transformation to an input flux
// sorry, A. Smirnov, I am using the word "oscillate" occasionally...
//  WARNING: assumes collective effects are unimportant
// K. Scholberg; some code from taken pinched.C (they should
//   really both read from a common library; todo
// 
// EARTH MATTER EFFECTS NOT WORKING YET
// Arguments:  
// mass ordering: "normal" or "inverted"
// th12
// pathlength through Earth in km
//  Output file named according to these; if Earth matter pathlength is zero, omits



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

const int maxpts = 1000;
double E_nu;
double E_nu_orig[maxpts];
double nue_unosc[maxpts];
double numu_unosc[maxpts];
double nutau_unosc[maxpts];
double nuebar_unosc[maxpts];
double numubar_unosc[maxpts];
double nutaubar_unosc[maxpts];

// [E]=GeV!!

double estep=0.0002;

double F[3];

void write(double a, double B[]);
void write_nh(double a, double B[], double th12);
void write_ih(double a, double B[], double th12);


//---------------------------------------------------------------------------------------------------------------------------------------------------------

void usage() {
  // print usage information
  printf("\n oscillate_flux  by K. Scholberg, \n");
  printf("Usage:\n");
  printf("./oscillate_flux fluxname MO(1=normal or -1=inverted) th12[rad] L in Earth[km]\n");
  printf("e.g. ./oscillate_flux livermore 0 0.588366 8000.\n");
}


int main(int argc, char **argv){


  if (argc != 5) {
  usage();
  exit(-1);
}

  // Get the input parameters
double th12=0.;
 
string influxname;
influxname = argv[1];
std::stringstream ss;
ss<<influxname<<".dat";
string infilename = ss.str();

cout << "infilename "<<infilename<<endl;
int mo = atof(argv[2]);
 
string ordering;
ordering = (mo>=0 ? "normal" : "inverted");
cout << "Assumed mass ordering: "<<mo <<" "<<ordering<<endl;
th12 = atof(argv[3]);
cout << "Assuming MSW with th12= "<<th12<<" radians"<<endl;
double pathlength=0.;
 string pathstring;
 pathstring = argv[4];
pathlength = atof(argv[4]);
cout << "Assuming pathlength in Earth=  "<<pathlength<<" km"<<endl;

// Read the input file
 
int i=0;

 infile.open(infilename);
 if (!infile.good()) {
   cout << "Can't open input file "<<infilename<<endl;
   exit(0);
 }
	// Loop over all lines of input files
 while ( infile.good() ){

   infile>>E_nu_orig[i]>>nue_unosc[i]>>numu_unosc[i]>>nutau_unosc[i]>>nuebar_unosc[i]>>numubar_unosc[i]>>nutaubar_unosc[i]; 
   
   cout << i<<" "<<E_nu_orig[i]<<" "<<nue_unosc[i]<<" "<<numu_unosc[i]<<" "<<nutau_unosc[i]<<" "<<nuebar_unosc[i]<<" "<<numubar_unosc[i]<<" "<<nutaubar_unosc[i]<<endl;


   if (!infile.good()) break;

   i++;
 }

 int numebin=i;
 cout << "Found "<<i<<" points"<<endl;

 infile.close();

   // E must be in GeV
   

 std::stringstream ss2;
 
 if (pathlength>0) {
   ss2<<influxname<<"_"<<ordering<<"-"<<pathstring<<"-km.dat";
 }
 else {
   ss2<<influxname<<"_"<<ordering<<".dat";
 }

 string outfilename;

 outfilename=ss2.str();

 cout << "Output file: "<<outfilename<<endl;


 
 cout << "--------------"<<endl;
   
 outfile.open(outfilename.c_str());
		
 if (!outfile.good()) {
     cout << "Outfile "<<outfilename<<" not opened..."<<endl;
     exit(0);
 }


		// File should have flux in the 0.0002 GeV bin, so multiply by bin size

   for(int i=0; i<=numebin; i++) {
		// Create fluxes F^0

     E_nu = E_nu_orig[i];
     // Assumes mu, tau, mubar, taubar already all equal
     F[0] = nue_unosc[i];
     F[1] = nuebar_unosc[i];
     F[2] = numu_unosc[i];

     // Write data to output file 
     //    // energies in output file need to be in GeV!
     if (mo>=0) {
       write_nh(E_nu,F,th12);
     } else {
       write_ih(E_nu,F,th12);
     }
   }	
   outfile.close();

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------



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

  outfile << setw(8) << a << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << setw(8) << (B[0]+B[2])/2. << "\t " ;
  outfile << setw(8) << (B[0]+B[2])/2. << "\t " ;
  outfile << setw(8) << c2th12*B[1]+s2th12*B[2] << "\t " ;
  outfile << setw(8) << ((1.-c2th12)*B[1]+(1+c2th12)*B[2])/2. << "\t " ;
  outfile << setw(8) << ((1.-c2th12)*B[1]+(1+c2th12)*B[2])/2. << "\t " ;
  outfile << endl;

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
  outfile << setw(8) << a << "\t " ;
  outfile << setw(8) << s2th12*B[0]+c2th12*B[2] << "\t " ;
  outfile << setw(8) << ((1-s2th12)*B[0]+(1+s2th12)*B[2])/2. << "\t " ;
  outfile << setw(8) << ((1-s2th12)*B[0]+(1+s2th12)*B[2])/2. << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << setw(8) << (B[1]+B[2])/2. << "\t " ;
  outfile << setw(8) << (B[1]+B[2])/2. << "\t " ;
  outfile << endl;


}






