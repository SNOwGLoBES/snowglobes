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
#include "supernova_mixing.h"

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

 infile.open(infilename.c_str());
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

     /// Write data to output file 
     //    // energies in output file need to be in GeV!
     if (mo>=0) {
       write_nh(E_nu,F,th12, outfile);
     } else {
       write_ih(E_nu,F,th12, outfile);
     }
   }	
   outfile.close();

}
