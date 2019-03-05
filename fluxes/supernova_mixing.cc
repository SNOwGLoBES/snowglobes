// See header for info on sources of these formulae

#include <iomanip>
#include <fstream>
#include <math.h>
#include "supernova_mixing.h"
using namespace std;

// Function to write energy and fluxes into fluxfile 
void write(double a, double B[], ofstream& outfile){

  // No flavor transitions assumed
  // Input:  B[0]: nue, B[1]: nuebar, B[2]: nux (any one of them)
  // Output flux file: first neutrinos, e, mu, tau, then antinus, ebar, mubar, taubar
  // This is basically writing out the flux before mixing

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
void write_nh(double a, double B[], double th12, ofstream& outfile){

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
void write_ih(double a, double B[], double th12, ofstream& outfile){


  // Input:  B[0]: nue, B[1]: nuebar, B[2]: nux (any one of them, all equal)
  // Output flux file: first neutrinos, e, mu, tau, then antinus, ebar, mubar, taubar

  // Inverted ordering:
  // nuebar=nuxbar0
  // nue = sin^2th12 nue0 + cos^2th12 nux0
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



