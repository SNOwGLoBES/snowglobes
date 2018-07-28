// supernova_mixing.h
// contains shared functions for computing neutrino mixing in the matter of a supernova,
// and outputting the results.  Note that the approximations used here are only valid
// for the early part of a supernova.  The formulae come from:
// Supernova Neutrinos: Production, Oscillations and Detection. (2015). 
// Mirizzi et al. arXiv:1508.00785  http://doi.org/10.1393/ncr/i2016-10120-8
// and
// Supernova Signatures of Neutrino Mass Ordering. (2017, July 20). 
// Scholberg, K. arXiv:1707.06384 https://doi.org/10.1088/1361-6471/aa97be

// (* don't just add in the scripts, check them against the values in the paper!! *)
// (* also go through the includes and remove the ones which are unnecessary! *)

#ifndef SNMIX_HEADER_FILE
#define SNMIX_HEADER_FILE

#include <iomanip>
#include <fstream>

// Declare shared functions for neutrino mixing and printing

void write(double a, double B[], std::ofstream& outfile);
void write_nh(double a, double B[], double th12, std::ofstream& outfile);
void write_ih(double a, double B[], double th12, std::ofstream& outfile);


#endif
