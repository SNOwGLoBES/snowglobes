#ifndef _xscns_
#define _xscns_

double diffxscnvec(double, double, double);
double diffxscnaxial(double, double, double);
double diffxscninterf(double, double, double);
double diffxscnmag(double, double);

double diffangdist(double, double);

void sm_vector_couplings(int, double*);
void sm_axial_couplings(int, int, double*);

double GV_SM(int,int, int);
double GA_SM(int,int, int, int, int, int);


double GV_nsi_nonuniv(int, int, int, double, double, double, double, double, double, double, double, double, double);

double GV_nsi_fc2(int, int, int, double, double, double, double, double, double, double, double, double, double);

double chgradcorr(int,int);
double mufactor(double);
double taufactor(double);

#endif
