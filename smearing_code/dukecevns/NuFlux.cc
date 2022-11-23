#include <cstdlib>
#include <ctime>
#include <cstring>
#include <sstream>
#include <iostream>
#include <math.h>
#include "NuFlux.h"


double PiDAR::fluxval(double Enu, int flavor, double ebinsize) 
{
    

  // Energies in MeV
  // 1 = e, 2 = mu
  //const double mmu = 105.6;
  const double mmu = 105.66837;
  //  const double Enumu = 29.9;
  const double Enumu = 29.792;

  const double a= 2/mmu;
  double flux = 0.;

  if (Enu>mmu/2.) {
    flux = 0.;
    return flux;
  }

  if (flavor == 1) {
    flux = norm*12*pow(a*Enu,2)*(1-a*Enu)*a*ebinsize;
   } else if (flavor == 2) {
    
    if (fabs(Enu-Enumu)<ebinsize/2.) {
      flux = norm;
    }

  } else if (flavor == -2) {
     flux = norm*2*pow(a*Enu,2)*(3-2*a*Enu)*a*ebinsize;
  } else {

    flux = 0.;
  }

  // Oscillate if requested

  //  std::cout << "1: "<< flux << std::endl;
  double sin22th;
  if (doosc==1) {
    
    if (abs(flavor) == 1) {
      sin22th = sin22thes;
    } else if (abs(flavor) == 2) {
      sin22th = sin22thmus;
    }  else if (abs(flavor) == 3) {
      sin22th = sin22thtaus;
    } else {
      sin22th = 0.;
    }

    // Simple sterile disappearance
    flux *= (1-sin22th*pow(sin(1.27*dm2*(baseline/100.)/Enu),2));


  }

  //  std::cout << 1.-sin22th<<" "<<dm2<<" "<<Enu<<" "<<" "<<baseline<<" "<<pow(sin(1.27*dm2*baseline/Enu),2)<<std::endl;
  //std::cout << flux << std::endl;
  
  return flux;

}

double PiDAR::maxEnu() 
{

  // Return the maximum energy in MeV

  // To compare with old
  //  double maxEnu = 105.6/2.;
  double maxEnu = 105.66837/2.;
  return maxEnu;

}

////

double Reactor::fluxval(double Enu, int flavor, double ebinsize) 
{
 
 // Polynomials from Mueller 2011.  Gives flux in per MeV per fission

  // parent array gives the relative contributions from each fissioning parent, which will actually vary with time

  double alpha235U[6] = {3.217, -3.111, 1.395, -3.690e-1, 4.445e-2, -2.053e-3};
  double alpha238U[6] = {4.833e-1, 1.927e-1, -1.283e-1, -6.762e-3, 2.233e-3, -1.536e-4};

  double alpha239Pu[6] = {6.413,-7.432,3.535,-0.8820,0.1025,-.004550};

  double alpha241Pu[6] = {3.251,-3.204,1.428,-.3675,0.04254,-0.001896};


  int numterms = 6;
  int p;

  // nuebar only
  if (flavor != -1) { 
    return 0;
  }

  double flux235 = 0.;
  double flux238 = 0.;
  double flux239 = 0.;
  double flux241 = 0.;
  for (p=1; p<=numterms; p++) {
    flux235 += alpha235U[p-1]*pow(Enu,p-1);
    flux238 += alpha238U[p-1]*pow(Enu,p-1);
    flux239 += alpha239Pu[p-1]*pow(Enu,p-1);
    flux241 += alpha241Pu[p-1]*pow(Enu,p-1);
  }

  flux235 = exp(flux235)*ebinsize;
  flux238 = exp(flux238)*ebinsize;
  flux239 = exp(flux239)*ebinsize;
  flux241 = exp(flux241)*ebinsize;

  double fluxtot = flux235*parentfrac[0]+flux238*parentfrac[1]+flux239*parentfrac[2]+flux241*parentfrac[3];

  return fluxtot*norm;


}


double Reactor::maxEnu() 
{

  // Return the maximum energy in MeV

  double maxEnu = 8.;

  return maxEnu;

}

void Reactor::SetParentFrac(double* frac) {

  int i;
  for (i=0;i<numparent;i++) {
    parentfrac[i] = frac[i];
  }

}

double* Reactor::GetParentFrac() {
  return parentfrac;
}


/////////

void NumericalFlux::ReadFluxFile()
{

  // Input file has energies in GeV
  double enu,nue,numu,nutau,nuebar,numubar,nutaubar;
  std::ifstream fluxfile;
  std::string fluxfilename = filename;
  fluxfile.open(fluxfilename.c_str());
  if (!fluxfile) {
    std::cout << "File "<<fluxfilename<<" does not exist!" <<std::endl;
    exit(-1);
  } else {
    while(! fluxfile.eof() ) 
      {
        fluxfile >> enu >> nue>>numu>>nutau>>nuebar>>numubar>>nutaubar;
        if (! fluxfile.eof()) {
	  _nuefluxmap[enu*1000.] = nue;
	  _numufluxmap[enu*1000.] = numu;
	  _nutaufluxmap[enu*1000.] = nutau;
	  _nuebarfluxmap[enu*1000.] = nuebar;
	  _numubarfluxmap[enu*1000.] = numubar;
	  _nutaubarfluxmap[enu*1000.] = nutaubar;

	}
      }
  }

  fluxfile.close();
}


double NumericalFlux::fluxval(double enu,int flavor, double ebinsize) 
{
  double flux = 1;

  // Output is flux per ebinsize, for ebinsize in MeV.  enu is in MeV.

    //http://www.bnikolic.co.uk/blog/cpp-map-interp.html
    // Interpolate from the map.  Must have been initalized for output to make sense

  typedef std::map<double, double>::const_iterator i_t;

  // Q is scaled by Rfac (see email from Chuck, Oct 17, 2017)

  std::map<double, double> _fluxmap;
  
  switch (flavor) {
  case 1: _fluxmap = _nuefluxmap;
    break;
  case 2: _fluxmap = _numufluxmap;
    break;
  case 3: _fluxmap = _nutaufluxmap;
    break;
  case -1: _fluxmap = _nuebarfluxmap;
    break;
  case -2: _fluxmap = _numubarfluxmap;
    break;
  case -3: _fluxmap = _nutaubarfluxmap;
    break;
    std::cout<< "Wrong flavor "<<std::endl;
  }



  i_t i=_fluxmap.upper_bound(enu);
  if(i==_fluxmap.end())
    {
      return (--i)->second;
    }
  if (i==_fluxmap.begin())
    {
      return i->second;
    }
  i_t l=i; --l;
  
  const double delta=(enu- l->first)/(i->first - l->first);
  flux= delta*i->second +(1-delta)*l->second;

  if (isnan(flux)) {flux=0.;}

  return flux*ebinsize*norm;

}

void NumericalFlux::SetFluxFilename(const char * fname) {
  strcpy(filename, fname);
}

const char * NumericalFlux::GetFluxFilename() {
  return filename;
}

double NumericalFlux::maxEnu() 
{

  // Return the maximum energy in MeV.  Map is such that
  //  this should be the same for all flavors

  typedef std::map<double, double>::const_reverse_iterator i_t;
  i_t it = _nuefluxmap.rbegin();
  double maxEnu = it->first;

  return maxEnu;

}


////

double PinchedThermal::fluxval(double Enu, int flavor, double ebinsize) 
{
 
  //  flavor index goes nue, numu, nutau, nuebar, numubar, nutaubar
  // -3, -2, -1, 1, 2, 3
  // Energy should be in MeV

  int j;

  switch (flavor) {
  case 1: j=0; 
    break;
  case 2: j=1;
    break;
  case 3: j=2;
    break;
  case -1: j=3;
    break;
  case -2: j=4;
    break;
  case -3: j=5;
    break;
  default: std::cout<< "Incorrect flavor "<<std::endl;
    exit(-1);
  }


  // Conversions for flux at 10 kpc

  double fluxtot= 0.;
  if (alpha[j]>0 && avgen[j]>0  && luminosity[j]>0) {
    const double dist=3.08568025e22; // [dist]=cm, 10 kpc

    double N=pow((alpha[j]+1.),(alpha[j]+1.))/(avgen[j]*tgamma(alpha[j]+1.));
    double phi=N*pow((Enu/avgen[j]),alpha[j])*exp((-1.)*(alpha[j]+1.)*Enu/avgen[j]); 

    fluxtot = 1./(4*M_PI*dist*dist)*luminosity[j]/avgen[j]*phi*ebinsize;

  }
  return fluxtot*norm;

}


double PinchedThermal::maxEnu() 
{

  // Return the maximum energy in MeV

  double maxEnu = 50.;

  return maxEnu;

}

void PinchedThermal::SetLuminosity(double* lumi) {

  // Input luminosity is erg/s, store as MeV/s
  const double mevpererg = 624150.;

  int i;
  for (i=0;i<6;i++) {
    luminosity[i] = lumi[i]*mevpererg;
  }

}

double* PinchedThermal::GetLuminosity() {
  return luminosity;
}

void PinchedThermal::SetAvgEn(double* en) {

  // In MeV
  int i;
  for (i=0;i<6;i++) {
    avgen[i] = en[i];
  }

}

double* PinchedThermal::GetAvgEn() {
  return avgen;
}

void PinchedThermal::SetAlpha(double* a) {

  int i;
  for (i=0;i<6;i++) {
    alpha[i] = a[i];
  }

}

double* PinchedThermal::GetAlpha() {
  return alpha;
}


////////

NuFlux::NuFlux(){
  norm = 1.;
}

NuFlux::NuFlux(const char * type)
{
   strcpy(fluxtype,type);
   norm = 1.;
}


void NuFlux::Setfluxtype(const char * type) {
  strcpy(fluxtype, type);
}

const char * NuFlux::Getfluxtype() {
  return fluxtype;
}

void NuFlux::SetNorm(double normval) {
  norm = normval;
}

double NuFlux::GetNorm() {
  return norm;
}


void NuFlux::SetOscParam(double* ua4, double m, double b) {
  doosc = 1;
  // Unitarity constraint... user must ensure not violated
  double us4_2 = 1.-pow(ua4[0],2)-pow(ua4[1],2)-pow(ua4[2],2);
  sin22thes = 4.*pow(ua4[0],2)*us4_2;
  sin22thmus = 4.*pow(ua4[0],2)*us4_2;
  sin22thtaus = 4.*pow(ua4[0],2)*us4_2;
  dm2 = m;
  baseline = b;
}

void NuFlux::GetOscParam(double* oscparam) {

  oscparam[0]=sin22thes;
  oscparam[1]=sin22thmus;
  oscparam[2]=sin22thtaus;
  oscparam[3]= dm2;

}

////
