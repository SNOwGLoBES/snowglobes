#ifndef _NuFlux_
#define _NuFlux_

#include <map>
#include <vector>
#include <fstream>
#include <math.h>
#include <string>

class NuFlux
{

 protected:

  char fluxtype[80];
  double norm;

  // Oscillation parameters... simple 2-param sterile osc but could expand to standard 3flav
  int doosc=0;
  double sin22thes;
  double sin22thmus;
  double sin22thtaus;
  double dm2; // eV^2
  double baseline; // in meters
    
 public: 
  NuFlux();
  NuFlux(const char *);
  ~NuFlux(){};

  // Arguments: energy in MeV, flavor, ebinsize in MeV.
  // Output is differential flux per ebinsize 
  virtual double fluxval(double, int, double) = 0;
  // Should make a flavor specific one too
  virtual double maxEnu() = 0;

  void Setfluxtype(const char *);
  const char * Getfluxtype();

  void SetNorm(double);
  double GetNorm();

  void SetOscParam(double*, double, double);
  void GetOscParam(double*);


};

class PiDAR: public NuFlux {

 protected:

 public:
  // PiDAR(){}
  PiDAR() : NuFlux("pidar") {}
  double fluxval(double,int, double);
  double maxEnu();

};

class Reactor: public NuFlux {

 protected:

 public:
  // Reactor() {}

 Reactor() : NuFlux("reactor") {}
 double fluxval(double, int, double);
 double maxEnu();

 const int numparent = 4;
 double parentfrac[4];
 void SetParentFrac(double*);
 double* GetParentFrac();

};

class NumericalFlux: public NuFlux {

  // Read generic flux from a file

 protected:
  std::map<double,double> _nuefluxmap;
  std::map<double,double> _numufluxmap;
  std::map<double,double> _nutaufluxmap;
  std::map<double,double> _nuebarfluxmap;
  std::map<double,double> _numubarfluxmap;
  std::map<double,double> _nutaubarfluxmap;
  char filename[80];

 public:

 NumericalFlux() : NuFlux("numericalflux") {}
 double fluxval(double, int, double);
 void SetFluxFilename(const char * filename);
 const char * GetFluxFilename();
 void ReadFluxFile();
 double maxEnu();


};

class PinchedThermal: public NuFlux {

  // Pinched-thermal flux
  // Indices: 0 nue, 1: nuebar, 2: nuex

 protected:

 public:

 PinchedThermal() : NuFlux("pinchedthermal") {}
 double fluxval(double, int, double);
 double maxEnu();

 const int numparams = 3;  // Number of parameters describing the flux

 double luminosity[6];
 double avgen[6];
 double alpha[6];

 void SetLuminosity(double*);
 double* GetLuminosity();
 void SetAvgEn(double*);
 double* GetAvgEn();
 void SetAlpha(double*);
 double* GetAlpha();

};


#endif
