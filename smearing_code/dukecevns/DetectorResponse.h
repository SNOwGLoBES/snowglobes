#ifndef _DetectorResponse_
#define _DetectorResponse_

#include <map>
#include <vector>
#include <fstream>
#include <math.h>
#include <string>

// Energies in MeV
class DetectorResponse
{

 protected:

  char detectortype[80];

  // For QF in numerical format
  std::map<double,double> _qfmap;
  char qffilename[80];

  // For QF in polynomial format
  char qfpolyfilename[80];
  std::vector<double> qfpolycoeff;
  double qfpolyrange[2]; // Range of validity for polynomial


  // For Gaussian smearing in polynomial format
  char gspolyfilename[80];
  std::vector<double> gspolycoeff;
  double gspolyrange[2]; // Range of validity for polynomial
  char gstype[0];

  int NEeeBin;
  double maxEee;
  double maxSmearEn;
  int NSmearBin;
  double** SmearingMatrix;

  int qcbinning=0; // Bin size for qc variable


  // Efficiency type:  recoil or electron-equivalent
  char efftype[80];

  // For a step-function threshold in MeVr

  double step_thresh=0.;

  // For a step-function upper threshold in MeVr

  double upper_thresh=0.;

  // For efficiency in numerical format
  std::map<double,double> _efficmap;
  char efficfilename[80];


 public: 
  DetectorResponse();
  DetectorResponse(const char *);
  ~DetectorResponse(){};

  void SetDetectorType(const char *);
  const char * GetDetectorType();

  // For QF in numerical format
  void SetQFFilename(const char * qffilename);
  const char * GetQFFilename();
  void ReadQFFile();
  double qfnum(double);
  double qfnumderiv(double);
  double maxErec();

  // For QF in polynomial format
  void SetQFPolyRange(double*);
  double* GetQFPolyRange();

  void SetQFPolyFilename(const char * qfpolyfilename);
  const char * GetQFPolyFilename();
  void ReadQFPolyFile();
  double qfpoly(double);
  double qfpolyderiv(double);


  // For Gaussian smearing in polynomial formats

  void SetGSType(const char *);
  const char * GetGSType();

  void SetGSPolyRange(double*);
  double* GetGSPolyRange();

  void SetGSPolyFilename(const char * gspolyfilename);
  const char * GetGSPolyFilename();
  void ReadGSPolyFile();
  double gspoly(double);
  double gspolysqrt(double);

  void SetNEeeBin(int);
  int GetNEeeBin();

  void SetQCBinning(int);
  int GetQCBinning();


  void SetMaxEee(double);
  double GetMaxEee();


  void SetMaxSmearEn(double);
  double GetMaxSmearEn();

  void SetNSmearBin(int);
  int GetNSmearBin();

  // Not bothering to clean this up with a delete method, I'm a bad person

  void SetGaussSmearingMatrix();
  void SetPoissonSmearingMatrix();

  std::map<double,double> Smear(std::map<double,double>);

  // For efficiency as a function of Erec, file in numerical format

  void SetEfficType(const char *);
  const char * GetEfficType();

  void SetEfficFilename(const char * efficfilename);
  const char * GetEfficFilename();
  void ReadEfficFile();
  double efficnum(double);
  double maxEfficErec();

  void SetStepThresh(double);
  double GetStepThresh();
  void SetUpperThresh(double);
  double GetUpperThresh();

};



#endif
