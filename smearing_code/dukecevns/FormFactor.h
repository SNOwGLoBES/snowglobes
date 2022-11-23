#ifndef _FormFactor_
#define _FormFactor_

#include <map>
#include <vector>
#include <fstream>
#include <math.h>
#include <string>

using namespace std;

class FormFactor
{

 protected:

  int A; 
  int Z;
  double Rfac = 1;

  char fftype[80];
 

 public: 
  FormFactor();
  FormFactor(const char *);
  ~FormFactor(){};

  virtual double FFval(double) = 0;

  void SetA(int);
  int GetA();

  void SetZ(int);
  int GetZ();

  // Variation of Rn (as fraction of nominal)

  void SetRfac(double);
  double GetRfac();


  void Setfftype(const char *);
  const char * Getfftype();

};

class Helm: public FormFactor {

 protected:
  double sval;

 public:
  Helm() : FormFactor("helm") {}
  double FFval(double);
  void Setsval(double);
  double Getsval(); 

};

class Klein: public FormFactor {

 protected:
  double akval;
  double skinfac=0;

 public:
  Klein() : FormFactor("klein") {}
  double FFval(double);
  void Setakval(double);
  double Getakval(); 

  // Skin factor (zero for protons), used for some form factors

  void Setskinfac(double);
  double Getskinfac();

};


class Horowitz: public FormFactor {

 protected:
  std::map<double,double> _ffmap;
  char filename[80];

 public:
  Horowitz() : FormFactor("horowitz") {}
  double FFval(double);
  void SetFFfilename(const char * filename);
  const char * GetFFfilename();
  void ReadFFfile();


};

#endif
