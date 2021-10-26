#include "FormFactor.h"

#include <cstdlib>
#include <ctime>
#include <cstring>
#include <sstream>
#include <iostream>
#include <math.h>


///////

void Helm::Setsval(double s) {

  sval = s;

}

double Helm::Getsval() { return sval;}


double Helm::FFval(double Q) 
{

  double ff = 1;
  //double R = 1.14*pow(A,1./3.);
  double R = 1.2*pow(A,1./3.);

  // This is for varying Rn but keeping sval fixed (will distort the shape
  // but in practice very small difference )

  // double Rnorig = sqrt(3./5.*pow(R,2)+3*sval*sval);
  //double Rmod = sqrt(5./3.*(pow(Rnorig*Rfac,2)-3*sval*sval));  
  //double qR = Q*Rmod;

  // This scales the radius
  Q *= Rfac;

   double qR = Q*R;
  
  ff= (3*(sin(qR)/(qR*qR)-cos(qR)/qR)/(qR))*exp(-1.*Q*Q*sval*sval/2.);
  //  ff2= pow(3*(sin(qR)/(qR*qR)-cos(qR)/qR)/(qR),2)*exp(-1.*Q*Q*sval*sval);
  if (isnan(ff)) {ff=1.;}
    
  return ff;

}

////


void Klein::Setakval(double ak) {

  akval = ak;

}

double Klein::Getskinfac() { return skinfac;}

void Klein::Setskinfac(double sf) {

  skinfac = sf;

}

double Klein::Getakval() { return akval;}



double Klein::FFval(double Q) 
{
 
  double ff = 1;

  // Gutlein... probably wrong
  // double R2 = 1.14*pow(A,1./3.);

  // Incorrect way to add skin
  //double R2 = 1.2*pow(A,1./3.)+skinfac*1.01*(double(A)-2.*Z)/double(A);

  // Adding a skin
  double skindelta = skinfac*1.01*(double(A)-2.*Z)/double(A);
  double R2 = 1.2*pow(A,1./3.);
  if (skindelta != 0) {
    
    R2 = sqrt(R2*R2+ 2*sqrt(15)/3*sqrt(R2*R2+10*akval*akval)*skindelta + 5*skindelta*skindelta/3);
  }

  //double Ravg = sqrt(3*R2*R2/5+6*akval*akval);
  //  std::cout << "A, Z, R2, skindelta, Ravg: "<<A<<" "<<Z<<" "<<R2<<" "<<skindelta<<" "<<Ravg<<std::endl;
  
      // This scales the radius by Rfac
  Q *= Rfac;

  double qR = Q*R2;


  ff= (3*(sin(qR)/(qR*qR)-cos(qR)/qR)/(qR))*(1./(1+akval*akval*Q*Q));
  //  ff2= pow(3*(sin(qR)/(qR*qR)-cos(qR)/qR)/(qR),2)*pow(1./(1+akval*akval*Q*Q),2);

    
  if (isnan(ff)) {ff=1.;}

  return ff;

}

/////////

void Horowitz::ReadFFfile()
{

  double q,ff;
  std::ifstream fffile;
  std::string fffilename = filename;
  fffile.open(fffilename.c_str());
  if (!fffile) {
    std::cout << "File "<<fffilename<<" does not exist!" <<std::endl;
    exit(-1);
  } else {
    while(! fffile.eof() ) 
      {
        fffile >> q >> ff;
	//	std::cout << q << " "<<ff << std::endl;
        _ffmap[q] = ff;
      }
  }

  fffile.close();
}

double Horowitz::FFval(double Q) 
{
  double ff = 1;

    //http://www.bnikolic.co.uk/blog/cpp-map-interp.html
    // Interpolate from the map.  Must have been initalize for output to make sense

  typedef std::map<double, double>::const_iterator i_t;

  // Q is scaled by Rfac (see email from Chuck, Oct 17, 2017)


  // Scale the radius
  Q *=Rfac;

  i_t i=_ffmap.upper_bound(Q);
  if(i==_ffmap.end())
    {
      return (--i)->second;
    }
  if (i==_ffmap.begin())
    {
      return i->second;
    }
  i_t l=i; --l;
  
  const double delta=(Q- l->first)/(i->first - l->first);
  ff= delta*i->second +(1-delta)*l->second;

  if (isnan(ff)) {ff=1.;}

  // Note not squared in the file
  return ff;

}


void Horowitz::SetFFfilename(const char * fname) {
  strcpy(filename, fname);
}

const char * Horowitz::GetFFfilename() {
  return filename;
}


////////

FormFactor::FormFactor()
{
  Rfac=1.;
}

FormFactor::FormFactor(const char * type)
{
  strcpy(fftype,type);
  Rfac=1.;

}


void FormFactor::SetA(int Aval) {

  A = Aval;

}

int FormFactor::GetA() { return A;}


void FormFactor::SetZ(int Zval) {

  Z = Zval;

}

int FormFactor::GetZ() { return Z;}



void FormFactor::SetRfac(double Rfacval) {

  Rfac = Rfacval;

}

double FormFactor::GetRfac() { return Rfac;}

void FormFactor::Setfftype(const char * type) {
  strcpy(fftype, type);
}

const char * FormFactor::Getfftype() {
  return fftype;
}
