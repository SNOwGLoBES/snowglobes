#include <cstdlib>
#include <ctime>
#include <cstring>
#include <sstream>
#include <iostream>
#include <math.h>
#include "DetectorResponse.h"

// Numerical QF file related methods
// Energies in MeV

void DetectorResponse::ReadQFFile()
{

  double erec;
  double qf;
  std::ifstream qffile;
  std::string filename = qffilename;
  qffile.open(filename.c_str());
  if (!qffile) {
    std::cout << "File "<<filename<<" does not exist!" <<std::endl;
    exit(-1);
  } else {
    while(! qffile.eof() ) 
      {
        qffile >> erec >> qf;

        if (! qffile.eof()) {
          _qfmap[erec] = qf;
        }
      }
  }

  qffile.close();
}


////

double DetectorResponse::qfnum(double erec) 
{
  double qf = 1;


    //http://www.bnikolic.co.uk/blog/cpp-map-interp.html
    // Interpolate from the map.  Must have been initalized for output to make sense

  typedef std::map<double, double>::const_iterator i_t;

  //  std::map<double, double> _qfmap;
  
  i_t i=_qfmap.upper_bound(erec);

  if(i==_qfmap.end())
    {
      return (--i)->second;
    }
  if (i==_qfmap.begin())
    {
      return i->second;
    }
  i_t l=i; --l;
  
  const double delta=(erec- l->first)/(i->first - l->first);
  qf= delta*i->second +(1-delta)*l->second;

  if (isnan(qf)) {qf=0.;}

  return qf;

}

double DetectorResponse::qfnumderiv(double erec)
{

  double qfnumderiv = 0.;

  // Compute the derivative of the numerical map at value erec... it's derivative of erec*qf

  typedef std::map<double, double>::const_iterator i_t;

  i_t i=_qfmap.upper_bound(erec);

  double delta;
  double er1,er2;
  double qf1,qf2;
  double rise;
  double run;

  i_t l;
  if(i==_qfmap.end()) {
    // Actually same as normal case
      i_t np = i; --np;
      er1 = np->first;
      er2 = i->first;
      qf1 = np->second;
      qf2 = i->second;

      //      rise = qf2 * er2- qf1 * er1;
      //run = er2-er1;
  }
  else if (i==_qfmap.begin()){
      i_t nl = i; nl++;
      er1 = i->first;
      er2 = nl->first;
      qf1 = i->second;
      qf2 = nl->second;

      
  } else {
    l=i; --l;
    er1 = l->first;
    er2 = i->first;
    
    qf1 = l->second;
    qf2 = i->second;

    
  }

  //double qferec = qf1+(qf2-qf1)/(er2-er1)*(erec-er1);
  // This gives problems if erec=er1
  //  rise = qferec * erec - qf1 * er1;
  rise = qf2 * er2 - qf1 * er1;
  run = er2-er1;

  //  std::cout << er1<<" "<<er2<<" "<<erec<<" "<<qf1<<" "<<qf2<<" "<<qferec<<std::endl;
  delta= rise/run;

  qfnumderiv= delta;

  if (isnan(qfnumderiv)) {qfnumderiv=0.;}

  return qfnumderiv;

}

void DetectorResponse::SetQFFilename(const char * fname) {
  strcpy(qffilename, fname);
}

const char * DetectorResponse::GetQFFilename() {
  return qffilename;
}

double DetectorResponse::maxErec() 
{

  // Return the maximum energy in MeV.   For the numerical file 

  typedef std::map<double, double>::const_reverse_iterator i_t;
  i_t it = _qfmap.rbegin();
  double maxErec = it->first;

  return maxErec;

}

/////// 

// Polynomial QF-related methods

void DetectorResponse::ReadQFPolyFile() {

  double coeff;
  std::ifstream qfpolyfile;
  std::string filename = qfpolyfilename;
  qfpolyfile.open(filename.c_str());
  if (!qfpolyfile) {
    std::cout << "File "<<filename<<" does not exist!" <<std::endl;
    exit(-1);
  } else {
    qfpolyfile >> qfpolyrange[0]>>qfpolyrange[1];

    while(! qfpolyfile.eof() ) 
      {
	qfpolyfile>>coeff;
        if (! qfpolyfile.eof()) {
	  qfpolycoeff.push_back(coeff);
        }
      }
  }

  qfpolyfile.close();

}


double DetectorResponse::qfpoly(double erec) {

  // erec in MeV
  double qf=0;
  for (size_t i=0; i<qfpolycoeff.size();i++) {
    if (erec>=qfpolyrange[0]  && erec<=qfpolyrange[1] ) {
      qf += qfpolycoeff[i]*pow(erec,i);
    } else if (erec<qfpolyrange[0]) {
      qf += qfpolycoeff[i]*pow(qfpolyrange[0],i);
    } else {
      qf += qfpolycoeff[i]*pow(qfpolyrange[1],i);

    }
  }

  return qf;
} 

double DetectorResponse::qfpolyderiv(double erec) {

  // Return the value of the derivative of the polynomial times Erec (assume Eee = qf(Erec)*Erec
  // This useful for binning quenched distributions, dN/dEee = dN/dEr*dEr/dEee

  // erec in MeV
  double qfderiv = 0;
  for (size_t i=0; i<qfpolycoeff.size();i++) {
      
    if (erec>=qfpolyrange[0]  && erec<=qfpolyrange[1] ) {

      qfderiv += (i+1)*qfpolycoeff[i]*pow(erec,i);
    } else if (erec<qfpolyrange[0]) {
      qfderiv += (i+1)*qfpolycoeff[i]*pow(qfpolyrange[0],i);
    } else {
      qfderiv += (i+1)*qfpolycoeff[i]*pow(qfpolyrange[1],i);

    }

  }

  return qfderiv;
} 


void DetectorResponse::SetQFPolyFilename(const char * fname) {
  strcpy(qfpolyfilename, fname);
}

const char * DetectorResponse::GetQFPolyFilename() {
  return qfpolyfilename;
}


void DetectorResponse::SetQFPolyRange(double* range) {

  qfpolyrange[0] = range[0];
  qfpolyrange[1] = range[1];

}

double* DetectorResponse::GetQFPolyRange() { return qfpolyrange;}

// For Gaussian smearing

// Polynomial GS-related methods


void DetectorResponse::SetGSType(const char * type) {
  strcpy(gstype, type);
}

const char * DetectorResponse::GetGSType() {
  return gstype;
}


void DetectorResponse::ReadGSPolyFile() {

  double coeff;
  std::ifstream gspolyfile;
  std::string filename = gspolyfilename;
  gspolyfile.open(filename.c_str());
  if (!gspolyfile) {
    std::cout << "File "<<filename<<" does not exist!" <<std::endl;
    exit(-1);
  } else {
    char smeartype[80];
    gspolyfile >> smeartype;

    strcpy(gstype, smeartype);
    std::cout << "gstype "<<gstype << std::endl;
    gspolyfile >> gspolyrange[0]>>gspolyrange[1];

    while(! gspolyfile.eof() ) 
      {
	gspolyfile>>coeff;
        if (! gspolyfile.eof()) {
	  gspolycoeff.push_back(coeff);
        }
      }
  }

  gspolyfile.close();

}


double DetectorResponse::gspoly(double en) {

  // returns the sigma as a function of energy (usually Eee)

  // erec in MeV
  double gs=0;
  for (size_t i=0; i<gspolycoeff.size();i++) {

    if (en>=gspolyrange[0]  && en<=gspolyrange[1] ) {
      gs += gspolycoeff[i]*pow(en,i);
    } else if (en<gspolyrange[0]) {
      // Use enpoint values if out of range
      gs += gspolycoeff[i]*pow(gspolyrange[0],i);
    } else {
      gs += gspolycoeff[i]*pow(gspolyrange[1],i);
    }


  }
    
  return gs;
} 


double DetectorResponse::gspolysqrt(double en) {

  // returns the sigma as a function of energy (usually Eee)

  // erec in MeV
  double gs=0;
  double entoeval = en;

  if (en<gspolyrange[0]) {
    // Use endpoint values if out of range
    entoeval = gspolyrange[0];
    
  } else if (en>gspolyrange[1]) {
    entoeval = gspolyrange[0];
  }
    

  gs=  gspolycoeff[0]*sqrt(entoeval)+gspolycoeff[1]*entoeval;
   
  return gs;
} 



void DetectorResponse::SetGSPolyFilename(const char * fname) {
  strcpy(gspolyfilename, fname);
}

const char * DetectorResponse::GetGSPolyFilename() {
  return gspolyfilename;
}


void DetectorResponse::SetGSPolyRange(double* range) {

  gspolyrange[0] = range[0];
  gspolyrange[1] = range[1];

}

double* DetectorResponse::GetGSPolyRange() { return gspolyrange;}

void DetectorResponse::SetMaxEee(double maxeee) {
  maxEee = maxeee;
}

double DetectorResponse::GetMaxEee() {
  return maxEee;
}

void DetectorResponse::SetNEeeBin(int neeebin) {
  NEeeBin = neeebin;
}

int DetectorResponse::GetNEeeBin() {
  return NEeeBin;
}

void DetectorResponse::SetQCBinning(int qcbin) {
  qcbinning = qcbin;
}

int DetectorResponse::GetQCBinning() {
  return qcbinning;
}



void DetectorResponse::SetMaxSmearEn(double maxsmearen) {
  maxSmearEn = maxsmearen;
}

double DetectorResponse::GetMaxSmearEn() {
  return maxSmearEn;
}


void DetectorResponse::SetNSmearBin(int nsmearbin) {
  NSmearBin = nsmearbin;
}

int DetectorResponse::GetNSmearBin() {
  return NSmearBin;
}


void DetectorResponse::SetGaussSmearingMatrix() {

	// Do a binned normalization for each column, or else it's subject to binning effects and edge effects, and won't be unitary


   SmearingMatrix = new double*[NSmearBin];
   for(int i = 0; i < NSmearBin; ++i) {
    SmearingMatrix[i] = new double[NSmearBin];
   }

  double eni=0.;
  double enj=0.;
  int ien, jen;

  double enstep;
  if (qcbinning>0) {
    enstep = qcbinning;
  } else {
    enstep  = maxSmearEn/NSmearBin;
  }
  std::cout << "Setting smearing of gstype: "<<gstype<<" "<<NSmearBin<<" "<<maxSmearEn<<std::endl;
     
  for (ien=0;ien<NSmearBin;ien++) {
    eni += enstep;

    enj = 0.;
    for (jen=0;jen<NSmearBin;jen++) {
      enj += enstep;
      double sigma;
      if (strncmp(gstype,"polyfrac",8)==0) {
	sigma = this->gspoly(enj)*enj;
      } else if (strncmp(gstype,"polysqrt",8)==0) {
	sigma = this->gspolysqrt(enj);
	//	std::cout << gstype<<" "<<ien<<" "<<jen<<" "<<eni<<" "<<enj<<" "<<sigma <<std::endl;
      } else if (strncmp(gstype,"sqrtofpoly",10)==0) {
	sigma = sqrt(this->gspoly(enj));
      } else if (strncmp(gstype,"poly",4)==0) {
	sigma = this->gspoly(enj);
      }
      else {
	std::cout << "Unknown smearing type; please set gstype "<<std::endl;
	exit(1);
      }

      //      std::cout << "sigma "<<eni<<" "<<enj<<" "<<sigma <<std::endl;

      SmearingMatrix[ien][jen] = exp(-pow(eni-enj,2)/(2*pow(sigma,2)))/(sqrt(2*M_PI)*sigma);

      
    }

    if (isnan(SmearingMatrix[ien][jen]) ) {SmearingMatrix[ien][jen] = 0.;}
	    
  } // End of loop over rows


  // Now normalize over rows in a given column

  for (jen=0;jen<NSmearBin;jen++) {
    
    double totincolumn = 0.;
    for (ien=0;ien<NSmearBin;ien++) {
      totincolumn += SmearingMatrix[ien][jen];
    }

    for (ien=0;ien<NSmearBin;ien++) {
      if (totincolumn>0) {
	SmearingMatrix[ien][jen] /= totincolumn;
      } else {
	SmearingMatrix[ien][jen] = 0.;
      }
    }
  } // End of normalizing

}



void DetectorResponse::SetPoissonSmearingMatrix() {

	// Do a binned normalization for each column, or else it's subject to binning effects and edge effects, and won't be unitary


   SmearingMatrix = new double*[NSmearBin];
   for(int i = 0; i < NSmearBin; ++i) {
    SmearingMatrix[i] = new double[NSmearBin];
   }

  double pei=0.;
  double pej=0.;
  int ipe, jpe;

  double pestep = maxSmearEn/NSmearBin; // to synch with Smear method 
     
  for (ipe=0;ipe<NSmearBin;ipe++) {
    
    pej = 0.;
    for (jpe=0;jpe<NSmearBin;jpe++) {

      // Get the Poisson prob for pei given mean pej

      double poissprob;
      if (pej==0) {
         poissprob=0;
      } else if (pej>500) { // gaussian approximation of poisson
         poissprob = exp(-pow((pei-pej)/pej, 2.0)/2)/(pej*sqrt(2*M_PI));
      } else {
         double lnfact = 0.;
         for (int j=0;j<=int(pei)-1;j++) lnfact += log(float(int(pei)-j));

         double temp1 = int(pei)*log(pej)-lnfact;
         double temp2 = exp(temp1);
         poissprob = temp2*exp(-1.*pej);
      }
      if (isnan(poissprob)) std::cout << pej<<" "<<pei<<" "<<poissprob<<std::endl;

      SmearingMatrix[ipe][jpe] = poissprob;
      pej += pestep;
    }

    if (isnan(SmearingMatrix[ipe][jpe]) ) {SmearingMatrix[ipe][jpe] = 0.;}

    pei += pestep;
	    
  } // End of loop over rows


  // Now normalize over rows in a given column

  for (jpe=0;jpe<NSmearBin;jpe++) {
    
    double totincolumn = 0.;
    for (ipe=0;ipe<NSmearBin;ipe++) {
      totincolumn += SmearingMatrix[ipe][jpe];
    }

    for (ipe=0;ipe<NSmearBin;ipe++) {
      if (totincolumn>0) {
	SmearingMatrix[ipe][jpe] /= totincolumn;
      } else {
	SmearingMatrix[ipe][jpe] = 0.;
      }
    }
  } // End of normalizing

}



std::map<double, double> DetectorResponse::Smear(std::map<double,double> _unsmeared) {

  std::map<double,double> _smeared = _unsmeared;
 
  // For norm check
  double totsmeared = 0;
  double totunsmeared = 0;

	// Do the smearing

  
  double enstep;
  if (qcbinning>0) {
    enstep = qcbinning;
  } else {
    enstep  = maxSmearEn/NSmearBin;
  }

  int ien, jen;
  double eni=0.;

  // Check unitarity-- don't expect for the Poisson smearing due to
  //   zero column
  //for (jen=0;jen<NSmearBin;jen++) {
  //  double totincol = 0.;
  // for (ien=0;ien<NSmearBin;ien++) {
  //    totincol += SmearingMatrix[ien][jen];
  //  }
    //    std::cout << "col "<<jen<<" "<<totincol<<std::endl;
  // }

  

  for (ien=0;ien<NSmearBin;ien++) {
    
    _smeared[eni] = 0.;
    double enj = 0;

    for (jen=0;jen<NSmearBin;jen++) {
	  
      _smeared[eni]  +=   _unsmeared[enj]*SmearingMatrix[ien][jen];

      //      if (SmearingMatrix[ien][jen]>0) { 
      //std::cout << maxSmearEn<<" "<<NSmearBin<<" "<<enstep<<" "<<ien<<" "<<jen<<" "<< eni<<" "<<enj<<" "<<_unsmeared[enj]<<" "<<SmearingMatrix[ien][jen]<<" "<<_smeared[eni]<<std::endl;
      //}
      enj += enstep;

    } // End of loop over columns for this row
	  
    totunsmeared += _unsmeared[eni];
    totsmeared += _smeared[eni];
    eni += enstep;


  } // End of loop over rows

  	std::cout << "Total smeared events: "<<totsmeared*enstep<<" unsmeared "<<totunsmeared*enstep<<std::endl;

  return _smeared;

}


// Generic methods

DetectorResponse::DetectorResponse(){}

DetectorResponse::DetectorResponse(const char * type)
{
   strcpy(detectortype,type);
}


void DetectorResponse::SetDetectorType(const char * type) {
  strcpy(detectortype, type);
}

const char * DetectorResponse::GetDetectorType() {
  return detectortype;
}

// Numerical efficiency file related methods

void DetectorResponse::ReadEfficFile()
{

  double erec;
  double effic;
  std::ifstream efficfile;
  std::string filename = efficfilename;
  efficfile.open(filename.c_str());
  if (!efficfile) {
    std::cout << "File "<<filename<<" does not exist!" <<std::endl;
    exit(-1);
  } else {
    while(! efficfile.eof() ) 
      {
        efficfile >> erec >> effic;

        if (! efficfile.eof()) {
          _efficmap[erec] = effic;
        }
      }
  }

  efficfile.close();
}


////

double DetectorResponse::efficnum(double en) 
{
  double effic = 1;


    //http://www.bnikolic.co.uk/blog/cpp-map-interp.html
    // Interpolate from the map.  Must have been initalized for output to make sense

  typedef std::map<double, double>::const_iterator i_t;

  //  std::map<double, double> _efficmap;
  
  i_t i=_efficmap.upper_bound(en);

  if(i==_efficmap.end())
    {
      return (--i)->second;
    }
  if (i==_efficmap.begin())
    {
      return i->second;
    }
  i_t l=i; --l;
  
  const double delta=(en- l->first)/(i->first - l->first);
  effic= delta*i->second +(1-delta)*l->second;

  if (isnan(effic)) {effic=0.;}

  return effic;

}

void DetectorResponse::SetEfficType(const char * etype) {
  strcpy(efftype, etype);
}

const char * DetectorResponse::GetEfficType() {
  return efftype;
}



void DetectorResponse::SetEfficFilename(const char * fname) {
  strcpy(efficfilename, fname);
}

const char * DetectorResponse::GetEfficFilename() {
  return efficfilename;
}

double DetectorResponse::maxEfficErec() 
{

  // Return the maximum energy in MeV.   For the numerical file 

  typedef std::map<double, double>::const_reverse_iterator i_t;
  i_t it = _efficmap.rbegin();
  double maxErec = it->first;

  return maxErec;

}


void DetectorResponse::SetStepThresh(double thresh) {

  step_thresh = thresh;
 
}

double DetectorResponse::GetStepThresh() { return step_thresh;}

void DetectorResponse::SetUpperThresh(double thresh) {

  upper_thresh = thresh;
 
}

double DetectorResponse::GetUpperThresh() { return upper_thresh;}


/////// 
