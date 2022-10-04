#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

#include "DetectorResponse.h"
#include "FormFactor.h"
#include "NuFlux.h"


#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <math.h>
#include <map>

#include "xscns.h"

void get_flavor_weight(int, double,double,double*, double*, double*);


int main(int argc, char * argv[] )
{

  if (argc<2) {
    std::cout << "Usage:  ./sns_rates [jsonfile]"<<std::endl;
    exit(0);
  }

  const char * jsonfile = argv[1];
  
  std::string jsonfilename = "jsonfiles/"+std::string(jsonfile)+".json";
  

// Read a JSON file with the parameters
    std::ifstream i(jsonfilename);
    json j;
    i >> j;

    // print values
    std::cout << j << '\n';

    std::cout << j["flux"]["nusperprotonperflavor"]<<std::endl;

    // Made these before adding json functionality
#include "isomaps.h"
  // Info for relevant mixtures
#include "mixtures.h"
	

  // Flavor weighting from file

  double wnumu=1.;
  double wnumubar=1.;
  double wnue=1.;

  int convolved=0;

  // Can't seem to do this with a nested key easily; make convolved not nested
  if (j.find("convolved") != j.end()) {
      convolved = j["convolved"];
  }
  std::cout << "Convolved "<<convolved<<std::endl;

  int perpot = 0;
  
 if (j.find("perpot") != j.end()) {
      perpot = j["perpot"];
  }

 // A generic normalization factor for the recoil spectrum
  double usernorm = 1.;

  if (j.find("usernorm") != j.end()) {
    usernorm = j["usernorm"];
    std::cout << "User normalization: "<<usernorm<<std::endl;
  }

 
  // Don't use flavor weights if using snsflux numerical flux; it that should take care of the weighting
  //    get_flavor_weight(1400.,7400.,&wnumu,&wnumubar,&wnue);

  double tw1 = j["timewindow"]["start"];
  double tw2 = j["timewindow"]["end"];
  get_flavor_weight(convolved,tw1,tw2,&wnumu,&wnumubar,&wnue);


  std::cout << "Flavor weights: "<< wnumu<<" "<<wnumubar<<" "<<wnue<<std::endl;

  // Zenodo for first microsecond
  // wnue=0.120603;
  // wnumubar = 0.120603;
  // wnumu =  0.974297;
   
  // Set up the form factor

  std::string ffname = j["formfactor"]["type"];

  int noff = 0.;
  if (ffname == "unity") {
    noff = 1.; 
  }
  // Array of pointers to form factors, for protons and neutrons separately, axial and vector separately
  // (although small differences)

  FormFactor** ffpv;
  ffpv = new FormFactor*[max_components];
  FormFactor** ffpa;
  ffpa = new FormFactor*[max_components];

  FormFactor** ffnv;
  ffnv = new FormFactor*[max_components];
  FormFactor** ffna;
  ffna = new FormFactor*[max_components];

  // Set up the flux

  //PiDAR* pidarflux = new PiDAR();
  //NumericalFlux* snsflux = new NumericalFlux();
  //snsflux->SetFluxFilename("time_integral.dat");
  //snsflux->ReadFluxFile();

    PiDAR* snsflux = new PiDAR();
  
  double kmax = snsflux->maxEnu();
  std::cout << "kmax "<<kmax << std::endl;
  // Normalize flux for 19.3 m from SNS, per cm^2 per s
  //snsflux->SetNorm(5.e14/(4*M_PI*1950.*1950.));

  double mevperproton = j["flux"]["mevperproton"];
  double jperproton= mevperproton*1.e6*1.6021e-19;
  double beampower = j["flux"]["power"];
  beampower*=1.e6; // in Joules/s
  double protonspersec = beampower/jperproton;
  double nusperprotonperflavor = j["flux"]["nusperprotonperflavor"];
  double nuspersecperflavor = nusperprotonperflavor*protonspersec;
  double dist = j["distance"];
  std::cout << "Nus per sec per flavor "<<nuspersecperflavor<<" "<<dist<<std::endl;
  snsflux->SetNorm(nuspersecperflavor/(4*M_PI*dist*dist));
  // Gives flux per pidk per energy bin per second, energy bin in MeV,  normalize for 5e14 decays/s 

  // Set up oscillations


  if (j.find("doosc") != j.end()) {
    int doosc = j["doosc"];
    if (doosc == 1) {
      double ua4[3];
      ua4[0] = j["osc"]["ue4"];
      ua4[1] = j["osc"]["umu4"];
      ua4[2] = j["osc"]["utau4"];
      double dm2 = j["osc"]["dm2"];
      snsflux->SetOscParam(ua4,dm2,dist);
    }
  }


  
  // Set up the detector response-- this is an overall detector response

  
  DetectorResponse* detresp = new DetectorResponse();

  double recoilthresh = 0.; //MeVr
  double eethresh = 0.; //MeVr
  double qcthresh = 0.; //Collected charge

  double recoilupperthresh = 0.; //MeVr
  double eeupperthresh = 0.; //MeVr
  double qcupperthresh = 0.; //Collected charge

  std::string eff_type = j["detectorresponse"]["efftype"];
  detresp->SetEfficType(eff_type.c_str());
  
  std::string effname = j["detectorresponse"]["efficiencyfile"];

  std::cout << "efficiency filename: "<<effname<<std::endl;
  if (effname != "none") {

    std::string eff_filename;
    eff_filename = "eff/"+std::string(effname);
    detresp->SetEfficFilename(eff_filename.c_str());
    detresp->ReadEfficFile();
  }
  if (eff_type == "erecoil") {
    recoilthresh = j["detectorresponse"]["stepthresh"];
    detresp->SetStepThresh(recoilthresh); // Not actually needed
    recoilupperthresh = j["detectorresponse"]["upperthresh"];
    detresp->SetUpperThresh(recoilupperthresh);
  } else if (eff_type == "eee") {
    eethresh = j["detectorresponse"]["stepthresh"];
    detresp->SetStepThresh(eethresh);
    eeupperthresh = j["detectorresponse"]["upperthresh"];
    detresp->SetUpperThresh(eeupperthresh);
  } else if (eff_type == "qc") {
    qcthresh = j["detectorresponse"]["stepthresh"];
    detresp->SetStepThresh(qcthresh);
    qcupperthresh = j["detectorresponse"]["upperthresh"];
    detresp->SetUpperThresh(qcupperthresh);
  }


    // Set up the material

    //    std::string material = "Ar40";

  std::string material = j["material"];
    

  std::ofstream outfile;
  std::string outfilename;
  //  outfilename = "sns_diff_rates-"+material+"-"+std::string(ffname)+".out";
  outfilename = "out/sns_diff_rates-"+std::string(jsonfile)+"-"+material+"-"+ffname+".out";
  outfile.open(outfilename);
  std::cout << outfilename <<std::endl;


   // Array for quenched total rates... could make this a stl vec 
   // but this is probably more efficient

  const int maxiq = 20000;
  double Er[maxiq];
  double Eee[max_components][maxiq];
  double dNdEee[max_components][maxiq];
  double dNdEr[max_components][maxiq];
  double dNdErall[maxiq]={0.};


  double M;
  double Delta;
  int Nn,Z,A;
  int Zdiff, Ndiff;

  std::string matname = material;

 // These are defined in mixtures include
  std::vector<double> fraction = molar_fraction[material];
  std::vector<std::string> isotope_component= isotopes[material];

  std::cout << "Material "<<matname<<std::endl;

  // Quenching filename info

  std::string qftype = j["detectorresponse"]["qftype"];
  std::string qfname = j["detectorresponse"]["qfname"];

  int is=0;
  std::vector<std::string>::iterator v = isotope_component.begin();
  std::string isotope;

  // First get the total mass.  Get also the maximum recoil values

  //double erecmaxvals[max_components];
  double Mtot = 0;
  v = isotope_component.begin();

  DetectorResponse** qffunc;
  qffunc = new DetectorResponse*[max_components];

  // First loop over materials, to set up material-specific arrays
  double minM = 1.e10;
  while( v != isotope_component.end()) {

    isotope = *v;
    std::cout << "isotope"<< isotope << std::endl;
    std::string isoname = std::string(isotope);

    Z = Zs[std::string(isotope)];
    Nn = Ns[std::string(isotope)];
    Delta = Deltas[std::string(isotope)];
    M = (Z+Nn)*amu - Z*me + Delta;
    if (M<minM) {minM=M;}
    Mtot += M*fraction[is];
    //erecmaxvals[is] = 2*kmax*kmax/(M+2*kmax);


    // Set up the form factor for this isotope

    //       std::cout << "Mass "<<M<<std::endl;


    double nvrfact = j["formfactor"]["nvrfact"];
    double narfact = j["formfactor"]["narfact"];
    double pvrfact = j["formfactor"]["pvrfact"];
    double parfact = j["formfactor"]["parfact"];


    if (ffname == "helm") {
          
      double nvsfact = j["formfactor"]["nvsfact"];
      double nasfact = j["formfactor"]["nasfact"];
      double pvsfact = j["formfactor"]["pvsfact"];
      double pasfact = j["formfactor"]["pasfact"];


      Helm* helmffnv= new Helm();
      ffnv[is] = helmffnv;
      helmffnv->Setsval(nvsfact);
      helmffnv->SetRfac(nvrfact);

      Helm* helmffna= new Helm();
      ffna[is] = helmffna;
      helmffna->Setsval(nasfact);
      helmffna->SetRfac(narfact);

      Helm* helmffpv= new Helm();
      ffpv[is] = helmffpv;
      helmffpv->Setsval(pvsfact);
      helmffpv->SetRfac(pvrfact);


      Helm* helmffpa= new Helm();
      ffpa[is] = helmffpa;
      helmffpa->Setsval(pasfact);
      helmffpa->SetRfac(parfact);

    }
    else if (ffname == "klein") {

      double nvak = j["formfactor"]["nvak"];
      double naak = j["formfactor"]["naak"];
      double pvak = j["formfactor"]["pvak"];
      double paak = j["formfactor"]["paak"];
      double nvskin = j["formfactor"]["nvskin"];
      double naskin = j["formfactor"]["naskin"];
      double pvskin = j["formfactor"]["pvskin"];
      double paskin = j["formfactor"]["paskin"];

      Klein* kleinffnv = new Klein();
      ffnv[is] = kleinffnv;
      kleinffnv->Setakval(nvak);
      kleinffnv->SetRfac(nvrfact);
      kleinffnv->Setskinfac(nvskin);

      Klein* kleinffna = new Klein();
      ffna[is] = kleinffna;
      kleinffna->Setakval(naak);
      kleinffna->SetRfac(narfact);
      kleinffna->Setskinfac(naskin);

      Klein* kleinffpv = new Klein();
      ffpv[is] = kleinffpv;
      kleinffpv->Setakval(pvak);
      kleinffpv->SetRfac(pvrfact);
      kleinffpv->Setskinfac(pvskin);


      Klein* kleinffpa = new Klein();
      ffpa[is] = kleinffpa;
      kleinffpa->Setakval(paak);
      kleinffpa->SetRfac(parfact);
      kleinffpv->Setskinfac(paskin);


    } 
    else if  (ffname =="horowitz"){
      Horowitz* horowitzffnv = new Horowitz();
      ffnv[is] = horowitzffnv;

      Horowitz* horowitzffna = new Horowitz();
      ffna[is] = horowitzffna;

      Horowitz* horowitzffpv = new Horowitz();
      ffpv[is] = horowitzffpv;

      Horowitz* horowitzffpa = new Horowitz();
      ffpa[is] = horowitzffpa;


      std::transform(isoname.begin(), isoname.end(),isoname.begin(), ::toupper);
      std::string horowitz_filename = isoname+".FF";
      std::cout << horowitz_filename << std::endl;
      horowitzffnv->SetFFfilename(horowitz_filename.c_str());
      horowitzffnv->ReadFFfile();
      horowitzffnv->SetRfac(nvrfact);

      horowitzffna->SetFFfilename(horowitz_filename.c_str());
      horowitzffna->ReadFFfile();
      horowitzffna->SetRfac(narfact);

      //      std::cout << horowitz_filename << std::endl;
      // Not really appropriate for protons, but using the structure
      horowitzffpv->SetFFfilename(horowitz_filename.c_str());
      horowitzffpv->ReadFFfile();
      horowitzffpv->SetRfac(pvrfact);


      horowitzffpa->SetFFfilename(horowitz_filename.c_str());
      horowitzffpa->ReadFFfile();
      horowitzffpa->SetRfac(parfact);


    }

    A = Nn + Z;
    if (noff == 0) {
      ffnv[is]->SetA(A);
      ffna[is]->SetA(A);
      ffpv[is]->SetA(A);
      ffpa[is]->SetA(A);
      
      ffnv[is]->SetZ(Z);
      ffna[is]->SetZ(Z);
      ffpv[is]->SetZ(Z);
      ffpa[is]->SetZ(Z);

    }
      
  // Set up detector quenching factors for each component

    
    if (qftype == "poly") {
      std::string qffilename;
      qffilename = "qf/"+std::string(qfname)+"_"+isoname+"_qf.txt";
      DetectorResponse* qf = new DetectorResponse();
	
      std::cout << "Quenching factor: "<<qffilename<<std::endl;
      qffunc[is] = qf;
      qffunc[is]->SetQFPolyFilename(qffilename.c_str());
      qffunc[is]->ReadQFPolyFile();
    } else if (qftype == "numerical") {
      std::string qffilename;
      qffilename = "qf/"+std::string(qfname)+"_"+isoname+"_qfnum.txt";
      DetectorResponse* qf = new DetectorResponse();
	
      std::cout << "Quenching factor: "<<qffilename<<std::endl;
      qffunc[is] = qf;
      qffunc[is]->SetQFFilename(qffilename.c_str());
      qffunc[is]->ReadQFFile();
    }
    v++; is++;
  }

    // Overall norm factor

  double detector_mass = j["mass"]; // tons

  double hours =j["flux"]["hours"];
  double exposure = 3600.*hours;

  double norm_factor = detector_mass*exposure*usernorm;

  if (perpot == 1) {
    double pot = protonspersec*exposure;
    norm_factor/= pot;
    std::cout << "Protons on target: "<<pot<<" ; output per pot "<<std::endl;
  }
  


  // Use the mass of the lightest component
  double erecmaxall = 2*kmax*kmax/(minM+2*kmax);
  
  //double recoilthresh = 0.013; //MeVr
  double erecstart = recoilthresh;
  double erecend = recoilupperthresh > recoilthresh ?
                   std::min(erecmaxall, recoilupperthresh) :
                   erecmaxall;
  //  double erecstep = 0.0001;
  double erecstep = 0.0001;

   // Now compute the differential recoil spectra

     
    double Erec;
    double knu;

    std::cout << "erecmaxall "<<erecmaxall<<std::endl;

    double knustep = 0.0001;
    if (j.find("knustep") != j.end()) {
      knustep = j["knustep"];
    }

   // The totals
   double toterecoil = 0.;
   double totevents = 0.;


   int iq=0;
   // Loop over recoil energy
   for (Erec=erecstart+erecstep;Erec<=erecend; Erec+=erecstep) {

     Er[iq] = Erec;

     // Contributions for each component
     double diffrate_e_vec[max_components]={0.};
     double diffrate_ebar_vec[max_components]={0.};
     double diffrate_mu_vec[max_components]={0.};
     double diffrate_mubar_vec[max_components]={0.};
     double diffrate_tau_vec[max_components]={0.};
     double diffrate_taubar_vec[max_components]={0.};
     
     double diffrate_e_axial[max_components]={0.};
     double diffrate_ebar_axial[max_components]={0.};
     double diffrate_mu_axial[max_components]={0.};
     double diffrate_mubar_axial[max_components]={0.};
     double diffrate_tau_axial[max_components]={0.};
     double diffrate_taubar_axial[max_components]={0.};

     double diffrate_e_interf[max_components]={0.};
     double diffrate_ebar_interf[max_components]={0.};
     double diffrate_mu_interf[max_components]={0.};
     double diffrate_mubar_interf[max_components]={0.};
     double diffrate_tau_interf[max_components]={0.};
     double diffrate_taubar_interf[max_components]={0.};

     double diffrate_e_mag[max_components]={0.};
     double diffrate_ebar_mag[max_components]={0.};
     double diffrate_mu_mag[max_components]={0.};
     double diffrate_mubar_mag[max_components]={0.};
     double diffrate_tau_mag[max_components]={0.};
     double diffrate_taubar_mag[max_components]={0.};



     // Sum for each component,  not quenched

     double sum_diffrate_e_vec=0;
     double sum_diffrate_ebar_vec=0;
     double sum_diffrate_mu_vec=0;
     double sum_diffrate_mubar_vec=0;
     double sum_diffrate_tau_vec=0;
     double sum_diffrate_taubar_vec=0;
     
     double sum_diffrate_e_axial=0;
     double sum_diffrate_ebar_axial=0;
     double sum_diffrate_mu_axial=0;
     double sum_diffrate_mubar_axial=0;
     double sum_diffrate_tau_axial=0;
     double sum_diffrate_taubar_axial=0;

     double sum_diffrate_e_interf=0;
     double sum_diffrate_ebar_interf=0;
     double sum_diffrate_mu_interf=0;
     double sum_diffrate_mubar_interf=0;
     double sum_diffrate_tau_interf=0;
     double sum_diffrate_taubar_interf=0;

     double sum_diffrate_e_mag=0;
     double sum_diffrate_ebar_mag=0;
     double sum_diffrate_mu_mag=0;
     double sum_diffrate_mubar_mag=0;
     double sum_diffrate_tau_mag=0;
     double sum_diffrate_taubar_mag=0;

     
     
     // With efficiency, which is a function of Erec in MeV in this formuation
     
     double recoil_eff_factor = 1;
     if (effname != "none" && eff_type == "erecoil") {
       recoil_eff_factor = detresp->efficnum(Erec);
     }
     
     //     double eff_factor = 1.;
     //     std::cout << "recoil eff factor "<<recoil_eff_factor<<std::endl;
     // Skip if too small contribution
     if (recoil_eff_factor>0) {
	  

       v = isotope_component.begin();
       // Now loop over components
       is=0;
       while( v != isotope_component.end()) {
	    
	 isotope = *v;
	 //	  std::cout << "isotope"<< isotope << std::endl;
	    
	 Z = Zs[std::string(isotope)];
	 Nn = Ns[std::string(isotope)];
	 Delta = Deltas[std::string(isotope)];
	 M = (Z+Nn)*amu - Z*me + Delta;
	    	    
	 Zdiff = Zdiffs[std::string(isotope)];
	 Ndiff = Ndiffs[std::string(isotope)];
	    
	 mass_fraction[is] = M/Mtot*fraction[is];
	    
	 A = Nn + Z;
	 //	  std::cout << " Z "<<Z<<" N "<<Nn<<" A "<<A<<" M "<<M << " "<<mass_fraction[is]<<std::endl;

	 // Loop over neutrino energy contributions
	  
	 // Minimum neutrino energy contributing to a given recoil energy

	 double knumin = 0.5*(Erec+sqrt(Erec*Erec+2*M*Erec));
	 double hbarc = 197.327; // MeV-fm, convert for Q in MeV for ff
	 double Q = sqrt(2*M*Erec+Erec*Erec); // MeV
	 //	 double Q = sqrt(2*M*Erec); // MeV

	 double qq = Q/hbarc;
	 //    double ff2 = helmff->FFval(qq);
	 double ffnvval;
	 double ffnaval;
	 double ffpvval;
	 double ffpaval;
	 if (noff == 0) {
	   ffnvval = ffnv[is]->FFval(qq);
	   ffnaval = ffna[is]->FFval(qq);
	   ffpvval = ffpv[is]->FFval(qq);
	   ffpaval = ffpa[is]->FFval(qq);
	 } else {
	   ffnvval = 1.;
	   ffnaval = 1.;
	   ffpvval = 1.;
	   ffpaval = 1.;
	 }
	   
	 //	 double ff2 = pow(ff[is]->FFval(qq),2);
	 //std::cout << "knumin, Erec, Q, ff2 "<<knumin<<" "<<Erec<<" "<<Q<<" "<<ff2<<" "<<M<<" "<<mass_fraction[is]<<std::endl;
	    
	 // SM Couplings

	 //	 double GV_sm = GV_SM(2015,Z,Nn);
	 //double GA_sm = GA_SM(2015,1,Z,Nn,Zdiff,Ndiff);
	 //double GA_sm_bar = GA_SM(2015,-1,Z,Nn,Zdiff,Ndiff);
	   
	 double gv[2], ga[2], gabar[2];
	 int pdgyr = j["couplings"]["pdgyear"];
	 sm_vector_couplings(pdgyr,gv);
	 sm_axial_couplings(pdgyr,1,ga);
	 sm_axial_couplings(pdgyr,-1,gabar);


	 // Bundle the form factor contributions with the SM couplings, separately for p and n
	 double GV_sm_wff = Z*gv[0]*ffpvval+Nn*gv[1]*ffnvval;
	 double GA_sm_wff = Zdiff*ga[0]*ffpaval+Ndiff*ga[1]*ffnaval;
	 double GA_sm_bar_wff = Zdiff*gabar[0]*ffpaval+Ndiff*gabar[1]*ffnaval;


	 // Charge radius correction
	 double mufact=1.;

	 if (j["couplings"]["chargeradiusfactor"] == "sehgal") {
	   mufact = mufactor(Q);
	 }	   

	double GV_sm_wff_e=GV_sm_wff;
	double GV_sm_wff_ebar=GV_sm_wff;
	double GV_sm_wff_mu=GV_sm_wff;
	double GV_sm_wff_mubar=GV_sm_wff;
	double GV_sm_wff_tau= GV_sm_wff;
	double GV_sm_wff_taubar= GV_sm_wff;

	double chgradcorr_e = 0.;
	double chgradcorr_ebar = 0.;
	double chgradcorr_mu = 0.;
	double chgradcorr_mubar = 0.;
	double chgradcorr_tau = 0.;
	double chgradcorr_taubar = 0.;

	if  (j["couplings"]["chargeradiusfactor"] == "erler") {

	  chgradcorr_e = chgradcorr(1,1);
	  chgradcorr_ebar = chgradcorr(-1,1);
	  chgradcorr_mu = chgradcorr(2,1);
	  chgradcorr_mubar = chgradcorr(-2,1);
	  chgradcorr_tau = chgradcorr(3,1);
	  chgradcorr_taubar = chgradcorr(-3,1);

	}

	
	if  (j["couplings"]["chargeradiusfactor"] == "giunti") {

	  chgradcorr_e = chgradcorr(1,2);
	  chgradcorr_ebar = chgradcorr(-1,2);
	  chgradcorr_mu = chgradcorr(2,2);
	  chgradcorr_mubar = chgradcorr(-2,2);
	  chgradcorr_tau = chgradcorr(3,2);
	  chgradcorr_taubar = chgradcorr(-3,2);

	}


	GV_sm_wff_e= Z*(gv[0]+chgradcorr_e)*ffpvval+Nn*gv[1]*ffnvval;
	GV_sm_wff_ebar= Z*(gv[0]+chgradcorr_ebar)*ffpvval+Nn*gv[1]*ffnvval;
	GV_sm_wff_mu= Z*(gv[0]+chgradcorr_mu)*ffpvval+Nn*gv[1]*ffnvval;
	GV_sm_wff_mubar= Z*(gv[0]+chgradcorr_mubar)*ffpvval+Nn*gv[1]*ffnvval;
	GV_sm_wff_tau= Z*(gv[0]+chgradcorr_tau)*ffpvval+Nn*gv[1]*ffnvval;
	GV_sm_wff_taubar= Z*(gv[0]+chgradcorr_taubar)*ffpvval+Nn*gv[1]*ffnvval;
	
	// Default is SM
	double GV_wff_e = GV_sm_wff_e;
	double GV_wff_ebar = GV_sm_wff_ebar;
	double GV_wff_mu = GV_sm_wff_mu;
	double GV_wff_mubar = GV_sm_wff_mubar;
	double GV_wff_tau = GV_sm_wff_tau;
	double GV_wff_taubar = GV_sm_wff_taubar;

	// Axial couplings are SM
	double GA_wff = GA_sm_wff;
	double GA_bar_wff = GA_sm_bar_wff;

	double eeeuV, eeedV;
	double eetauuV,  eetaudV;
	double eemuuV, eemudV;
	double emumuuV,  emumudV;
	double emutauuV,  emutaudV;
	double etautauuV,  etautaudV;

	// NSI vector couplings
	int donsi = 0;
	if (j.find("nsi") != j.end()) {
	  donsi = j["nsi"];
	}
	if (donsi != 0) {

	  eeeuV = j["couplings"]["eeeuV"];
	  eeedV = j["couplings"]["eeedV"];
	  eemuuV = j["couplings"]["eemuuV"];
	  eemudV = j["couplings"]["eemudV"];
	  eetauuV = j["couplings"]["eetauuV"];
	  eetaudV = j["couplings"]["eetaudV"];
	  emumuuV = j["couplings"]["emumuuV"];
	  emumudV = j["couplings"]["emumudV"];
	  emutauuV = j["couplings"]["emutauuV"];
	  emutaudV = j["couplings"]["emutaudV"];
	  etautauuV = j["couplings"]["etautauuV"];
	  etautaudV = j["couplings"]["etautaudV"];

	  // Not completely sure I am handling form factors properly for flavor-changing, but should be OK for equal form factors
	  double gV_nsi_wff_e= sqrt(pow(Z*(gv[0]+chgradcorr_e+2*eeeuV+ eeedV)*ffpvval
					+Nn*(gv[1]+eeeuV+2*eeedV)*ffnvval,2)
				    +pow(Z*(2*eemuuV+eemudV)*ffpvval
					 +Nn*(eemuuV+2*eemudV)*ffnvval,2)
				   +pow(Z*(2*eetauuV+eetaudV)*ffpvval
					+Nn*(eetauuV+2*eetaudV)*ffnvval,2));

	  double gV_nsi_wff_ebar= sqrt(pow(Z*(gv[0]+chgradcorr_ebar+2*eeeuV+ eeedV)*ffpvval
					+Nn*(gv[1]+eeeuV+2*eeedV)*ffnvval,2)
				    +pow(Z*(2*eemuuV+eemudV)*ffpvval
					 +Nn*(eemuuV+2*eemudV)*ffnvval,2)
				   +pow(Z*(2*eetauuV+eetaudV)*ffpvval
					+Nn*(eetauuV+2*eetaudV)*ffnvval,2));


	  double gV_nsi_wff_mu= sqrt(pow(Z*(gv[0]+chgradcorr_mu+2*emumuuV+ emumudV)*ffpvval
					+Nn*(gv[1]+emumuuV+2*emumudV)*ffnvval,2)
				    +pow(Z*(2*eemuuV+eemudV)*ffpvval
					 +Nn*(eemuuV+2*eemudV)*ffnvval,2)
				   +pow(Z*(2*emutauuV+emutaudV)*ffpvval
					+Nn*(emutauuV+2*emutaudV)*ffnvval,2));

	  double gV_nsi_wff_mubar= sqrt(pow(Z*(gv[0]+chgradcorr_mubar
					       +2*emumuuV+ emumudV)*ffpvval
					+Nn*(gv[1]+emumuuV+2*emumudV)*ffnvval,2)
				    +pow(Z*(2*eemuuV+eemudV)*ffpvval
					 +Nn*(eemuuV+2*eemudV)*ffnvval,2)
				   +pow(Z*(2*emutauuV+emutaudV)*ffpvval
					+Nn*(emutauuV+2*emutaudV)*ffnvval,2));

	  double gV_nsi_wff_tau= sqrt(pow(Z*(gv[0]+chgradcorr_tau+2*etautauuV+ etautaudV)*ffpvval
					+Nn*(gv[1]+etautauuV+2*etautaudV)*ffnvval,2)
				    +pow(Z*(2*eetauuV+eetaudV)*ffpvval
					 +Nn*(eetauuV+2*eetaudV)*ffnvval,2)
				   +pow(Z*(2*emutauuV+emutaudV)*ffpvval
					+Nn*(emutauuV+2*emutaudV)*ffnvval,2));

	  double gV_nsi_wff_taubar= sqrt(pow(Z*(gv[0]+chgradcorr_taubar
					       +2*etautauuV+ etautaudV)*ffpvval
					+Nn*(gv[1]+etautauuV+2*etautaudV)*ffnvval,2)
				    +pow(Z*(2*eetauuV+eetaudV)*ffpvval
					 +Nn*(eetauuV+2*eetaudV)*ffnvval,2)
				   +pow(Z*(2*emutauuV+emutaudV)*ffpvval
					+Nn*(emutauuV+2*emutaudV)*ffnvval,2));


	  GV_wff_e = gV_nsi_wff_e;
	  GV_wff_ebar = gV_nsi_wff_ebar;
	  GV_wff_mu = gV_nsi_wff_mu;
	  GV_wff_mubar = gV_nsi_wff_mubar;
	  GV_wff_tau = gV_nsi_wff_tau;
	  GV_wff_taubar = gV_nsi_wff_taubar;
	  
	} // End of donsi case
	
	
	// Magnetic moment..by flavor not quite right but OK for now

	double munu_e=0.;
	double munu_ebar=0.;
	double munu_mu=0.;
	double munu_mubar=0.;
	double munu_tau=0.;
	double munu_taubar=0.;

	if (j.find("magmom") != j.end()) {

	  munu_e = j["magmom"]["nue"];
	  munu_ebar = j["magmom"]["nuebar"];
	  munu_mu = j["magmom"]["numu"];
	  munu_mubar = j["magmom"]["numubar"];
	  munu_tau = j["magmom"]["nutau"];
	  munu_taubar = j["magmom"]["nutaubar"];
	  munu_e *= 1.e-10;
	  munu_ebar *= 1.e-10;
	  munu_mu *= 1.e-10;
	  munu_mubar *= 1.e-10;
	  munu_tau *= 1.e-10;
	  munu_taubar *= 1.e-10;
	  
	}

	 // Normalize for one ton of material
	 // Will be weighted by mass fraction
	    
	 double Nt = 1.e6/(M/amu)*6.0221409e23;
	    
	 // A2: G^2/(2Pi) * hbarcinmeters^-4 
	 double ntfac = Nt;
	    
	  // Quenching factor for this component and Eee for this Erec

	 double qfderiv=1;
	 if (qftype == "poly") {
	  Eee[is][iq] = qffunc[is]->qfpoly(Erec)*Erec;
	  qfderiv = abs(qffunc[is]->qfpolyderiv(Erec));
	 } else if (qftype == "numerical") {
	  Eee[is][iq] = qffunc[is]->qfnum(Erec)*Erec;
	  qfderiv = abs(qffunc[is]->qfnumderiv(Erec));
	 }


	 double drate_e_vec=0;
	 double drate_ebar_vec=0;
	 double drate_mu_vec=0;
	 double drate_mubar_vec=0;
	 double drate_tau_vec=0;
	 double drate_taubar_vec=0;
	    

	 double drate_e_axial=0;
	 double drate_ebar_axial=0;
	 double drate_mu_axial=0;
	 double drate_mubar_axial=0;
	 double drate_tau_axial=0;
	 double drate_taubar_axial=0;
	    

	 double drate_e_interf=0;
	 double drate_ebar_interf=0;
	 double drate_mu_interf=0;
	 double drate_mubar_interf=0;
	 double drate_tau_interf=0;
	 double drate_taubar_interf=0;


	 double drate_e_mag=0;
	 double drate_ebar_mag=0;
	 double drate_mu_mag=0;
	 double drate_mubar_mag=0;
	 double drate_tau_mag=0;
	 double drate_taubar_mag=0;


	 // Dumb integral, could be more clever to make it faster
	 for (knu=knumin;knu<=kmax;knu+=knustep) {
	  	  
	   drate_e_vec += diffxscnvec(knu,M,Erec)*snsflux->fluxval(knu,1,knustep);
	   drate_ebar_vec += diffxscnvec(knu,M,Erec)*snsflux->fluxval(knu,-1,knustep);
	   drate_mu_vec += diffxscnvec(knu,M,Erec)*snsflux->fluxval(knu,2,knustep);
	   drate_mubar_vec += diffxscnvec(knu,M,Erec)*snsflux->fluxval(knu,-2,knustep);
	   drate_tau_vec +=  diffxscnvec(knu,M,Erec)*snsflux->fluxval(knu,3,knustep);
	   drate_taubar_vec +=  diffxscnvec(knu,M,Erec)*snsflux->fluxval(knu,-3,knustep);


	   drate_e_axial += diffxscnaxial(knu,M,Erec)*snsflux->fluxval(knu,1,knustep);
	   drate_ebar_axial += diffxscnaxial(knu,M,Erec)*snsflux->fluxval(knu,-1,knustep);
	   drate_mu_axial += diffxscnaxial(knu,M,Erec)*snsflux->fluxval(knu,2,knustep);
	   drate_mubar_axial += diffxscnaxial(knu,M,Erec)*snsflux->fluxval(knu,-2,knustep);
	   drate_tau_axial +=  diffxscnaxial(knu,M,Erec)*snsflux->fluxval(knu,3,knustep);
	   drate_taubar_axial +=  diffxscnaxial(knu,M,Erec)*snsflux->fluxval(knu,-3,knustep);


	   drate_e_interf += diffxscninterf(knu,M,Erec)*snsflux->fluxval(knu,1,knustep);
	   drate_ebar_interf += diffxscninterf(knu,M,Erec)*snsflux->fluxval(knu,-1,knustep);
	   drate_mu_interf += diffxscninterf(knu,M,Erec)*snsflux->fluxval(knu,2,knustep);
	   drate_mubar_interf += diffxscninterf(knu,M,Erec)*snsflux->fluxval(knu,-2,knustep);
	   drate_tau_interf +=  diffxscninterf(knu,M,Erec)*snsflux->fluxval(knu,3,knustep);
	   drate_taubar_interf += diffxscninterf(knu,M,Erec)*snsflux->fluxval(knu,-3,knustep);

	   drate_e_mag += diffxscnmag(knu,Erec)*snsflux->fluxval(knu,1,knustep);
	   drate_ebar_mag += diffxscnmag(knu,Erec)*snsflux->fluxval(knu,-1,knustep);
	   drate_mu_mag += diffxscnmag(knu,Erec)*snsflux->fluxval(knu,2,knustep);
	   drate_mubar_mag += diffxscnmag(knu,Erec)*snsflux->fluxval(knu,-2,knustep);
	   drate_tau_mag +=  diffxscnmag(knu,Erec)*snsflux->fluxval(knu,3,knustep);
	   drate_taubar_mag += diffxscnmag(knu,Erec)*snsflux->fluxval(knu,-3,knustep);


	      
	   //    std::cout << Erec << " "<<knu<<" "<<snsflux->fluxval(knu,1,knustep)<<" "<<diffrate_e_vec<<std::endl;

	 } // End of loop over neutrino energy contributions
	   
	 // Now multiply by target-dependent factors and add up this recoil energy bin


	 diffrate_e_vec[is] += ntfac*pow(GV_wff_e,2)*mass_fraction[is]*drate_e_vec*wnue;
	 diffrate_ebar_vec[is] += ntfac*pow(GV_wff_ebar,2)*mass_fraction[is]*drate_ebar_vec;
	 diffrate_mu_vec[is] += ntfac*pow(GV_wff_mu,2)*mass_fraction[is]*drate_mu_vec*mufact*wnumu;
	 diffrate_mubar_vec[is] += ntfac*pow(GV_wff_mubar,2)*mass_fraction[is]*drate_mubar_vec*mufact*wnumubar;
	 diffrate_tau_vec[is] +=  ntfac*pow(GV_wff_tau,2)*mass_fraction[is]*drate_tau_vec;
	 diffrate_taubar_vec[is] += ntfac*pow(GV_wff_taubar,2)*mass_fraction[is]*drate_taubar_vec;


	 diffrate_e_axial[is] += ntfac*pow(GA_wff,2)*mass_fraction[is]*drate_e_axial*wnue;
	 diffrate_ebar_axial[is] += ntfac*pow(GA_bar_wff,2)*mass_fraction[is]*drate_ebar_axial;
	 diffrate_mu_axial[is] += ntfac*pow(GA_wff,2)*mass_fraction[is]*drate_mu_axial*mufact*wnumu;
	 diffrate_mubar_axial[is] += ntfac*pow(GA_bar_wff,2)*mass_fraction[is]*drate_mubar_axial*mufact*wnumubar;
	 diffrate_tau_axial[is] +=  ntfac*pow(GA_wff,2)*mass_fraction[is]*drate_tau_axial;
	 diffrate_taubar_axial[is] +=  ntfac*pow(GA_bar_wff,2)*mass_fraction[is]*drate_taubar_axial;


	 diffrate_e_interf[is] += ntfac*GV_wff_e*GA_wff*mass_fraction[is]*drate_e_interf*wnue;
	 diffrate_ebar_interf[is] += ntfac*GV_wff_ebar*GA_bar_wff*mass_fraction[is]*drate_ebar_interf;
	 diffrate_mu_interf[is] += ntfac*GV_wff_mu*GA_wff*mass_fraction[is]*drate_mu_interf*mufact*wnumu;
	 diffrate_mubar_interf[is] += ntfac*GV_wff_mubar*GA_bar_wff*mass_fraction[is]*drate_mubar_interf*mufact*wnumubar;
	 diffrate_tau_interf[is] +=  ntfac*GV_wff_tau*GA_wff*mass_fraction[is]*drate_tau_interf;
	 diffrate_taubar_interf[is] +=  ntfac*GV_wff_taubar*GA_bar_wff*mass_fraction[is]*drate_taubar_interf;

	 diffrate_e_mag[is] += ntfac*pow(munu_e,2)*pow(Z,2)*mass_fraction[is]*drate_e_mag*wnue;
	 diffrate_ebar_mag[is] += ntfac*pow(munu_ebar,2)*pow(Z,2)*mass_fraction[is]*drate_ebar_mag;
	 diffrate_mu_mag[is] += ntfac*pow(munu_mu,2)*pow(Z,2)*mass_fraction[is]*drate_mu_mag*wnumu;
	 diffrate_mubar_mag[is] += ntfac*pow(munu_mubar,2)*pow(Z,2)*mass_fraction[is]*drate_mubar_mag*wnumubar;
	 diffrate_tau_mag[is] +=  ntfac*pow(munu_tau,2)*pow(Z,2)*mass_fraction[is]*drate_tau_mag;
	 diffrate_taubar_mag[is] +=  ntfac*pow(munu_taubar,2)*pow(Z,2)*mass_fraction[is]*drate_taubar_mag;


	  // Now add the contribution from this isotope to the sum
	  

	 sum_diffrate_e_vec += diffrate_e_vec[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_ebar_vec += diffrate_ebar_vec[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mu_vec += diffrate_mu_vec[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mubar_vec += diffrate_mubar_vec[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_tau_vec += diffrate_tau_vec[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_taubar_vec += diffrate_taubar_vec[is]*norm_factor*recoil_eff_factor;
	 
	 sum_diffrate_e_axial += diffrate_e_axial[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_ebar_axial += diffrate_ebar_axial[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mu_axial += diffrate_mu_axial[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mubar_axial += diffrate_mubar_axial[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_tau_axial += diffrate_tau_axial[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_taubar_axial += diffrate_taubar_axial[is]*norm_factor*recoil_eff_factor;
	  
	 sum_diffrate_e_interf += diffrate_e_interf[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_ebar_interf += diffrate_ebar_interf[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mu_interf += diffrate_mu_interf[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mubar_interf += diffrate_mubar_interf[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_tau_interf += diffrate_tau_interf[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_taubar_interf += diffrate_taubar_interf[is]*norm_factor*recoil_eff_factor;


	 sum_diffrate_e_mag += diffrate_e_mag[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_ebar_mag += diffrate_ebar_mag[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mu_mag += diffrate_mu_mag[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mubar_mag += diffrate_mubar_mag[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_tau_mag += diffrate_tau_mag[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_taubar_mag += diffrate_taubar_mag[is]*norm_factor*recoil_eff_factor;



	  // Sum for this Erec and isotope
	  double sum_events_iso = 0;
	  sum_events_iso = diffrate_e_vec[is] + diffrate_ebar_vec[is] + diffrate_mu_vec[is]+ diffrate_mubar_vec[is]+ diffrate_tau_vec[is] + diffrate_taubar_vec[is];
	  sum_events_iso += diffrate_e_axial[is] + diffrate_ebar_axial[is] + diffrate_mu_axial[is]+ diffrate_mubar_axial[is]+ diffrate_tau_axial[is] + diffrate_taubar_axial[is];

	  sum_events_iso+= diffrate_e_interf[is] + diffrate_ebar_interf[is] + diffrate_mu_interf[is]+ diffrate_mubar_interf[is]+ diffrate_tau_interf[is] + diffrate_taubar_interf[is];


	  sum_events_iso+= diffrate_e_mag[is] + diffrate_ebar_mag[is] + diffrate_mu_mag[is]+ diffrate_mubar_mag[is]+ diffrate_tau_mag[is] + diffrate_taubar_mag[is];



	  sum_events_iso *= norm_factor*recoil_eff_factor;

	  // Now apply the quenching for this Ee and isotope component
	// sum_events_iso is dNderec

	  dNdEr[is][iq] = sum_events_iso;
	  dNdErall[iq] += sum_events_iso;
	    
	  if (qfderiv>0) {
	    dNdEee[is][iq] = sum_events_iso/qfderiv;
	  } else {
	    dNdEee[is][iq] = 0.;
	  }


	  //	 std::cout << is<<" "<<"Erec "<<Erec<<" mass frac "<<mass_fraction[is]<<" "<<" ntfac "<<ntfac<<" GV "<<GV_wff<<" GA "<<GA_wff<<std::endl;

	  //	 cout << "is, rate "<<is<<" "<<diffrate_e_vec[is]<<endl;

	 v++;is++;


       } // End of loop over material components

     } // End of efficiency factor check

     //     std::cout <<Erec<<scientific<<" "<<sum_diffrate_e_vec<<" "<<sum_diffrate_ebar_vec<<" "<<sum_diffrate_mu_vec<<" "<<sum_diffrate_mubar_vec<<" "<<sum_diffrate_tau_vec<<" "<<sum_diffrate_taubar_vec<<" "<<sum_diffrate_e_axial<<" "<<sum_diffrate_ebar_axial<<" "<<sum_diffrate_mu_axial<<" "<<sum_diffrate_mubar_axial<<" "<<sum_diffrate_tau_axial<<" "<<sum_diffrate_taubar_axial<<" "<<sum_diffrate_e_interf<<" "<<sum_diffrate_ebar_interf<<" "<<sum_diffrate_mu_interf<<" "<<sum_diffrate_mubar_interf<<" "<<sum_diffrate_tau_interf<<" "<<sum_diffrate_taubar_interf <<std::endl;

     // Only want diff values in scientific format
     std::cout.unsetf(ios::fixed | ios::scientific);


//      sum_diffrate_e_vec *= norm_factor*wnue;
//      sum_diffrate_ebar_vec *= norm_factor;
//      sum_diffrate_mu_vec *= norm_factor*wnumu;
//      sum_diffrate_mubar_vec *= norm_factor*wnumubar;
//      sum_diffrate_tau_vec *= norm_factor;
//      sum_diffrate_taubar_vec *= norm_factor;
     
//      sum_diffrate_e_axial *= norm_factor*wnue;
//      sum_diffrate_ebar_axial *= norm_factor;
//      sum_diffrate_mu_axial *= norm_factor*wnumu;
//      sum_diffrate_mubar_axial *= norm_factor*wnumubar;
//      sum_diffrate_tau_axial *= norm_factor;
//      sum_diffrate_taubar_axial *= norm_factor;
	
//      sum_diffrate_e_interf *= norm_factor*wnue;
//      sum_diffrate_ebar_interf *= norm_factor;
//      sum_diffrate_mu_interf *= norm_factor*wnumu;
//      sum_diffrate_mubar_interf *= norm_factor*wnumubar;
//      sum_diffrate_tau_interf *= norm_factor;
//      sum_diffrate_taubar_interf *= norm_factor;



     outfile << Erec<<scientific<<" "<<sum_diffrate_e_vec<<" "<<sum_diffrate_ebar_vec<<" "<<sum_diffrate_mu_vec<<" "<<sum_diffrate_mubar_vec<<" "<<sum_diffrate_tau_vec<<" "<<sum_diffrate_taubar_vec<<" "<<sum_diffrate_e_axial<<" "<<sum_diffrate_ebar_axial<<" "<<sum_diffrate_mu_axial<<" "<<sum_diffrate_mubar_axial<<" "<<sum_diffrate_tau_axial<<" "<<sum_diffrate_taubar_axial<<" "<<sum_diffrate_e_interf<<" "<<sum_diffrate_ebar_interf<<" "<<sum_diffrate_mu_interf<<" "<<sum_diffrate_mubar_interf<<" "<<sum_diffrate_tau_interf<<" "<<sum_diffrate_taubar_interf <<" "<<sum_diffrate_e_mag<<" "<<sum_diffrate_ebar_mag<<" "<<sum_diffrate_mu_mag<<" "<<sum_diffrate_mubar_mag<<" "<<sum_diffrate_tau_mag<<" "<<sum_diffrate_taubar_mag<<std::endl;
	// Reset the format
     std::cout.unsetf(ios::fixed | ios::scientific);
     
	double events=0;
	events = sum_diffrate_e_vec + sum_diffrate_ebar_vec + sum_diffrate_mu_vec+ sum_diffrate_mubar_vec+ sum_diffrate_tau_vec + sum_diffrate_taubar_vec;
	events += sum_diffrate_e_axial + sum_diffrate_ebar_axial + sum_diffrate_mu_axial+ sum_diffrate_mubar_axial+ sum_diffrate_tau_axial + sum_diffrate_taubar_axial;
        events += sum_diffrate_e_interf + sum_diffrate_ebar_interf + sum_diffrate_mu_interf+ sum_diffrate_mubar_interf+ sum_diffrate_tau_interf + sum_diffrate_taubar_interf;

        events += sum_diffrate_e_mag + sum_diffrate_ebar_mag + sum_diffrate_mu_mag+ sum_diffrate_mubar_mag+ sum_diffrate_tau_mag + sum_diffrate_taubar_mag;

	totevents+=events*erecstep;

	toterecoil += events*Erec*erecstep;

	// Increment bin for quenching
	
	iq++;

  } // End of loop over Erec

   if (recoilupperthresh > recoilthresh) {
      std::cout << "Total events over "<< recoilthresh*1000.<<" keVr and under "<<recoilupperthresh*1000<<" keVr: "<<totevents<< std::endl;
   }
   else {
      std::cout << "Total events over "<< recoilthresh*1000.<<" keVr: "<<totevents<< std::endl;
   }
   std::cout << "Total recoil energy deposited:  "<< toterecoil<< std::endl;

   outfile.close();


  std::ofstream integraloutfile;
  outfilename = "out/sns_diff_rates-"+std::string(jsonfile)+"-"+material+"-"+ffname+"-integral.out";

  std::cout << outfilename << std::endl;
  integraloutfile.open(outfilename);

  
  integraloutfile << j << '\n';
  if (recoilupperthresh > recoilthresh) {
    integraloutfile << "Total events over "<< recoilthresh*1000.<<" keVr and under "<<recoilupperthresh*1000<<" keVr: "<<totevents<< std::endl;
  }
  else {
    integraloutfile << "Total events over "<< recoilthresh*1000.<<" keVr: "<<totevents<< std::endl;
  }
  integraloutfile.close();

  // Output by isotope, integrated over flavor.
  //Can also output quenched stuff here
   // Loop over isotopes


  v = isotope_component.begin();
  // Now loop over components
  is=0;
  while( v != isotope_component.end()) {
    
    isotope = *v;
    //	  std::cout << "isotope"<< isotope << std::endl;
    std::string isoname = std::string(isotope);
    
    std::ofstream isooutfile;
    outfilename = "out/sns_diff_rates-"+std::string(jsonfile)+"-"+material+"-"+ffname+"-"+isoname+".out";

    std::cout << outfilename << std::endl;
    isooutfile.open(outfilename);

    int ie;
    for (ie=0;ie<iq;ie++) {

      isooutfile << Er[ie]<< "  "<<dNdEr[is][ie]<<endl;
    
    }
    isooutfile.close();
    v++;is++;
  }

  // Integrated over flavor and isotope

    std::ofstream allisooutfile;
    outfilename = "out/sns_diff_rates_alliso-"+std::string(jsonfile)+"-"+material+"-"+ffname+".out";

    std::cout << outfilename << std::endl;
    allisooutfile.open(outfilename);

    int ie;
    for (ie=0;ie<iq;ie++) {

      allisooutfile << Er[ie]<< "  "<<dNdErall[ie]<<endl;
    
    }
    allisooutfile.close();

    ////////////  QUENCHED  RESPONSE ///////////////

    // Now dump the quenched output, by isotope.  This has non-uniform Eee energy bins.  At the same time fill some maps to be interpolated for a sum, and get the maximum quenched energy for use for that/

    // Don't have this broken down by flavor and interaction... need to do that

    if (qfname != "none") {

      double maxeee = 0;
      // One of these per component
      std::map<double, double> _quenchedmap[max_components];

      // The total response
      std::map<double, double> _quenchedtot;

      v = isotope_component.begin();
      // Now loop over components
      is=0;

      while( v != isotope_component.end()) {
	
	isotope = *v;
	//	  std::cout << "isotope"<< isotope << std::endl;
	std::string isoname = std::string(isotope);
	
	std::ofstream qisooutfile;
	outfilename = "out/sns_diff_rates_quenched-"+std::string(jsonfile)+"-"+material+"-"+ffname+"-"+isoname+".out";
	
	std::cout << outfilename << std::endl;
	qisooutfile.open(outfilename);
	
	int ie;
	for (ie=0;ie<iq;ie++) {
	  
	  if (Eee[is][ie]>maxeee) {maxeee = Eee[is][ie];}
	  qisooutfile << Eee[is][ie]<< "  "<<dNdEee[is][ie]<<std::endl;
	  _quenchedmap[is][Eee[is][ie]] = dNdEee[is][ie];

	}
	qisooutfile.close();
	v++;is++;
      }

      // Retrieve the smearing function, which should have Gaussian sigma as a function of Eee
      std::string gsname = j["detectorresponse"]["gsname"];

      // Now interpolated rates for quenched, summed over components

      std::ofstream qoutfile;
      outfilename = "out/sns_diff_rates_quenched-alliso-"+std::string(jsonfile)+"-"+material+"-"+ffname+".out";
	
      std::cout << outfilename << std::endl;
      qoutfile.open(outfilename);

      // Fixed number of steps in quenched energy
      int  ieee = 0;

      //      int neeebin = j["detectorresponse"]["neeebin"];
      int neeebin = iq;

      double eeestep = maxeee/neeebin;
      double eee = 0;
      double dndeee=0.;

      double nquenchedtot=0;
      
      for (ieee=0;ieee<neeebin;ieee++) {
	
	eee += eeestep;
	_quenchedtot[eee] = 0.;

	// Now loop over components
	v = isotope_component.begin();

	is=0;
	while( v != isotope_component.end()) {
	  
	  isotope = *v;
	  //	  std::cout << "isotope"<< isotope << std::endl;
	  std::string isoname = std::string(isotope);
	  
	  // Interpolate dNdEee value for isotope is, at this eee
	  // Should encapsulate this in an interpolation routine

	  typedef std::map<double, double>::const_iterator i_t;

	  i_t i=_quenchedmap[is].upper_bound(eee);
	  if(i==_quenchedmap[is].end())
	    {
	      dndeee = (--i)->second;
	    } else if (i==_quenchedmap[is].begin()) 
	    {
	      dndeee =  i->second;
	    } else {
	    i_t l=i; --l;
	  
	    const double delta=(eee- l->first)/(i->first - l->first);
	    dndeee= delta*i->second +(1-delta)*l->second;
	  }	    
	  if (::isnan(dndeee)) {dndeee=0.;}
	  
	  _quenchedtot[eee] += dndeee;

	  v++;is++;

	} // End of loop over components


	  // Apply the Eee efficiency here, if requested and not smeared
	double eee_eff_factor = 0.;
	if (gsname == "none"  && eff_type == "eee"){
	  
	  if (eee>=eethresh &&
              (eeupperthresh > eethresh ? eee <= eeupperthresh : true)) {
	    if (effname != "none") {
	      eee_eff_factor = detresp->efficnum(eee);
	    } else {
	      eee_eff_factor = 1;
	    }
	  }

	} else {
	  eee_eff_factor = 1;
	}
	  
	
      // Now output the total quenched output, per MeVee
	qoutfile <<eee<<" "<<_quenchedtot[eee]*eee_eff_factor<<std::endl;
	std::cout.unsetf(ios::fixed | ios::scientific);

	nquenchedtot += _quenchedtot[eee]*eee_eff_factor*eeestep;



      } // End of loop over Eee

      qoutfile.close();


      std::cout << "Total quenched: "<<nquenchedtot<<std::endl;


      // Now pad the end of the quenched array to allow for smearing

      for (ieee=neeebin;ieee<neeebin*2;ieee++) {
	eee += eeestep;
	_quenchedtot[eee] = 0.;
      }


    // Now do Gaussian smearing, if requested.  Quenching must be requested also if this is to be invoked.  Apply also efficiency here


      if (gsname != "none") {
	

	// Read the smearing parameters from the file and set them
	std::string gsfilename;
	gsfilename = "gs/"+std::string(gsname)+"_gs.txt";

	DetectorResponse* gs = new DetectorResponse();


	gs->SetGSPolyFilename(gsfilename.c_str());
	gs->ReadGSPolyFile();

	gs->SetMaxSmearEn(maxeee*2);
	gs->SetNSmearBin(neeebin*2);
	gs->SetGaussSmearingMatrix();

	// Do the smearing
	std::map<double,double> _smearedmap= gs->Smear(_quenchedtot);

      // Output the smeared output file

	std::ofstream smoutfile;
	outfilename = "out/sns_diff_rates_smeared-alliso-"+std::string(jsonfile)+"-"+material+"-"+ffname+".out";
	
	std::cout << outfilename << std::endl;
	smoutfile.open(outfilename);

	double eeei=0.;
	double totweff = 0.;
	for (ieee=0;ieee<neeebin*2;ieee++) {

	  // Apply the Eee efficiency here, if requested
	  double eee_eff_factor = 1.;
	  if (eeei>=eethresh &&
              (eeupperthresh > eethresh ? eeei <= eeupperthresh : true)) {
	    if (effname != "none" && eff_type == "eee"){
	      eee_eff_factor = detresp->efficnum(eeei);	    
	    }
	    eeei += eeestep;
	    smoutfile << eeei<<" "<<_smearedmap[eeei]*eee_eff_factor<<std::endl;
	    totweff += _smearedmap[eeei]*eee_eff_factor*eeestep;
	  } else {
	    eeei += eeestep;
	  }
	}

	smoutfile.close();
	std::cout << "Total smeared with efficiency: "<<totweff<<std::endl;

      } // End of do-smearing case

    ////////////  QC-LEVEL RESPONSE ///////////////
    // For this, require quenching also
      // Collected charge, either pe or ADC, etc.


      // In principle this can be a nonlinear response function... linear for now
      double qcperkeVee = j["detectorresponse"]["qcperkeVee"];

      // The total qc distribution (all components... could break it out by isotope)
      std::map<double, double> _qcmapall;
      //Eee in MeVee

      double qcperMeVee = qcperkeVee*1000.;
      double maxmeanqc = maxeee*qcperMeVee;
      double maxqc = maxeee*qcperMeVee*2; // for smearing matrix

      //      std::cout << "maxqc "<<maxqc<<" maxmeanqc "<<maxmeanqc<<std::endl;
      // Loop over qc's, as integers

      // Do the Poisson qc smearing if requested

      std::string qcsmearing = j["detectorresponse"]["qcsmearing"];

      if (qcsmearing != "none") {
      
	int iqc;
	double qc;

	double totinqc = 0.;
	int qcbinning = j["detectorresponse"]["qcbinning"];
      
	for (iqc=0;iqc<=int(maxqc);iqc+=qcbinning) {

	  // Interpolate dNdqc from the quenchedmap

	  qc = double(iqc);
	  if (qc<=maxmeanqc+0.01) {
	  

	    typedef std::map<double, double>::const_iterator i_t;
	  
	    //	  double mevee = qc/qcperMeVee;
	    
	  // Do a more fine-grained interpolation and integrate over qc bin,
	  // to reduce binned integration error
	  // Not always really necessary

	    double fracqc;
	    double qcstep = 0.1*qcbinning;
	  
	    double dndqcinbin=0.;
	    double dndqcinterp;
	    double dndqc;
	  
	    // if <qc+0.5 includes last point
	    for (fracqc=qc-0.5*qcbinning;fracqc<qc+0.49*qcbinning;fracqc+=qcstep) {
	      if (fracqc>0) {
		
	      double mevee2 = fracqc/qcperMeVee;
	      
	      i_t i=_quenchedtot.upper_bound(mevee2);
	      
	      if(i==_quenchedtot.end())
		{
		  dndqc = (--i)->second;
		}
	      else if (i==_quenchedtot.begin())
		{
		  dndqc =  i->second;
		  
		} else {
		
		i_t l=i; --l;
		
		const double delta=(mevee2- l->first)/(i->first - l->first);
		dndqc = delta*i->second +(1-delta)*l->second;
	      }
	      dndqcinbin += dndqc*qcstep; 
	      
	      //	      std::cout <<qc <<" "<<fracqc<<" "<<mevee2<<" "<<dndqc*qcstep<<" "<<dndqcinbin<<std::endl;
	      
	      } // End of >0 qc case
	      
	  } // End of loop over fractional qc integration
	    //	std::cout <<qc <<" "<<mevee<<" "<<dndqcinbin<<std::endl;
	  
	  
	    dndqcinterp = dndqcinbin/qcperMeVee;
	    
	    if (::isnan(dndqcinterp)) {dndqcinterp=0.;}
	  
	  //	if (qc>0) {
	    totinqc += dndqcinterp;	
	    //	}
	    _qcmapall[qc] = dndqcinterp;

	    //	    std::cout << qc<<" "<<dndqcinterp<<" "<<totinqc<<std::endl;
	    
	  } else {
	    
	    // Pad the end with zeroes to allow for smearing
	    _qcmapall[qc] = 0.;
	    
	  } // end of qc<=maxmeanqc check
	  
	  
	}
	
	std::cout << "Integral of qc dist, including zero bin: "<<totinqc<<std::endl;
	
	
	// Poisson smear includes the zero bin
	DetectorResponse* qcsmear = new DetectorResponse();
	
	qcsmear->SetQCBinning(qcbinning);
	qcsmear->SetNSmearBin(int(maxqc/qcbinning)+1);
	qcsmear->SetMaxSmearEn(double(int(maxqc)+1));
	
	if (qcsmearing == "poisson") {
	  qcsmear->SetPoissonSmearingMatrix();
	  
	} else {
	  std::string qcgsname = j["detectorresponse"]["qcgsname"];

	  std::string qcsmearfilename;
	  qcsmearfilename = "gs/"+std::string(qcgsname)+"_qcsmear.txt";
	  qcsmear->SetGSPolyFilename(qcsmearfilename.c_str());
	  qcsmear->ReadGSPolyFile();
	  qcsmear->SetGaussSmearingMatrix();
	}


	std::cout << "do qc smearing" <<std::endl;
	// Do the smearing
	std::map<double,double> _smearedqcmap = qcsmear->Smear(_qcmapall);
	
      
	// Output the qc distribution, applying efficiency if requested
	// (should not also have recoil or quenched efficiency)
	
	std::ofstream qcoutfile;
	outfilename = "out/sns_diff_rates_qc-alliso-"+std::string(jsonfile)+"-"+material+"-"+ffname+".out";
	
	std::cout << outfilename << std::endl;
	qcoutfile.open(outfilename);
	
	double totev = 0.;
	double totevunsmeared = 0.;
	
	for (iqc=0;iqc<=int(maxqc);iqc+=qcbinning) {
	  
	// Apply the qc efficiency here, if requested
	
	  qc = double(iqc);
	  double qc_eff_factor = 1.;
	
	  if (iqc>=qcthresh &&
              (qcupperthresh > qcthresh ? iqc <= qcupperthresh : true)) {
	    if (effname != "none" && eff_type == "qc"){
	      qc_eff_factor = detresp->efficnum(qc);	    
	    }
	  
	    qcoutfile << iqc <<" "<<_smearedqcmap[qc]<<" "<<_smearedqcmap[qc]*qc_eff_factor<<" "<<_qcmapall[qc]<<" "<<_qcmapall[qc]*qc_eff_factor<<std::endl;

	    //	    std::cout << iqc <<" "<<_smearedqcmap[qc]<<" "<<_smearedqcmap[qc]*qc_eff_factor<<" "<<_qcmapall[qc]<<" qc "<<qc<<" eff "<<qc_eff_factor<<" "<<_qcmapall[qc]*qc_eff_factor<<std::endl;

	    // It's events per qc bin
	    totev += _smearedqcmap[qc]*qc_eff_factor;
	    totevunsmeared += _qcmapall[qc]*qc_eff_factor;
	  }
	}
      

	qcoutfile.close();

      
	cout << "Total qc events: "<<totev<<" unsmeared "<<totevunsmeared<<endl;

      } // End of do qc smearing case

    }  // End of do-quenching case


  return 0;
  
}


