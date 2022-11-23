#include <iostream>
#include <fstream>

#include "DetectorResponse.h"
#include "FormFactor.h"
#include "NuFlux.h"
// #include "TFile.h"
// #include "TGraph.h"
// #include "TObjArray.h"
// #include "TString.h"
// #include "TMath.h"
// #include "TH2D.h"

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <math.h>
#include <map>

#include "xscns.h"

void get_flavor_weight(double,double,double*, double*, double*);
// Check the CsI calculation

int main(int argc, char * argv[] )
{

#include "isomaps.h"
  // Info for relevant mixtures
#include "mixtures.h"


  // Flavor weighting from file

  double wnumu=1.;
  double wnumubar=1.;
  double wnue=1.;


  // Don't use flavor weights if using snsflux numerical flux; it that should take care of the weighting
  //    get_flavor_weight(1400.,7400.,&wnumu,&wnumubar,&wnue);
  get_flavor_weight(0.,6000.,&wnumu,&wnumubar,&wnue);

  std::cout << "Flavor weights: "<< wnumu<<" "<<wnumubar<<" "<<wnue<<std::endl;

  // Set up the form factor

  const char * ffname = "helm";

  // Array of pointers to form factors, for protons and neutrons separately, axial and vector separately
  // (although small differences

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
    snsflux->SetNorm(5.2e14/(4*M_PI*2900.*2900.));
  // Gives flux per pidk per energy bin per second, energy bin in MeV,  normalize for 5e14 decays/s 

  // Set up the detector response

  
    //   DetectorResponse* csitl = new DetectorResponse();

    //csitl->SetEfficFilename("bjorn-eff.dat");
    //csitl->ReadEfficFile();

    // Set up the material

    //    std::string material = "Na23";
    std::string material = "Ar";

  std::ofstream outfile;
  std::string outfilename;
  outfilename = "lar_sns_diff_rates-"+material+"-"+std::string(ffname)+".out";
  outfile.open(outfilename);
 

  double M;
  double Delta;
  int Nn,Z,A;
  int Zdiff, Ndiff;

  std::string matname = material;

 // These are defined in mixtures include
  std::vector<double> fraction = molar_fraction[material];
  std::vector<std::string> isotope_component= isotopes[material];

  std::cout << "Material "<<matname<<std::endl;


  int is=0;
  std::vector<std::string>::iterator v = isotope_component.begin();
  std::string isotope;

  // First get the total mass.  Get also the maximum recoil values

  double erecmaxvals[max_components];
  double Mtot = 0;
  v = isotope_component.begin();

  double minM = 1.e10;
  while( v != isotope_component.end()) {

    isotope = *v;
        std::cout << "isotope"<< isotope << std::endl;

    Z = Zs[std::string(isotope)];
    Nn = Ns[std::string(isotope)];
    Delta = Deltas[std::string(isotope)];
    M = (Z+Nn)*amu - Z*me + Delta;
    if (M<minM) {minM=M;}
    Mtot += M*fraction[is];
    erecmaxvals[is] = 2*kmax*kmax/(M+2*kmax);


    // Set up the form factor for this isotope

    //       std::cout << "Mass "<<M<<std::endl;


    double rfac=1.;
    if (strcmp(ffname, "helm")==0) {
      Helm* helmffnv= new Helm();
      ffnv[is] = helmffnv;
      helmffnv->Setsval(0.9);

      Helm* helmffna= new Helm();
      ffna[is] = helmffna;
      helmffna->Setsval(0.9);

      Helm* helmffpv= new Helm();
      ffpv[is] = helmffpv;
      helmffpv->Setsval(0.9);


      Helm* helmffpa= new Helm();
      ffpa[is] = helmffpa;
      helmffpa->Setsval(0.9);


    }
    else if (strcmp(ffname, "klein")==0) {

      Klein* kleinffnv = new Klein();
      ffnv[is] = kleinffnv;
      kleinffnv->Setakval(0.7);
      kleinffnv->SetRfac(1.);
      kleinffnv->Setskinfac(1.);

      Klein* kleinffna = new Klein();
      ffna[is] = kleinffna;
      kleinffna->Setakval(0.7);
      kleinffna->SetRfac(1.);
      kleinffna->Setskinfac(1.);

      Klein* kleinffpv = new Klein();
      ffpv[is] = kleinffpv;
      kleinffpv->Setakval(0.7);
      kleinffpv->SetRfac(1.);
      kleinffpv->Setskinfac(0.);


      Klein* kleinffpa = new Klein();
      ffpa[is] = kleinffpa;
      kleinffpa->Setakval(0.7);
      kleinffpa->SetRfac(1.);
      kleinffpv->Setskinfac(0.);


    } 
    else if  (strcmp(ffname, "horowitz")==0){
      Horowitz* horowitzffnv = new Horowitz();
      ffnv[is] = horowitzffnv;

      Horowitz* horowitzffna = new Horowitz();
      ffna[is] = horowitzffna;

      Horowitz* horowitzffpv = new Horowitz();
      ffnv[is] = horowitzffpv;

      Horowitz* horowitzffpa = new Horowitz();
      ffna[is] = horowitzffpa;


      std::string isoname = std::string(isotope);
      std::transform(isoname.begin(), isoname.end(),isoname.begin(), ::toupper);
      std::string horowitz_filename = isoname+".FF";
      //      std::cout << horowitz_filename << std::endl;
      horowitzffnv->SetFFfilename(horowitz_filename.c_str());
      horowitzffnv->ReadFFfile();
      horowitzffnv->SetRfac(rfac);

      horowitzffna->SetFFfilename(horowitz_filename.c_str());
      horowitzffna->ReadFFfile();
      horowitzffna->SetRfac(rfac);

      //      std::cout << horowitz_filename << std::endl;
      // Not really appropriate for protons, but using the structure
      horowitzffpv->SetFFfilename(horowitz_filename.c_str());
      horowitzffpv->ReadFFfile();
      horowitzffpv->SetRfac(1.);

      horowitzffpa->SetFFfilename(horowitz_filename.c_str());
      horowitzffpa->ReadFFfile();
      horowitzffpa->SetRfac(1.);

    }

    A = Nn + Z;
    ffnv[is]->SetA(A);
    ffna[is]->SetA(A);
    ffpv[is]->SetA(A);
    ffpa[is]->SetA(A);

    ffnv[is]->SetZ(Z);
    ffna[is]->SetZ(Z);
    ffpv[is]->SetZ(Z);
    ffpa[is]->SetZ(Z);
 

    v++; is++;
  }



  // Use the mass of the lightest component
  double erecmaxall = 2*kmax*kmax/(minM+2*kmax);
  
  //double recoilthresh = 0.013; //MeVr
  double recoilthresh = 0.; //MeVr
  double erecstart = recoilthresh;
  double erecend = erecmaxall;
  //  double erecstep = 0.0001;
  double erecstep = 0.0001;

   // Now compute the differential recoil spectra

     
    double Erec;
    double knu;

    std::cout << "erecmaxall "<<erecmaxall<<std::endl;

    //double knustep = 1.;
    double knustep = 0.0001;

   // The totals
   double toterecoil = 0.;
   double totevents = 0.;
   
   // Loop over recoil energy
   for (Erec=erecstart+erecstep;Erec<=erecend; Erec+=erecstep) {
     
     double diffrate_e_vec = 0;
     double diffrate_ebar_vec = 0;
     double diffrate_mu_vec = 0;
     double diffrate_mubar_vec = 0;
     double diffrate_tau_vec = 0;
     double diffrate_taubar_vec = 0;
     
     double diffrate_e_axial = 0;
     double diffrate_ebar_axial = 0;
     double diffrate_mu_axial = 0;
     double diffrate_mubar_axial = 0;
     double diffrate_tau_axial = 0;
     double diffrate_taubar_axial = 0;
     
     double diffrate_e_interf = 0;
     double diffrate_ebar_interf = 0;
     double diffrate_mu_interf = 0;
     double diffrate_mubar_interf = 0;
     double diffrate_tau_interf = 0;
     double diffrate_taubar_interf = 0;
     
     // With efficiency, which is a function of Erec in MeV in this formuation
     
     //     double eff_factor = csitl->efficnum(Erec);

     
     double eff_factor = 1.;
     std::cout << "eff factor "<<eff_factor<<std::endl;
     // Skip if too small contribution
     if (eff_factor>0.0001) {
	  

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

	 double qq = Q/hbarc;
	 //    double ff2 = helmff->FFval(qq);
	 double ffnvval = ffnv[is]->FFval(qq);
	 double ffnaval = ffna[is]->FFval(qq);
	 double ffpvval = ffpv[is]->FFval(qq);
	 double ffpaval = ffpa[is]->FFval(qq);
 
	 //	 double ff2 = pow(ff[is]->FFval(qq),2);
	 //std::cout << "knumin, Erec, Q, ff2 "<<knumin<<" "<<Erec<<" "<<Q<<" "<<ff2<<" "<<M<<" "<<mass_fraction[is]<<std::endl;
	    
	 // SM Couplings

	 //	 double GV_sm = GV_SM(2015,Z,Nn);
	 //double GA_sm = GA_SM(2015,1,Z,Nn,Zdiff,Ndiff);
	 //double GA_sm_bar = GA_SM(2015,-1,Z,Nn,Zdiff,Ndiff);
	   
	 double gv[2], ga[2], gabar[2];
	 int pdgyr = 2015;
	 sm_vector_couplings(pdgyr,gv);
	 sm_axial_couplings(pdgyr,1,ga);
	 sm_axial_couplings(pdgyr,-1,gabar);

	 // Bundle the form factor contributions with the SM couplings, separately for p and n
	 double GV_sm_wff = Z*gv[0]*ffpvval+Nn*gv[1]*ffnvval;
	 double GA_sm_wff = Zdiff*ga[0]*ffpaval+Ndiff*ga[1]*ffnaval;
	 double GA_sm_bar_wff = Zdiff*gabar[0]*ffpaval+Ndiff*gabar[1]*ffnaval;
 
	 // Normalize for one ton of material
	 // Will be weighted by mass fraction
	    
	 double Nt = 1.e6/(M/amu)*6.0221409e23;
	    
	 // A2: G^2/(2Pi) * hbarcinmeters^-4 
	 double norm = Nt;
	    
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

	      
	   //    std::cout << Erec << " "<<knu<<" "<<snsflux->fluxval(knu,1,knustep)<<" "<<diffrate_e_vec<<std::endl;

	 } // End of loop over neutrino energy contributions
	   
	 // Now multiply by target-dependent factors and add up this recoil energy bin

	 //	 double mufact = 1.037; too large!
	 double mufact = mufactor(Q);
	 diffrate_e_vec += norm*pow(GV_sm_wff,2)*mass_fraction[is]*drate_e_vec;
	 diffrate_ebar_vec += norm*pow(GV_sm_wff,2)*mass_fraction[is]*drate_ebar_vec;
	 diffrate_mu_vec += norm*pow(GV_sm_wff,2)*mass_fraction[is]*drate_mu_vec*mufact;
	 diffrate_mubar_vec += norm*pow(GV_sm_wff,2)*mass_fraction[is]*drate_mubar_vec*mufact;
	 diffrate_tau_vec +=  norm*pow(GV_sm_wff,2)*mass_fraction[is]*drate_tau_vec;
	 diffrate_taubar_vec += norm*pow(GV_sm_wff,2)*mass_fraction[is]*drate_taubar_vec;


	 diffrate_e_axial += norm*pow(GA_sm_wff,2)*mass_fraction[is]*drate_e_axial;
	 diffrate_ebar_axial += norm*pow(GA_sm_bar_wff,2)*mass_fraction[is]*drate_ebar_axial;
	 diffrate_mu_axial += norm*pow(GA_sm_wff,2)*mass_fraction[is]*drate_mu_axial*mufact;
	 diffrate_mubar_axial += norm*pow(GA_sm_bar_wff,2)*mass_fraction[is]*drate_mubar_axial*mufact;
	 diffrate_tau_axial +=  norm*pow(GA_sm_wff,2)*mass_fraction[is]*drate_tau_axial;
	 diffrate_taubar_axial +=  norm*pow(GA_sm_bar_wff,2)*mass_fraction[is]*drate_taubar_axial;


	 diffrate_e_interf += norm*GV_sm_wff*GA_sm_wff*mass_fraction[is]*drate_e_interf;
	 diffrate_ebar_interf += norm*GV_sm_wff*GA_sm_bar_wff*mass_fraction[is]*drate_ebar_interf;
	 diffrate_mu_interf += norm*GV_sm_wff*GA_sm_wff*mass_fraction[is]*drate_mu_interf*mufact;
	 diffrate_mubar_interf += norm*GV_sm_wff*GA_sm_bar_wff*mass_fraction[is]*drate_mubar_interf*mufact;
	 diffrate_tau_interf +=  norm*GV_sm_wff*GA_sm_wff*mass_fraction[is]*drate_tau_interf;
	 diffrate_taubar_interf +=  norm*GV_sm_wff*GA_sm_bar_wff*mass_fraction[is]*drate_taubar_interf;

	 std::cout << is<<" "<<"Erec "<<Erec<<" mass frac "<<mass_fraction[is]<<" "<<" norm "<<norm<<" GV "<<GV_sm_wff<<" GA "<<GA_sm_wff<<std::endl;

	 cout << "is, rate "<<is<<" "<<diffrate_e_vec<<endl;

	 v++;is++;


       } // End of loop over material components

     } // End of efficiency factor check

     std::cout << Erec<<" "<<diffrate_e_vec<<" "<<diffrate_ebar_vec<<" "<<diffrate_mu_vec<<" "<<diffrate_mubar_vec<<" "<<diffrate_tau_vec<<" "<<diffrate_taubar_vec<<" "<<diffrate_e_axial<<" "<<diffrate_ebar_axial<<" "<<diffrate_mu_axial<<" "<<diffrate_mubar_axial<<" "<<diffrate_tau_axial<<" "<<diffrate_taubar_axial<<" "<<diffrate_e_interf<<" "<<diffrate_ebar_interf<<" "<<diffrate_mu_interf<<" "<<diffrate_mubar_interf<<" "<<diffrate_tau_interf<<" "<<diffrate_taubar_interf <<std::endl;



     //	double detector_mass = 0.01457; // tons
     double detector_mass = 0.750; // tons
	double exposure = 365.25*24.*3600.; //  1 year * 0.892

	// Overall norm factor

	double norm_factor = eff_factor*detector_mass*exposure;

	diffrate_e_vec *= norm_factor*wnue;
	diffrate_ebar_vec *= norm_factor;
	diffrate_mu_vec *= norm_factor*wnumu;
	diffrate_mubar_vec *= norm_factor*wnumubar;
	diffrate_tau_vec *= norm_factor;
	diffrate_taubar_vec *= norm_factor;

	diffrate_e_axial *= norm_factor*wnue;
	diffrate_ebar_axial *= norm_factor;
	diffrate_mu_axial *= norm_factor*wnumu;
	diffrate_mubar_axial *= norm_factor*wnumubar;
	diffrate_tau_axial *= norm_factor;
	diffrate_taubar_axial *= norm_factor;

	diffrate_e_interf *= norm_factor*wnue;
	diffrate_ebar_interf *= norm_factor;
	diffrate_mu_interf *= norm_factor*wnumu;
	diffrate_mubar_interf *= norm_factor*wnumubar;
	diffrate_tau_interf *= norm_factor;
	diffrate_taubar_interf *= norm_factor;



	outfile << Erec<<" "<<diffrate_e_vec<<" "<<diffrate_ebar_vec<<" "<<diffrate_mu_vec<<" "<<diffrate_mubar_vec<<" "<<diffrate_tau_vec<<" "<<diffrate_taubar_vec<<" "<<diffrate_e_axial<<" "<<diffrate_ebar_axial<<" "<<diffrate_mu_axial<<" "<<diffrate_mubar_axial<<" "<<diffrate_tau_axial<<" "<<diffrate_taubar_axial<<" "<<diffrate_e_interf<<" "<<diffrate_ebar_interf<<" "<<diffrate_mu_interf<<" "<<diffrate_mubar_interf<<" "<<diffrate_tau_interf<<" "<<diffrate_taubar_interf <<std::endl;
   
	double events=0;
	events = diffrate_e_vec + diffrate_ebar_vec + diffrate_mu_vec+ diffrate_mubar_vec+ diffrate_tau_vec + diffrate_taubar_vec;
	events += diffrate_e_axial + diffrate_ebar_axial + diffrate_mu_axial+ diffrate_mubar_axial+ diffrate_tau_axial + diffrate_taubar_axial;
        events += diffrate_e_interf + diffrate_ebar_interf + diffrate_mu_interf+ diffrate_mubar_interf+ diffrate_tau_interf + diffrate_taubar_interf;

	totevents+=events*erecstep;

	toterecoil += events*Erec*erecstep;

  } // End of loop over Erec

   std::cout << "Total events over "<< recoilthresh*1000.<<" keVr: "<<totevents<< std::endl;
   std::cout << "Total recoil energy deposited:  "<< toterecoil<< std::endl;

   return 0;
   outfile.close();
}


