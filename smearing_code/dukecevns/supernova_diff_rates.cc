#include <iostream>
#include <fstream>

#include "FormFactor.h"
#include "NuFlux.h"
#include "DetectorResponse.h"
#include "TFile.h"
#include "TGraph.h"
#include "TObjArray.h"
#include "TString.h"
#include "TMath.h"
#include "TH2D.h"

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <math.h>
#include <map>

#include "xscns.h"

// Energies in MeV

int main(int argc, char * argv[] )
{

#include "isomaps.h"
  // Info for relevant mixtures
#include "mixtures.h"

  
  if (argc<4) {
    std::cout << "Usage:  ./supernova_diff_rates [target material] [form factor] [qfname] [Rfac] [photons per MeVee]"<<std::endl;
    exit(0);
  }


  // Set up the form factor

  const char * ffname = argv[2];

  const char * qfname = argv[3];

  double rfac=1.;
  if (argc >=5) {
     rfac = (double)atof(argv[4]);
  }


  // Photons per Eee in MeV

  double lightyield=24000.;
  if (argc >=6) {
     lightyield = (double)atof(argv[5]);
  }



  // Array of pointers to form factors
  FormFactor** ff;
  ff = new FormFactor*[max_components];

  // Array of pointers to quenching factors

  DetectorResponse** qffunc;
  qffunc = new DetectorResponse*[max_components];

    // Set up the flux  (here, fluence)

 //  NumericalFlux* livermore = new NumericalFlux();
//   livermore->SetFluxFilename("livermore.dat");
//   livermore->ReadFluxFile();
//   double kmax = livermore->maxEnu();
//   std::cout << "Max neutrino energy: "<<kmax<<std::endl;

  PinchedThermal* pinched = new PinchedThermal();

  double alpha[6]={3.,2.5,2.5,3.,2.5,2.5};
  double avgen[6]={10.,15.,15.,14.,15.,15.};
  double lumi[6]={1.6e52,1.6e52,1.6e52,1.6e52,1.6e52,1.6e52};
  pinched->SetAlpha(alpha);
  pinched->SetLuminosity(lumi);
  pinched->SetAvgEn(avgen);
  double kmax = pinched->maxEnu();

    // Set up the material

  std::string material = argv[1];

  std::ofstream outfile;
  std::string outfilename;
  outfilename = "out/supernova_diff_rates-"+material+"-"+std::string(ffname)+".out";
  outfile.open(outfilename);


  std::ofstream phoutfile;
  std::string phoutfilename;
  phoutfilename = "out/supernova_diff_rates-"+material+"-"+std::string(ffname)+"-photons.out";
  phoutfile.open(phoutfilename);
  
   // Array for quenched total rates

   const int maxiq = 10000;
   double Eee[max_components][maxiq];
   double dNdEee[max_components][maxiq];
   double dNdEr[max_components][maxiq];

  //  std::cout << "Material "<<material << std::endl;

  //  std::string material = "Ar";

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

  // First loop over isotope components
  // First get the total mass.  Get also the maximum recoil values
  // Set up form factors and quenching factors for each

  double erecmaxvals[max_components];
  double Mtot = 0;
  v = isotope_component.begin();

  double minM = 1.e10;
  int num_components;

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

    std::string isoname = std::string(isotope);

    if (strcmp(ffname, "helm")==0) {
      Helm* helmff= new Helm();
      ff[is] = helmff;
      helmff->Setsval(0.9);
    }
    else if (strcmp(ffname, "klein")==0) {
      Klein* kleinff = new Klein();
      ff[is] = kleinff;
      kleinff->Setakval(0.7);
    } 
    else if  (strcmp(ffname, "horowitz")==0){
      Horowitz* horowitzff = new Horowitz();
      ff[is] = horowitzff;
      std::transform(isoname.begin(), isoname.end(),isoname.begin(), ::toupper);
      std::string horowitz_filename = "ff/"+isoname+".FF";
      //      std::cout << horowitz_filename << std::endl;
      horowitzff->SetFFfilename(horowitz_filename.c_str());
      horowitzff->ReadFFfile();
      horowitzff->SetRfac(rfac);
    }

    A = Nn + Z;
    ff[is]->SetA(A);
    ff[is]->SetZ(Z);
 


  // Set up detector quenching factors for each component
    
    std::string qffilename;
    qffilename = "qf/"+std::string(qfname)+"_"+isoname+"_qf.txt";

    std::cout << "Quenching factor: "<<qffilename<<std::endl;
    DetectorResponse* qf = new DetectorResponse();
    qffunc[is] = qf;
    qffunc[is]->SetQFPolyFilename(qffilename.c_str());
    qffunc[is]->ReadQFPolyFile();


    v++; is++;
  }  // End of first loop over isotope components


  num_components = is;

  // Use the mass of the lightest component
  double erecmaxall = 2*kmax*kmax/(minM+2*kmax);
  
  double erecstart = 0.;
  double erecend = erecmaxall;
  //  double erecstep = 0.0001;
  double erecstep = 0.001;

   // Now compute the differential recoil spectra

    double Erec;
    double knu;

    std::cout << "erecmaxall "<<erecmaxall<<std::endl;

    //   double knustep = 0.0001;
    double knustep = 0.001;

   // The totals
   double toterecoil = 0.;
   double totevents = 0.;
   

   int iq=0;
   // Loop over recoil energy
   for (Erec=erecstart+erecstep;Erec<=erecend; Erec+=erecstep) {

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
	  Double_t hbarc = 197.327; // MeV-fm, convert for Q in MeV for ff
	  double Q = sqrt(2*M*Erec+Erec*Erec); // MeV
	  double qq = Q/hbarc;
	  //    double ff2 = helmff->FFval(qq);

	  double ff2 = pow(ff[is]->FFval(qq),2);

	  // SM Couplings

	  double GV_sm = GV_SM(2015,Z,Nn);
	  double GA_sm = GA_SM(2015,1,Z,Nn,Zdiff,Ndiff);
	  double GA_sm_bar = GA_SM(2015,-1,Z,Nn,Zdiff,Ndiff);

	// Normalize for one ton of material
	// Weight by mass fraction
	
	  double Nt = 1.e6/(M/amu)*6.022e23;

	// A2: G^2/(2Pi) * hbarcinmeters^-4 
	  double norm = Nt;

	  // Quenching factor for this component and Eee for this Erec
	  
	  Eee[is][iq] = qffunc[is]->qfpoly(Erec)*Erec;
	  double qfderiv = abs(qffunc[is]->qfpolyderiv(Erec));
	
	  // For the sum over neutrino energy
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

	  for (knu=knumin;knu<=kmax;knu+=knustep) {


	    drate_e_vec += diffxscnvec(knu,M,Erec)*pinched->fluxval(knu,1,knustep);
	    drate_ebar_vec += diffxscnvec(knu,M,Erec)*pinched->fluxval(knu,-1,knustep);
	    drate_mu_vec += diffxscnvec(knu,M,Erec)*pinched->fluxval(knu,2,knustep);
	    drate_mubar_vec += diffxscnvec(knu,M,Erec)*pinched->fluxval(knu,-2,knustep);
	    drate_tau_vec +=  diffxscnvec(knu,M,Erec)*pinched->fluxval(knu,3,knustep);
	    drate_taubar_vec +=  diffxscnvec(knu,M,Erec)*pinched->fluxval(knu,-3,knustep);


	    drate_e_axial += diffxscnaxial(knu,M,Erec)*pinched->fluxval(knu,1,knustep);
	    drate_ebar_axial += diffxscnaxial(knu,M,Erec)*pinched->fluxval(knu,-1,knustep);
	    drate_mu_axial += diffxscnaxial(knu,M,Erec)*pinched->fluxval(knu,2,knustep);
	    drate_mubar_axial += diffxscnaxial(knu,M,Erec)*pinched->fluxval(knu,-2,knustep);
	    drate_tau_axial +=  diffxscnaxial(knu,M,Erec)*pinched->fluxval(knu,3,knustep);
	    drate_taubar_axial +=  diffxscnaxial(knu,M,Erec)*pinched->fluxval(knu,-3,knustep);


	    drate_e_interf += diffxscninterf(knu,M,Erec)*pinched->fluxval(knu,1,knustep);
	    drate_ebar_interf += diffxscninterf(knu,M,Erec)*pinched->fluxval(knu,-1,knustep);
	    drate_mu_interf += diffxscninterf(knu,M,Erec)*pinched->fluxval(knu,2,knustep);
	    drate_mubar_interf += diffxscninterf(knu,M,Erec)*pinched->fluxval(knu,-2,knustep);
	    drate_tau_interf +=  diffxscninterf(knu,M,Erec)*pinched->fluxval(knu,3,knustep);
	    drate_taubar_interf += diffxscninterf(knu,M,Erec)*pinched->fluxval(knu,-3,knustep);

	      
	    //	    std::cout << Erec << " "<<knu<<" "<<std::endl;

	  } // End of loop over neutrino energy contributions
	    

	  // Now multiply by target-dependent factors and add up this recoil energy bin

	  diffrate_e_vec[is] = norm*pow(GV_sm,2)*ff2*mass_fraction[is]*drate_e_vec;
	  diffrate_ebar_vec[is] = norm*pow(GV_sm,2)*ff2*mass_fraction[is]*drate_ebar_vec;
	  diffrate_mu_vec[is] = norm*pow(GV_sm,2)*ff2*mass_fraction[is]*drate_mu_vec;
	  diffrate_mubar_vec[is] = norm*pow(GV_sm,2)*ff2*mass_fraction[is]*drate_mubar_vec;
	  diffrate_tau_vec[is] =  norm*pow(GV_sm,2)*ff2*mass_fraction[is]*drate_tau_vec;
	  diffrate_taubar_vec[is] = norm*pow(GV_sm,2)*ff2*mass_fraction[is]*drate_taubar_vec;


	  diffrate_e_axial[is] = norm*pow(GA_sm,2)*ff2*mass_fraction[is]*drate_e_axial;
	  diffrate_ebar_axial[is] = norm*pow(GA_sm_bar,2)*ff2*mass_fraction[is]*drate_ebar_axial;
	  diffrate_mu_axial[is] = norm*pow(GA_sm,2)*ff2*mass_fraction[is]*drate_mu_axial;
	  diffrate_mubar_axial[is] = norm*pow(GA_sm_bar,2)*ff2*mass_fraction[is]*drate_mubar_axial;
	  diffrate_tau_axial[is] =  norm*pow(GA_sm,2)*ff2*mass_fraction[is]*drate_tau_axial;
	  diffrate_taubar_axial[is] =  norm*pow(GA_sm_bar,2)*ff2*mass_fraction[is]*drate_taubar_axial;


	  diffrate_e_interf[is] = norm*GV_sm*GA_sm*ff2*mass_fraction[is]*drate_e_interf;
	  diffrate_ebar_interf[is] = norm*GV_sm*GA_sm_bar*ff2*mass_fraction[is]*drate_ebar_interf;
	  diffrate_mu_interf[is] = norm*GV_sm*GA_sm*ff2*mass_fraction[is]*drate_mu_interf;
	  diffrate_mubar_interf[is] = norm*GV_sm*GA_sm_bar*ff2*mass_fraction[is]*drate_mubar_interf;
	  diffrate_tau_interf[is] =  norm*GV_sm*GA_sm*ff2*mass_fraction[is]*drate_tau_interf;
	  diffrate_taubar_interf[is] =  norm*GV_sm*GA_sm_bar*ff2*mass_fraction[is]*drate_taubar_interf;

	    //	    std::cout << is<<" "<<"Erec "<<Erec<<" mass frac "<<mass_fraction[is]<<" "<<ff2<<" "<<diffrate_e<<std::endl;

	  // Now add the contribution from this isotope to the sum
	  
	  sum_diffrate_e_vec += diffrate_e_vec[is] ;
	  sum_diffrate_ebar_vec += diffrate_ebar_vec[is];
	  sum_diffrate_mu_vec += diffrate_mu_vec[is];
	  sum_diffrate_mubar_vec += diffrate_mubar_vec[is];
	  sum_diffrate_tau_vec += diffrate_tau_vec[is];
	  sum_diffrate_taubar_vec += diffrate_taubar_vec[is];
	  
	  sum_diffrate_e_axial += diffrate_e_axial[is];
	  sum_diffrate_ebar_axial += diffrate_ebar_axial[is];
	  sum_diffrate_mu_axial += diffrate_mu_axial[is];
	  sum_diffrate_mubar_axial += diffrate_mubar_axial[is];
	  sum_diffrate_tau_axial += diffrate_tau_axial[is];
	  sum_diffrate_taubar_axial += diffrate_taubar_axial[is];
	  
	  sum_diffrate_e_interf += diffrate_e_interf[is];
	  sum_diffrate_ebar_interf += diffrate_ebar_interf[is];
	  sum_diffrate_mu_interf += diffrate_mu_interf[is];
	  sum_diffrate_mubar_interf += diffrate_mubar_interf[is];
	  sum_diffrate_tau_interf += diffrate_tau_interf[is];
	  sum_diffrate_taubar_interf += diffrate_taubar_interf[is];

	  // Sum for this Erec and isotope
	  double sum_events_iso = 0;
	  sum_events_iso = diffrate_e_vec[is] + diffrate_ebar_vec[is] + diffrate_mu_vec[is]+ diffrate_mubar_vec[is]+ diffrate_tau_vec[is] + diffrate_taubar_vec[is];
	  sum_events_iso += diffrate_e_axial[is] + diffrate_ebar_axial[is] + diffrate_mu_axial[is]+ diffrate_mubar_axial[is]+ diffrate_tau_axial[is] + diffrate_taubar_axial[is];

	  sum_events_iso+= diffrate_e_interf[is] + diffrate_ebar_interf[is] + diffrate_mu_interf[is]+ diffrate_mubar_interf[is]+ diffrate_tau_interf[is] + diffrate_taubar_interf[is];

	  // Sum for this Erec, Eee and isotope

	  // Now apply the quenching for this Ee and isotope component
	// sum_events_iso is dNderec

	  dNdEr[is][iq] = sum_events_iso;
	    
	  if (qfderiv>0) {
	    dNdEee[is][iq] = sum_events_iso/qfderiv;
	  } else {
	    dNdEee[is][iq] = 0.;
	  }

	  v++;is++;

	  
	} // End of loop over material components

	// This is events per MeV

     // This is the total for all components 
	double events=0;
	events = sum_diffrate_e_vec + sum_diffrate_ebar_vec + sum_diffrate_mu_vec+ sum_diffrate_mubar_vec+ sum_diffrate_tau_vec + sum_diffrate_taubar_vec;
	events += sum_diffrate_e_axial + sum_diffrate_ebar_axial + sum_diffrate_mu_axial+ sum_diffrate_mubar_axial+ sum_diffrate_tau_axial + sum_diffrate_taubar_axial;
        events += sum_diffrate_e_interf + sum_diffrate_ebar_interf + sum_diffrate_mu_interf+ sum_diffrate_mubar_interf+ sum_diffrate_tau_interf + sum_diffrate_taubar_interf;

	std::cout << Erec<<" "<<events<<" "<<sum_diffrate_e_vec<<" "<<sum_diffrate_ebar_vec<<" "<<sum_diffrate_mu_vec<<" "<<sum_diffrate_mubar_vec<<" "<<sum_diffrate_tau_vec<<" "<<sum_diffrate_taubar_vec<<" "<<sum_diffrate_e_axial<<" "<<sum_diffrate_ebar_axial<<" "<<sum_diffrate_mu_axial<<" "<<sum_diffrate_mubar_axial<<" "<<sum_diffrate_tau_axial<<" "<<sum_diffrate_taubar_axial<<" "<<sum_diffrate_e_interf<<" "<<sum_diffrate_ebar_interf<<" "<<sum_diffrate_mu_interf<<" "<<sum_diffrate_mubar_interf<<" "<<sum_diffrate_tau_interf<<" "<<sum_diffrate_taubar_interf <<std::endl;
	outfile << Erec<<" "<<events<<" "<<sum_diffrate_e_vec<<" "<<sum_diffrate_ebar_vec<<" "<<sum_diffrate_mu_vec<<" "<<sum_diffrate_mubar_vec<<" "<<sum_diffrate_tau_vec<<" "<<sum_diffrate_taubar_vec<<" "<<sum_diffrate_e_axial<<" "<<sum_diffrate_ebar_axial<<" "<<sum_diffrate_mu_axial<<" "<<sum_diffrate_mubar_axial<<" "<<sum_diffrate_tau_axial<<" "<<sum_diffrate_taubar_axial<<" "<<sum_diffrate_e_interf<<" "<<sum_diffrate_ebar_interf<<" "<<sum_diffrate_mu_interf<<" "<<sum_diffrate_mubar_interf<<" "<<sum_diffrate_tau_interf<<" "<<sum_diffrate_taubar_interf <<std::endl;
   
	// Now write the quenched output
	// Note this will be per MeVee, but uneven bins and different x scale for each component

	// Loop over components
	int is2=0;
	phoutfile  << Erec<< "  "<<events<<" ";
	for (is2=0;is2<num_components;is2++) {
	  phoutfile <<dNdEr[is2][iq]<<" "<<Eee[is2][iq]<<" "<<dNdEee[is2][iq]<<" ";
	}
	phoutfile <<std::endl;

	iq++;

	totevents+=events*erecstep;

	toterecoil += events*Erec*erecstep;

  } // End of loop over Erec

   double time_interval = 10.; // Assume flux is over 10 seconds
   std::cout << "Total events:  "<< totevents*time_interval<< std::endl;
   std::cout << "Total recoil energy deposited:  "<< toterecoil*time_interval<< std::endl;

   outfile.close();
   phoutfile.close();


   // Interpolate TGraphs of quenched differential spectrum to get evenly spaced Eee bins

   std::ofstream phoutfile2;
   std::string phoutfilename2;
   phoutfilename2 = "out/supernova_diff_rates-"+material+"-"+std::string(ffname)+"-photons2.out";
    phoutfile2.open(phoutfilename2);

   int nquenched = iq;  // Number of quenched points
   TGraph** quenchedspecs;
   quenchedspecs = new TGraph*[max_components];
   
   Double_t Energy_ee[maxiq];
   Double_t quenched_diffspec[maxiq];

   Double_t max_enee=0;
   for (is=0;is<num_components;is++) {
     // Probably a more elegant way to do this
     
     Int_t iq2;

     for (iq2=0;iq2<nquenched;iq2++) {
       Energy_ee[iq2] = Eee[is][iq2];

       if (Energy_ee[iq2]>max_enee) {
	 max_enee = Energy_ee[iq2];
       }
       quenched_diffspec[iq2] = dNdEee[is][iq2];

     }
     quenchedspecs[is]= new TGraph(nquenched,Energy_ee,quenched_diffspec);


   } // End of loop over components


   // Now write interpolated quenched spectra to file

   double enee;
   double eneestep=0.0001;
   double sumquencheden=0;
   double sumquenchedevents=0;

   std::cout << "Max enee "<<max_enee<<std::endl;
   for (enee=0;enee<=max_enee;enee+=eneestep) {

     std::cout << "Quenched: "<<enee << " ";
     phoutfile2 << enee << " ";
     
     Double_t totquencheden=0;
     Double_t totquenchedevents=0;
     // Loop over components
     for (is=0;is<num_components;is++) {

       double qspecval = quenchedspecs[is]->Eval(enee);
       std::cout <<qspecval<<" ";
       phoutfile2 <<qspecval<<" ";
       totquenchedevents += qspecval;
       totquencheden +=qspecval*enee;
     }
     std::cout<<totquencheden<<std::endl;
     sumquenchedevents += totquenchedevents;
     sumquencheden += totquencheden*eneestep;

     std::cout<<totquenchedevents<<std::endl;

   } // End of loop over quenched energies



   std::cout<< "Total quenched energy deposited in MeV "<<sumquencheden*time_interval<<" events "<<sumquenchedevents*time_interval<<" photons: "<<sumquencheden*lightyield*time_interval<<std::endl;
    phoutfile2.close();
   

   return 0;

}


