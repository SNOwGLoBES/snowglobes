//#include "Math/Polynomial.h"
#include "Math/Interpolator.h"

void stpi_flux()
{

  using namespace std;

    // SNS flux parameters

  Double_t MeVperproton = 1010.;

    Double_t Jperproton = MeVperproton*1e6*1.6021e-19;
    Double_t Beampower = 1.4e6;
    Double_t Protonspersec = Beampower/Jperproton;
    Double_t Nusperprotonperflavor = 0.09;
    Double_t Nuspersecperflavor = Nusperprotonperflavor*Protonspersec;
    Double_t flux_per_s_percm2_at_20m = Nuspersecperflavor/(4*TMath::Pi()*2000*2000);
    cout << "Flux per second per cm^2 at 20 m per flavor: "<< flux_per_s_percm2_at_20m <<endl;

  // Read in the data from the file

  const Int_t maxpoints = 1000;
  Double_t en[maxpoints],nue[maxpoints],nuebar[maxpoints],nux[maxpoints];
  Double_t numu[maxpoints],nutau[maxpoints],numubar[maxpoints],nutaubar[maxpoints];


  const Int_t numipoints=501;
  
  Double_t ien[numipoints];
  Double_t inue[numipoints],inuebar[numipoints];
  Double_t inumu[numipoints],inutau[numipoints],inumubar[numipoints],inutaubar[numipoints];
  
  Int_t i;

  Double_t ie;

  Double_t minen=0.5;
  Double_t maxen=100.;

  // IMPORTANT TO USE THIS BINNING
  //  to match smearing matrices, or else Globes does not handle the 
  // monochromatic numu properly

  Double_t numipoints2 = 200;
  Double_t step=(Double_t)(maxen-minen)/(numipoints2);
  const Double_t mmu = 105.66837;

  cout << "Step  "<<step<<endl;
  Double_t extra = 0.;
  ie=minen+step/2.+extra;

  ofstream out;
  out.open("stpi2.dat");

  Double_t sumnue= 0;
  for (i=0;i<numipoints;i++) {

    inuebar[i]=0.;
    inutau[i]=0.;
    inutaubar[i]=0.;

    if (ie<=mmu/2.) {
      
      inumubar[i] = 2*pow(2./mmu*ie,2)*(3-2*((2./mmu)*ie))*(2./mmu);
      inue[i] = 12*pow(2./mmu*ie,2)*(1-(2./mmu*ie))*(2./mmu);

    } else {

      inumubar[i]=0;
      inue[i]=0;

    }


    ien[i]=ie;
    cout << i<<" "<<ie<<" "<<inue[i]<<" "<<inumu[i]<<" "<<
      inutau[i]<<" "<<inuebar[i]<<" "<<inumubar[i]<<" "<<inutaubar[i]<<endl;    
    // Output for Globes needs to be in GeV

    // Normalization factor

    // Scale for 1.5 km since Globes will multiply by 1/L^2

    //    Double_t distfac=(1.5e5)*(1.5e5);
    //    Double_t distfac= pow(20.,2)/pow(1000.,2);

    // Make it for 20 m

    Double_t distfac = 1.;

    // Time bin in the file, file is per sec
    //    Double_t tbinsize=0.1;

    // For integrated flux over time, factor is 1

    Double_t tbinsize = 1.;
    // File has per eV in bins of 125000 eV 
    //    Double_t ebinfac=125000;  //
    //    Double_t ebinfac=125000.;  //

    // This is for flux in the bin
    //    Double_t ebinfac=step;  //

    // This one is anomalous in binning due to the monochromatic numu issue
    //  Want flux per 0.2 MeV
    Double_t ebinfac=0.2;  

    // Norm: want flux in the bin at 1 km, assumed by Globes
    //  Aug 6, 2018:  this is for 0.08 protons per 1 MW, 


    Double_t fac=flux_per_s_percm2_at_20m*distfac*tbinsize*ebinfac;

    // All the numu goes in one bin
    //    Double_t numuen = 29.792;
    Double_t numuen = 29.6;

    if (fabs(ie-numuen)<step/2.) {
      inumu[i] = 1./step;
    }

    sumnue += inue[i]*fac;

    
    out <<TString::Format("%12.10f",ie/1000) <<" "<<inue[i]*fac<<" "<<inumu[i]*fac<<" "<<
      inutau[i]*fac<<" "<<inuebar[i]*fac<<" "<<inumubar[i]*fac<<" "<<inutaubar[i]*fac<<endl;    

    ie+=step;
  }

  out.close();

  cout << "Sumnue "<<sumnue<<endl;
  Double_t xmin = 0.;
  Double_t xmax = 55.;
  Double_t ymin = 0.;
  Double_t ymax = 0.1;

  TCanvas* canv = new TCanvas("c1"," ",800,700);

  TString plot_title = "";
  TH2F *hr = new TH2F("hr",plot_title,2,xmin,xmax,2,ymin,ymax);
  hr->SetXTitle("");
  hr->SetYTitle("   ");
  hr->GetXaxis()->SetTitleSize(0.03);
  hr->GetXaxis()->SetTitleOffset(0.85);
  hr->GetXaxis()->SetLabelSize(0.03);
  hr->GetYaxis()->SetTitleSize(0.03);
  hr->GetYaxis()->SetTitleOffset(0.7);
  hr->GetYaxis()->SetLabelSize(0.03);


  gStyle->SetStatH(0);

  hr->Draw();


  TGraph* gr0= new TGraph(numipoints,ien,inumubar);
  gr0->SetLineColor(1);
  gr0->Draw("same");


  TGraph* gr1= new TGraph(numipoints,ien,inue);
  gr1->SetLineColor(2);
  gr1->Draw("same");

  TGraph* gr2= new TGraph(numipoints,ien,inumu);
  gr2->SetLineColor(3);
  gr2->Draw("same");
  //  gr2->Print();

  cout << gr0->Integral() << endl;
  cout << gr1->Integral() << endl;
  cout << gr2->Integral() << endl;
}
