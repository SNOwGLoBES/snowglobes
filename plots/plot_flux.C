// Plot fluxes from the GLoBES-formatted flux files

void plot_flux(TString fluxname, Int_t titleopt)
{

  using namespace std;


  // Read in the data from the file

  ifstream in;
  //  in.open("../fluxes/kneller.dat");
  TString fluxfilename = "../fluxes/"+fluxname+".dat";
  in.open(fluxfilename);
  Int_t nlines = 0;

  const Int_t maxpoints = 1000;
  Double_t en[maxpoints],nue[maxpoints],nuebar[maxpoints],nux[maxpoints];
  Double_t numu[maxpoints],nutau[maxpoints],numubar[maxpoints],nutaubar[maxpoints];

  Int_t i=0;
  while (1) {
    //    in >> en[i]>>nue[i]>>nuebar[i]>>numu[i]>>numubar[i]>>nutau[i]>>nutaubar[i];
    in >> en[i]>>nue[i]>>numu[i]>>nutau[i]>>nuebar[i]>>numubar[i]>>nutaubar[i];

    // Energy in file is in eV
    en[i] *= 1000.;  // in MeV for the plots

    // File has actually flux per bin

    Double_t binfact=1;
    nue[i] /=  binfact;
    nuebar[i] /=  binfact;
    numu[i] /=  binfact;
    numubar[i] /=  binfact;
    nutau[i] /=  binfact;
    nutaubar[i] /=  binfact;

    nux[i] = numu[i]+numubar[i]+nutau[i]+nutaubar[i];
    
    cout <<i<<" "<<en[i]<<" "<<nue[i]<<" "<<nuebar[i]<<" "<<nux[i]<<endl;


    if (!in.good()) break;
    i++;

  }

  Int_t numpoints=i;

  printf(" found %d points\n",numpoints);

  in.close();

  Double_t xmin = 0.;
  Double_t xmax = 100.;
  Double_t ymin = 0.;

  // Because graph method not working
  Int_t p;
  Double_t ymax=0;
  for (p=0;p<numpoints;p++) {
    if (nux[p]>ymax) {
      ymax = nux[p];
    }
  }

  ymax*=1.1;

  TCanvas* canv = new TCanvas("c1"," ",800,700);

  TString plot_title;
  if (titleopt == 0) {
    plot_title = " ";
  } else {
    plot_title= "Flux: "+fluxname;
  }

  gStyle->SetTitleFontSize(.05);

  TH2F *hr = new TH2F("hr",plot_title,2,xmin,xmax,2,ymin,ymax);
  hr->SetXTitle(" Neutrino Energy (MeV) ");

  // Note it may be fluence, if time-integrated
  hr->SetYTitle(" Fluence (neutrinos per 0.2 MeV per cm^{2})");
  hr->GetXaxis()->SetTitleSize(0.04);
  hr->GetXaxis()->SetTitleOffset(0.85);
  hr->GetXaxis()->SetLabelSize(0.04);
  hr->GetYaxis()->SetTitleSize(0.04);
  hr->GetYaxis()->SetTitleOffset(1.25);
  hr->GetYaxis()->SetLabelSize(0.032);


  gStyle->SetOptStat(0);

  hr->Draw();


  gr0= new TGraph(numpoints-1,en,nuebar);
  gr0->SetLineWidth(4);

  gr0->Draw("same");
  //  gr0->Print();

  gr1= new TGraph(numpoints-1,en,nue);
  gr1->SetLineColor(2);
  gr1->SetLineWidth(4);

  gr1->Draw("same");

  gr2= new TGraph(numpoints-1,en,nux);
  gr2->SetLineColor(3);
  gr2->SetLineWidth(4);

   gr2->Draw("same");

   leg = new TLegend(0.45,0.7,0.82,0.85);
   leg->AddEntry(gr2,"#nu_{x} (#nu_{#mu}+#bar{#nu}_{#mu}+#nu_{#tau}+#bar{#nu}_{#tau})","l");
   leg->AddEntry(gr1,"#nu_{e}","l");
   leg->AddEntry(gr0,"#bar{#nu}_{e}","l");
   leg->SetFillColor(10);
   leg->SetTextFont((Font_t) 62);
   leg->SetTextSize(0.05);
   leg->SetBorderSize(0);
                      
   leg->Draw();

   
   TString printfilename = "flux_"+fluxname+".gif";
   canv->Print(printfilename);
   printfilename = "flux_"+fluxname+".eps";
   canv->Print(printfilename);
   printfilename = "flux_"+fluxname+".pdf";
   canv->Print(printfilename);

}
