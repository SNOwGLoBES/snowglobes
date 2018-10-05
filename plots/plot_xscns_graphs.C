void plot_xscns_graphs(TString chanfilename)
{
  // Plot graphs, having made them first with make_xscn_graphs.C
  // This is specific for water, scint, argon channels files

  // First read in the channels to plot

  const Int_t maxchan = 32;
  TString channame[maxchan];
  Int_t channum[maxchan];
  Int_t numchans;
  TString cpstate[maxchan];
  TString flav[maxchan];
  Double_t num_target_factor[maxchan];


  ifstream chans;
  TString inchanfilename = "../channels/channels_"+chanfilename+".dat";

  chans.open(inchanfilename);
  Int_t j=0;
  while(1) {
    chans >> channame[j]>> channum[j] ;
    if (!chans.good()) break; // check in the middle of the line read to avoid endline shenanigans
    chans >>cpstate[j]>>flav[j]>>num_target_factor[j];
    cout << "Channel: "<<j<<" "<<channum[j]<<" "<<channame[j]<<" "<<num_target_factor[j]<<endl;
    j++;
  }
  numchans=j;

  chans.close();

  cout << "Number of channels: "<< numchans<<endl;

  TString graphfilename = "graphs_"+chanfilename+".root";
  TFile f(graphfilename);

  TObjArray xscngraphlist(0);
  TGraph* gr;

  Int_t i;
  for (i=0;i<numchans;i++) {
    gr = (TGraph*)f.Get(channame[channum[i]]);
    xscngraphlist.Add(gr);
  }

  TCanvas* canv = new TCanvas("c1"," ",800,700);
  canv->SetLogy(1);
  canv->SetGridx(1);
  canv->SetGridy(1);

  //  TString plot_title = "Cross-sections "+chanfilename;
  TString plot_title = " ";

  gStyle->SetTitleFontSize(.04);

  Double_t xmin = 5.;
  Double_t xmax = 100.;
  Double_t ymin = 0.0000001;

  // Because graph method not working
  Int_t p;
  Double_t ymax=100.;
  //  for (p=0;p<numpoints;p++) {
  // if (nuescan[p]>ymax) {
  //   ymax = nux[p];
  // }
  //}

  ymax*=1.1;

  TH2F *hr = new TH2F("hr",plot_title,2,xmin,xmax,2,ymin,ymax);
  hr->SetXTitle(" Neutrino Energy (MeV) ");
  hr->SetYTitle(" Cross section (10^{-38} cm^{2})");
  hr->GetXaxis()->SetTitleSize(0.04);
  hr->GetXaxis()->SetTitleOffset(0.9);
  hr->GetXaxis()->SetLabelSize(0.04);
  hr->GetYaxis()->SetTitleSize(0.04);
  hr->GetYaxis()->SetTitleOffset(1.2);
  hr->GetYaxis()->SetLabelSize(0.04);

  gStyle->SetOptStat(0);

  hr->Draw();


   TLegend * leg = new TLegend(0.18,0.62,0.4,0.88);

   for (i=0;i<numchans;i++) {

    gr = (TGraph*)xscngraphlist.At(i);

    //    gr->Print();
    gr->Draw("same");
    
    if( chanfilename.Contains("water") == 1){

      if (i<3 || i==7 || i==8) {
	//	leg->AddEntry(gr,channame[channum[i]],"l");
	if (i==0) leg->AddEntry(gr,"IBD","l");
	if (i==1) leg->AddEntry(gr,"#nu_{e} e","l");
	if (i==2) leg->AddEntry(gr,"#bar{#nu}_{e} e","l");
	if (i==7) leg->AddEntry(gr,"#nu_{e}-^{16}O","l");
	if (i==8) leg->AddEntry(gr,"#bar{#nu}_{e}-^{16}O","l");

      } else {
	if (i==3) {
	  leg->AddEntry(gr,"#nu_{x}-e","l");
	}
	if (i==4) {
	  leg->AddEntry(gr,"#bar{#nu}_{x}-e","l");
	}
	if (i==9) {
	  leg->AddEntry(gr,"NC ^{16}O","l");
	}
      }
    } else  {
      if( chanfilename.Contains("argon") == 1){


      // argon
      if (i<2|| i==6 || i==7 ) {
	//	leg->AddEntry(gr,channame[channum[i]],"l");
	if (i==0) leg->AddEntry(gr,"#nu_{e} e","l");
	if (i==1) leg->AddEntry(gr,"#bar{#nu}_{e} e","l");

	if (i==6) leg->AddEntry(gr,"#nu_{e} ^{40}Ar","l");
	if (i==7) leg->AddEntry(gr,"#bar{#nu}_{e} ^{40}Ar","l");

      } else {
	if (i==2) {
	  leg->AddEntry(gr,"#nu_{x} e","l");
	}
	if (i==3) {
	  leg->AddEntry(gr,"#bar{#nu}_{x} e","l");
	}
	if (i==8) {
	  leg->AddEntry(gr,"coh #nu-A","l");

	}
        if (i==9) {
          leg->AddEntry(gr,"NC inl #nu ^{40}Ar", "l");
          // if we put coherent scattering in, we need to move the legend
          leg->SetX1(0.12); leg->SetX2(0.60); leg->SetY1(0.70); leg->SetY2(0.89);
          leg->SetNColumns(2);
        }

      }
      } 
      else {
      // scint
	if (i==0) leg->AddEntry(gr,"IBD","l");
	if (i==1) leg->AddEntry(gr,"#nu_{e}-e","l");
	if (i==2) leg->AddEntry(gr,"#bar{#nu}_{e}-e","l");
	if (i==3) leg->AddEntry(gr,"#nu_{x}-e","l");
	if (i==4) leg->AddEntry(gr,"#bar{#nu}_{x}-e","l");

	if (i==7) leg->AddEntry(gr,"#nu_{e}-^{12}C","l");
	if (i==8) leg->AddEntry(gr,"#bar{#nu}_{e}-^{12}C","l");
	if (i==9) leg->AddEntry(gr,"NC ^{12}C","l");
        if (i==15) leg->AddEntry(gr,"#nu_{e}-^{13}C", "l");
        if (i==16) leg->AddEntry(gr,"NC ^{13}C", "l");
        leg->SetX1(0.12); leg->SetX2(0.60); leg->SetY1(0.70); leg->SetY2(0.89);
        leg->SetNColumns(2);

      }
      
    } 
   
  //  gr = (TGraph*)xscngraphlist.At(11);
  //gr->Draw();
   }

   leg->SetFillColor(10);
   leg->SetTextFont((Font_t) 62);
   leg->SetTextSize(0.037);
   leg->SetBorderSize(0);
                      
   leg->Draw();

  TString printfilename = "xscns_"+chanfilename+".gif";
  canv->Print(printfilename);

  printfilename = "xscns_"+chanfilename+".eps";
  canv->Print(printfilename);

  printfilename = "xscns_"+chanfilename+".pdf";
  canv->Print(printfilename);

  //  f.Close();
}
