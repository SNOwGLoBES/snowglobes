// Plot event rates
void plot_event_rates(TString fluxname, TString chanfile, TString exptconfig, Int_t logopt, Int_t titleopt)
{

  using namespace std;


  gStyle->SetOptStat(0);

  // First read in the channels to plot

  const Int_t maxchan = 32;
  TString channame[maxchan];
  Int_t channum[maxchan];
  Int_t numchans;
  TString cpstate[maxchan];
  TString flav[maxchan];
  Double_t num_target_factor[maxchan];


  ifstream chans;
  TString chanfilename = "../channels/channels_"+chanfile+".dat";
  chans.open(chanfilename);
  Int_t j=0;
  while(1) {
    chans >> channame[j]>> channum[j] >>cpstate[j]>>flav[j]>>num_target_factor[j];
    cout << "Channel: "<<j<<" "<<channum[j]<<" "<<channame[j]<<" "<<num_target_factor[j]<<endl;
    j++;
    if (!chans.good()) break;
  }
  numchans=j;

  chans.close();

  cout << "Number of channels: "<< numchans<<endl;

  TString graphfilename = "graphs_"+fluxname+"_"+chanfile+"_"+exptconfig+".root";
  TFile f(graphfilename);

  TGraph* gr;
  TGraph* gr_smeared;

  // Get the maximum  
  TString graphname=fluxname+"_tot_"+exptconfig+"_smeared";
  gr = (TGraph*)f.Get(graphname);

  Double_t x1,y1,x2,y2;
  gr->ComputeRange(x1,y1,x2,y2);
  //  maxy *=1.6;

  Double_t xmin;
  if (exptconfig.CompareTo("ar40kt") == 0 || exptconfig.CompareTo("ar40ktres1")==0|| exptconfig.CompareTo("ar40ktres2") == 0) {
    xmin = 5.;
  } else {
    xmin = 5.;
  }

  Double_t xmax = 100.;
  Double_t ymin = 0.1;
  Double_t ymax;

  //  cout << "Max y "<<y2<<endl;
  if (logopt != 0) {
    ymax = y2*1.6;
  } else {
    ymax = y2*1.1;
  }
  //  Double_t ymax = 1100.;

  TCanvas* canv = new TCanvas("c1"," ",800,700);

  TString plotexpt;
  if (exptconfig == "100kt30pc") {
    plotexpt = "100kt 30%";
  } else {
    if (exptconfig == "100kt15pc") {
      plotexpt = "100kt 15%";
    } else {
      plotexpt = exptconfig;
    }
  }

  TString plot_title;
  if (titleopt == 0) {
    plot_title = " "; 
  } else {
    plot_title = "Flux: "+fluxname+"      Detector: "+plotexpt;
  }
  TH2F *hr = new TH2F("hr",plot_title,2,xmin,xmax,2,ymin,ymax);
  gStyle->SetTitleFontSize(.04);
  hr->SetXTitle("Energy (MeV) ");
  hr->SetYTitle("Events per 0.5 MeV");
  hr->GetXaxis()->SetTitleSize(0.04);
  hr->GetXaxis()->SetTitleOffset(1.3);
  hr->GetXaxis()->SetLabelSize(0.04);
  hr->GetYaxis()->SetTitleSize(0.04);
  hr->GetYaxis()->SetTitleOffset(1.);
  hr->GetYaxis()->SetLabelSize(0.04);

  canv->SetLogy(logopt);

  canv->cd();
  hr->Draw();
  TLegend* leg = new TLegend(0.6,0.7,0.8,0.75);
  leg->SetFillColor(10);
  leg->SetTextFont((Font_t) 62);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);

  graphname=fluxname+"_tot_"+exptconfig;
  gr = (TGraph*)f.Get(graphname);
  gr->SetLineWidth(4);
  gr->SetLineColor(6);

  gr->Print("all");

  graphname=fluxname+"_tot_"+exptconfig+"_smeared";
  gr_smeared = (TGraph*)f.Get(graphname);
  gr_smeared->SetLineWidth(4);
  gr_smeared->SetLineColor(6);
  gr_smeared->SetLineStyle(2);

  gr_smeared->Print("all");

  // Graph with error bars from the smeared graph-- doesn't seem to be a builtin method

  const Int_t maxpoints=1000;
  Int_t numgrpts = gr_smeared->GetN();
  Double_t x[maxpoints],y[maxpoints];
  Double_t ex[maxpoints],ey[maxpoints];
 
  Int_t k;
  for (k=0;k<numgrpts;k++) {
    gr_smeared->GetPoint(k,x[k],y[k]);

    if (y[k]>=0) {
      ey[k] = sqrt(y[k]);
    } else {
      ey[k]=0.;
    }

    ex[k]=0.;

  }

 TGraphErrors* gr_smeared_errors = new TGraphErrors(numgrpts,x,y,ex,ey);

  gr->Draw();
  leg->AddEntry(gr,"Total","l");

  gr_smeared->Draw("same");
  //gr_smeared_errors->Draw();
  leg->AddEntry(gr_smeared,"Total, smeared","l");

  leg->Draw();

  TString printfilename = "tot_rates_"+fluxname+"_"+exptconfig+".pdf";

  canv->Print(printfilename);


//-------------------------------------------------------------------- interaction rates -------------------------------------------------------------------

  TCanvas* canv2 = new TCanvas("c2"," ",800,700);


   canv2->cd();
   hr->Draw();
   canv2->SetLogy(logopt);
   TLegend* leg2 = new TLegend(0.6,0.6,0.88,0.85);
   leg2->SetFillColor(10);
   leg2->SetTextFont((Font_t) 62);
   leg2->SetTextSize(0.04);
   leg2->SetBorderSize(0);

   
   graphname=fluxname+"_tot_"+exptconfig;
   gr = (TGraph*)f.Get(graphname);
   gr->SetLineWidth(3);
   gr->SetLineColor(6);
   gr->Draw("same");
   leg2->AddEntry(gr,"Total","l");
   
   graphname=fluxname+"_es_"+exptconfig;
   gr = (TGraph*)f.Get(graphname);
   gr->SetLineWidth(3);
   gr->SetLineColor(2);
   gr->Draw("same");
   leg2->AddEntry(gr,"ES","l");
     
      if( exptconfig.CompareTo("ar40kt") == 0 || exptconfig.CompareTo("ar40ktres1") == 0 || exptconfig.CompareTo("ar40ktres2")==0){

     graphname=fluxname+"_nue_Ar40_"+exptconfig;
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineWidth(3);
     gr->SetLineColor(3);
     gr->Draw("same");
     leg2->AddEntry(gr,"#nu_{e} ^{40}Ar","l");

     graphname=fluxname+"_nuebar_Ar40_"+exptconfig;
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineWidth(3);
     gr->SetLineColor(4);
     gr->Draw("same");
     leg2->AddEntry(gr,"#bar{#nu}_{e} ^{40}Ar","l");


   } 
   else  {

     if( exptconfig.CompareTo("scint50kt") != 0) {

     graphname=fluxname+"_ibd_"+exptconfig;
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineWidth(3);
     gr->SetLineColor(1);
     gr->Draw("same");
     leg2->AddEntry(gr,"IBD","l");

     graphname=fluxname+"_nue_O16_"+exptconfig;
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineWidth(3);
     gr->SetLineColor(3);
     gr->Draw("same");
     leg2->AddEntry(gr,"#nu_{e}-^{16}O","l");

     graphname=fluxname+"_nuebar_O16_"+exptconfig;
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineWidth(3);
     gr->SetLineColor(4);
     gr->Draw("same");
     leg2->AddEntry(gr,"#bar{#nu}_{e}-^{16}O","l");

     graphname=fluxname+"_ncO16_"+exptconfig;
     gr = (TGraph*)f.Get(graphname);
          gr->SetLineWidth(3);
     gr->SetLineColor(7);
     gr->Draw("same");
     leg2->AddEntry(gr,"NC ^{16}O","l");

     }

     else {

     graphname=fluxname+"_ibd_"+exptconfig;
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineWidth(3);
     gr->SetLineColor(1);
     gr->Draw("same");
     leg2->AddEntry(gr,"IBD","l");

     graphname=fluxname+"_nue_C12_"+exptconfig;
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineWidth(3);
     gr->SetLineColor(3);
     gr->Draw("same");
     leg2->AddEntry(gr,"#nu_{e}-^{12}C","l");

     graphname=fluxname+"_nuebar_C12_"+exptconfig;
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineWidth(3);
     gr->SetLineColor(4);
     gr->Draw("same");
     leg2->AddEntry(gr,"#bar{#nu}_{e}-^{12}C","l");

     graphname=fluxname+"_ncC12_"+exptconfig;
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineWidth(3);
     gr->SetLineColor(7);
     gr->Draw("same");
     leg2->AddEntry(gr,"NC ^{12}C","l");

     }     

   }

   leg2->Draw();

   printfilename = "interaction_rates_"+fluxname+"_"+exptconfig+".pdf";

   canv2->Print(printfilename);


//--------------------------------------------------------------- smeared rates channel contribution -------------------------------------------------------

   TCanvas* canv3 = new TCanvas("c3"," ",800,700);

   canv3->cd();
   hr->Draw();
   canv3->SetLogy(logopt);

   TLegend* leg3 = new TLegend(0.6,0.6,0.88,0.85);
   leg3->SetFillColor(10);
   leg3->SetTextFont((Font_t) 62);
   leg3->SetTextSize(0.04);
   leg3->SetBorderSize(0);

   
   graphname=fluxname+"_tot_"+exptconfig+"_smeared";
   gr = (TGraph*)f.Get(graphname);
   gr->SetLineWidth(3);
   gr->SetLineColor(6);
   gr->SetLineStyle(2);
   gr->Draw("same");
   leg3->AddEntry(gr,"Total","l");
   
   graphname=fluxname+"_es_"+exptconfig+"_smeared";
   gr = (TGraph*)f.Get(graphname);
   gr->SetLineWidth(3);
   gr->SetLineColor(2);
   gr->SetLineStyle(2);
   gr->Draw("same");
   leg3->AddEntry(gr,"ES","l");
     
      if( exptconfig.CompareTo("ar40kt") == 0 || exptconfig.CompareTo("ar40ktres1") == 0 || exptconfig.CompareTo("ar40ktres2")==0){

     graphname=fluxname+"_nue_Ar40_"+exptconfig+"_smeared";
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineStyle(2);
     gr->SetLineWidth(3);
     gr->SetLineColor(3);
     gr->Draw("same");
     leg3->AddEntry(gr,"#nu_{e} ^{40}Ar","l");

     graphname=fluxname+"_nuebar_Ar40_"+exptconfig+"_smeared";
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineStyle(2);
     gr->SetLineWidth(3);
     gr->SetLineColor(4);
     gr->Draw("same");
     leg3->AddEntry(gr,"#bar{#nu}_{e} ^{40}Ar","l");


   }

   if( exptconfig.CompareTo("wc100kt30prct") == 0 || exptconfig.CompareTo("wc100kt15prct") == 0) {

     graphname=fluxname+"_ibd_"+exptconfig+"_smeared";
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineStyle(2);
     gr->SetLineWidth(3);
     gr->SetLineColor(1);
     gr->Draw("same");
     leg3->AddEntry(gr,"IBD","l");

     graphname=fluxname+"_nue_O16_"+exptconfig+"_smeared";
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineStyle(2);
     gr->SetLineWidth(3);
     gr->SetLineColor(3);
     gr->Draw("same");
     leg3->AddEntry(gr,"#nu_{e}-^{16}O","l");

     graphname=fluxname+"_nuebar_O16_"+exptconfig+"_smeared";
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineStyle(2);
     gr->SetLineWidth(3);
     gr->SetLineColor(4);
     gr->Draw("same");
     leg3->AddEntry(gr,"#bar{#nu}_{e}-^{16}O","l");

     graphname=fluxname+"_ncO16_"+exptconfig+"_smeared";
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineStyle(2);
     gr->SetLineWidth(3);
     gr->SetLineColor(7);
     gr->Draw("same");
     leg3->AddEntry(gr,"NC ^{16}O","l");
   }

   if( exptconfig.CompareTo("scint50kt") == 0) {

     graphname=fluxname+"_ibd_"+exptconfig+"_smeared";
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineStyle(2);
     gr->SetLineWidth(3);
     gr->SetLineColor(1);
     gr->Draw("same");
     leg3->AddEntry(gr,"IBD","l");

     graphname=fluxname+"_nue_C12_"+exptconfig+"_smeared";
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineStyle(2);
     gr->SetLineWidth(3);
     gr->SetLineColor(3);
     gr->Draw("same");
     leg3->AddEntry(gr,"#nu_{e}-^{12}C","l");

     graphname=fluxname+"_nuebar_C12_"+exptconfig+"_smeared";
     gr = (TGraph*)f.Get(graphname);  
     gr->SetLineStyle(2);
     gr->SetLineWidth(3);
     gr->SetLineColor(4);
     gr->Draw("same");
     leg3->AddEntry(gr,"#bar{#nu}_{e}-^{12}C","l");

     graphname=fluxname+"_ncC12_"+exptconfig+"_smeared";
     gr = (TGraph*)f.Get(graphname);
     gr->SetLineStyle(2);
     gr->SetLineWidth(3);
     gr->SetLineColor(7);
     gr->Draw("same");
     leg3->AddEntry(gr,"NC ^{12}C","l");
   }


   leg3->Draw();

   printfilename = "smeared_rates_"+fluxname+"_"+exptconfig+".pdf";

   canv3->Print(printfilename);

//-------------------------------------------------------------------true flavor content ES-----------------------------------------------------------------

  TCanvas* canv4 = new TCanvas("c4"," ",800,700);

   canv4->cd();
   hr->Draw();
   canv4->SetLogy(logopt);

   TLegend* leg4 = new TLegend(0.6,0.6,0.88,0.85);
   leg4->SetFillColor(10);
   leg4->SetTextFont((Font_t) 62);
   leg4->SetTextSize(0.04);
   leg4->SetBorderSize(0);

   
   graphname=fluxname+"_es_"+exptconfig+"_smeared";
   gr = (TGraph*)f.Get(graphname);
   gr->SetLineWidth(3);
   gr->SetLineColor(2);
   gr->SetLineStyle(2);
   gr->Draw("same");
   leg4->AddEntry(gr,"ES total","l");
   
   graphname=fluxname+"_nue_es_"+exptconfig+"_smeared";
   gr = (TGraph*)f.Get(graphname);
   gr->SetLineWidth(3);
   gr->SetLineColor(3);
   gr->SetLineStyle(2);
   gr->Draw("same");
   leg4->AddEntry(gr,"#nu_{e}","l");

   graphname=fluxname+"_nuebar_es_"+exptconfig+"_smeared";
   gr = (TGraph*)f.Get(graphname);
   gr->SetLineWidth(3);
   gr->SetLineColor(4);
   gr->SetLineStyle(2);
   gr->Draw("same");
   leg4->AddEntry(gr,"#bar{#nu}_{e}","l");

   graphname=fluxname+"_nux_es_"+exptconfig+"_smeared";
   gr = (TGraph*)f.Get(graphname);
   gr->SetLineWidth(3);
   gr->SetLineColor(7);
   gr->SetLineStyle(2);
   gr->Draw("same");
   leg4->AddEntry(gr,"#nu_{x}","l");

   leg4->Draw();

   printfilename = "true_flavor_content_ES_"+fluxname+"_"+exptconfig+".pdf";

   canv4->Print(printfilename);

}
