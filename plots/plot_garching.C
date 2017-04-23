// Pretty Garching plots of fluxes 

void plot_garching()
{
  gStyle->SetPadTopMargin(0.01);

  // Get the events per bin from the file made by garching_time_plot.C

  TFile f("eventsperbin.root");

  TH1D* eventsperbin = (TH1D*)f.Get("eventsperbin");

  TH1D* rateinbin = (TH1D*)f.Get("rateinbin");

  // Return to home directory
  gROOT->cd();
  TH1D* eventsperbin2 = (TH1D*)eventsperbin->Clone("eventsperbin2");
  TH1D* rateinbin2 = (TH1D*)rateinbin->Clone("rateinbin2");

  f.Close();

  TString quantity = "luminosity";
  gStyle->SetOptStat(0);

  ifstream infile;
  // Two timescales, 6 flavors, three types of quantities (lum, eavg, alpha)
  // One graph array per quantity, indexed by flavor.  Glom together both timescales, read from two files.

  const Int_t numflavor = 3;

  TString flavname[numflavor] = {"nue","nuebar","numu"};

  TGraph** lumgraphs = new TGraph*[numflavor];
  TGraph** avgengraphs = new TGraph*[numflavor];
  TGraph** alphagraphs = new TGraph*[numflavor];
  
  const Int_t maxpoints=2000;
  Double_t time[maxpoints];
  Double_t yval[maxpoints];

  Int_t i=0;
  Int_t j=0;

  quantity = "luminosity";

  for (i=0;i<numflavor;i++){

    ifstream in;

    TString infilename2 = quantity+"_"+flavname[i]+"_garching.txt";

    cout << infilename2<<endl;

    in.open(infilename2);
    j=0;
    Int_t numpoints;
    
    while(1) {
      in >> time[j]>>yval[j];
      //      cout << "time "<<time[j]<<" "<<yval[j]<<endl;
      
      // For log x scale
      time[j]+=0.02;
      if (!in.good()) break;
      j++;

    } // End of file reading loop

    numpoints = j;
    in.close();
    
    lumgraphs[i]= new TGraph(numpoints,time,yval);
    lumgraphs[i]->SetLineWidth(3);


  } // End of loop over flavors


  lumgraphs[1]->SetLineColor(2);
  lumgraphs[1]->SetLineStyle(2);

  lumgraphs[2]->SetLineColor(3);
  lumgraphs[2]->SetLineStyle(3);

  lumgraphs[2]->Print();

  quantity = "avgen";

  for (i=0;i<numflavor;i++){
    cout << "i "<<i<< endl;

    ifstream in;

    TString infilename2 = quantity+"_"+flavname[i]+"_garching.txt";

    cout << infilename2<<endl;

    in.open(infilename2);
    j=0;
    Int_t numpoints;
    
    while(1) {
      in >> time[j]>>yval[j];
      //      cout << "time "<<time[j]<<" "<<yval[j]<<endl;
      
      // For log x scale
      time[j]+=0.02;
      if (!in.good()) break;
      j++;

    } // End of file reading loop

    numpoints = j;
    in.close();
    
    avgengraphs[i]= new TGraph(numpoints,time,yval);
    avgengraphs[i]->SetLineWidth(3);


  } // End of loop over flavors


  avgengraphs[1]->SetLineColor(2);
  avgengraphs[1]->SetLineStyle(2);

  avgengraphs[2]->SetLineColor(3);
  avgengraphs[2]->SetLineStyle(3);

  //////

  quantity = "alpha";

  for (i=0;i<numflavor;i++){

    ifstream in;

    TString infilename2 = quantity+"_"+flavname[i]+"_garching.txt";

    cout << infilename2<<endl;

    in.open(infilename2);
    j=0;
    Int_t numpoints;
    
    while(1) {
      in >> time[j]>>yval[j];
      //      cout << "time "<<time[j]<<" "<<yval[j]<<endl;
      
      // For log x scale
      time[j]+=0.02;
      if (!in.good()) break;
      j++;

    } // End of file reading loop

    numpoints = j;
    in.close();
    
    alphagraphs[i]= new TGraph(numpoints,time,yval);
    alphagraphs[i]->SetLineWidth(3);


  } // End of loop over flavors


  alphagraphs[1]->SetLineColor(2);
  alphagraphs[1]->SetLineStyle(2);

  alphagraphs[2]->SetLineColor(3);
  alphagraphs[2]->SetLineStyle(3);


  ///////////////////////////////


  TCanvas* canv4 = new TCanvas("c4"," ",800,700);

  canv4->Divide(1,4,0.001,0.001);

  TLegend* leg = new TLegend(0.75,0.5,0.9,0.88);
  leg->AddEntry(lumgraphs[0],"#nu_{e}","l");
  leg->AddEntry(lumgraphs[1],"#bar{#nu}_{e}","l");
  leg->AddEntry(lumgraphs[2],"#nu_{x}","l");

  leg->SetFillColor(10);
  leg->SetTextFont((Font_t) 62);
  leg->SetTextSize(0.15);
  leg->SetBorderSize(0);
           


  Double_t xmin = 0.005;
  Double_t xmax = 9.;

  TString title = " ";

  Double_t ymin = 0.01;
  Double_t ymax = 200.;

  TH2F *hr2 = new TH2F("hr2","",2,xmin,xmax,2,ymin,ymax);
  //  hr2->SetXTitle(" Time (seconds) ");
  title ="L (10^{52} ergs/s)   ";

  hr2->SetYTitle(title);
  hr2->GetXaxis()->SetTitleSize(0.1);
  hr2->GetXaxis()->SetTitleOffset(1.1);
  hr2->GetXaxis()->SetLabelSize(0.1);
  hr2->GetYaxis()->SetTitleSize(0.1);
  hr2->GetYaxis()->SetTitleOffset(.35);
  hr2->GetYaxis()->SetLabelSize(0.1);


  ymin = 5.;
  ymax = 15.;

  TH2F *hr3 = new TH2F("hr3","",2,xmin,xmax,2,ymin,ymax);
  //  hr3->SetXTitle(" Time (seconds) ");

  title = "<E> (MeV)   ";
  hr3->SetYTitle(title);
  hr3->GetXaxis()->SetTitleSize(0.1);
  hr3->GetXaxis()->SetTitleOffset(1.1);
  //  hr3->GetXaxis()->SetLabelSize(0.1);
  hr3->GetXaxis()->SetLabelSize(0.);
  hr3->GetYaxis()->SetTitleSize(0.1);
  hr3->GetYaxis()->SetTitleOffset(.3);
  hr3->GetYaxis()->SetLabelSize(0.1);

  ymin = 2.;
  ymax = 5;

  TH2F *hr4 = new TH2F("hr4","",2,xmin,xmax,2,ymin,ymax);
  //  hr4->SetXTitle(" Time (seconds) ");
  title = "Alpha     ";
  hr4->SetYTitle(title);
  hr4->GetXaxis()->SetTitleSize(0.1);
  hr4->GetXaxis()->SetTitleOffset(1.1);
  //  hr4->GetXaxis()->SetLabelSize(0.1);
  hr4->GetXaxis()->SetLabelSize(0.);
  hr4->GetYaxis()->SetTitleSize(0.1);
  hr4->GetYaxis()->SetTitleOffset(.3);
  hr4->GetYaxis()->SetLabelSize(0.1);


  ymin = 0.1;
  ymax = 75.;
  //  ymax = 1000.;
  TH2F *hr5 = new TH2F("hr5","",2,xmin,xmax,2,ymin,ymax);
  //  hr5->SetXTitle(" Time (seconds) ");
  title = "Events per bin";
  hr5->SetYTitle(title);
  hr5->GetXaxis()->SetTitleSize(0.1);
  hr5->GetXaxis()->SetTitleOffset(1.1);
  hr5->GetXaxis()->SetLabelSize(0.1);
  hr5->GetYaxis()->SetTitleSize(0.1);
  hr5->GetYaxis()->SetTitleOffset(.3);
  hr5->GetYaxis()->SetLabelSize(0.1);



  ////

  // Plot the stuff

  gStyle->SetOptStat(0);

  canv4->cd(1);
  //  TPad *current_pad = (TPad*)gROOT->GetSelectedPad();
  //current_pad->SetTopMargin(0.);
  //current_pad->SetLogx(1.);
  //current_pad->SetLogy(1.);

  gPad->SetLogx(1.);
  gPad->SetLogy(1.);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0.01);
  hr2->Draw();
  leg->Draw("p");

  lumgraphs[0]->Draw("same");
  lumgraphs[1]->Draw("same");
  lumgraphs[2]->Draw("same");

    Double_t x1,x2,y1,y2;
    x1 = 0.02;
    x2 = 0.02;

    y1 = 0;
    y2 = 200.;

  TLine* bouncetime1 = new TLine(x1,y1,x2,y2);
  bouncetime1->SetLineWidth(2);
  bouncetime1->SetLineStyle(2);
  bouncetime1->SetLineColor(1);
  bouncetime1->Draw("same");

  x1= 0.02+0.05;
  x2= 0.02+0.05;


  TLine* accretiontime1 = new TLine(x1,y1,x2,y2);
  accretiontime1->SetLineWidth(2);
  accretiontime1->SetLineStyle(3);
  accretiontime1->SetLineColor(1);
  accretiontime1->Draw("same");

  x1= 0.02+0.2;
  x2= 0.02+0.2;


  TLine* coolingtime1 = new TLine(x1,y1,x2,y2);
  coolingtime1->SetLineWidth(2);
  coolingtime1->SetLineStyle(3);
  coolingtime1->SetLineColor(1);
  coolingtime1->Draw("same");

  TText tinfall;
  tinfall.SetTextAlign(32);
  tinfall.SetTextAngle(0);
  tinfall.SetTextSize(0.08);
  tinfall.DrawText(0.012,100,"Infall");

  TText tneutronization;
  tneutronization.SetTextAlign(32);
  tneutronization.SetTextAngle(0);
  tneutronization.SetTextSize(0.08);
  tneutronization.DrawText(0.06,100,"Neutronization");

  TText taccretion;
  taccretion.SetTextAlign(32);
  taccretion.SetTextAngle(0);
  taccretion.SetTextSize(0.08);
  taccretion.DrawText(0.17,100,"Accretion");

  TText tcooling;
  tcooling.SetTextAlign(32);
  tcooling.SetTextAngle(0);
  tcooling.SetTextSize(0.08);
  tcooling.DrawText(.7,100,"Cooling");


  // Replace the y-labels, otherwise they get cut off
  //http://root.cern.ch/root/roottalk/roottalk01/2812.html
  // Ugh, Root is tedious

  hr2->GetYaxis()->SetLabelOffset(99);
  const Int_t ny=3;
    char *hr2ylab[ny] = {"0.1","1","10"};
  //char *hr2ylab[ny] = {"A","B","C"};
  Float_t hr2ylabpos[ny] = {0.1,1,10};

  Float_t x,y;
  //x= gPad->GetUxmin() - 10.0*hr2->GetXaxis()->GetBinWidth(1);
  x= 0.0049;
  TText t;
  t.SetTextAlign(32);
  t.SetTextAngle(0);
  t.SetTextSize(0.11);
  for (i=0;i<ny;i++) {
  //   y = hr2->GetYaxis()->GetBinCenter(i+1);
    cout << "i "<<i<<" "<<x<<" "<<y<<endl;
    y=hr2ylabpos[i];
    t.DrawText(x,y,hr2ylab[i]);
  }

  /////

  canv4->cd(2);
  //  TPad *current_pad = (TPad*)gROOT->GetSelectedPad();
  //current_pad->SetTopMargin(0.);

  gPad->SetLogx(1.);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0.01);
  hr3->Draw();

  avgengraphs[0]->Draw("same");
  avgengraphs[1]->Draw("same");
  avgengraphs[2]->Draw("same");

  x1 = 0.02;
  x2 = 0.02;

  y1 = 0;
  y2 = 100.;

  TLine* bouncetime2 = new TLine(x1,y1,x2,y2);
  bouncetime2->SetLineWidth(2);
  bouncetime2->SetLineStyle(2);
  bouncetime2->SetLineColor(1);
  bouncetime2->Draw("same");

  x1= 0.02+0.05;
  x2= 0.02+0.05;


  TLine* accretiontime2 = new TLine(x1,y1,x2,y2);
  accretiontime2->SetLineWidth(2);
  accretiontime2->SetLineStyle(3);
  accretiontime2->SetLineColor(1);
  accretiontime2->Draw("same");

  x1= 0.02+0.2;
  x2= 0.02+0.2;


  TLine* coolingtime2 = new TLine(x1,y1,x2,y2);
  coolingtime2->SetLineWidth(2);
  coolingtime2->SetLineStyle(3);
  coolingtime2->SetLineColor(1);
  coolingtime2->Draw("same");


  //////


  hr3->GetYaxis()->SetLabelOffset(99);
  const Int_t ny3=5;
  char *hr3ylab[ny3] = {"6","8","10","12","14"};
  Float_t hr3ylabpos[ny3] = {6.,8.,10.,12.,14.};

  t.SetTextAlign(32);
  t.SetTextAngle(0);
  t.SetTextSize(0.11);
  for (i=0;i<ny3;i++) {
  //   y = hr2->GetYaxis()->GetBinCenter(i+1);
    cout << "i "<<i<<" "<<x<<" "<<y<<endl;
    y=hr3ylabpos[i];
    t.DrawText(x,y,hr3ylab[i]);
  }



  canv4->cd(3);
  gPad->SetLogx(1.);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0.01);
  hr4->Draw();


  alphagraphs[0]->Draw("same");
  alphagraphs[1]->Draw("same");
  alphagraphs[2]->Draw("same");

    x1 = 0.02;
    x2 = 0.02;

    y1 = 0;
    y2 = 100.;

  TLine* bouncetime3 = new TLine(x1,y1,x2,y2);
  bouncetime3->SetLineWidth(2);
  bouncetime3->SetLineStyle(2);
  bouncetime3->SetLineColor(1);
  bouncetime3->Draw("same");

  x1= 0.02+0.05;
  x2= 0.02+0.05;


  TLine* accretiontime3 = new TLine(x1,y1,x2,y2);
  accretiontime3->SetLineWidth(2);
  accretiontime3->SetLineStyle(3);
  accretiontime3->SetLineColor(1);
  accretiontime3->Draw("same");

  x1= 0.02+0.2;
  x2= 0.02+0.2;


  TLine* coolingtime3 = new TLine(x1,y1,x2,y2);
  coolingtime3->SetLineWidth(2);
  coolingtime3->SetLineStyle(3);
  coolingtime3->SetLineColor(1);
  coolingtime3->Draw("same");


  hr4->GetYaxis()->SetLabelOffset(99);
  const Int_t ny4=5;
  char *hr4ylab[ny4] = {"2.5","3","3.5","4","4.5"};
  Float_t hr4ylabpos[ny4] = {2.5,3,3.5,4,4.5};

  t.SetTextAlign(32);
  t.SetTextAngle(0);
  t.SetTextSize(0.11);
  for (i=0;i<ny4;i++) {
  //   y = hr2->GetYaxis()->GetBinCenter(i+1);
    cout << "i "<<i<<" "<<x<<" "<<y<<endl;
    y=hr4ylabpos[i];
    t.DrawText(x,y,hr4ylab[i]);
  }





  canv4->cd(4);
  gPad->SetLogx(1.);
  gPad->SetTopMargin(0.015);
  gPad->SetBottomMargin(0.15);
  //  gPad->SetLogy(1.);
  hr5->Draw();
  hr5->SetXTitle(" Time (seconds) ");
  hr5->GetXaxis()->SetTitleSize(0.1);
  hr5->GetXaxis()->SetTitleOffset(0.5);
  hr5->GetXaxis()->SetLabelSize(0.1);

  eventsperbin2->Draw("pesame");

  // Note first (0) and last bins are under and overflow
  cout << " Events in infall: " <<eventsperbin2->Integral(1,9,"")<<endl;
  cout << " Events in neutronization: " <<eventsperbin2->Integral(10,17,"")<<endl;
  cout << " Events in accretion: " <<eventsperbin2->Integral(18,25,"")<<endl;
  cout << " Events in cooling: " <<eventsperbin2->Integral(26,eventsperbin2->GetNbinsX(),"")<<endl;
  cout << "Total events: "<<eventsperbin2->Integral()<<endl;
  cout << "Nbins "<< eventsperbin2->GetNbinsX()<<endl;
  eventsperbin2->Print("all");

  //  rateinbin2->Print("all");
  //rateinbin2->Draw("pesame");
  
  //  canv4->Update();

    x1 = 0.02;
    x2 = 0.02;

    y1 = 0;
    y2 = 100.;

  TLine* bouncetime4 = new TLine(x1,y1,x2,y2);
  bouncetime4->SetLineWidth(2);
  bouncetime4->SetLineStyle(2);
  bouncetime4->SetLineColor(1);
  bouncetime4->Draw("same");

  x1= 0.02+0.05;
  x2= 0.02+0.05;


  TLine* accretiontime4 = new TLine(x1,y1,x2,y2);
  accretiontime4->SetLineWidth(2);
  accretiontime4->SetLineStyle(3);
  accretiontime4->SetLineColor(1);
  accretiontime4->Draw("same");

  x1= 0.02+0.2;
  x2= 0.02+0.2;


  TLine* coolingtime4 = new TLine(x1,y1,x2,y2);
  coolingtime4->SetLineWidth(2);
  coolingtime4->SetLineStyle(3);
  coolingtime4->SetLineColor(1);
  coolingtime4->Draw("same");


  hr5->GetYaxis()->SetLabelOffset(99);
  const Int_t ny5=7;
  char *hr5ylab[ny5] = {"10","20","30","40","50","60","70"};
  Float_t hr5ylabpos[ny5] = {10.,20.,30.,40.,50.,60.,70.};

  t.SetTextAlign(32);
  t.SetTextAngle(0);
  t.SetTextSize(0.11);
  for (i=0;i<ny5;i++) {
  //   y = hr2->GetYaxis()->GetBinCenter(i+1);
    //    cout << "i "<<i<<" "<<x<<" "<<y<<endl;
    y=hr5ylabpos[i];
    t.DrawText(x,y,hr5ylab[i]);
  }


}

