void plot_xscns_lead()
{
  // Plot lead cross-sections; graphs must first have been made by make_xscn_graphs.C

  TObjArray xscngraphlist(0);

  TString graphfilename4 = "graphs_lead.root";
  TFile f4(graphfilename4);


  TGraph* gr;

  TCanvas* canv = new TCanvas("c1"," ",800,700);
  canv->SetLogy(1);

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
  hr->SetYTitle(" Cross-section (10^{-38} cm^{2})");
  hr->GetXaxis()->SetTitleSize(0.03);
  hr->GetXaxis()->SetTitleOffset(0.85);
  hr->GetXaxis()->SetLabelSize(0.03);
  hr->GetYaxis()->SetTitleSize(0.03);
  hr->GetYaxis()->SetTitleOffset(1.7);
  hr->GetYaxis()->SetLabelSize(0.03);

  gStyle->SetOptStat(0);

  hr->Draw();


   leg = new TLegend(0.12,0.7,0.3,0.89);
   leg2 = new TLegend(0.3,0.75,0.45,0.88);


   gr = (TGraph*)f4.Get("nue_Pb208_1n");
   gr->SetLineColor(3);
   gr->SetLineStyle(1);
   gr->Draw();

   // Assume all graphs have same x axis, I think true
   Double_t xval[1000],ytot[1000],y1ntot[1000],y2ntot[1000];
   Int_t numpoint = gr->GetN();

   Int_t i;
   for (i=0;i<numpoint;i++) 
     {
       Double_t x,y;
       Int_t ret=gr->GetPoint(i,x,y);
       xval[i]=x;
       ytot[i]=y;
       y1ntot[i]=y;
       //       cout <<i<<" "<< x<<" "<<y<<" "<<endl;
     }

   leg2->AddEntry(gr,"#nu_{e}-^{208}Pb 1n","l");

   gr = (TGraph*)f4.Get("nue_Pb208_2n");
   gr->SetLineColor(3);
   gr->SetLineStyle(2);
   gr->Draw();

   for (i=0;i<numpoint;i++) 
     {
       Double_t x,y;
       Int_t ret=gr->GetPoint(i,x,y);
       xval[i]=x;
       ytot[i]+=y;
       y2ntot[i]=y;
     }

   leg2->AddEntry(gr,"#nu_{e}-^{208}Pb 2n","l");

   gr = (TGraph*)f4.Get("nc_nue_Pb208_1n");
   gr->SetLineColor(7);
   gr->SetLineStyle(3);
   gr->Draw();
   for (i=0;i<numpoint;i++) 
     {
       Double_t x,y;
       Int_t ret=gr->GetPoint(i,x,y);
       xval[i]=x;
       ytot[i]+=y;
       y1ntot[i]+=y;

     }

   leg2->AddEntry(gr,"NC-#nu ^{208}Pb 1n","l");

   gr = (TGraph*)f4.Get("nc_nue_Pb208_2n");
   gr->SetLineColor(7);
   gr->SetLineStyle(4);
   gr->Draw();
   for (i=0;i<numpoint;i++) 
     {
       Double_t x,y;
       Int_t ret=gr->GetPoint(i,x,y);
       xval[i]=x;
       ytot[i]+=y;
       y2ntot[i]+=y;
     }
   leg2->AddEntry(gr,"NC-#nu ^{208}Pb 2n","l");

   gr = (TGraph*)f4.Get("nc_nuebar_Pb208_1n");
   gr->SetLineColor(7);
   gr->SetLineStyle(5);
   gr->Draw();
   for (i=0;i<numpoint;i++) 
     {
       Double_t x,y;
       Int_t ret=gr->GetPoint(i,x,y);
       xval[i]=x;
       ytot[i]+=y;
       y1ntot[i]+=y;
     }
   leg2->AddEntry(gr,"NC-#bar{#nu} ^{208}Pb 1n","l");

   gr = (TGraph*)f4.Get("nc_nuebar_Pb208_2n");
   gr->SetLineColor(7);
   gr->SetLineStyle(6);
   gr->Draw();
   for (i=0;i<numpoint;i++) 
     {
       Double_t x,y;
       Int_t ret=gr->GetPoint(i,x,y);
       xval[i]=x;
       ytot[i]+=y;
       y2ntot[i]+=y;

     }
   leg2->AddEntry(gr,"NC-#bar{#nu} ^{208}Pb 2n","l");


   leg->SetFillColor(10);
   leg->SetTextFont((Font_t) 62);
   leg->SetTextSize(0.02);
   leg->SetBorderSize(0);
                      
   leg->Draw();



   leg2->SetFillColor(10);
   leg2->SetTextFont((Font_t) 62);
   leg2->SetTextSize(0.02);
   leg2->SetBorderSize(0);
                      
   leg2->Draw();


   TString printfilename = "xscns_lead.gif";
   canv->Print(printfilename);

   printfilename = "xscns_lead.eps";
   canv->Print(printfilename);


   printfilename = "xscns_lead.pdf";
   canv->Print(printfilename);

   //   The sums: 
   TGraph* grtot = new TGraph(numpoint,xval,ytot);

   TGraph* grtot1n = new TGraph(numpoint,xval,y1ntot);
   TGraph* grtot2n = new TGraph(numpoint,xval,y2ntot);
   grtot->Print();

   ymin = 0.;
   ymax = 4.5;
   TH2F *hr2 = new TH2F("hr",plot_title,2,xmin,xmax,2,ymin,ymax);
   hr2->SetXTitle(" Neutrino Energy (MeV) ");
   hr2->SetYTitle(" Cross-section (10^{-38} cm^{2})");
   hr2->GetXaxis()->SetTitleSize(0.03);
   hr2->GetXaxis()->SetTitleOffset(0.85);
   hr2->GetXaxis()->SetLabelSize(0.03);
   hr2->GetYaxis()->SetTitleSize(0.03);
   hr2->GetYaxis()->SetTitleOffset(1.7);
   hr2->GetYaxis()->SetLabelSize(0.03);

  TCanvas* canv2 = new TCanvas("c2"," ",800,700);
  canv2->SetLogy(0);
  canv2->cd();
  hr2->Draw();
  grtot->Draw("same");
  grtot->SetLineColor(1);
  grtot->SetLineStyle(1);
  grtot->SetLineWidth(2);


  grtot1n->Draw("same");
  grtot1n->SetLineColor(2);
  grtot1n->SetLineStyle(2);
  grtot1n->SetLineWidth(2);

  grtot2n->Draw("same");
  grtot2n->SetLineColor(3);
  grtot2n->SetLineStyle(3);
  grtot2n->SetLineWidth(2);



  printfilename = "xscns_lead_tot.gif";
  canv2->Print(printfilename);

  printfilename = "xscns_lead_tot.eps";
  canv2->Print(printfilename);

  printfilename = "xscns_lead_tot.pdf";
  canv2->Print(printfilename);

  //  f.Close();
}
