// Example of time-dependent plotting
// Plot Garching rate as a function of time
// Assume each file is fluence-- i.e. has actual number of events in the interval (dt already already applied)
// Multiply by 2.4 here for 40 kton

void garching_time_plot(TString indir, Double_t mass_factor)
{

  // Histogram of events per time bin, in logarithmic time bins, for plotting 

  const Int_t ntimebins = 50;
  Double_t timebins[ntimebins+1];


  Double_t t_offset = 0.02;
  Double_t t_start = -0.02+t_offset+0.005 ;
  Double_t t_end = 10.+t_offset;
  Double_t t_step = TMath::Log10(t_end/t_start)/double(ntimebins);

  Int_t it;
  for (it=0;it<ntimebins;it++) {
    timebins[it] = t_start*TMath::Power(10.0, double(it)*t_step);
    cout << it<<" "<<timebins[it]<<endl;
  }

  // This one for log x plot

  TH1D* eventsperbin = new TH1D("eventsperbin"," ",ntimebins-1,timebins);

  // This one for linear x plot
  TH1D* eventsperbinlin = new TH1D("eventsperbinlin"," ",100,0.,9.);

  TH1D* rateinbin = new TH1D("rateinbin"," ",ntimebins-1,timebins);


  // Now get the stuff from the files
  const Int_t maxpoints=1000;
  Int_t filenum[maxpoints];
  Double_t time[maxpoints];
  Double_t events[maxpoints];
  Double_t rate[maxpoints];
  Double_t dt[maxpoints];

  // Read the key file

  ifstream in;
  TString infilename = "garching_pinched_info_key.dat";

  in.open(infilename);
  Int_t j=0;
  Int_t numpoints;
    while(1) {
      in >> filenum[j]>>time[j]>>dt[j];
      cout << filenum[j]<< " "<<time[j]<<" "<< dt[j] <<endl;

      // Offset for log x scale in time
      time[j]+=0.02;
      
      if (!in.good()) break;
      j++;

    } // End of file reading loop
    
    numpoints = j;
    //    cout << "Read "<<numpoints<<" points from "<<infilename<<endl;
    
    in.close();

    Int_t i;

    
    Double_t allevents = 0; // Total events for the supernova

    for (i=0;i<numpoints;i++) {
      
      infilename = indir+"pinched_"+TString::Format("%d",filenum[i])+"_smeared_sum.dat";

      in.open(infilename);

      Double_t totevents = 0.;
      Double_t ev;
      Double_t en;
      Int_t k=0;


      // totevents here is integrating events over energy in the file
      // Apply the factor of 2.4 here for 40 kton (Snowglobes does 17 kt)

      while(1) {
	in >>en>>ev;
	totevents += ev*2.4;

	if (!in.good()) break;
	k++;


      }
      //      cout << "Read "<<k<<" points from "<<infilename<<endl;
      in.close();

      events[i] = totevents;
      // Divide by dt for rate
      rate[i] = totevents/dt[i];
      allevents += totevents;
      eventsperbin->Fill(time[i],events[i]);
      Int_t binnum = eventsperbin->FindBin(time[i]);

      // Root 6 default histo error behavior does not do the
      // right thing.  I want the sqrt of contents.  Can't seem
      // to find a Root method that does this


      eventsperbin->SetBinError(binnum,TMath::Sqrt(eventsperbin->GetBinContent(binnum)));

      //      cout << "Bin number: "<<binnum<<" events "<<events[i]<<" Content "<<eventsperbin->GetBinContent(binnum)<<" error "<<eventsperbin->GetBinError(binnum)<<endl;

      eventsperbinlin->Fill(time[i],events[i]);
      rateinbin->Fill(time[i],rate[i]);

      cout << i<< " time "<<time[i]<<" events "<<events[i]<<endl;

    }


    TGraph* rategraph = new TGraph(numpoints,time,rate);
    TGraph* timegraph = new TGraph(numpoints,time,events);

    TCanvas* canv1 = new TCanvas("c1"," ",800,700);

    Double_t xmin = 0.001;
    Double_t xmax = 9.;
    Double_t ymin = 0.1;
    Double_t ymax = 20000.;
    
    canv1->SetLogy(1.);
    canv1->SetLogx(1.);

    TH2F *hr = new TH2F("hr","",2,xmin,xmax,2,ymin,ymax);
    hr->SetXTitle(" Time (seconds) ");
    hr->SetYTitle(" ");

    hr->GetXaxis()->SetTitleSize(0.04);
    hr->GetXaxis()->SetTitleOffset(1.1);
    hr->GetXaxis()->SetLabelSize(0.04);
    hr->GetYaxis()->SetTitleSize(0.04);
    hr->GetYaxis()->SetTitleOffset(1.);
    hr->GetYaxis()->SetLabelSize(0.04);

    gStyle->SetOptStat(0);

    hr->GetYaxis()->SetRange(0.02,50);
    hr->Draw();

    rategraph->SetMarkerStyle(22);
    rategraph->SetMarkerSize(1.);
    rategraph->SetLineWidth(3.);
    //    rategraph->Print();
    rategraph->Draw("same");

    //////////

    TCanvas* canv2 = new TCanvas("c2"," ",800,700);

    xmin = 0.001;
    xmax = 9.;
    ymin = 0.001;
    ymax = 10.;
    
    //    canv2->SetLogy(1.);
    canv2->SetLogx(1.);

    TH2F *hr2 = new TH2F("hr2","",2,xmin,xmax,2,ymin,ymax);
    hr2->SetXTitle(" Time (seconds) ");
    hr2->SetYTitle(" ");

    hr2->GetXaxis()->SetTitleSize(0.04);
    hr2->GetXaxis()->SetTitleOffset(1.1);
    hr2->GetXaxis()->SetLabelSize(0.04);
    hr2->GetYaxis()->SetTitleSize(0.04);
    hr2->GetYaxis()->SetTitleOffset(1.);
    hr2->GetYaxis()->SetLabelSize(0.04);


    hr2->Draw();

    timegraph->SetMarkerStyle(22);
    timegraph->SetMarkerSize(1.);
    timegraph->SetLineWidth(3.);
    //    timegraph->Print();
    timegraph->Draw("psame");

    cout << "Rate graph integral: "<<rategraph->Integral(-1,-1)<<endl;
    cout << "All events: "<<allevents<<endl;


    //////

    TCanvas* canv3 = new TCanvas("c3"," ",800,700);

    canv3->SetLogx(1.);
    eventsperbin->SetXTitle("Time (seconds)");
    eventsperbin->SetYTitle("Events per bin");
    eventsperbin->SetMarkerStyle(22);
    eventsperbin->SetMarkerColor(4);
    eventsperbin->GetXaxis()->SetTitleSize(0.04);
    eventsperbin->GetXaxis()->SetTitleOffset(1.1);
    eventsperbin->GetXaxis()->SetLabelSize(0.04);
    eventsperbin->GetYaxis()->SetTitleSize(0.04);
    eventsperbin->GetYaxis()->SetTitleOffset(1.);
    eventsperbin->GetYaxis()->SetLabelSize(0.04);

//    Double_t xbins[ntimebins];
//    eventsperbin->Rebin(ntimebins-1,"eventsperbin",xbins);

    eventsperbin->Draw("pe");

    Double_t x1,x2,y1,y2;
    x1 = 0.02;
    x2 = 0.02;

    //    y1 = eventsperbin->GetMinimum();
    //y2 = eventsperbin->GetMaximum();
    y1 = -3.5;
    y2 = 66.5;


    cout << y1<< " "<<y2<<endl;
    TLine* collapsetime = new TLine(x1,y1,x2,y2);
    collapsetime->SetLineWidth(3);
    collapsetime->SetLineColor(2);
    collapsetime->Draw("same");


    //////

    //    TCanvas* canv4 = new TCanvas("c4"," ",800,700);

    //    canv4->SetLogy(1.);
    eventsperbinlin->SetMarkerStyle(22);
    eventsperbinlin->SetMarkerColor(3);
    //eventsperbinlin->Draw("pe");

    
    TFile f("eventsperbin.root","recreate");

    eventsperbin->Write();
    rateinbin->Write();
    f.Close();



}




