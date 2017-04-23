// Tool to turn digitized Garching files into pinched flux files for time slices
// Reads the digitized files, makes garching_pinched_info.dat for use with
// Snowglobes pinched.C file
// Also makes garching_key.dat to associate flux number with time slide
// Note: produces luminosities appropriate for *fluences* in Snowglobes-- luminosities are in ergs/s, then integrated over the given time bin.
// K. Scholberg December 2013


void make_garching_fluxes()
{

  gStyle->SetOptStat(0);

  ifstream infile;
  // Two timescales, 6 flavors, three types of quantities (lum, eavg, alpha)
  // One graph array per quantity, indexed by flavor. 

  const Int_t numflavor = 3;
  //  TString flavname[6] = {"nue","nuebar","numu","numubar","nutau","nutaubar"};
  TString flavname[numflavor] = {"nue","nuebar","numu"};

  // numu is same as nux

  TGraph** luminosity_graphs = new TGraph*[numflavor];
  TGraph** avgen_graphs = new TGraph*[numflavor];
  TGraph** alpha_graphs = new TGraph*[numflavor];
  

  const Int_t maxpoints=2000;
  Double_t time[maxpoints];
  Double_t yval[maxpoints];

  Int_t i=0;
  Int_t j=0;

  for (i;i<numflavor;i++){

    // Luminosity
    ifstream in;

    TString infilename2 = "luminosity_"+flavname[i]+"_garching.txt";

    cout << "Reading "<< infilename2<<endl;

    in.open(infilename2);
    j=0;
    Int_t numpoints;
    
    while(1) {
      in >> time[j]>>yval[j];
      //      cout << "time "<<time[j]<<" "<<yval[j]<<endl;
      if (yval[j]<0) {yval[j]=0;}
      
      if (!in.good()) break;
      j++;

    } // End of file reading loop

    numpoints = j;
    in.close();
    
    luminosity_graphs[i]= new TGraph(numpoints,time,yval);

    // Average energy

    TString infilename3 = "avgen_"+flavname[i]+"_garching.txt";

    cout <<"Reading " << infilename3<<endl;

    in.open(infilename3);
    j=0;
    
    while(1) {
      in >> time[j]>>yval[j];
      //     cout << "time "<<time[j]<<" "<<yval[j]<<endl;
      
      if (yval[j]<0) {yval[j]=0;}

      if (!in.good()) break;
      j++;

    } // End of file reading loop

    numpoints = j;
    in.close();
    
    avgen_graphs[i]= new TGraph(numpoints,time,yval);

    // Alpha

    TString infilename4 = "alpha_"+flavname[i]+"_garching.txt";

    cout << "Reading "<<infilename4<<endl;

    in.open(infilename4);
    j=0;
    
    while(1) {
      in >> time[j]>>yval[j];
      //      cout << "time "<<time[j]<<" "<<yval[j]<<endl;
      
      if (yval[j]<0) {yval[j]=0;}

      if (!in.good()) break;
      j++;

    } // End of file reading loop

    numpoints = j;
    in.close();
    
    alpha_graphs[i]= new TGraph(numpoints,time,yval);

    //    graphs[i]->Print();

  } // End of loop over flavors

  // Now choose time intervals in log steps

  const Int_t ntimebins = 300;
  Double_t timebins[ntimebins+1];

  // Want log steps but starts negative, so put an offset

  Double_t t_offset = 0.02;
  Double_t t_start = -0.02+t_offset+0.001 ;
  Double_t t_end = 9.+t_offset;
  Double_t t_step = TMath::Log10(t_end/t_start)/double(ntimebins);

  //  cout << "t_step "<<t_step<<endl;

  Double_t t;

  Double_t tot_dt=0;

  // Get the ranges of each of the graphs

  Int_t k;
  Double_t xx,yy;
  Double_t alpha_first[numflavor];
  Double_t avgen_first[numflavor];
  Double_t luminosity_first[numflavor];
  Double_t alpha_last[numflavor];
  Double_t avgen_last[numflavor];
  Double_t luminosity_last[numflavor];

  //Saved for nue flux
  Double_t alpha_val_first[numflavor];
  Int_t ret;

  for (k=0;k<numflavor;k++) {
    Int_t numpt;
    numpt = alpha_graphs[k]->GetN();
    if (numpt>0) {
      ret = alpha_graphs[k]->GetPoint(0,xx,yy);
      alpha_first[k] = xx;
      alpha_val_first[k] = yy;
      ret = alpha_graphs[k]->GetPoint(numpt-1,xx,yy);
      alpha_last[k] = xx;

    }

    numpt = avgen_graphs[k]->GetN();
    if (numpt>0) {
      ret = avgen_graphs[k]->GetPoint(0,xx,yy);
      avgen_first[k] = xx;
      ret = avgen_graphs[k]->GetPoint(numpt-1,xx,yy);
      avgen_last[k] = xx;
    }

    numpt = luminosity_graphs[k]->GetN();
    if (numpt>0) {
      ret = luminosity_graphs[k]->GetPoint(0,xx,yy);
      luminosity_first[k] = xx;
      ret = luminosity_graphs[k]->GetPoint(numpt-1,xx,yy);
      luminosity_last[k] = xx;
    }


    //    cout << "flav "<<k<<" alpha "<< alpha_first[k]<<" "<<alpha_last[k]<<endl;
    //cout << "flav "<<k<<" avgen "<< avgen_first[k]<<" "<<avgen_last[k]<<endl;
    //cout << "flav "<<k<<" luminosity "<< luminosity_first[k]<<" "<<luminosity_last[k]<<endl;
  }


  // Output files

  ofstream outfile;
  outfile.open("garching_pinched_info.dat");


  ofstream outfile2;
  outfile2.open("garching_pinched_info_key.dat");


  for (i=0;i<ntimebins;i++)   {
      t = t_start*TMath::Power(10.0, double(i)*t_step);
      timebins[i]= t-t_offset;

      //      cout <<timebins[i]<<" "<<luminosity_graphs[0]->Eval(timebins[i],0,"")<<endl;

      // Approximate dt
      Double_t dt = t_start*TMath::Power(10.0,double(i+0.5)*t_step)-
                      t_start*TMath::Power(10.0,double(i-0.5)*t_step);
      //      cout << i<<" "<<t_start*TMath::Power(10.0,double(i-0.5)*t_step)-t_offset<<" "<<timebins[i]<<" "<<t_start*TMath::Power(10.0,double(i+0.5)*t_step)-t_offset<<" "<<dt<<endl;

      tot_dt += dt;  // for check


      Double_t alpha_nue, avgen_nue, luminosity_nue;
      Double_t alpha_nuebar, avgen_nuebar, luminosity_nuebar;
      Double_t alpha_nux, avgen_nux, luminosity_nux;

      // Get the values from the plots.  But first check we are within the meaningful range for each graph for each flavor.


      // Fine-tuning for nue for this flux:  nue luminosity and energies are non-zero at negative times; assume alpha is the same as at the beginning for the range.

      if (timebins[i]<alpha_last[0]
	  && timebins[i]>avgen_first[0] && timebins[i]<avgen_last[0] 
          && timebins[i]>luminosity_first[0] && timebins[i]<luminosity_last[0]) {

	if (timebins[i]>alpha_first[0]) {
	  alpha_nue = alpha_graphs[0]->Eval(timebins[i],0,""); 
	  avgen_nue = avgen_graphs[0]->Eval(timebins[i],0,""); 
	  luminosity_nue = luminosity_graphs[0]->Eval(timebins[i],0,""); 
	} else {

	  // Use the saved value
	  alpha_nue = alpha_val_first[0];
	  avgen_nue = avgen_graphs[0]->Eval(timebins[i],0,""); 
	  luminosity_nue = luminosity_graphs[0]->Eval(timebins[i],0,""); 

	}

	
      } else {
	alpha_nue=0.;
	avgen_nue = 0.;
	luminosity_nue=0.;
      }


      // Nuebar
      if (timebins[i]>alpha_first[1] && timebins[i]<alpha_last[1]
	  && timebins[i]>avgen_first[1] && timebins[i]<avgen_last[1] 
          && timebins[i]>luminosity_first[1] && timebins[i]<luminosity_last[1]) {

	alpha_nuebar = alpha_graphs[1]->Eval(timebins[i],0,""); 
	avgen_nuebar = avgen_graphs[1]->Eval(timebins[i],0,""); 
	luminosity_nuebar = luminosity_graphs[1]->Eval(timebins[i],0,""); 
	
      } else {
	alpha_nuebar=0.;
	avgen_nuebar = 0.;
	luminosity_nuebar=0.;
      }

      //Numu

      if (timebins[i]>alpha_first[2] && timebins[i]<alpha_last[2]
	  && timebins[i]>avgen_first[2] && timebins[i]<avgen_last[2] 
          && timebins[i]>luminosity_first[2] && timebins[i]<luminosity_last[2]) {

	alpha_nux = alpha_graphs[2]->Eval(timebins[i],0,""); 
	avgen_nux = avgen_graphs[2]->Eval(timebins[i],0,""); 
	luminosity_nux = luminosity_graphs[2]->Eval(timebins[i],0,""); 
	
      } else {
	alpha_nux=0.;
	avgen_nux = 0.;
	luminosity_nux=0.;
      }

  // Save the time and dt for each flux file number in the key file
      outfile2 << i<< " "<<timebins[i]<<" "<<dt<<endl;

      // Correct luminosity units
      //      cout <<"timebin"<<timebins[i]<<" "<<dt<<endl;
      dt*= 1.e52;

  outfile << i<< " "
       <<alpha_nue<<" "
       <<alpha_nuebar<<" "
       <<alpha_nux<<" "
       <<avgen_nue<<" "
       <<avgen_nuebar<<" "
       <<avgen_nux<<" "
       <<luminosity_nue*dt<<" "
       <<luminosity_nuebar*dt<<" "
       <<luminosity_nux*dt<<endl;



  //  cout << i<< " "<<timebins[i]<<" "
  //     <<alpha_nue<<" "
  //     <<alpha_nuebar<<" "
  //     <<alpha_nux<<" "
  //     <<avgen_nue<<" "
  //     <<avgen_nuebar<<" "
  //     <<avgen_nux<<" "
  //     <<luminosity_nue*dt<<" "
  //     <<luminosity_nuebar*dt<<" "
  //     <<luminosity_nux*dt<<endl;



    }

  //  cout << "tot_dt "<< tot_dt <<" "<<t_end-t_start<<endl;

  outfile.close();
}

