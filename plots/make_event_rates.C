// Output graphs with event rates

void make_event_rates(TString fluxname, TString chanfile, TString exptconfig)
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
  chanfilename = "../channels/channels_"+chanfile+".dat";
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

  // Now get the info from the files.  Assumption is that all have the same binning.

  const Int_t maxpoints = 1000;
  Double_t total_events[maxpoints];
  Double_t total_events_smeared[maxpoints];
  Double_t en[maxpoints], eng[maxpoints];
  Double_t events[maxchan][maxpoints];
  Double_t events_smeared[maxchan][maxpoints];
  Double_t nue_es_smeared[maxpoints], nuebar_es_smeared[maxpoints], nux_es_smeared[maxpoints];

  Int_t ifile;
  TString filename[maxchan];
  TString filename_smeared[maxchan];

  Int_t k;
  for (k=0;k<maxpoints;k++) {
    total_events[k]=0;
    total_events_smeared[k]=0;
    nue_es_smeared[k]=0;
    nuebar_es_smeared[k]=0;
    nux_es_smeared[k]=0;
  }

  Int_t numpoints, numpoints_smeared;

  for (ifile=0; ifile<numchans; ifile++) {

    filename[ifile] = "../out/"+fluxname+"_"+channame[ifile]+"_"+exptconfig+"_events.dat";	// Build filename
    cout << ifile<<" "<<filename[ifile] <<endl;

    ifstream in;
    in.open(filename[ifile]);

    Int_t i=0;

    while (1) {

      in >> en[i]>>events[ifile][i];			// Energy in file is in GeV

      en[i]*= 1000;  					// in MeV for the plots

      eng[i]=en[i];
    
      //cout <<ifile<<" "<<i<<" "<<en[i]<<" "<<events[ifile][i]<<endl;

      //      events[ifile][i] *= num_target_factor[ifile];	//Account for the number of targets relative to reference target-- for unweighted file only


      total_events[i] += events[ifile][i];

      if (!in.good()) break;
      i++;

    }

    numpoints=i;

    printf(" found %d points\n",numpoints);

    in.close();

    // Now the smeared file

    filename_smeared[ifile] = "../out/"+fluxname+"_"+channame[ifile]+"_"+exptconfig+"_events_smeared.dat";
    cout << ifile<<" "<<filename_smeared[ifile] <<endl;

    in.open(filename_smeared[ifile]);
    
    i=0;
    while (1) {

      in >> en[i]>>events_smeared[ifile][i];			// Energy in file is in GeV
    

      en[i]*= 1000;  						// in MeV for the plots

      //      events_smeared[ifile][i] *= num_target_factor[ifile];	//Account for the number of targets relative to reference target-- for unweighted file only


      if(channame[ifile].CompareTo("nue_e") == 0) {			//
	nue_es_smeared[i] += events_smeared[ifile][i];			//	nue_e
      }									//

      if(channame[ifile].CompareTo("nuebar_e") == 0) {			//
	nuebar_es_smeared[i] += events_smeared[ifile][i];		//	nuebar_e
      }									//

      if(channame[ifile].CompareTo("numu_e") == 0 || channame[ifile].CompareTo("numubar_e") == 0 || channame[ifile].CompareTo("nutau_e") == 0 || channame[ifile].CompareTo("nutaubar_e") == 0) {					//
	nux_es_smeared[i] += events_smeared[ifile][i];			//	nux_e
      }									//

      total_events_smeared[i] += events_smeared[ifile][i];

      if (!in.good()) break;
      i++;

    }

    numpoints_smeared=i;

    printf(" found %d points\n",numpoints_smeared);

    in.close();


  }

  // Plots for total events

  TObjArray changraphlist(0);
  TObjArray changraph_smearedlist(0);
  TGraph* gr;
  TGraph* gr_smeared;

  gr= new TGraph(numpoints,en,total_events);
  gr->SetLineWidth(3);
  gr->SetLineColor(6);
  TString graphname=fluxname+"_tot_"+exptconfig;
  gr->SetName(graphname);

  changraphlist.Add(gr);

  gr_smeared= new TGraph(numpoints,en,total_events_smeared);
  gr_smeared->SetLineWidth(3);
  gr_smeared->SetLineColor(6);
  gr_smeared->SetLineStyle(2);
  graphname=fluxname+"_tot_"+exptconfig+"_smeared";
  gr_smeared->SetName(graphname);

  changraph_smearedlist.Add(gr_smeared);

  gr_smeared= new TGraph(numpoints,en,nue_es_smeared);
  gr_smeared->SetLineWidth(3);
  gr_smeared->SetLineColor(6);
  gr_smeared->SetLineStyle(2);
  graphname=fluxname+"_nue_es_"+exptconfig+"_smeared";
  gr_smeared->SetName(graphname);

  changraph_smearedlist.Add(gr_smeared);

  gr_smeared= new TGraph(numpoints,en,nuebar_es_smeared);
  gr_smeared->SetLineWidth(3);
  gr_smeared->SetLineColor(6);
  gr_smeared->SetLineStyle(2);
  graphname=fluxname+"_nuebar_es_"+exptconfig+"_smeared";
  gr_smeared->SetName(graphname);

  changraph_smearedlist.Add(gr_smeared);

  gr_smeared= new TGraph(numpoints,en,nux_es_smeared);
  gr_smeared->SetLineWidth(3);
  gr_smeared->SetLineColor(6);
  gr_smeared->SetLineStyle(2);
  graphname=fluxname+"_nux_es_"+exptconfig+"_smeared";
  gr_smeared->SetName(graphname);

  changraph_smearedlist.Add(gr_smeared);

  // Make the individual channel graphs

  Int_t water_colors[maxchan]={1,2,2,2,2,2,2,3,4,7,7,7,7,7,7,0,
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t argon_colors[maxchan]={2,2,2,2,2,2,3,4,7,7,7,7,7,7,0,
			       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t colors[maxchan];
  Int_t ic=0;

  if( exptconfig.CompareTo("ar17kt") == 0){
    for (ic=0;ic<numchans;ic++) {
      colors[ic] = argon_colors[ic];
    } 
  }
  else {
    // works for scint too
    for (ic=0;ic<numchans;ic++) {
      colors[ic] = water_colors[ic];
    } 
  }
 

  Int_t m,n;

  for (m=0;m<numchans;m++) {

    gr= new TGraph(numpoints,en,events[m]);
    gr->SetLineWidth(3); 
    gr->SetLineColor(colors[m]);
    graphname=fluxname+"_"+channame[channum[m]]+"_"+exptconfig;
    gr->SetName(graphname);

    changraphlist.Add(gr);

    gr_smeared = new TGraph(numpoints,en,events_smeared[m]);
    gr_smeared->SetLineWidth(3);
    gr_smeared->SetLineColor(colors[m]);
    gr_smeared->SetLineStyle(2);
    graphname=fluxname+"_"+channame[channum[m]]+"_"+exptconfig+"_smeared";
    gr_smeared->SetName(graphname);

    changraph_smearedlist.Add(gr_smeared);

  }

  // Sum the NC and ES channels

  Int_t ichan;
  Double_t es_events[maxpoints],es_events_smeared[maxpoints];
  Double_t ncO16_events[maxpoints],ncO16_events_smeared[maxpoints],ncC12_events[maxpoints],ncC12_events_smeared[maxpoints];

  Double_t nc_nu_Pb_1n_events[maxpoints],nc_nu_Pb_2n_events[maxpoints],nc_nubar_Pb_1n_events[maxpoints],nc_nubar_Pb_2n_events[maxpoints];
  Double_t tot_n_events[maxpoints];

  for (j=0;j<numpoints;j++) {
    es_events[j]=0;
    es_events_smeared[j]=0;
    ncO16_events[j]=0;
    ncO16_events_smeared[j]=0;
    ncC12_events[j]=0;
    ncC12_events_smeared[j]=0;

    tot_n_events[j]=0;
    nc_nu_Pb_1n_events[j]=0;
    nc_nu_Pb_2n_events[j]=0;
    nc_nubar_Pb_1n_events[j]=0;
    nc_nubar_Pb_2n_events[j]=0;

    Int_t firsteschan, lasteschan;
    if( exptconfig.CompareTo("ar17kt") == 0 || exptconfig.CompareTo("halo1")==0 || exptconfig.CompareTo("halo2") == 0){
      firsteschan = 0;
      lasteschan = 5;
    } 
    else {
      firsteschan = 1;
      lasteschan = 6;
    }


    for (ichan=firsteschan;ichan<=lasteschan;ichan++) {

      es_events[j] += events[ichan][j];
      es_events_smeared[j] += events_smeared[ichan][j];

    }
    if( exptconfig.CompareTo("wc100kt30prct") == 0 ||exptconfig.CompareTo("wc100kt15prct") == 0 ) {
      for (ichan=9;ichan<15;ichan++) {
	ncO16_events[j] += events[ichan][j];
	ncO16_events_smeared[j] += events_smeared[ichan][j];
      }
      //	cout << "O16  " << ncO16_events_smeared[j] << endl;
    }

    if( exptconfig.CompareTo("scint50kt") == 0) {
      for (ichan=9;ichan<15;ichan++) {
	ncC12_events[j] += events[ichan][j];
	ncC12_events_smeared[j] += events_smeared[ichan][j];
      }

    }

    if( exptconfig.CompareTo("halo1") == 0 || exptconfig.CompareTo("halo2") == 0) {
      nc_nu_Pb_1n_events[j] = events[8][j]+events[10][j]+events[12][j];
      nc_nubar_Pb_1n_events[j] = events[9][j]+events[11][j]+events[13][j];
      nc_nu_Pb_2n_events[j] = events[14][j]+events[16][j]+events[18][j];
      nc_nubar_Pb_2n_events[j] = events[15][j]+events[17][j]+events[19][j];
      tot_n_events[j] = nc_nu_Pb_1n_events[j]+nc_nubar_Pb_1n_events[j]+
	nc_nu_Pb_2n_events[j]+nc_nubar_Pb_2n_events[j]+
	events[6][j]+events[7][j];
      //	cout << "C12  " << ncC12_events_smeared[j] << endl;
    }
  

    
  }


  gr = new TGraph(numpoints,en,es_events);
  gr->SetLineWidth(3); 
  gr->SetLineColor(colors[1]);
  graphname=fluxname+"_es_"+exptconfig;
  gr->SetName(graphname);

  changraphlist.Add(gr);
  
  gr_smeared = new TGraph(numpoints,en,es_events_smeared);
  gr_smeared->SetLineWidth(3); 
  gr_smeared->SetLineColor(colors[1]);
  gr_smeared->SetLineStyle(2); 
  graphname=fluxname+"_es_"+exptconfig+"_smeared";
  gr_smeared->SetName(graphname);

  changraph_smearedlist.Add(gr_smeared);


    // Oxygen

  if( exptconfig.CompareTo("wc100kt30prct") == 0 ||exptconfig.CompareTo("wc100kt15prct") == 0 ) {
    gr = new TGraph(numpoints,en,ncO16_events);
    gr->SetLineWidth(3); 
    gr->SetLineColor(colors[9]);
    graphname=fluxname+"_ncO16_"+exptconfig;
    gr->SetName(graphname);

    changraphlist.Add(gr);

    gr_smeared = new TGraph(numpoints,en,ncO16_events_smeared);
    gr_smeared->SetLineWidth(3); 
    gr_smeared->SetLineStyle(2); 
    gr_smeared->SetLineColor(colors[9]);
    graphname=fluxname+"_ncO16_"+exptconfig+"_smeared";
    gr_smeared->SetName(graphname);

    changraph_smearedlist.Add(gr_smeared);

   }

     // C12
  if( exptconfig.CompareTo("scint50kt") == 0) {
    gr = new TGraph(numpoints,en,ncC12_events);
    gr->SetLineWidth(3); 
    gr->SetLineColor(colors[9]);
    graphname=fluxname+"_ncC12_"+exptconfig;
    gr->SetName(graphname);

    changraphlist.Add(gr);

    gr_smeared = new TGraph(numpoints,en,ncC12_events_smeared);
    gr_smeared->SetLineWidth(3); 
    gr_smeared->SetLineStyle(2); 
    gr_smeared->SetLineColor(colors[9]);
    graphname=fluxname+"_ncC12_"+exptconfig+"_smeared";
    gr_smeared->SetName(graphname);

    changraph_smearedlist.Add(gr_smeared);
  }

  if( exptconfig.CompareTo("halo1") == 0 || exptconfig.CompareTo("halo2") == 0) {

    gr = new TGraph(numpoints,en,tot_n_events);
    gr->SetLineWidth(3); 
    gr->SetLineColor(6);
    graphname=fluxname+"_tot_n_"+exptconfig;
    gr->SetName(graphname);

    changraphlist.Add(gr);

    gr = new TGraph(numpoints,en,nc_nu_Pb_1n_events);
    gr->SetLineWidth(3); 
    gr->SetLineColor(colors[9]);
    graphname=fluxname+"_nc_nu_Pb_1n_"+exptconfig;
    gr->SetName(graphname);

    changraphlist.Add(gr);

    gr = new TGraph(numpoints,en,nc_nubar_Pb_1n_events);
    gr->SetLineWidth(3); 
    gr->SetLineColor(colors[9]);
    graphname=fluxname+"_nc_nubar_Pb_1n_"+exptconfig;
    gr->SetName(graphname);

    changraphlist.Add(gr);

    gr = new TGraph(numpoints,en,nc_nu_Pb_2n_events);
    gr->SetLineWidth(3); 
    gr->SetLineColor(colors[9]);
    graphname=fluxname+"_nc_nu_Pb_2n_"+exptconfig;
    gr->SetName(graphname);

    changraphlist.Add(gr);

    gr = new TGraph(numpoints,en,nc_nubar_Pb_2n_events);
    gr->SetLineWidth(3); 
    gr->SetLineColor(colors[9]);
    graphname=fluxname+"_nc_nubar_Pb_2n_"+exptconfig;
    gr->SetName(graphname);

    changraphlist.Add(gr);


  }
  // Output the graphs to a file

  TString graphfilename = "graphs_"+fluxname+"_"+chanfile+"_"+exptconfig+".root";
  TFile f(graphfilename,"recreate");
  changraphlist.Write();
  changraph_smearedlist.Write();
  f.Close();

//---------------------------------------------------------------------------------------------------------------------------------------------------------


}
