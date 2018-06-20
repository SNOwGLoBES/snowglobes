// Make cross-section graphs
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


void make_xscn_graphs(TString chanfilename)
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

  // Loop through channels and read the x-scn file

  TString xscnfilename[maxchan];
  //  TGraph *xscngraphs[maxchan];



  TObjArray xscngraphlist(0);
  TGraph* gr;


  Int_t water_colors[maxchan]={1,2,2,2,2,2,2,3,4,7,7,7,7,7,7,0,
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t argon_colors[maxchan]={2,2,2,2,2,2,3,4,7,7,7,7,7,7,0,
			       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t colors[maxchan];
  Int_t ic=0;
  if( chanfilename.Contains("argon") == 1){
    for (ic=0;ic<maxchan;ic++) {
      colors[ic] = argon_colors[ic];
    } 
  }
  else {
    for (ic=0;ic<maxchan;ic++) {
      colors[ic] = water_colors[ic];
    } 
  }

  const Int_t max_xscn_points = 1001;
    
  // Loop over channels and fill the histos
  for (j=0;j<numchans;j++) {
    xscnfilename[j]="../xscns/xs_"+channame[j]+".dat";

    cout << j<<" " <<xscnfilename[j]<<endl;
    ifstream xscnin;
    xscnin.open(xscnfilename[j]);
    Int_t nlines = 0;


    Double_t energy[max_xscn_points];
    Double_t nuexscn[max_xscn_points],numuxscn[max_xscn_points],nutauxscn[max_xscn_points],nuebarxscn[max_xscn_points],numubarxscn[max_xscn_points],nutaubarxscn[max_xscn_points];

    Double_t xscn[max_xscn_points];

    Int_t k=0;
    for (k=0;k<max_xscn_points;k++) {

      energy[k]=0;
      nuexscn[k]=0;
      numuxscn[k]=0;
      nutauxscn[k]=0;
      nuebarxscn[k]=0;
      numubarxscn[k]=0;
      nutaubarxscn[k]=0;
      xscn[k]=0;
    }

    Int_t i=0;
    Float_t logen=0;
    while (1) {
      // Ignore comments or null lines
      string line;
      string comchar1 = "#";
      getline(xscnin,line);
      string firstchar = line.substr(0,1);
      if (firstchar.compare(comchar1) == 0 || line.length()==1) {
	continue;
      }

      istringstream instream;
      instream.str(line);
      instream >> logen>>nuexscn[i]>>numuxscn[i]>>nutauxscn[i]>>nuebarxscn[i]
      	     >>numubarxscn[i]>>nutaubarxscn[i];
      // Skip zero energy (or blank)
      if (logen==0) {continue;}

      energy[i] = pow(10,logen);  //Energy in GeV

      // x-scn is divided by energy in GeV in file
      nuexscn[i] *= energy[i];
      nuebarxscn[i] *= energy[i];
      numuxscn[i] *= energy[i];
      numubarxscn[i] *= energy[i];
      nutauxscn[i] *= energy[i];
      nutaubarxscn[i] *= energy[i];
    

      // Convert to MeV
      energy[i]*=1000;

      cout <<i<<" "<<logen<<" "<<energy[i]<<" "<<" "<<nuexscn[i]<<endl;


      if (!xscnin.good()||i>=1000) break;
      i++;

    }

    Int_t numxscnpoints=i;

    printf(" Found %d points\n",numxscnpoints);

    xscnin.close();

    // Assume equidistant points

    Int_t linestyle;

    // Column to use depends on flavor, CP state

    cout << "Flavor, CP state: "<<j<<" "<<flav[j]<<" "<<cpstate[j]<<endl;
    if (flav[j].CompareTo("e")==0&&cpstate[j].CompareTo("+")==0) {
      for (k=0;k<numxscnpoints;k++) {
	xscn[k] = nuexscn[k];
      }
      linestyle = 2;

    } else { 
      if (flav[j].CompareTo("e")==0&&cpstate[j].CompareTo("-")==0) {
	for (k=0;k<numxscnpoints;k++) {
	  xscn[k] = nuebarxscn[k];
	}
	linestyle = 1;


      } else {  
	if (flav[j].CompareTo("m")==0&&cpstate[j].CompareTo("+")==0) {
	  for (k=0;k<numxscnpoints;k++) {
	    xscn[k] = numuxscn[k];
	  }
	  linestyle = 3;


	} else {
	  if (flav[j].CompareTo("m")==0&&cpstate[j].CompareTo("-")==0) {
	    for (k=0;k<numxscnpoints;k++) {
	      xscn[k] = numubarxscn[k];
	    }
	    linestyle = 3;



	  } else {
	    if (flav[j].CompareTo("t")==0&&cpstate[j].CompareTo("+")==0) {

	      for (k=0;k<numxscnpoints;k++) {
		xscn[k] = nutauxscn[k];
	      }
	      linestyle = 3;


		} else { 
	      if (flav[j].CompareTo("t")==0&&cpstate[j].CompareTo("-")==0) {

		for (k=0;k<numxscnpoints;k++) {
		  xscn[k] = nutaubarxscn[k];
		}
		linestyle = 3;

	      } else {
		cout << "Error! No flavor, cpstate for this channel" <<endl;
		exit(0);
	      }
	    }
	  }
	}
      }
    }


    // Fill the graphs


    gr = new TGraph(numxscnpoints,energy,xscn);

    cout << j <<" "<<channum[j]<<" "<<channame[channum[j]]<<endl;
    gr->Print("all");

    xscngraphlist.Add(gr);

    gr->SetName(channame[channum[j]]);

    if( channame[channum[j]].Contains("nc") == 1){ linestyle=1;}

    gr->SetLineStyle(linestyle);


    gr->SetLineWidth(3);
    gr->SetLineColor(colors[channum[j]]);

    //    cout << channame[j]<<" "<<channum[j]<<" setting color:"<<colors[channum[j]]<<" style "<<linestyle<<endl;
    //cout<< "In graph:"<<gr->GetLineColor()<<endl;
    

  }

  // Output the graphs to a file

  TString graphfilename = "graphs_"+chanfilename+".root";
  TFile f(graphfilename,"recreate");
  xscngraphlist.Write();
  f.Close();



  Double_t xmin = 0.;
  Double_t xmax = 100.;
  Double_t ymin = 0.00001;

  // Because graph method not working
  Int_t p;
  Double_t ymax=0.1;
  //  for (p=0;p<numpoints;p++) {
  // if (nuescan[p]>ymax) {
  //   ymax = nux[p];
  // }
  //}

  ymax*=1.1;


  TCanvas* canv = new TCanvas("c1"," ",800,700);
  canv->SetLogy(1);

  TString plot_title = "Cross-sections";

  gStyle->SetTitleFontSize(.04);

  TH2F *hr = new TH2F("hr",plot_title,2,xmin,xmax,2,ymin,ymax);
  hr->SetXTitle(" Neutrino Energy (MeV) ");
  hr->SetYTitle(" Cross-section ");
  hr->GetXaxis()->SetTitleSize(0.03);
  hr->GetXaxis()->SetTitleOffset(0.85);
  hr->GetXaxis()->SetLabelSize(0.03);
  hr->GetYaxis()->SetTitleSize(0.03);
  hr->GetYaxis()->SetTitleOffset(1.7);
  hr->GetYaxis()->SetLabelSize(0.03);

  gStyle->SetOptStat(0);

  hr->Draw();


  for (int i=0;i<numchans;i++) {

    gr = (TGraph*)xscngraphlist.At(i);

    gr->Draw();


  }


   
   //   TString printfilename = "xscns.gif";
   // canv->Print(printfilename);


}
