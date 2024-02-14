{
   TSystemDirectory folder(".", ".");
   TList *files = folder.GetListOfFiles();
   if (!files) {
      cout<<"the folder is empty."<<endl;
      return;
   }

   // open output pdf
   gROOT->SetStyle("Plain"); // pick up a good default drawing style to modify
   gStyle->SetTitleFont(132,"XY");
   gStyle->SetLabelFont(132,"XY");
   gStyle->SetLabelSize(0.05,"XY");
   gStyle->SetTitleSize(0.05,"XY");
   gStyle->SetPadRightMargin(0.01);
   TCanvas *can = new TCanvas;
   can->SetLogy();
   can->SetGridx(); can->SetGridy();
   can->Print("plots.pdf[");

   // draw plots
   const int n = files->GetSize();
   TGraph *gr[100] = {0};
   TSystemFile *file;
   TIter next(files);
   int i=0;
   while((file = (TSystemFile*) next())) {
      if (file->IsDirectory()) continue;
      TString name = file->GetName();
      if (!name.EndsWith(".out")) continue;
      cout<<"plot "<<file->GetName()<<endl;

      gr[i] = new TGraph(file->GetName());
      if (gr[i]->GetN()==0) {
         cout<<"no data to plot in "<<name<<endl;
         continue;
      }
      gr[i]->Draw("al");
      gr[i]->GetYaxis()->SetTitle("Event rate [/MeV]");
      gr[i]->GetXaxis()->SetTitle("Energy [MeV]");
      can->Print("plots.pdf");
      i++;
   }

   // close output pdf
   can->Print("plots.pdf]");
}
