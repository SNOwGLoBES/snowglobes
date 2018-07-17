
// Create smearing matrices by multiplying interaction and resolution files
void read_matrix(TString, TMatrixD*,Int_t);

void write_smear(TString,TString,TMatrixD*, Int_t);
void smear(TString config, TString channel, TString resolution)
{


  const Int_t matsize = 200;

  TMatrixD intmat(matsize,matsize);
  TMatrixD resmat(matsize,matsize);
  TMatrixD smearmat(matsize,matsize);

  // Read in the interaction matrix

  TString interaction = "interaction/interaction_"+channel+".ssv";
  cout << "Reading file: "<<interaction<<endl;

  read_matrix(interaction,&intmat,matsize);

  TString resolution_file = "resolution/resolution_"+resolution+".ssv";
  cout << "Reading file: "<<resolution_file<<endl;


  read_matrix(resolution_file,&resmat,matsize);

  smearmat = resmat*intmat;

  write_smear(channel,config,&smearmat,matsize);
  //  smearmat.Print();

}

// Read the matrix from the file, formatted by rows
// Assumes a square matrix of size matsize
void read_matrix(TString filename,TMatrixD* matrix, Int_t matsize) {

  ifstream in;
  in.open(filename);
  Int_t nlines = 0;

  Int_t irow=0;
  Int_t icol=0;
  while (in.good()&& irow<matsize) {
    for (icol=0;icol<matsize;icol++) {

      Double_t entry;
      in >> entry;
      //      cout << "irow "<<irow<<" icol "<<icol<<" "<<entry<<endl;

      matrix(irow,icol) = entry;

      //      if (!in.good() && (irow !=matsize-1 || icol != matsize-1)) {
      if (!in.good()) {


      //  cout << "Error reading "<<filename<<" at "<<irow<<" "<<icol<<endl;
      //  exit();
      //	}
      //else {
      	break;
      }
    }
    irow++;
  }

  in.close();

}



void write_smear(TString channel, TString config,TMatrixD* smearmat,Int_t matsize) {
  // Write the smearing matrix to the file

  TString outfilename = "out/smear_"+channel+"_"+config+".dat";
  cout << "Creating file: "<<outfilename<<endl;

  ofstream out;
  out.open(outfilename);

  Int_t irow=0;
  Int_t icol=0;


  TString firstline = "energy(#"+channel+"_smear)<";

  out << firstline <<endl;
  out << "@energy = ";
  for (irow=0;irow<matsize-1;irow++) {
    out << "{0,"<<matsize-1<<",";
    for (icol=0;icol<matsize-1;icol++) {
      out <<smearmat(irow,icol)<<",";
      
    }
    out << smearmat(irow,matsize-1)<<"}:"<<endl;

  }
  // Last line
  out << "{0,"<<matsize-1<<",";
  for (icol=0;icol<matsize-1;icol++) {
    out <<smearmat(matsize-1,icol)<<",";
    
  }
  out << smearmat(matsize-1,matsize-1)<<"};"<<endl;
  out << ">"<<endl;

  out.close();
}
