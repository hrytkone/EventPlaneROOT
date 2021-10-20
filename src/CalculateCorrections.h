const int nCorrections = 8;
const int ndet = 3;
const int nq = 2;
const int nbin = 400;
const double binmin = -0.2;
const double binmax = 0.2;

TFile *fin;
TTree *inTree;

TString detName[ndet] = {"FV0", "FT0A", "FT0C"};
TString saveFileName[ndet];
float corrections[nCorrections];
float qvec[ndet][nq];
TH2D *hQvec[ndet];

std::ofstream outputFile;

void InitOutput(TString outfile);
int LoadInput(TString infile);
void FillHistos();
void SaveCorrections();
void Print();
