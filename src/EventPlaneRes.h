const int ndet = 3;
const int nq = 2;
const TString detName[ndet] = {"FV0", "FT0A", "FT0C"};

TFile *fin, *fout;
TTree *inTree;

float epA[ndet];
float epB, epC;
float qvecA[ndet][nq];
float qvecB[nq];
float qvecC[nq];

// Resolution components
TH1D *hRAB[ndet];
TH1D *hRAC[ndet];
TH1D *hRBC;

// Event plane distributions
TH1D *hEPA[ndet];
TH1D *hEPB;
TH1D *hEPC;

int LoadInput(TString infile);
void InitOutput(TString outfile);
float GetEventPlane(float qx, float qy);
void FillHistos();
