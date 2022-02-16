// Calculate resolution & v2obs from tracks with 2 sub events
// This is to get "true" v2 value for comparion

const int nq = 2;

TFile *fin, *fout;
TTree *inTree;

float epB, epC, epFull;
float qvecB[nq];
float qvecC[nq];
float qvecFull[nq];
std::vector<float> *tpcPhi = 0;

// Resolution component
TH1D *hRsub;

// Event plane distributions
TH1D *hEPB;
TH1D *hEPC;
TH1D *hEPFull;

// vnobs
TH1D *hVnObs;

int LoadInput(TString infile);
void InitOutput(TString outfile);
float GetEventPlane(float qx, float qy);
float GetVnObs(float ep, int n);
void FillHistos();
