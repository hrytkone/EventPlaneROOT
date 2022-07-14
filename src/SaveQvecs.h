const double bmin = 0.; // To change impact parameter range if data is min bias
const double bmax = 20.;
const int nq = 2;
const int ndet = 3;
const int ncorr = 8;
const TString detName[ndet] = {"FV0", "FT0A", "FT0C"};
//const float asideOffsetX = 0.314;
const float asideOffsetX = 0.;
//const float asideOffsetY = 0.718;
const float asideOffsetY = 0.;

bool bDoCorrections = true;
bool bUseTPCeff = true;
TFile *fout;
TTree *outTree;
TFile *finTPCeff;           // efficiency for TPC
TH1D* hCoeff;
std::ifstream corrInFile;

float corr[ndet][ncorr];
float qvecFV0[nq];          // subevent A for FV0 (digits)
float qvecFT0A[nq];         // subevent A for FT0A (digits)
float qvecFT0C[nq];         // subevent A for FT0C (digits)
float qvecA[nq];            // subevent A for 'ideal' calculation (MCTracks)
float qvecB[nq];            // subevent B (0.1 < eta < 0.8) (MCTracks)
float qvecC[nq];            // subevent C (-0.8 < eta < -0.1) (MCTracks)
float qvecFull[nq];         // subevent B + subevent C
std::vector<float> * tpcPhi = 0;

TFile *finKine;
TFile *finFV0Digit;
TFile *finFT0Digit;

TTree *fKineTree;
TTree *fFV0DigitTree;
TTree *fFT0DigitTree;

std::vector<o2::MCTrack>* mctrack = nullptr;
o2::dataformats::MCEventHeader* mcheader = nullptr;

std::vector<o2::fv0::BCData> fv0BCData, *fv0BCDataPtr = &fv0BCData;
std::vector<o2::fv0::ChannelData> fv0ChData, *fv0ChDataPtr = &fv0ChData;
o2::dataformats::MCTruthContainer<o2::fv0::MCLabel> fv0Labels, *fv0LabelsPtr = &fv0Labels;

std::vector<o2::ft0::Digit> ft0BCData, *ft0BCDataPtr = &ft0BCData;
std::vector<o2::ft0::ChannelData> ft0ChData, *ft0ChDataPtr = &ft0ChData;
o2::dataformats::MCTruthContainer<o2::ft0::MCLabel> ft0labels, *ft0labelsPtr = &ft0labels;

void InitOutput(TString outfile);
int LoadInput(TString nameKineFile, TString nameFV0DigitFile, TString nameFT0DigitFile);
void LoadCorrections(TString cent);
void CloseFiles();

void FillQvectors();
void FillQvecBC(UInt_t ient);
void FillQvecA(UInt_t ient);
int FillQvecFV0(UInt_t ient, int bcbegin);
int FillQvecFT0(UInt_t ient, int bcbegin);
void SumQvec(TComplex &Qvec, double nch, int chno, TString det);
double GetFV0Phi(int chno);
double GetFT0APhi(int chno);
double GetFT0CPhi(int chno);
double GetEffFromHisto(TH1D* h, double pt);
bool GoodEvent();
bool IsDeadChannel(int ich);
