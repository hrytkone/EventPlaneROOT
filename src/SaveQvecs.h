const double bmin = 0.; // To change impact parameter range if data is min bias
const double bmax = 20.;
const int nq = 2;
const int ndet = 3;
const int ncorr = 8;
const bool bDoCorrections = true;
const TString detName[ndet] = {"FV0", "FT0A", "FT0C"};

TFile *fout;
TTree *outTree;
std::ifstream corrInFile;

float corr[ndet][ncorr];
float qvecFV0[nq];
float qvecFT0A[nq];
float qvecFT0C[nq];
float qvecB[nq];
float qvecC[nq];

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
int FillQvecFV0(UInt_t ient, int bcbegin);
int FillQvecFT0(UInt_t ient, int bcbegin);
void SumQvec(TComplex &Qvec, double nch, int chno, TString det);
double GetFV0Phi(int chno);
double GetFT0APhi(int chno);
double GetFT0CPhi(int chno);