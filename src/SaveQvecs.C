#include "SaveQvecs.h"
#include "Corrections.h"
#include "Const.h"

void SaveQvecs(TString kinefile="", TString fv0digitfile="", TString ft0digitfile="",
                TString outfile="output.root", TString cent="20-30")
{
    int isOpen = LoadInput(kinefile, fv0digitfile, ft0digitfile);
    if (!isOpen) return;
    if (bDoCorrections) LoadCorrections(cent);

    InitOutput(outfile);
    FillQvectors();
    fout->Write("", TObject::kOverwrite);
    CloseFiles();
}

//_____________________________________________________________
void InitOutput(TString outfile)
{
    if (!gSystem->AccessPathName(outfile.Data())) {
        fout = TFile::Open(outfile.Data(), "UPDATE");
        outTree = (TTree*)fout->Get("Qvecs");
        outTree->SetBranchAddress("qvecFV0", &qvecFV0);
        outTree->SetBranchAddress("qvecFT0A", &qvecFT0A);
        outTree->SetBranchAddress("qvecFT0C", &qvecFT0C);
        outTree->SetBranchAddress("qvecB", &qvecB);
        outTree->SetBranchAddress("qvecC", &qvecC);
    } else {
    	fout = TFile::Open(outfile.Data(), "RECREATE");
        outTree = new TTree("Qvecs", "Q-vectors for FIT detectors");
        outTree->Branch("qvecFV0", &qvecFV0, Form("qvecFV0[%d]/F", nq));
        outTree->Branch("qvecFT0A", &qvecFT0A, Form("qvecFT0A[%d]/F", nq));
        outTree->Branch("qvecFT0C", &qvecFT0C, Form("qvecFT0C[%d]/F", nq));
        outTree->Branch("qvecB", &qvecB, Form("QvecB[%d]/F", nq));
        outTree->Branch("qvecC", &qvecC, Form("QvecC[%d]/F", nq));
    }
}

int LoadInput(TString nameKineFile, TString nameFV0DigitFile, TString nameFT0DigitFile)
{
    finKine = TFile::Open(nameKineFile, "READ");
    finFV0Digit = TFile::Open(nameFV0DigitFile, "READ");
    finFT0Digit = TFile::Open(nameFT0DigitFile, "READ");

    if (!finKine || !finFV0Digit || !finFT0Digit) {
        std::cout << "Could not open the files, skipping them.." << std::endl;
        return 0;
    }

    fKineTree = (TTree*)finKine->Get("o2sim");
    fFV0DigitTree = (TTree*)finFV0Digit->Get("o2sim");
    fFT0DigitTree = (TTree*)finFT0Digit->Get("o2sim");

    fFV0DigitTree->SetBranchAddress("FV0DigitBC", &fv0BCDataPtr);
    fFV0DigitTree->SetBranchAddress("FV0DigitCh", &fv0ChDataPtr);
    fFV0DigitTree->SetBranchAddress("FV0DigitLabels", &fv0LabelsPtr);

    fFT0DigitTree->SetBranchAddress("FT0DIGITSBC", &ft0BCDataPtr);
    fFT0DigitTree->SetBranchAddress("FT0DIGITSCH", &ft0ChDataPtr);
    fFT0DigitTree->SetBranchAddress("FT0DIGITSMCTR", &ft0labelsPtr);

    return 1;
}

void LoadCorrections(TString cent)
{
    for (int idet=0; idet<ndet; idet++) {
        corrInFile.open(Form("corrections/corr_%s_%s.txt", detName[idet].Data(), cent.Data()));
        std::string line;
        int icorr = 0;
        while (std::getline(corrInFile, line)) {
            corr[idet][icorr] = atof(line.c_str());
            icorr++;
        }
        corrInFile.close();
    }
}

void CloseFiles()
{
    finKine->Close();
    finFV0Digit->Close();
    finFT0Digit->Close();
}

void FillQvectors()
{
    auto mcbr = fKineTree->GetBranch("MCTrack");
    auto hdrbr = fKineTree->GetBranch("MCEventHeader.");
    mcbr->SetAddress(&mctrack);
    hdrbr->SetAddress(&mcheader);

    UInt_t nEntries = fKineTree->GetEntries();
    int fv0bcbegin = 0;
    int ft0bcbegin = 0;
    for (UInt_t ient = 0; ient < nEntries; ient++) {

        double b = mcheader->GetB();
        if (b < bmin || b > bmax) continue;

        FillQvecBC(ient);
        fv0bcbegin = FillQvecFV0(ient, fv0bcbegin);
        ft0bcbegin = FillQvecFT0(ient, ft0bcbegin);
        if (bDoCorrections) {
            DoCorrections(qvecFV0[0], qvecFV0[1], corr[0]);
            DoCorrections(qvecFT0A[0], qvecFT0A[1], corr[1]);
            DoCorrections(qvecFT0C[0], qvecFT0C[1], corr[2]);
        }
        outTree->Fill();
    }
}

void FillQvecBC(UInt_t ient)
{
    fKineTree->GetEntry(ient);

    TComplex QvecB(0), QvecC(0);
    int nTracksB = 0, nTracksC = 0;
    for (auto &t : *mctrack) {

        if (t.GetPt() < 0.2) continue;

        Int_t pid = t.GetPdgCode();
        //if (TMath::Abs(pid)>1000000000) continue;
        //if (TDatabasePDG::Instance()->GetParticle(pid)==NULL) continue;
        Double_t charge = TDatabasePDG::Instance()->GetParticle(pid)->Charge();
        if (charge==0.0) continue;
        //double len = TMath::Sqrt(t.Vx()*t.Vx() + t.Vy()*t.Vy() + t.Vz()*t.Vz());
        //if (len > 0.5) continue;
        //if (!isHadron(pid)) continue;

        double eta = t.GetEta();
        double phi = TMath::ATan2(t.Py(), t.Px());

        if ( eta>-0.8 && eta<-0.1 ) {
            QvecB += TComplex(TMath::Cos(2.0 * phi), TMath::Sin(2.0 * phi));
            nTracksB++;
        }

        if ( eta>0.1 && eta<0.8 ) {
            QvecC += TComplex(TMath::Cos(2.0 * phi), TMath::Sin(2.0 * phi));
            nTracksC++;
        }
    }

    QvecB /= (double)nTracksB;
    QvecC /= (double)nTracksC;
    qvecB[0] = QvecB.Re();
    qvecB[1] = QvecB.Im();
    qvecC[0] = QvecC.Re();
    qvecC[1] = QvecC.Im();
}

int FillQvecFV0(UInt_t ient, int bcbegin)
{
    fFV0DigitTree->GetEntry(ient);

    TComplex Qvec(0);
    double signalSum = 0;
    int prevLabelEvent = -1;
    int nbc = fv0BCData.size();
    for (int ibc = bcbegin; ibc < nbc; ibc++) {
        const auto lb = fv0Labels.getLabels(ibc);
        int labelEvent = 0;
        if (lb.size() > 0) {
            labelEvent = lb[0].getEventID();
        } else {
            //std::cout << "Problem with labels" << std::endl;
            continue;
        }

        if (prevLabelEvent >= labelEvent) continue;
        if (labelEvent != ient) return ibc;
        prevLabelEvent = labelEvent;

        double charge;
        int channel;
        const auto &bcd = fv0BCData[ibc];
        int chEnt = bcd.ref.getFirstEntry();
        for (int ich = 0; ich < bcd.ref.getEntries(); ich++) { // Go through FV0 channels
            const auto &chd = fv0ChData[chEnt++];
            charge = chd.chargeAdc;
            channel = chd.pmtNumber;
            SumQvec(Qvec, charge, channel, "FV0");
            signalSum += charge;
        }

        if (signalSum != 0) {
            Qvec /= signalSum;
            qvecFV0[0] = Qvec.Re();
            qvecFV0[1] = Qvec.Im();
        } else {
            qvecFV0[0] = 0.;
            qvecFV0[1] = 0.;
        }
    }
    return -1;
}

int FillQvecFT0(UInt_t ient, int bcbegin)
{
    fFT0DigitTree->GetEntry(ient);

    TComplex QvecFT0A(0), QvecFT0C(0);
    double signalSumFT0A = 0, signalSumFT0C = 0;
    int prevLabelEvent = -1;
    int nbc = ft0BCData.size();
    for (int ibc = bcbegin; ibc < nbc; ibc++) {

        const auto lb = ft0labels.getLabels(ibc);
        int labelEvent = 0;
        if (lb.size() > 0) {
            labelEvent = lb[0].getEventID();
        } else {
            std::cout << "Problem with labels" << std::endl;
            continue;
        }

        if (labelEvent == prevLabelEvent) continue;
        if (labelEvent != ient) return ibc;
        prevLabelEvent = labelEvent;

        double charge;
        int channel;
        const auto &bcd = ft0BCData[ibc];
        int chEnt = bcd.ref.getFirstEntry();
        for (int ich = 0; ich < bcd.ref.getEntries(); ich++) { // Go through FT0 channels
            const auto &chd = ft0ChData[chEnt++];
            charge = chd.QTCAmpl;
            channel = chd.ChId;
            if (channel < 96) {
                SumQvec(QvecFT0A, charge, channel, "FT0A");
                signalSumFT0A += charge;
            } else {
                SumQvec(QvecFT0C, charge, channel, "FT0C");
                signalSumFT0C += charge;
            }
        }

        if (signalSumFT0A != 0) {
            QvecFT0A /= signalSumFT0A;
            qvecFT0A[0] = QvecFT0A.Re();
            qvecFT0A[1] = QvecFT0A.Im();
        } else {
            qvecFT0A[0] = 0.;
            qvecFT0A[1] = 0.;
        }

        if (signalSumFT0C != 0) {
            QvecFT0C /= signalSumFT0C;
            qvecFT0C[0] = QvecFT0C.Re();
            qvecFT0C[1] = QvecFT0C.Im();
        } else {
            qvecFT0C[0] = 0.;
            qvecFT0C[1] = 0.;
        }
    }
    return -1;
}


void SumQvec(TComplex &Qvec, double nch, int chno, TString det)
{
    double phi = 0.0;
    if (det == "FV0") {
        phi = GetFV0Phi(chno);
    } else if (det == "FT0A") {
        phi = GetFT0APhi(chno);
    } else if (det == "FT0C") {
        phi = GetFT0CPhi(chno-96);
    } else {
        std::cout << "Eventplane::SumQvec : warning, phi not taken from any detector" << std::endl;
    }
    Qvec += TComplex(nch*TMath::Cos(2.0 * phi), nch*TMath::Sin(2.0 * phi));
}

double GetFV0Phi(int chno)
{
    const double pi = TMath::Pi();
    int pmtMapSmallCh[8] = {5,4,3,2,6,7,0,1};
    std::map<int, int> pmtMapBigCh = {
        {46, 0}, {38, 1}, {47, 2}, {39, 3}, {43, 4}, {35, 5}, {42, 6}, {34, 7},
        {41, 8}, {33, 9}, {40, 10}, {32, 11}, {44, 12}, {36, 13}, {45, 14}, {37, 15}
    };

    if (chno>31) { // 5th ring has 16 segments
        return pi/16.0 + (pmtMapBigCh[chno])*pi/8.0 - pi;
    } else {
        return pi/8.0 + (pmtMapSmallCh[chno%8])*pi/4.0 - pi;
    }
}

double GetFT0APhi(int chno)
{
    return TMath::ATan2(ft0ay[chno], ft0ax[chno]);
}

double GetFT0CPhi(int chno)
{
    return TMath::ATan2(ft0cy[chno], ft0cx[chno]);
}