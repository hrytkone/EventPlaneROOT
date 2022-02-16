#include "EventPlaneResTrue.h"

void EventPlaneResTrue(TString infile, TString outfile)
{
    gROOT->ProcessLine("#include <vector>");

    int isOpen = LoadInput(infile);
    if (!isOpen) return;
    InitOutput(outfile);

    int nev = inTree->GetEntriesFast();
    for (int iev=0; iev<nev; iev++) {

        inTree->GetEntry(iev);

        epB    = GetEventPlane(qvecB[0], qvecB[1]);
        epC    = GetEventPlane(qvecC[0], qvecC[1]);
        epFull = GetEventPlane(qvecFull[0], qvecFull[1]);

        FillHistos();
    }
    fout->Write("", TObject::kOverwrite);
}

//_____________________________________________________________
int LoadInput(TString infile)
{
    fin = TFile::Open(infile, "READ");
    if (!fin) {
        std::cout << "Could not open input file, return" << std::endl;
        return 0;
    }
    inTree = (TTree*)fin->Get("Qvecs");
    inTree->SetBranchAddress("qvecB", &qvecB);
    inTree->SetBranchAddress("qvecC", &qvecC);
    inTree->SetBranchAddress("qvecFull", &qvecFull);
    inTree->SetBranchAddress("tpcPhi", &tpcPhi);

    return 1;
}

void InitOutput(TString outfile)
{
    fout    = TFile::Open(outfile, "RECREATE");
    hEPB    = new TH1D("hEPB", "hEPB", 200, -TMath::Pi()/2.0, TMath::Pi()/2.0);
    hEPC    = new TH1D("hEPC", "hEPC", 200, -TMath::Pi()/2.0, TMath::Pi()/2.0);
    hEPFull = new TH1D("hEPFull", "hEPFull", 200, -TMath::Pi()/2.0, TMath::Pi()/2.0);
    hRsub   = new TH1D("hRsub", "hRsub", 400, -1.0, 1.0);
    hVnObs  = new TH1D(Form("hVnObs_%s", detName[idet].Data()), Form("hVnObs_%s", detName[idet].Data()), 400, -1.0, 1.0);
}

float GetEventPlane(float qx, float qy)
{
    return TMath::ATan2(qy, qx)/2.0;
}

float GetVnObs(float ep, int n)
{
    int ntrack = tpcPhi->size();
    float vn = 0.;
    for (int itrack=0; itrack<ntrack; itrack++) {
        vn += TMath::Cos(n*(tpcPhi->at(itrack) - ep));
    }
    return vn/ntrack;
}

void FillHistos()
{
    // Calculate & save v2obs for validating the correction with EP res
    hVnObs->Fill(GetVnObs(epFull, 2.));
    hRsub->Fill(TMath::Cos(2*(epB - epC)));
    hEPB->Fill(epB);
    hEPC->Fill(epC);
    hEPFull->Fill(epFull);
}
