#include "EventPlaneRes.h"

void EventPlaneRes(TString infile, TString outfile)
{
    int isOpen = LoadInput(infile);
    if (!isOpen) return;
    InitOutput(outfile);

    int nev = inTree->GetEntriesFast();
    for (int iev=0; iev<nev; iev++) {

        inTree->GetEntry(iev);

        epB = GetEventPlane(qvecB[0], qvecB[1]);
        epC = GetEventPlane(qvecC[0], qvecC[1]);
        for (int idet=0; idet<ndet; idet++) {
            hQvec[idet]->Fill(qvecA[idet][0], qvecA[idet][1]);
            epA[idet] = GetEventPlane(qvecA[idet][0], qvecA[idet][1]);
        }

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
    for (int idet=0; idet<ndet; idet++)
        inTree->SetBranchAddress(Form("qvec%s",detName[idet].Data()), &qvecA[idet]);
    inTree->SetBranchAddress("qvecB", &qvecB);
    inTree->SetBranchAddress("qvecC", &qvecC);

    return 1;
}

void InitOutput(TString outfile)
{
    fout = TFile::Open(outfile, "RECREATE");
    hEPB = new TH1D("hEPB", "hEPB", 200, -TMath::Pi()/2.0, TMath::Pi()/2.0);
    hEPC = new TH1D("hEPC", "hEPC", 200, -TMath::Pi()/2.0, TMath::Pi()/2.0);
    hRBC = new TH1D("hRBC", "hRBC", 400, -1.0, 1.0);
    for (int idet=0; idet<ndet; idet++) {
        hEPA[idet] = new TH1D(Form("hEPA_%s", detName[idet].Data()), Form("hEPA_%s", detName[idet].Data()), 200, -TMath::Pi()/2.0, TMath::Pi()/2.0);
        hRAB[idet] = new TH1D(Form("hRAB_%s", detName[idet].Data()), Form("hRAB_%s", detName[idet].Data()), 400, -1.0, 1.0);
        hRAC[idet] = new TH1D(Form("hRAC_%s", detName[idet].Data()), Form("hRAC_%s", detName[idet].Data()), 400, -1.0, 1.0);
        hQvec[idet] = new TH2D(Form("hQvec_%s", detName[idet].Data()), Form("hQvec_%s", detName[idet].Data()), 400, -0.02, 0.02, 400, -0.02, 0.02);
    }
}

float GetEventPlane(float qx, float qy)
{
    return TMath::ATan2(qy, qx)/2.0;
}

void FillHistos()
{
    for (int idet=0; idet<ndet; idet++) {
        hEPA[idet]->Fill(epA[idet]);
        hRAB[idet]->Fill(TMath::Cos(2*(epA[idet] - epB)));
        hRAC[idet]->Fill(TMath::Cos(2*(epA[idet] - epC)));
    }
    hRBC->Fill(TMath::Cos(2*(epB - epC)));
    hEPB->Fill(epB);
    hEPC->Fill(epC);
}
