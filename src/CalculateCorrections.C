#include "CalculateCorrections.h"
#include "Corrections.h"

void CalculateCorrections(TString infile, TString cent)
{
    int isOpen = LoadInput(infile);
    if (!isOpen) return;
    InitOutput(cent);
    FillHistos();
    SaveCorrections();
}

//_____________________________________________________________
void InitOutput(TString cent)
{
    for (int idet=0; idet<ndet; idet++)
        saveFileName[idet] = Form("corr_%s_%s.txt", detName[idet].Data(), cent.Data());
}

int LoadInput(TString infile)
{
    fin = TFile::Open(infile, "READ");
    if (!fin) {
        std::cout << "Could not open input file, return" << std::endl;
        return 0;
    }
    inTree = (TTree*)fin->Get("Qvecs");
    for (int idet=0; idet<ndet; idet++)
        inTree->SetBranchAddress(Form("qvec%s",detName[idet].Data()), &qvec[idet]);
    return 1;
}

void FillHistos()
{
    for (int idet=0; idet<ndet; idet++) {
        hQvec[idet] = new TH2D(Form("hQvec_%s", detName[idet].Data()), Form("hQvec_%s", detName[idet].Data()), nbin, binmin, binmax, nbin, binmin, binmax);
    }

    int nev = inTree->GetEntriesFast();
    for (int iev=0; iev<nev; iev++) {
        inTree->GetEntry(iev);
        for (int idet=0; idet<ndet; idet++)
            hQvec[idet]->Fill(qvec[idet][0], qvec[idet][1]);
    }
}

void SaveCorrections()
{
    for (int idet=0; idet<ndet; idet++) {
        std::cout << "\nCalculating corrections.." << std::endl;
        corr = CalculateCorrections(hQvec[idet]);
        outputFile.open(saveFileName[idet].Data());
        for (int i = 0; i < nCorrections; i++) {
            outputFile << corr[i] << "\n";
        }
        outputFile.close();
        std::cout << "Corrections saved to file " << saveFileName[idet] << std::endl;
    }
}

void Print()
{
    std::cout << "<Q>      : " << corr[0] << "  " << corr[1] << std::endl;
    std::cout << "a+-      : " << corr[4] << "  " << corr[5] << std::endl;
    std::cout << "lambda+- : " << corr[6] << "  " << corr[7] << std::endl;
}
