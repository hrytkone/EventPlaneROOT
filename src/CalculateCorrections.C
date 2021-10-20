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

void GetCorrections(TH2D *hQ)
{
    float xmean, ymean, xdev, ydev, aplus, aminus, lambdaplus, lambdaminus;
    GetRecenteringCorrection(hQ, xmean, ymean);
    GetWidthCorrection(hQ, xdev, ydev);
    GetTwistAndRescaleCorrection(hQ, aplus, aminus, lambdaplus, lambdaminus);
    corrections[0] = xmean;
    corrections[1] = ymean;
    corrections[2] = xdev;
    corrections[3] = ydev;
    corrections[4] = aplus;
    corrections[5] = aminus;
    corrections[6] = lambdaplus;
    corrections[7] = lambdaminus;
}

void SaveCorrections()
{
    for (int idet=0; idet<ndet; idet++) {
        std::cout << "\nCalculating corrections.." << std::endl;
        GetCorrections(hQvec[idet]);
        outputFile.open(saveFileName[idet].Data());
        for (int i = 0; i < nCorrections; i++) {
            outputFile << corrections[i] << "\n";
        }
        outputFile.close();
        std::cout << "Corrections saved to file " << saveFileName[idet] << std::endl;
    }
}

void Print()
{
    std::cout << "<Q>      : " << corrections[0] << "  " << corrections[1] << std::endl;
    std::cout << "a+-      : " << corrections[4] << "  " << corrections[5] << std::endl;
    std::cout << "lambda+- : " << corrections[6] << "  " << corrections[7] << std::endl;
}
