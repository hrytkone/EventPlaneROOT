const int ncent = 1;
const int ndet = 3;
const int nx = 4, ny = 2;

TString files[ncent] = {"../cent20-30.root"};
TString cent[ncent] = {"20-30 %"};
TString detname[ndet] = {"FV0", "FT0A", "FT0C"};

TFile *fin[ncent];
TH1D *hEPA[ncent][ndet];
TCanvas *cEP[ndet];

void LoadData();
void ConfigPlots();
void ConfigCanvas();

void PlotEventPlanes()
{
    gStyle->SetOptStat(0);
    gStyle->SetLineScalePS(1);

    LoadData();
    ConfigPlots();
    ConfigCanvas();

    for (int idet=0; idet<ndet; idet++) {
        cEP[idet]->cd();
        for (int icent=0; icent<ncent; icent++) {
            hEPA[icent][idet]->Draw();
        }
    }
}

void LoadData()
{
    for (int icent=0; icent<ncent; icent++) {
        fin[icent] = TFile::Open(files[icent]);

        for (int idet=0; idet<ndet; idet++) {
            hEPA[icent][idet] = (TH1D*)fin[icent]->Get(Form("hEPA_%s", detname[idet].Data()));
            hEPA[icent][idet]->Scale(1./hEPA[icent][idet]->GetEntries(), "width");
        }
    }
}

void ConfigPlots()
{
    for (int icent=0; icent<ncent; icent++) {
        for (int idet=0; idet<ndet; idet++) {
            hEPA[icent][idet]->GetXaxis()->SetLabelFont(53);
            hEPA[icent][idet]->GetXaxis()->SetLabelSize(10);
            //hEPA[icent][idet]->Rebin(2);
            //hEPA[icent][idet]->GetXaxis()->SetRangeUser(-10., 100.);
            //hEPA[icent][idet]->GetYaxis()->SetRangeUser(0.1, 1000000.);
            hEPA[icent][idet]->GetYaxis()->SetLabelFont(53);
            hEPA[icent][idet]->GetYaxis()->SetLabelSize(10);
            hEPA[icent][idet]->SetLineColor(kBlack);
            hEPA[icent][idet]->SetTitle(Form("%s;;", detname[idet].Data()));
        }
    }
}

void ConfigCanvas()
{
    for (int idet=0; idet<ndet; idet++) {
        cEP[idet] = new TCanvas(Form("cEP%d", idet), Form("cEP%d", idet), 800, 800);
        //cEP[idet]->Divide(nx, ny, 0, 0);
    }
}
