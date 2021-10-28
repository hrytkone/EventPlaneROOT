const int nset = 1;
const int ncent = 1;
const int ndet = 3;
const int nx = 4, ny = 2;

//TString files[nset][ncent] = {{"../cent20-30_res.root"}, {"../cent20-30_recenter_res.root"}, {"../cent20-30_corr_res.root"}};
TString files[ncent][ncent] = {{"../ep_cent20-30_276TeV_corrected.root"}};
TString cent[ncent] = {"20-30 %"};
TString detname[ndet] = {"FV0", "FT0A", "FT0C"};
TString legEntry[nset] = {"no corrections"};//, "recentering", "all corrections"};
EColor mColor[nset] = {kBlue};//, kRed, kBlack};

TFile *fin[nset][ncent];
TH1D *hEPA[nset][ncent][ndet];
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
            for (int iset=0; iset<nset; iset++) {
                hEPA[iset][icent][idet]->Draw("HIST SAME");
            }
        }
    }
}

void LoadData()
{
    for (int iset=0; iset<nset; iset++) {
        for (int icent=0; icent<ncent; icent++) {
            fin[iset][icent] = TFile::Open(files[iset][icent]);
            for (int idet=0; idet<ndet; idet++) {
                hEPA[iset][icent][idet] = (TH1D*)fin[iset][icent]->Get(Form("hEPA_%s", detname[idet].Data()));
                hEPA[iset][icent][idet]->Rebin(4);
                hEPA[iset][icent][idet]->Scale(1./hEPA[iset][icent][idet]->GetEntries(), "width");
            }
        }
    }
}

void ConfigPlots()
{
    for (int iset=0; iset<nset; iset++) {
        for (int icent=0; icent<ncent; icent++) {
            for (int idet=0; idet<ndet; idet++) {
                //hEPA[iset][icent][idet]->GetXaxis()->SetLabelFont(53);
                hEPA[iset][icent][idet]->GetXaxis()->SetLabelSize(0.025);
                //hEPA[icent][idet]->GetXaxis()->SetRangeUser(-10., 100.);
                //hEPA[iset][iset][icent][idet]->Rebin(2);
                hEPA[iset][icent][idet]->GetYaxis()->SetRangeUser(0., .75);
                //hEPA[iset][icent][idet]->GetYaxis()->SetLabelFont(53);
                hEPA[iset][icent][idet]->GetYaxis()->SetLabelSize(0.025);
                hEPA[iset][icent][idet]->SetLineColor(mColor[iset]);
                hEPA[iset][icent][idet]->SetTitle(Form("%s;#Psi_{2}; 1/N dN/d#Psi_{2}", detname[idet].Data()));
            }
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
