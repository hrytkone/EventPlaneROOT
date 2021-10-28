const int ncent = 1;
const int ndet = 3;

TString colsys = "Pb#font[122]{-}Pb";
TString energy = "5.5 TeV";

float cent[ncent+1] = {20., 30.};
//TString files[ncent] = {"../res_cent20-30_tpc-eff_corr.root"};
TString files[ncent] = {"../cent20-30_corr_res.root"};
TString detname[ndet] = {"FT0C", "FV0", "FT0A"};
TString legentry[ndet] = {"FT0-C", "FV0", "FT0-A"};
EColor mColor[ndet] = {kRed, kBlue, kBlack};

TFile *fin[ncent];
TH1D *hRAB[ncent][ndet];
TH1D *hRAC[ncent][ndet];
TH1D *hRBC[ncent];
TGraphErrors *gRes[ndet];
TLegend *leg;
TCanvas *c1;

float res[ndet][ncent];
float resErr[ndet][ncent];

void LoadData();
void SetStyle(Bool_t graypalette);
void GetResAndErr(int icent, int idet);
void ConfigPlots();
void PlotToCanvas();

void PlotRes()
{
    SetStyle(0);
    LoadData();


    for (int idet=0; idet<ndet; idet++) {
        gRes[idet] = new TGraphErrors();
        for (int icent=0; icent<ncent; icent++) {
            GetResAndErr(icent, idet);
            gRes[idet]->SetPoint(icent, cent[icent] + (cent[icent+1]-cent[icent])/2., res[idet][icent]);
            gRes[idet]->SetPointError(icent, (cent[icent+1]-cent[icent])/2., resErr[idet][icent]);
        }
    }
    ConfigPlots();
    PlotToCanvas();
}

//________________________________
void LoadData()
{
    for (int icent=0; icent<ncent; icent++) {
        fin[icent] = TFile::Open(files[icent]);

        hRBC[icent] = (TH1D*)fin[icent]->Get("hRBC");
        for (int idet=0; idet<ndet; idet++) {
            hRAB[icent][idet] = (TH1D*)fin[icent]->Get(Form("hRAB_%s", detname[idet].Data()));
            hRAC[icent][idet] = (TH1D*)fin[icent]->Get(Form("hRAC_%s", detname[idet].Data()));
        }
    }
}

void SetStyle(Bool_t graypalette)
{
    cout << "Setting style!" << endl;

    gStyle->Reset("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //gStyle->SetLineScalePS(1);
    if(graypalette) gStyle->SetPalette(8,0);
    else gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    gStyle->SetPadTickX(0);
    gStyle->SetPadTickY(0);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetHistLineWidth(1);
    gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(2);
    gStyle->SetFuncColor(kGreen);
    gStyle->SetLineWidth(1);
    gStyle->SetLabelSize(0.035,"xyz");
    gStyle->SetLabelOffset(0.01,"y");
    gStyle->SetLabelOffset(0.01,"x");
    gStyle->SetLabelColor(kBlack,"xyz");
    gStyle->SetTitleSize(0.035,"xyz");
    gStyle->SetTitleOffset(1.25,"y");
    gStyle->SetTitleOffset(1.2,"x");
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetTextSizePixels(26);
    gStyle->SetTextFont(42);
    //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    //  gStyle->SetFillColor(kWhite);
    gStyle->SetLegendFont(42);
}

void GetResAndErr(int icent, int idet)
{
    float rab = hRAB[icent][idet]->GetMean();
    float rac = hRAC[icent][idet]->GetMean();
    float rbc = hRBC[icent]->GetMean();
    float rabErr = hRAB[icent][idet]->GetMeanError();
    float racErr = hRAC[icent][idet]->GetMeanError();
    float rbcErr = hRBC[icent]->GetMeanError();

    res[idet][icent] = TMath::Sqrt((rab * rac)/rbc);
    resErr[idet][icent] = 0.5*res[idet][icent]*TMath::Sqrt(TMath::Power(rabErr/rab, 2) + TMath::Power(racErr/rac, 2) + TMath::Power(rbcErr/rbc, 2));

    cout << "\tRes : " << res[idet][icent] << " +- " << resErr[idet][icent] << endl;
}

void ConfigPlots()
{
    leg = new TLegend(0.6, 0.63, 0.72, 0.79);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(gStyle->GetTextSize()*0.5);
    leg->SetHeader(Form("%s #sqrt{#it{s}_{NN}} = %s", colsys.Data(), energy.Data()));
    for (int idet=0; idet<ndet; idet++) {
        leg->AddEntry(gRes[idet], legentry[idet], "p");
        gRes[idet]->SetTitle("; centrality (%); R_{2}");
        //gRes[idet]->GetXaxis()->SetLabelSize(0.025);
        //gRes[idet]->GetXaxis()->SetTitleSize(0.025);
        //gRes[idet]->GetXaxis()->CenterTitle();
        gRes[idet]->GetXaxis()->SetLimits(0., 80.);
        //gRes[idet]->GetYaxis()->SetLabelSize(0.025);
        //gRes[idet]->GetYaxis()->SetTitleSize(0.025);
        //gRes[idet]->GetYaxis()->CenterTitle();
        gRes[idet]->GetYaxis()->SetRangeUser(0., 1.);
        gRes[idet]->SetMarkerColor(mColor[idet]);
        gRes[idet]->SetMarkerStyle(20);
    }
}

void PlotToCanvas()
{
    c1 = new TCanvas("c1", "c1", 800, 800);
    for (int idet=0; idet<ndet; idet++) {
        if (idet==0)
            gRes[idet]->Draw("AP");
        else
            gRes[idet]->Draw("P SAME");
    }
    leg->Draw("SAME");

    TLatex * tex = new TLatex(0.1,0.1, "ALICE Preliminary");
    tex->SetNDC();
    tex->SetTextFont(42);
    //tex->Draw();

    TLatex * texsim = new TLatex(0.602,0.8, "Simulation");
    texsim->SetNDC();
    texsim->SetTextSize(gStyle->GetTextSize()*0.5);
    texsim->Draw();
}
