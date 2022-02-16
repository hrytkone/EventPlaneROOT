const int ncent = 8;
const int ndet = 3;

TString colsys = "Pb#font[122]{-}Pb";
TString energy = "5.5 TeV";

float cent[ncent+1] = {0., 5., 10., 20., 30., 40., 50., 60., 80.};
TString dir = "data-v2_2022-02-09/res_qvecs-corr";
TString files[ncent] = {"cent00-05.root", "cent05-10.root", "cent10-20.root",
                        "cent20-30.root", "cent30-40.root", "cent40-50.root",
                        "cent50-60.root", "cent60-80.root"};
TString detname[ndet] = {"FT0C", "FV0", "FT0A"};
TString legentry[ndet] = {"FT0-C", "FV0", "FT0-A"};
EColor mColor[ndet] = {kRed, kBlue, kBlack};

TFile *fin[ncent];
TH1D *hVnObs[ncent][ndet];
TH1D *hRAB[ncent][ndet];
TH1D *hRAC[ncent][ndet];
TH1D *hRBC[ncent];
TGraphErrors *gV2[ndet];
TLegend *leg;
TCanvas *c1;

float res[ndet][ncent];
float resErr[ndet][ncent];
float v2[ndet][ncent];
float v2Err[ndet][ncent];

void LoadData();
void SetStyle(Bool_t graypalette);
void GetResAndErr(int icent, int idet);
void GetV2AndErr(int icent, int idet);
void ConfigPlots();
void PlotToCanvas();

void PlotV2()
{
    SetStyle(0);
    LoadData();

    for (int idet=0; idet<ndet; idet++) {
        gV2[idet] = new TGraphErrors();
        for (int icent=0; icent<ncent; icent++) {
            GetResAndErr(icent, idet);
            GetV2AndErr(icent, idet);
            gV2[idet]->SetPoint(icent, cent[icent] + (cent[icent+1]-cent[icent])/2., v2[idet][icent]);
            gV2[idet]->SetPointError(icent, (cent[icent+1]-cent[icent])/2., v2Err[idet][icent]);
        }
    }

    ConfigPlots();
    PlotToCanvas();
}

//________________________________
void LoadData()
{
    for (int icent=0; icent<ncent; icent++) {
        fin[icent] = TFile::Open(Form("%s/%s", dir.Data(), files[icent].Data()));

        hRBC[icent] = (TH1D*)fin[icent]->Get("hRBC");
        for (int idet=0; idet<ndet; idet++) {
            hRAB[icent][idet] = (TH1D*)fin[icent]->Get(Form("hRAB_%s", detname[idet].Data()));
            hRAC[icent][idet] = (TH1D*)fin[icent]->Get(Form("hRAC_%s", detname[idet].Data()));
            hVnObs[icent][idet] = (TH1D*)fin[icent]->Get(Form("hVnObs_%s", detname[idet].Data()));
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

void GetV2AndErr(int icent, int idet)
{
    float vnObs = hVnObs[icent][idet]->GetMean();
    float vnObsErr = hVnObs[icent][idet]->GetMeanError();

    v2[idet][icent] = vnObs/res[idet][icent];
    v2Err[idet][icent] = v2[idet][icent]*TMath::Sqrt((resErr[idet][icent]/res[idet][icent])*(resErr[idet][icent]/res[idet][icent])
                        + (vnObsErr/vnObs)*(vnObsErr/vnObs));
}

void ConfigPlots()
{
    leg = new TLegend(0.4, 0.23, 0.52, 0.38);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.035);
    //leg->SetHeader(Form("%s #sqrt{#it{s}_{NN}} = %s simulation", colsys.Data(), energy.Data()));
    for (int idet=0; idet<ndet; idet++) {
        leg->AddEntry(gV2[idet], legentry[idet], "p");
        gV2[idet]->SetTitle("; centrality (%); v_{2}");
        gV2[idet]->GetXaxis()->SetLabelSize(0.035);
        gV2[idet]->GetXaxis()->SetTitleSize(0.04);
        gV2[idet]->GetXaxis()->SetLimits(0., 80.);
        gV2[idet]->GetYaxis()->SetMaxDigits(1);
        gV2[idet]->GetYaxis()->SetLabelSize(0.035);
        gV2[idet]->GetYaxis()->SetTitleSize(0.04);
        //gV2[idet]->GetYaxis()->SetRangeUser(0., 1.);
        gV2[idet]->SetMarkerColor(mColor[idet]);
        gV2[idet]->SetMarkerStyle(20);
    }
}

void PlotToCanvas()
{
    c1 = new TCanvas("c1", "c1", 800, 800);
    for (int idet=0; idet<ndet; idet++) {
        if (idet==0)
            gV2[idet]->Draw("AP");
        else
            gV2[idet]->Draw("P SAME");
    }
    leg->Draw("SAME");

    //gV2V0A->Draw("P SAME");

    TLatex * tex = new TLatex(0.1,0.1, "ALICE Preliminary");
    tex->SetNDC();
    tex->SetTextFont(42);
    //tex->Draw();

    TLatex * texsim1 = new TLatex(0.402,0.47, Form("%s #sqrt{#it{s}_{NN}} = %s simulation", colsys.Data(), energy.Data()));
    texsim1->SetNDC();
    texsim1->SetTextSize(0.035);
    texsim1->Draw();

    TLatex * texsim2 = new TLatex(0.402,0.43, "-0.8 < #eta < 0.8");
    texsim2->SetNDC();
    texsim2->SetTextSize(0.035);
    texsim2->Draw();

    TLatex * texsim3 = new TLatex(0.402,0.39, "0.2 GeV/c < p_{T} < 5.0 GeV/c");
    texsim3->SetNDC();
    texsim3->SetTextSize(0.035);
    texsim3->Draw();
}
