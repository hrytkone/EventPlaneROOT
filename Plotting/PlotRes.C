const int ncent = 8;
const int ndet = 3;

TString colsys = "Pb#font[122]{-}Pb";
TString energy = "5.5 TeV";

float cent[ncent+1] = {0., 5., 10., 20., 30., 40., 50., 60., 80.};
//float cent[ncent+1] = {20., 30.};
//TString files[ncent] = {"../res_cent20-30_tpc-eff_corr.root"};
//TString dir = "2022-03-10_deadch-50/res_qvecs-corr";
//TString dir = "2022-05-10_v2-true-corrected/res_qvecs-corr";
//TString dir = "data-v2_2022-02-09/res_qvecs-corr";

TString dir = "20220516_data/res_qvecs-corr";
//TString dirnocorr = "20220517_data-mctracks/res_qvecs-mctrack";
TString dirnocorr = "20220516_data/res_qvecs-no-corr";
TString files[ncent] = {"cent00-05.root", "cent05-10.root", "cent10-20.root",
                        "cent20-30.root", "cent30-40.root", "cent40-50.root",
                        "cent50-60.root", "cent60-80.root"};
//TString files[ncent] = {"cent20-30.root"};
TString detname[ndet] = {"FT0C", "FV0", "FT0A"};
TString legentry[ndet] = {"FT0-C", "FV0", "FT0-A"};
int mColor[ndet] = {kRed+1, kBlack, kBlue};
int mStyle[ndet] = {kOpenCircle, kFullCircle, kOpenSquare};
int mFillStyle[ndet] = {3001, 3002, 3354};

TFile *fin[ncent];
TFile *fin2[ncent];
TH1D *hRAB[ncent][ndet];
TH1D *hRAC[ncent][ndet];
TH1D *hRBC[ncent];
TH1D *hRABnocorr[ncent][ndet];
TH1D *hRACnocorr[ncent][ndet];
TH1D *hRBCnocorr[ncent];
TH1D *hRsub[ncent];
TGraphErrors *gRes[ndet];
TGraphErrors *gResNocorr[ndet];
TGraphErrors *gResV0A;
TLegend *leg1;
TLegend *leg2;
TCanvas *c1;

float res[ndet][ncent];
float resErr[ndet][ncent];

void LoadData();
void SetStyle(Bool_t graypalette);
void GetResAndErr(int icent, int idet);
void GetResAndErrNocorr(int icent, int idet);
void ConfigPlots();
void PlotToCanvas();
void RedrawBorder();

float offset[ncent][2] = {
    {-.25,.25},
    {-.25,.25},
    {-.25,.25},
    {-.25,.25},
    {-.25,.25},
    {-.25,.25},
    {-.25,.25},
    {-.25,.25}
};

void PlotRes()
{
    SetStyle(0);
    LoadData();

    for (int idet=0; idet<ndet; idet++) {
        gRes[idet] = new TGraphErrors();
        cout << detname[idet] << endl;
        for (int icent=0; icent<ncent; icent++) {
            GetResAndErr(icent, idet);
            gRes[idet]->SetPoint(icent, cent[icent] + (cent[icent+1]-cent[icent])/2., res[idet][icent]);
            gRes[idet]->SetPointError(icent, (cent[icent+1]-cent[icent])/2., resErr[idet][icent]);
            //gRes[idet]->SetPointError(icent, 0., resErr[idet][icent]);
        }

        cout << "\nNo corrections" << endl;
        gResNocorr[idet] = new TGraphErrors();
        for (int icent=0; icent<ncent; icent++) {
            GetResAndErrNocorr(icent, idet);
            gResNocorr[idet]->SetPoint(icent, cent[icent] + (cent[icent+1]-cent[icent])/2., res[idet][icent]);
            gResNocorr[idet]->SetPointError(icent, (cent[icent+1]-cent[icent])/2., resErr[idet][icent]);
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
        }
        //hRsub[icent]      = (TH1D*)fin[icent]->Get("hRsub");
        //hRAB[icent][ndet] = (TH1D*)fin[icent]->Get("hRABtrue");
        //hRAC[icent][ndet] = (TH1D*)fin[icent]->Get("hRACtrue");

        fin2[icent]    = TFile::Open(Form("%s/%s", dirnocorr.Data(), files[icent].Data()));
        hRBCnocorr[icent] = (TH1D*)fin2[icent]->Get("hRBC");
        for (int idet=0; idet<ndet; idet++) {
            hRABnocorr[icent][idet] = (TH1D*)fin2[icent]->Get(Form("hRAB_%s", detname[idet].Data()));
            hRACnocorr[icent][idet] = (TH1D*)fin2[icent]->Get(Form("hRAC_%s", detname[idet].Data()));
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
    gStyle->SetFrameLineWidth(2);
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
    gStyle->SetLineWidth(2);
    gStyle->SetLabelSize(0.042,"xyz");
    gStyle->SetLabelOffset(0.01,"y");
    gStyle->SetLabelOffset(0.01,"x");
    gStyle->SetLabelColor(kBlack,"xyz");
    gStyle->SetTitleSize(0.042,"xyz");
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
    gStyle->SetHatchesSpacing(0.5);
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

void GetResAndErrNocorr(int icent, int idet)
{
    float rab = hRABnocorr[icent][idet]->GetMean();
    float rac = hRACnocorr[icent][idet]->GetMean();
    float rbc = hRBCnocorr[icent]->GetMean();
    float rabErr = hRABnocorr[icent][idet]->GetMeanError();
    float racErr = hRACnocorr[icent][idet]->GetMeanError();
    float rbcErr = hRBCnocorr[icent]->GetMeanError();

    res[idet][icent] = TMath::Sqrt((rab * rac)/rbc);
    resErr[idet][icent] = 0.5*res[idet][icent]*TMath::Sqrt(TMath::Power(rabErr/rab, 2) + TMath::Power(racErr/rac, 2) + TMath::Power(rbcErr/rbc, 2));

    cout << "\tRes : " << res[idet][icent] << " +- " << resErr[idet][icent] << endl;
}

void ConfigPlots()
{
    leg1 = new TLegend(0.23, 0.23, 0.32, 0.39);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(gStyle->GetTextSize()*0.7);
    //leg1->SetHeader("Digits");
    leg1->SetHeader("Corrected");
    for (int idet=0; idet<ndet; idet++) {
        leg1->AddEntry(gRes[idet], legentry[idet], "pe");
        gRes[idet]->SetTitle("; centrality (%); R_{2}");
        gRes[idet]->GetXaxis()->SetLimits(0., 80.);
        gRes[idet]->GetYaxis()->SetRangeUser(0., 1.);
        gRes[idet]->SetMarkerColor(mColor[idet]);
        gRes[idet]->SetLineColor(mColor[idet]);
        //gRes[idet]->SetLineWidth(1.);
        gRes[idet]->SetMarkerStyle(mStyle[idet]);
    }

    leg2 = new TLegend(0.43, 0.23, 0.52, 0.39);
    leg2->SetFillStyle(0); leg2->SetBorderSize(0); leg2->SetTextSize(gStyle->GetTextSize()*0.7);
    //leg2->SetHeader("MC tracks");
    leg2->SetHeader("Raw");
    for (int idet=0; idet<ndet; idet++) {
        leg2->AddEntry(gResNocorr[idet], legentry[idet], "f");
        gResNocorr[idet]->SetTitle("; centrality (%); R_{2}");
        gResNocorr[idet]->GetXaxis()->SetLimits(0., 80.);
        gResNocorr[idet]->GetYaxis()->SetRangeUser(0., 1.);
        gResNocorr[idet]->SetLineWidth(1.);
        //gResNocorr[idet]->SetFillStyle(mFillStyle[idet]);
        gResNocorr[idet]->SetLineColor(mColor[idet]-10);
        if (idet==1) {
            gResNocorr[idet]->SetFillColor(kGray);
            gResNocorr[idet]->SetMarkerColor(kGray+1);
        } else {
            gResNocorr[idet]->SetFillColor(mColor[idet]-10);
            gResNocorr[idet]->SetMarkerColor(mColor[idet]);
        }
        gResNocorr[idet]->SetMarkerStyle(kDot);
    }
    //gResV0A->SetMarkerColor(kBlue);
    //gResV0A->SetMarkerStyle(24);
    //leg->AddEntry(gResV0A, "VZERO-A", "p");
}

void PlotToCanvas()
{
    c1 = new TCanvas("c1", "c1", 800, 800);
    for (int idet=0; idet<ndet; idet++) {
        if (idet==0) {
            gResNocorr[idet]->Draw("A2");
            //gResNocorr[idet]->Draw("P");
        } else {
            gResNocorr[idet]->Draw("2");
            //gResNocorr[idet]->Draw("P SAME");
        }
    }

    //gRes[0]->Draw("AP");
    for (int idet=0; idet<ndet; idet++)
        gRes[idet]->Draw("P SAME");

    leg1->Draw("SAME");
    leg2->Draw("SAME");

    //gResV0A->Draw("P SAME");

    TLatex * texsim = new TLatex(0.202,0.41, Form("Simulation, %s #sqrt{#it{s}_{NN}} = %s", colsys.Data(), energy.Data()));
    texsim->SetNDC();
    texsim->SetTextSize(gStyle->GetTextSize()*0.7);
    texsim->Draw();

    RedrawBorder();
}

void RedrawBorder()
{
   gPad->Update();
   gPad->RedrawAxis();
   TLine l;
   l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
   l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());
}
