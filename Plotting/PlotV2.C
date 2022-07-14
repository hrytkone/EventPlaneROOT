const int ncent = 8;
const int ndet = 4; // 4th = true

TString colsys = "Pb#font[122]{-}Pb";
TString energy = "5.5 TeV";

float cent[ncent+1] = {0., 5., 10., 20., 30., 40., 50., 60., 80.};
//TString dir     = "2022-05-11_v2-true_no-corr_eta-1.2-7/res_qvecs-corr";
//TString dir     = "2022-05-10_v2-true-corrected/res_qvecs-no-corr";
TString dir     = "20220516_data/res_qvecs-corr";
//TString dir     = "20220517_data-mctracks/res_qvecs-mctrack";
//TString dir     = "data-v2_2022-02-09/res_qvecs-corr";
//TString dirtrue = "20220516_data/res_qvecs-corr";
//TString dirtrue = "data-true-res_2022-02-23/res_qvecs-corr";
TString files[ncent] = {"cent00-05.root", "cent05-10.root", "cent10-20.root",
                        "cent20-30.root", "cent30-40.root", "cent40-50.root",
                        "cent50-60.root", "cent60-80.root"};
TString detname[ndet] = {"Ideal", "FT0C", "FV0", "FT0A"};
TString legentry[ndet] = {"Ideal", "FT0-C", "FV0", "FT0-A"};
int mColor[ndet] = {kBlue, kRed, kBlue, kBlack};

float offset[ncent][ndet] = {
    {0.,-.6,0.,.6},
    {0.,-.6,0.,.6},
    {0.,-.5,0.,.5},
    {0.,-.5,0.,.5},
    {0.,-.5,0.,.5},
    {0.,-.5,0.,.5},
    {0.,-.5,0.,.5},
    {0.,-.5,0.,.5}
};

TFile *fin[ncent];
TFile *fintrue[ncent];
TH1D *hVnObs[ncent][ndet];
TH1D *hRAB[ncent][ndet];
TH1D *hRAC[ncent][ndet];
TH1D *hRABtrue[ncent];
TH1D *hRACtrue[ncent];
TH1D *hRBC[ncent];
TH1D *hRsub[ncent];
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
double R1(double khi);
double func(double *x, double *p);
double RIter(double x0, double R0, double err);
double CalculateRerror(double khi, double khiErr);
void ConfigPlots();
void PlotToCanvas();
void RedrawBorder();

void PlotV2()
{
    SetStyle(0);
    LoadData();

    for (int idet=0; idet<ndet; idet++) {
        cout << legentry[idet] << " : " << endl;
        gV2[idet] = new TGraphErrors();
        for (int icent=0; icent<ncent; icent++) {
            GetResAndErr(icent, idet);
            GetV2AndErr(icent, idet);
            if (idet==0) {
                gV2[idet]->SetPoint(icent, cent[icent] + (cent[icent+1]-cent[icent])/2., v2[idet][icent]);
                gV2[idet]->SetPointError(icent, (cent[icent+1]-cent[icent])/2., v2Err[idet][icent]);
            } else {
                gV2[idet]->SetPoint(icent, cent[icent] + (1.+offset[icent][idet])*(cent[icent+1]-cent[icent])/2., v2[idet][icent]);
                gV2[idet]->SetPointError(icent, 0., v2Err[idet][icent]);
            }
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
        for (int idet=1; idet<ndet; idet++) {
            hRAB[icent][idet]   = (TH1D*)fin[icent]->Get(Form("hRAB_%s", detname[idet].Data()));
            hRAC[icent][idet]   = (TH1D*)fin[icent]->Get(Form("hRAC_%s", detname[idet].Data()));
            hVnObs[icent][idet] = (TH1D*)fin[icent]->Get(Form("hVnObs_%s", detname[idet].Data()));
        }

        // Get "true" values
        //fintrue[icent]    = TFile::Open(Form("%s/%s", dirtrue.Data(), files[icent].Data()));
        //hRsub[icent]      = (TH1D*)fin[icent]->Get("hRsub");
        //hVnObs[icent][ndet-1] = (TH1D*)fin[icent]->Get("hVnObs");
        hRABtrue[icent]      = (TH1D*)fin[icent]->Get("hRABtrue");
        hRACtrue[icent]      = (TH1D*)fin[icent]->Get("hRACtrue");
        hVnObs[icent][0] = (TH1D*)fin[icent]->Get("hVnObsTrue");
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
}

void GetResAndErr(int icent, int idet)
{
    if (idet > 0) {
        float rab = hRAB[icent][idet]->GetMean();
        float rac = hRAC[icent][idet]->GetMean();
        float rbc = hRBC[icent]->GetMean();
        float rabErr = hRAB[icent][idet]->GetMeanError();
        float racErr = hRAC[icent][idet]->GetMeanError();
        float rbcErr = hRBC[icent]->GetMeanError();

        res[idet][icent] = TMath::Sqrt((rab * rac)/rbc);
        resErr[idet][icent] = 0.5*res[idet][icent]*TMath::Sqrt(TMath::Power(rabErr/rab, 2) + TMath::Power(racErr/rac, 2) + TMath::Power(rbcErr/rbc, 2));
    } else {
        float rab = hRABtrue[icent]->GetMean();
        float rac = hRACtrue[icent]->GetMean();
        float rbc = hRBC[icent]->GetMean();
        float rabErr = hRABtrue[icent]->GetMeanError();
        float racErr = hRACtrue[icent]->GetMeanError();
        float rbcErr = hRBC[icent]->GetMeanError();

        res[idet][icent] = TMath::Sqrt((rab * rac)/rbc);
        resErr[idet][icent] = 0.5*res[idet][icent]*TMath::Sqrt(TMath::Power(rabErr/rab, 2) + TMath::Power(racErr/rac, 2) + TMath::Power(rbcErr/rbc, 2));

        double khi0 = 0.5;
        double err = 0.0001;

        //res[idet][icent] = TMath::Sqrt(hRsub[icent]->GetMean());
        //double res0 = TMath::Sqrt(hRsub[icent]->GetMean());
        //double khi = RIter(khi0, res0, err);
        //cout << "\tKhi : " << khi << endl;
        //res[idet][icent] = R1(TMath::Sqrt(2)*khi);

        //resErr[idet][icent] = hRsub[icent]->GetMeanError()/(2.*TMath::Sqrt(hRsub[icent]->GetMean()));
        //resErr[idet][icent] = CalculateRerror(khi, err);
    }
    cout << "\tRes : " << res[idet][icent] << " +- " << resErr[idet][icent] << endl;
}

void GetV2AndErr(int icent, int idet)
{
    float vnObs = hVnObs[icent][idet]->GetMean();
    float vnObsErr = hVnObs[icent][idet]->GetMeanError();

    v2[idet][icent] = vnObs/res[idet][icent];
    //v2[idet][icent] = vnObs/(TMath::Sqrt(2.)*res[idet][icent]);
    v2Err[idet][icent] = v2[idet][icent]*TMath::Sqrt((resErr[idet][icent]/res[idet][icent])*(resErr[idet][icent]/res[idet][icent])
                        + (vnObsErr/vnObs)*(vnObsErr/vnObs));
}

double R1(double khi) {
    return (TMath::Sqrt(TMath::Pi())/2)*khi*TMath::Exp(-khi*khi/2)*(TMath::BesselI0(khi*khi/2) + TMath::BesselI1(khi*khi/2));
}

double func(double *x, double *p) {
    double khi = x[0];
    double Rk = p[0];
    return (TMath::Sqrt(TMath::Pi())/2)*khi*TMath::Exp(-khi*khi/2)*(TMath::BesselI0(khi*khi/2) + TMath::BesselI1(khi*khi/2)) - Rk;
}

double RIter(double x0, double R0, double err) {
    double x = 0;
    TF1 *fRes = new TF1("fRes", func, 0, 50.0, 1);
    fRes->SetParameter(0, R0);
    while (TMath::Abs(R1(x) - R0) > err) {
        x = x0 - fRes->Eval(x0)/fRes->Derivative(x0);
        x0 = x;
    }
    return x;
}

// Virheen yleisellä etenemisellä R(khi):n lausekkeesta
double CalculateRerror(double khi, double khiErr) {
    double bessel0 = -(khi*khi-2)*TMath::BesselI0(khi*khi/2);
    double bessel1 = 2*TMath::BesselI1(khi*khi/2);
    double bessel2 = khi*khi*TMath::BesselI(2,khi*khi/2);
    return TMath::Sqrt(TMath::Pi()/4)*TMath::Exp(-khi*khi)*(bessel0 + bessel1 + bessel2)*khiErr;
}

void ConfigPlots()
{
    leg = new TLegend(0.4, 0.2, 0.72, 0.3);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(gStyle->GetTextSize()*0.7);
    leg->SetNColumns(2);
    //leg->SetHeader(Form("%s #sqrt{#it{s}_{NN}} = %s simulation", colsys.Data(), energy.Data()));
    for (int idet=0; idet<ndet; idet++) {
        gV2[idet]->SetTitle("; centrality (%); v_{2}");
        gV2[idet]->GetYaxis()->SetRangeUser(0.032, 0.1299);
        gV2[idet]->GetXaxis()->SetLabelSize(0.035);
        gV2[idet]->GetXaxis()->SetTitleSize(0.04);
        gV2[idet]->GetXaxis()->SetLimits(0., 80.);
        gV2[idet]->GetYaxis()->SetMaxDigits(4);
        gV2[idet]->GetYaxis()->SetLabelSize(0.035);
        gV2[idet]->GetYaxis()->SetTitleSize(0.05);
        gV2[idet]->GetYaxis()->SetTitleOffset(0.8);
        //gV2[idet]->GetYaxis()->SetRangeUser(0., 1.);
        gV2[idet]->SetMarkerColor(mColor[idet]);
        gV2[idet]->SetLineColor(mColor[idet]);
        if (idet==0) {
            leg->AddEntry(gV2[idet], legentry[idet], "lef");
            gV2[idet]->SetFillColor(mColor[idet]-10);
            gV2[idet]->SetMarkerColor(mColor[idet]-8);
            gV2[idet]->SetMarkerStyle(kDot);
            gV2[idet]->SetLineColor(mColor[idet]-8);
            gV2[idet]->SetLineWidth(1);
        } else {
            leg->AddEntry(gV2[idet], legentry[idet], "pe");
            gV2[idet]->SetMarkerStyle(kOpenCircle);

        }
    }
}

void PlotToCanvas()
{
    c1 = new TCanvas("c1", "c1", 800, 800);
    for (int idet=0; idet<ndet; idet++) {
        if (idet==0) {
            gV2[idet]->Draw("A2");
            gV2[idet]->Draw("P");
        } else {
            gV2[idet]->Draw("P SAME");
        }
    }
    leg->Draw("SAME");

    //gV2V0A->Draw("P SAME");

    TLatex * tex = new TLatex(0.1,0.1, "ALICE Preliminary");
    tex->SetNDC();
    tex->SetTextFont(42);
    //tex->Draw();

    TLatex * texsim1 = new TLatex(0.402,0.42, Form("%s #sqrt{#it{s}_{NN}} = %s simulation", colsys.Data(), energy.Data()));
    texsim1->SetNDC();
    texsim1->SetTextSize(0.035);
    texsim1->Draw();

    TLatex * texsim2 = new TLatex(0.402,0.37, "-0.8 < #eta < 0.8");
    texsim2->SetNDC();
    texsim2->SetTextSize(0.035);
    texsim2->Draw();

    TLatex * texsim3 = new TLatex(0.402,0.32, "0.2 GeV/c < p_{T} < 5.0 GeV/c");
    texsim3->SetNDC();
    texsim3->SetTextSize(0.035);
    texsim3->Draw();

    RedrawBorder();
    c1->SaveAs("v2-true.eps");
}

void RedrawBorder()
{
   gPad->Update();
   gPad->RedrawAxis();
   TLine l;
   l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
   l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());
}
