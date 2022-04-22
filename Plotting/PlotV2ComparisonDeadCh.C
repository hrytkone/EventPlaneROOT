const int ncent = 8;
const int ndet = 2;
const int nset = 4;

TString colsys = "Pb#font[122]{-}Pb";
TString energy = "5.5 TeV";

float cent[ncent+1] = {0., 5., 10., 20., 30., 40., 50., 60., 80.};
//float cent[ncent+1] = {20., 30.};
//TString files[ncent] = {"../res_cent20-30_tpc-eff_corr.root"};
//TString dir[nset] = {"data-v2_2022-02-09/res_qvecs-corr", "2022-03-10_deadch-20/res_qvecs-corr", "2022-03-10_deadch-50/res_qvecs-corr", "2022-03-10_deadch-100/res_qvecs-corr"};
TString dir[nset] = {"data-v2_2022-02-09/res_qvecs-corr", "2022-04-19_deadch-10percent/res_qvecs-corr", "2022-04-19_deadch-25percent/res_qvecs-corr", "2022-04-19_deadch-50percent/res_qvecs-corr"};

TString files[ncent] = {"cent00-05.root", "cent05-10.root", "cent10-20.root",
                        "cent20-30.root", "cent30-40.root", "cent40-50.root",
                        "cent50-60.root", "cent60-80.root"};
//TString files[ncent] = {"cent20-30.root"};
TString detname[ndet] = {"FT0C", "FT0A"};
TString legentry[nset] = {"no dead ch", "10 dead ch", "25 dead ch", "50 dead ch"};
EColor mColor[nset] = {kBlack, kMagenta, kRed, kBlue};

float offset[nset] = {-0.1,-0.3,0.1,0.3};

TFile *fin[nset][ncent];
TH1D *hVnObs[nset][ncent][ndet];
TH1D *hRAB[nset][ncent][ndet];
TH1D *hRAC[nset][ncent][ndet];
TH1D *hRBC[nset][ncent];
TGraphErrors *gV2[nset][ndet];
TLegend *leg[ndet];
TCanvas *c1[ndet];

float res[nset][ndet][ncent];
float resErr[nset][ndet][ncent];
float v2[nset][ndet][ncent];
float v2Err[nset][ndet][ncent];

void LoadData();
void SetStyle(Bool_t graypalette);
void GetResAndErr(int iset, int icent, int idet);
void GetV2AndErr(int iset, int icent, int idet);
void ConfigPlots();
void PlotToCanvas();

void PlotV2ComparisonDeadCh()
{
    SetStyle(0);
    LoadData();

    for (int iset=0; iset<nset; iset++) {
        for (int idet=0; idet<ndet; idet++) {
            gV2[iset][idet] = new TGraphErrors();
            cout << detname[idet] << endl;
            for (int icent=0; icent<ncent; icent++) {
                GetResAndErr(iset, icent, idet);
                GetV2AndErr(iset, icent, idet);
                //gV2[iset][idet]->SetPoint(icent, cent[icent] + (cent[icent+1]-cent[icent])/2., v2[iset][idet][icent]);
                //gV2[iset][idet]->SetPointError(icent, (cent[icent+1]-cent[icent])/2., v2Err[iset][idet][icent]);
                gV2[iset][idet]->SetPoint(icent, cent[icent] + (1.+offset[iset])*(cent[icent+1]-cent[icent])/2., v2[iset][idet][icent]);
                gV2[iset][idet]->SetPointError(icent, 0., v2Err[iset][idet][icent]);
            }
        }
    }

    ConfigPlots();
    PlotToCanvas();
}

//________________________________
void LoadData()
{
    for (int iset=0; iset<nset; iset++) {
        for (int icent=0; icent<ncent; icent++) {
            fin[iset][icent] = TFile::Open(Form("%s/%s", dir[iset].Data(), files[icent].Data()));

            hRBC[iset][icent] = (TH1D*)fin[iset][icent]->Get("hRBC");
            for (int idet=0; idet<ndet; idet++) {
                hRAB[iset][icent][idet] = (TH1D*)fin[iset][icent]->Get(Form("hRAB_%s", detname[idet].Data()));
                hRAC[iset][icent][idet] = (TH1D*)fin[iset][icent]->Get(Form("hRAC_%s", detname[idet].Data()));
                hVnObs[iset][icent][idet] = (TH1D*)fin[iset][icent]->Get(Form("hVnObs_%s", detname[idet].Data()));
            }
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

void GetResAndErr(int iset, int icent, int idet)
{
    float rab = hRAB[iset][icent][idet]->GetMean();
    float rac = hRAC[iset][icent][idet]->GetMean();
    float rbc = hRBC[iset][icent]->GetMean();
    float rabErr = hRAB[iset][icent][idet]->GetMeanError();
    float racErr = hRAC[iset][icent][idet]->GetMeanError();
    float rbcErr = hRBC[iset][icent]->GetMeanError();

    res[iset][idet][icent] = TMath::Sqrt((rab * rac)/rbc);
    resErr[iset][idet][icent] = 0.5*res[iset][idet][icent]*TMath::Sqrt(TMath::Power(rabErr/rab, 2) + TMath::Power(racErr/rac, 2) + TMath::Power(rbcErr/rbc, 2));

    cout << "\tRes : " << res[iset][idet][icent] << " +- " << resErr[iset][idet][icent] << endl;
}

void GetV2AndErr(int iset, int icent, int idet)
{
    float vnObs = hVnObs[iset][icent][idet]->GetMean();
    float vnObsErr = hVnObs[iset][icent][idet]->GetMeanError();

    v2[iset][idet][icent] = vnObs/res[iset][idet][icent];
    //v2[idet][icent] = vnObs/(TMath::Sqrt(2.)*res[idet][icent]);
    v2Err[iset][idet][icent] = v2[iset][idet][icent]*TMath::Sqrt((resErr[iset][idet][icent]/res[iset][idet][icent])*(resErr[iset][idet][icent]/res[iset][idet][icent])
                        + (vnObsErr/vnObs)*(vnObsErr/vnObs));
}

void ConfigPlots()
{
    for (int idet=0; idet<ndet; idet++) {
        leg[idet] = new TLegend(0.42, 0.2, 0.65, 0.42);
        leg[idet]->SetFillStyle(0); leg[idet]->SetBorderSize(0); leg[idet]->SetTextSize(gStyle->GetTextSize()*0.65);
        leg[idet]->SetHeader(Form("%s #sqrt{#it{s}_{NN}} = %s", colsys.Data(), energy.Data()));
        for (int iset=0; iset<nset; iset++) {
            leg[idet]->AddEntry(gV2[iset][idet], legentry[iset], "p");
            gV2[iset][idet]->SetTitle("; centrality (%); v_{2}");
            gV2[iset][idet]->GetXaxis()->SetLabelSize(0.035);
            gV2[iset][idet]->GetXaxis()->SetTitleSize(0.04);
            gV2[iset][idet]->GetXaxis()->SetLimits(0., 80.);
            gV2[iset][idet]->GetYaxis()->SetMaxDigits(4);
            gV2[iset][idet]->GetYaxis()->SetLabelSize(0.035);
            gV2[iset][idet]->GetYaxis()->SetTitleSize(0.04);
            //gV2[iset][idet]->GetYaxis()->SetRangeUser(0., 1.);
            gV2[iset][idet]->SetMarkerColor(mColor[iset]);
            gV2[iset][idet]->SetMarkerStyle(kOpenCircle);
        }
    }
}

void PlotToCanvas()
{
    for (int idet=0; idet<ndet; idet++) {
        c1[idet] = new TCanvas(Form("c1_%d", idet), "c1", 800, 800);
        for (int iset=0; iset<nset; iset++) {
            if (iset==0)
                gV2[iset][idet]->Draw("AP");
            else
                gV2[iset][idet]->Draw("P SAME");
        }
        leg[idet]->Draw("SAME");

        TLatex * texsim = new TLatex(0.42,0.43, Form("Simulation for %s", detname[idet].Data()));
        texsim->SetNDC();
        texsim->SetTextSize(gStyle->GetTextSize()*0.65);
        texsim->Draw();
    }
}
