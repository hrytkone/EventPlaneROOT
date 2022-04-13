const int ncent = 8;
const int ndet = 2;
const int nset = 4;

TString colsys = "Pb#font[122]{-}Pb";
TString energy = "5.5 TeV";

float cent[ncent+1] = {0., 5., 10., 20., 30., 40., 50., 60., 80.};
//float cent[ncent+1] = {20., 30.};
//TString files[ncent] = {"../res_cent20-30_tpc-eff_corr.root"};
TString dir[nset] = {"data-v2_2022-02-09/res_qvecs-corr", "2022-03-10_deadch-20/res_qvecs-corr", "2022-03-10_deadch-50/res_qvecs-corr", "2022-03-10_deadch-100/res_qvecs-corr"};
TString files[ncent] = {"cent00-05.root", "cent05-10.root", "cent10-20.root",
                        "cent20-30.root", "cent30-40.root", "cent40-50.root",
                        "cent50-60.root", "cent60-80.root"};
//TString files[ncent] = {"cent20-30.root"};
TString detname[ndet] = {"FT0C", "FT0A"};
TString legentry[nset] = {"no dead ch", "10 dead ch", "25 dead ch", "50 dead ch"};
EColor mColor[nset] = {kBlack, kMagenta, kRed, kBlue};

// Leg coordinates
double xi[ndet] = {0.20, 0.58};
double yi[ndet] = {0.20, 0.61};
double xf[ndet] = {0.45, 0.72};
double yf[ndet] = {0.42, 0.82};

TFile *fin[nset][ncent];
TH1D *hRAB[nset][ncent][ndet];
TH1D *hRAC[nset][ncent][ndet];
TH1D *hRBC[nset][ncent];
TGraphErrors *gRes[nset][ndet];
TLegend *leg[ndet];
TCanvas *c1[ndet];

float res[nset][ndet][ncent];
float resErr[nset][ndet][ncent];

void LoadData();
void SetStyle(Bool_t graypalette);
void GetResAndErr(int iset, int icent, int idet);
void ConfigPlots();
void PlotToCanvas();

void PlotResComparisonDeadCh()
{
    SetStyle(0);
    LoadData();

    for (int iset=0; iset<nset; iset++) {
        for (int idet=0; idet<ndet; idet++) {
            gRes[iset][idet] = new TGraphErrors();
            cout << detname[idet] << endl;
            for (int icent=0; icent<ncent; icent++) {
                GetResAndErr(iset, icent, idet);
                gRes[iset][idet]->SetPoint(icent, cent[icent] + (cent[icent+1]-cent[icent])/2., res[iset][idet][icent]);
                gRes[iset][idet]->SetPointError(icent, (cent[icent+1]-cent[icent])/2., resErr[iset][idet][icent]);
            }
        }
    }

    ConfigPlots();
    PlotToCanvas();
}

//______________________________________________________________________________
void LoadData()
{
    for (int iset=0; iset<nset; iset++) {
        for (int icent=0; icent<ncent; icent++) {
            fin[iset][icent] = TFile::Open(Form("%s/%s", dir[iset].Data(), files[icent].Data()));

            hRBC[iset][icent] = (TH1D*)fin[iset][icent]->Get("hRBC");
            for (int idet=0; idet<ndet; idet++) {
                hRAB[iset][icent][idet] = (TH1D*)fin[iset][icent]->Get(Form("hRAB_%s", detname[idet].Data()));
                hRAC[iset][icent][idet] = (TH1D*)fin[iset][icent]->Get(Form("hRAC_%s", detname[idet].Data()));
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

void ConfigPlots()
{
    for (int idet=0; idet<ndet; idet++) {
        leg[idet] = new TLegend(xi[idet], yi[idet], xf[idet], yf[idet]);
        leg[idet]->SetFillStyle(0); leg[idet]->SetBorderSize(0); leg[idet]->SetTextSize(gStyle->GetTextSize()*0.65);
        leg[idet]->SetHeader(Form("%s #sqrt{#it{s}_{NN}} = %s", colsys.Data(), energy.Data()));
        for (int iset=0; iset<nset; iset++) {
            leg[idet]->AddEntry(gRes[iset][idet], legentry[iset], "p");
            gRes[iset][idet]->SetTitle("; centrality (%); R_{2}");
            //gRes[idet]->GetXaxis()->SetLabelSize(0.025);
            //gRes[idet]->GetXaxis()->SetTitleSize(0.025);
            //gRes[idet]->GetXaxis()->CenterTitle();
            gRes[iset][idet]->GetXaxis()->SetLimits(0., 80.);
            //gRes[idet]->GetYaxis()->SetLabelSize(0.025);
            //gRes[idet]->GetYaxis()->SetTitleSize(0.025);
            //gRes[idet]->GetYaxis()->CenterTitle();
            gRes[iset][idet]->GetYaxis()->SetRangeUser(0., 1.);
            gRes[iset][idet]->SetMarkerColor(mColor[iset]);
            gRes[iset][idet]->SetMarkerStyle(kCircle);
        }
    }
}

void PlotToCanvas()
{
    for (int idet=0; idet<ndet; idet++) {
        c1[idet] = new TCanvas(Form("c1_%d", idet), "c1", 800, 800);
        for (int iset=0; iset<nset; iset++) {
            if (iset==0)
                gRes[iset][idet]->Draw("AP");
            else
                gRes[iset][idet]->Draw("P SAME");
        }
        leg[idet]->Draw("SAME");

        TLatex * texsim = new TLatex(xi[idet],yf[idet]+0.01, Form("Simulation for %s", detname[idet].Data()));
        texsim->SetNDC();
        texsim->SetTextSize(gStyle->GetTextSize()*0.65);
        texsim->Draw();
    }
}
