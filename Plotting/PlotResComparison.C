const int nset = 3;
const int ndet = 2;

TString colsys = "Pb#font[122]{-}Pb";
TString energy = "5.5 TeV";

TString legentry[ndet] = {"FV0", "FT0-A"};
EColor mColor[ndet] = {kRed, kBlue};
TString leglabel[nset] = {"default geometry", "new location", "FV0 fibre + new location"};

TGraphErrors *gRes[ndet];
TLegend *leg;
TCanvas *c1;

float res[ndet][nset] = {{0.770082,0.760733,0.760607},
                        {0.756099,0.737033,0.738304}};
float resErr[ndet][nset] = {{0.00419979,0.00625369,0.00596725},
                            {0.00422471,0.00632708,0.00608569}};

void SetStyle(Bool_t graypalette);
void ConfigPlots();
void PlotToCanvas();

void PlotResComparison()
{
    SetStyle(0);

    for (int idet=0; idet<ndet; idet++) {
        gRes[idet] = new TGraphErrors();
        for (int iset=0; iset<nset; iset++) {
            gRes[idet]->SetPoint(iset, (double)iset+.5, res[idet][iset]);
            gRes[idet]->SetPointError(iset, 0., resErr[idet][iset]);
        }
        TAxis *ax = gRes[idet]->GetHistogram()->GetXaxis();
        Double_t x1 = ax->GetBinLowEdge(1);
        Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
        gRes[idet]->GetHistogram()->GetXaxis()->Set(3,x1,x2);

        for (Int_t k=0;k<nset;k++) {
            gRes[idet]->GetHistogram()->GetXaxis()->SetBinLabel(k+1, leglabel[k].Data());
        }
    }

    ConfigPlots();
    PlotToCanvas();
}

//________________________________
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

void ConfigPlots()
{
    leg = new TLegend(0.6, 0.63, 0.72, 0.79);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(gStyle->GetTextSize()*0.5);
    leg->SetHeader(Form("%s #sqrt{#it{s}_{NN}} = %s", colsys.Data(), energy.Data()));
    for (int idet=0; idet<ndet; idet++) {
        leg->AddEntry(gRes[idet], legentry[idet], "p");
        gRes[idet]->SetTitle("; ; R_{2}");
        //gRes[idet]->GetXaxis()->SetLabelSize(0.025);
        //gRes[idet]->GetXaxis()->SetTitleSize(0.025);
        //gRes[idet]->GetXaxis()->CenterTitle();
        //gRes[idet]->GetXaxis()->SetLimits(0., 80.);
        //gRes[idet]->GetYaxis()->SetLabelSize(0.025);
        //gRes[idet]->GetYaxis()->SetTitleSize(0.025);
        //gRes[idet]->GetYaxis()->CenterTitle();
        gRes[idet]->GetYaxis()->SetRangeUser(0.7, 0.8);
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

    TLatex * texsim = new TLatex(0.602,0.8, "Simulation");
    texsim->SetNDC();
    texsim->SetTextSize(gStyle->GetTextSize()*0.5);
    texsim->Draw();
}
