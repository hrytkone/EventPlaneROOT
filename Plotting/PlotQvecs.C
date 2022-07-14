const int nset = 1;
const int ncorr = 2;
const int nq = 2;

const double xi = 0.60;
const double yi = 0.75;
const double xf = 0.90;
const double yf = 0.95;

TString colsys = "Pb#font[122]{-}Pb";
TString energy = "5.5 TeV";

TString dir[nset]     = {"2022-05-10_v2-true-corrected"};
//TString dir[nset] = {"data-v2_2022-02-09", "2022-03-10_deadch-20", "2022-03-10_deadch-50", "2022-03-10_deadch-100"};
TString files[ncorr] = {"qvecs-no-corr/cent20-30.root", "qvecs-corr/cent20-30.root"};
TString legentry[nset] = {"no dead ch"};//, "10 dead ch", "25 dead ch", "50 dead ch"};

double binMin[ncorr] = {-0.49, -9.99};
double binMax[ncorr] = {0.49, 9.99};

TFile *fin[nset][ncorr];
TTree *fTree[nset][ncorr];
float qvecFT0A[nset][ncorr][nq];
float qvecFT0C[nset][ncorr][nq];
float qvecA[nset][ncorr][nq];
TH2D *hQvecFT0A[nset][ncorr];
TH2D *hQvecFT0C[nset][ncorr];
TH2D *hQvecA[nset][ncorr];
TH1D *hEPFT0A[nset][ncorr];
TH1D *hEPFT0C[nset][ncorr];
TCanvas *cFT0A;
TCanvas *cFT0C;
TCanvas *cA;
TCanvas *cFT0Aep;
TCanvas *cFT0Cep;
TLegend *leg[nset][ncorr];

void LoadData();
void CreateHistos();
void SetStyle(Bool_t graypalette);
void GetResAndErr(int iset, int icent, int idet);
void ConfigPlots();
void PlotToCanvas();
float GetEventPlane(float qx, float qy);

void PlotQvecs()
{
    SetStyle(0);
    LoadData();
    CreateHistos();
    ConfigPlots();
    PlotToCanvas();
}

//________________________________
void LoadData()
{
    for (int iset=0; iset<nset; iset++) {
        for (int icorr=0; icorr<ncorr; icorr++) {
            fin[iset][icorr] = TFile::Open(Form("%s/%s", dir[iset].Data(), files[icorr].Data()));
            fTree[iset][icorr] = (TTree*)fin[iset][icorr]->Get("Qvecs");
            //fTree[iset][icorr]->Print();
            fTree[iset][icorr]->SetBranchAddress("qvecFT0A", &qvecFT0A[iset][icorr]);
            fTree[iset][icorr]->SetBranchAddress("qvecFT0C", &qvecFT0C[iset][icorr]);
            fTree[iset][icorr]->SetBranchAddress("qvecA", &qvecA[iset][icorr]);
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

void CreateHistos()
{
    for (int iset=0; iset<nset; iset++) {
        for (int icorr=0; icorr<ncorr; icorr++) {
            hQvecFT0A[iset][icorr] = new TH2D(Form("hQvecFT0A_%d_%d", iset, icorr), "", 200, binMin[icorr], binMax[icorr], 200, binMin[icorr], binMax[icorr]);
            hQvecFT0C[iset][icorr] = new TH2D(Form("hQvecFT0C_%d_%d", iset, icorr), "", 200, binMin[icorr], binMax[icorr], 200, binMin[icorr], binMax[icorr]);
            hQvecA[iset][icorr] = new TH2D(Form("hQvecA_%d_%d", iset, icorr), "", 200, binMin[icorr], binMax[icorr], 200, binMin[icorr], binMax[icorr]);
            hEPFT0A[iset][icorr] = new TH1D(Form("hEPFT0A_%d_%d", iset, icorr), "", (int)TMath::Pi()/0.1, -TMath::Pi()/2., TMath::Pi()/2.);
            hEPFT0A[iset][icorr]->Sumw2();
            hEPFT0C[iset][icorr] = new TH1D(Form("hEPFT0C_%d_%d", iset, icorr), "", (int)TMath::Pi()/0.1, -TMath::Pi()/2., TMath::Pi()/2.);
            hEPFT0C[iset][icorr]->Sumw2();
            int nent = fTree[iset][icorr]->GetEntries();
            for (int ient=0; ient<nent; ient++) {
                fTree[iset][icorr]->GetEntry(ient);
                hQvecFT0A[iset][icorr]->Fill(qvecFT0A[iset][icorr][0], qvecFT0A[iset][icorr][1]);
                hQvecFT0C[iset][icorr]->Fill(qvecFT0C[iset][icorr][0], qvecFT0C[iset][icorr][1]);
                hQvecA[iset][icorr]->Fill(qvecA[iset][icorr][0], qvecA[iset][icorr][1]);
                hEPFT0A[iset][icorr]->Fill(GetEventPlane(qvecFT0A[iset][icorr][0], qvecFT0A[iset][icorr][1]));
                hEPFT0C[iset][icorr]->Fill(GetEventPlane(qvecFT0C[iset][icorr][0], qvecFT0C[iset][icorr][1]));
            }
        }
    }
}

void ConfigPlots()
{
    for (int iset=0; iset<nset; iset++) {
        for (int icorr=0; icorr<ncorr; icorr++) {
            hEPFT0A[iset][icorr]->Scale(1./hEPFT0A[iset][icorr]->GetEntries(), "width");
            int imaxFT0A = hEPFT0A[iset][icorr]->GetMaximumBin();
            hEPFT0A[iset][icorr]->GetYaxis()->SetRangeUser(0, hEPFT0A[iset][icorr]->GetBinContent(imaxFT0A)+0.2);
            hEPFT0A[iset][icorr]->SetTitle(";#Psi_{2};1/N dN/d#Psi_{2}");

            hEPFT0C[iset][icorr]->Scale(1./hEPFT0C[iset][icorr]->GetEntries(), "width");
            int imaxFT0C = hEPFT0C[iset][icorr]->GetMaximumBin();
            hEPFT0C[iset][icorr]->GetYaxis()->SetRangeUser(0, hEPFT0C[iset][icorr]->GetBinContent(imaxFT0C)+0.2);
            hEPFT0C[iset][icorr]->SetTitle(";#Psi_{2};1/N dN/d#Psi_{2}");

            hQvecFT0A[iset][icorr]->SetTitle(";Q_{x};Q_{y}");
            hQvecFT0C[iset][icorr]->SetTitle(";Q_{x};Q_{y}");
            hQvecA[iset][icorr]->SetTitle(";Q_{x};Q_{y}");

            leg[iset][icorr] = new TLegend(xi, yi, xf, yf);
            leg[iset][icorr]->SetFillStyle(0); leg[iset][icorr]->SetBorderSize(0); leg[iset][icorr]->SetTextSize(gStyle->GetTextSize()*1.2);
            if (icorr==0)
                leg[iset][icorr]->SetHeader(Form("#splitline{%s}{before corr.}", legentry[iset].Data()));
            else
                leg[iset][icorr]->SetHeader(Form("#splitline{%s}{after corr.}", legentry[iset].Data()));
        }
    }
}

void PlotToCanvas()
{
    cFT0A = new TCanvas("cFT0A", "cFT0A", 1600, 800);
    cFT0A->Divide(4,2,0.001, 0.0001);
    for (int icorr=0; icorr<ncorr; icorr++) {
        for (int iset=0; iset<nset; iset++) {
            cFT0A->cd(icorr*nset + iset + 1);
            gPad->SetLeftMargin(0.1);
            gPad->SetBottomMargin(0.1);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.05);
            hQvecFT0A[iset][icorr]->Draw("COL");
            leg[iset][icorr]->Draw("SAME");
        }
    }

    cFT0C = new TCanvas("cFT0C", "cFT0C", 1600, 800);
    cFT0C->Divide(4,2,0.001,0.001);
    for (int icorr=0; icorr<ncorr; icorr++) {
        for (int iset=0; iset<nset; iset++) {
            cFT0C->cd(icorr*nset + iset + 1);
            gPad->SetLeftMargin(0.1);
            gPad->SetBottomMargin(0.1);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.05);
            hQvecFT0C[iset][icorr]->Draw("COL");
            leg[iset][icorr]->Draw("SAME");
        }
    }

    cA = new TCanvas("cA", "cA", 1600, 800);
    cA->Divide(4,2,0.001,0.001);
    for (int icorr=0; icorr<ncorr; icorr++) {
        for (int iset=0; iset<nset; iset++) {
            cA->cd(icorr*nset + iset + 1);
            gPad->SetLeftMargin(0.1);
            gPad->SetBottomMargin(0.1);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.05);
            hQvecA[iset][icorr]->Draw("COL");
            leg[iset][icorr]->Draw("SAME");
        }
    }

    cFT0Aep = new TCanvas("cFT0Aep", "cFT0Aep", 1600, 800);
    cFT0Aep->Divide(4,2,0.001,0.001);
    for (int icorr=0; icorr<ncorr; icorr++) {
        for (int iset=0; iset<nset; iset++) {
            cFT0Aep->cd(icorr*nset + iset + 1);
            gPad->SetLeftMargin(0.1);
            gPad->SetBottomMargin(0.1);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.05);
            hEPFT0A[iset][icorr]->Draw("HIST E");
            leg[iset][icorr]->Draw("SAME");
        }
    }

    cFT0Cep = new TCanvas("cFT0Cep", "cFT0Cep", 1600, 800);
    cFT0Cep->Divide(4,2,0.001,0.001);
    for (int icorr=0; icorr<ncorr; icorr++) {
        for (int iset=0; iset<nset; iset++) {
            cFT0Cep->cd(icorr*nset + iset + 1);
            gPad->SetLeftMargin(0.1);
            gPad->SetBottomMargin(0.1);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.05);
            hEPFT0C[iset][icorr]->Draw("HIST E");
            leg[iset][icorr]->Draw("SAME");
        }
    }

    cFT0A->SaveAs("Qvecs-ft0a.pdf");
    cFT0C->SaveAs("Qvecs-ft0c.pdf");
    cFT0Aep->SaveAs("eps-ft0a.pdf");
    cFT0Cep->SaveAs("eps-ft0c.pdf");
}

float GetEventPlane(float qx, float qy)
{
    return TMath::ATan2(qy, qx)/2.0;
}
