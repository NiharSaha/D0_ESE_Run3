// plot_SPdist_qbins.C
//
// Standalone macro: draws every SP-method v2/v3 distribution histogram
// (histograms_SP/h_v2_hist_<cent>_q2bin<Q>_<pT>  and
//  histograms_SP/h_v3_hist_<cent>_q3bin<Q>_<pT>)
// from D0_Flow_12NUQbin_diffq2q3_out_combined_Aug24.root as an individual
// single-pad PDF, then tiles all 12 q bins for a given (cent, pT) onto one
// page (4 cols x 3 rows) via tile_SPdist_pdfs.py.
//
//   v2: 5 cent x 9 pT x 12 q = 540 plots -> histogram_SP/SPdist_v2_AllQ2Bins.pdf (45 pages)
//   v3: 5 cent x 6 pT x 12 q = 360 plots -> histogram_SP/SPdist_v3_AllQ3Bins.pdf (30 pages)
//
// Requirements: python3 with PyMuPDF (pip install pymupdf --break-system-packages)
//
// Usage:  root -l -q plot_SPdist_qbins.C

#include "TFile.h"
#include "TDirectoryFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

static void setSPStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    TGaxis::SetMaxDigits(2);
}

static void DrawCMSHeader(double leftMargin, double rightMargin, double y)
{
    TLatex cmsLat;
    cmsLat.SetNDC();
    cmsLat.SetTextFont(42);
    cmsLat.SetTextSize(0.048);
    cmsLat.SetTextAlign(11);
    cmsLat.DrawLatex(leftMargin + 0.08, y, "#font[61]{CMS} #font[52]{Preliminary}");

    TLatex enLat;
    enLat.SetNDC();
    enLat.SetTextFont(42);
    enLat.SetTextSize(0.048);
    enLat.SetTextAlign(31);
    enLat.DrawLatex(1.0 - rightMargin, y, "PbPb #sqrt{s_{NN}} = 5.36 TeV");
}

void plot_SPdist_qbins()
{
    setSPStyle();
    gROOT->ForceStyle();
    gROOT->SetBatch(kTRUE);

    // ── open file ────────────────────────────────────────────────────────
    TFile* f = TFile::Open("D0_Flow_12NUQbin_diffq2q3_out_combined_Aug24.root");
    if (!f || f->IsZombie()) {
        Printf("ERROR: cannot open input ROOT file");
        return;
    }
    TDirectoryFile* d = (TDirectoryFile*)f->Get("histograms_SP");
    if (!d) {
        Printf("ERROR: histograms_SP directory not found");
        f->Close(); return;
    }

    // ── bin definitions (must match the ROOT file's naming convention) ────
    const int   NCENT            = 5;
    const char* centTag  [NCENT] = { "cent0to10", "cent10to20", "cent20to30",
                                      "cent30to40", "cent40to50" };
    const char* centLabel[NCENT] = { "0-10%", "10-20%", "20-30%",
                                      "30-40%", "40-50%" };

    const int   NPT2            = 9;
    const char* ptTag2   [NPT2] = { "pT2to3",  "pT3to4",  "pT4to5",
                                     "pT5to6",  "pT6to8",  "pT8to10",
                                     "pT10to15","pT15to30","pT30to100" };
    const char* ptLo2    [NPT2] = { "2", "3", "4", "5", "6", "8", "10", "15", "30" };
    const char* ptHi2    [NPT2] = { "3", "4", "5", "6", "8", "10", "15", "30", "100" };

    const int   NPT3            = 6;
    const char* ptTag3   [NPT3] = { "pT2to4", "pT4to6", "pT6to8",
                                     "pT8to10", "pT10to20", "pT20to50" };
    const char* ptLo3    [NPT3] = { "2", "4", "6", "8", "10", "20" };
    const char* ptHi3    [NPT3] = { "4", "6", "8", "10", "20", "50" };

    const int NQ = 12;

    TString outBase = "histogram_SP";
    gSystem->mkdir(outBase, kTRUE);
    TString v2TmpDir = outBase + "/v2_pdfs_tmp";
    TString v3TmpDir = outBase + "/v3_pdfs_tmp";
    gSystem->mkdir(v2TmpDir, kTRUE);
    gSystem->mkdir(v3TmpDir, kTRUE);

    TCanvas* cv = new TCanvas("cv_single", "", 600, 500);

    // ========================================================================
    // Generic per-harmonic loop
    // ========================================================================
    auto drawHarmonic = [&](int harmonic, int NPT, const char* const* ptTag,
                             const char* const* ptLo, const char* const* ptHi,
                             const TString& tmpDir)
    {
        for (int ic = 0; ic < NCENT; ic++) {
            for (int ip = 0; ip < NPT; ip++) {
                for (int iq = 0; iq < NQ; iq++) {

                    cv->Clear();
                    cv->SetLeftMargin  (0.16);
                    cv->SetRightMargin (0.05);
                    cv->SetBottomMargin(0.14);
                    cv->SetTopMargin   (0.10);

                    TString hname = Form("h_v%d_hist_%s_q%dbin%d_%s",
                                          harmonic, centTag[ic], harmonic, iq, ptTag[ip]);
                    TH1D* h = (TH1D*)d->Get(hname);

                    if (!h) {
                        Printf("WARNING: %s not found", hname.Data());
                        TLatex miss;
                        miss.SetNDC();
                        miss.SetTextFont(42);
                        miss.SetTextSize(0.07);
                        miss.SetTextAlign(22);
                        miss.DrawLatex(0.5, 0.5, "no data");
                        cv->SaveAs(Form("%s/sp_c%d_p%d_q%d.pdf", tmpDir.Data(), ic, ip, iq));
                        continue;
                    }

                    TH1D* hc = (TH1D*)h->Clone(Form("hc_sp_%d_%d_%d_%d", harmonic, ic, ip, iq));
                    hc->SetDirectory(nullptr);

                    hc->GetXaxis()->SetTitle(Form("SP v_{%d}", harmonic));
                    hc->GetYaxis()->SetTitle("Entries");
                    hc->GetXaxis()->SetTitleSize (0.055);
                    hc->GetYaxis()->SetTitleSize (0.055);
                    hc->GetXaxis()->SetLabelSize (0.050);
                    hc->GetYaxis()->SetLabelSize (0.050);
                    hc->GetXaxis()->SetTitleOffset(0.95);
                    hc->GetYaxis()->SetTitleOffset(1.15);
                    hc->GetXaxis()->SetNdivisions(505);

                    hc->SetLineColor(kBlack);
                    hc->SetLineWidth(2);
                    hc->SetFillColorAlpha(kAzure - 4, 0.35);

                    hc->Draw("HIST");

                    // cent / q-bin / pT label, top-left INSIDE the frame
                    TLatex lat1, lat2;
                    lat1.SetNDC(); lat2.SetNDC();
                    lat1.SetTextFont(42); lat2.SetTextFont(42);
                    lat1.SetTextSize(0.050); lat2.SetTextSize(0.050);
                    lat1.SetTextAlign(12); lat2.SetTextAlign(12);
                    lat1.DrawLatex(0.20, 0.87,
                        Form("Cent : %s (q_{%d} bin =%d)", centLabel[ic], harmonic, iq));
                    lat2.DrawLatex(0.20, 0.80,
                        Form("%s <p_{T}< %s GeV", ptLo[ip], ptHi[ip]));

                    DrawCMSHeader(0.16, 0.05, 0.925);

                    cv->Modified();
                    cv->Update();
                    cv->SaveAs(Form("%s/sp_c%d_p%d_q%d.pdf", tmpDir.Data(), ic, ip, iq));
                    delete hc;
                }
            }
            Printf("Saved v%d SP PDFs: cent %s", harmonic, centLabel[ic]);
        }
    };

    drawHarmonic(2, NPT2, ptTag2, ptLo2, ptHi2, v2TmpDir);
    drawHarmonic(3, NPT3, ptTag3, ptLo3, ptHi3, v3TmpDir);

    delete cv;
    f->Close();

    // ========================================================================
    // Tile both sets into their final multi-page PDFs
    // ========================================================================
    Printf("Tiling v2 SP PDFs (45 pages, 12/page) ...");
    gSystem->Exec("/usr/local/bin/python3.11 tile_SPdist_pdfs.py v2");

    Printf("Tiling v3 SP PDFs (30 pages, 12/page) ...");
    gSystem->Exec("/usr/local/bin/python3.11 tile_SPdist_pdfs.py v3");

    Printf("Done.");
}
