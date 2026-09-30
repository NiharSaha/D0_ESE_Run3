// plot_mass_qbins.C
//
// Step 1: save every DATA (h_mass_cent<C>_q2bin<Q>_pT<PT>, 540 histos) and
//         MC (h_mc_sig_v2_<ic>_<ip>_<ip>_cent<C>, 45 histos) mass histogram,
//         together with its embedded fit, as an individual single-pad PDF.
//         Single-pad canvases avoid the TPaveText/TLegend/TLatex
//         repositioning issues that come with TCanvas::Divide() multi-pad
//         layouts -- each plot is saved exactly as ROOT renders it.
//
// Step 2: tile_mass_pdfs.py (PyMuPDF) merges:
//           data_pdfs_tmp/  (540 files) -> DataMass_AllQ2Bins.pdf
//                                          (45 pages, 4 cols x 3 rows = 12/page,
//                                           one page per cent/pT combination)
//           mc_pdfs_tmp/    (45  files) -> MCMass_Sig_AllPtBins.pdf
//                                          (5 pages, 3 cols x 3 rows = 9/page,
//                                           one page per centrality)
//
// Requirements: python3 with PyMuPDF  (pip install pymupdf --break-system-packages)
//
// Usage:  root -l -q plot_mass_qbins.C

#include "TFile.h"
#include "TDirectoryFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

void setMassStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);

    // Force axis labels to a maximum of 2 significant digits; beyond that
    // ROOT switches to the common "xN x10^k" exponent notation instead of
    // printing e.g. "100,200,...,800".
    TGaxis::SetMaxDigits(2);
}

// Standard CMS-style header: "CMS Preliminary" (or "CMS Simulation" for MC)
// left-aligned at top-left, and the collision-system/energy label
// right-aligned at top-right, both sitting just above the frame.
//
// The CMS text starts noticeably right of the pad's left margin (rather
// than exactly at it) because ROOT draws the axis's "x10^k" exponent
// multiplier in that same top-left corner, just above the y-axis -- flush
// alignment causes the two to overlap.
static void DrawCMSHeader(double leftMargin, double rightMargin, double y, bool isMC)
{
    TLatex cmsLat;
    cmsLat.SetNDC();
    cmsLat.SetTextFont(42);
    cmsLat.SetTextSize(0.048);
    cmsLat.SetTextAlign(11); // left, bottom
    TString cmsText = isMC ? "#font[61]{CMS} #font[52]{Simulation}"
                            : "#font[61]{CMS} #font[52]{Preliminary}";
    cmsLat.DrawLatex(leftMargin + 0.08, y, cmsText);

    TLatex enLat;
    enLat.SetNDC();
    enLat.SetTextFont(42);
    enLat.SetTextSize(0.048);
    enLat.SetTextAlign(31); // right, bottom
    enLat.DrawLatex(1.0 - rightMargin, y, "PbPb #sqrt{s_{NN}} = 5.36 TeV");
}

void plot_mass_qbins()
{
    setMassStyle();
    gROOT->ForceStyle();
    gROOT->SetBatch(kTRUE);

    // ── open file ────────────────────────────────────────────────────────
    TFile* f = TFile::Open("D0_Flow_12NUQbin_diffq2q3_out_combined_Aug24.root");
    if (!f || f->IsZombie()) {
        Printf("ERROR: cannot open input ROOT file");
        return;
    }
    TDirectoryFile* d = (TDirectoryFile*)f->Get("mass_fitted");
    if (!d) {
        Printf("ERROR: mass_fitted directory not found");
        f->Close(); return;
    }

    // ── bin definitions (must match the ROOT file's naming convention) ────
    const int   NCENT            = 5;
    const char* centTag  [NCENT] = { "cent0to10", "cent10to20", "cent20to30",
                                      "cent30to40", "cent40to50" };
    const char* centLabel[NCENT] = { "0-10%", "10-20%", "20-30%",
                                      "30-40%", "40-50%" };

    const int   NPT             = 9;
    const char* ptTag    [NPT]  = { "pT2to3",  "pT3to4",  "pT4to5",
                                     "pT5to6",  "pT6to8",  "pT8to10",
                                     "pT10to15","pT15to30","pT30to100" };
    const char* ptLabel  [NPT]  = { "2-3",  "3-4",  "4-5",
                                     "5-6",  "6-8",  "8-10",
                                     "10-15","15-30","30-100" };
    const char* ptLo     [NPT]  = { "2", "3", "4", "5", "6", "8", "10", "15", "30" };
    const char* ptHi     [NPT]  = { "3", "4", "5", "6", "8", "10", "15", "30", "100" };

    const int NQ = 12;

    TString dataTmpDir = "data_pdfs_tmp";
    TString mcTmpDir   = "mc_pdfs_tmp";
    gSystem->mkdir(dataTmpDir, kTRUE);
    gSystem->mkdir(mcTmpDir,   kTRUE);

    TCanvas* cv = new TCanvas("cv_single", "", 600, 500);

    // ========================================================================
    // DATA: h_mass_cent<CENT>_q2bin<Q>_pT<PT>   (540 plots)
    // Style: black points+error bars, fit curves as saved in the file
    //        (total/signal/swap/bkg/KK-reflection each already have their own
    //        colors). "Cent / q bin / pT" label sits top-left INSIDE the
    //        frame. The Yield/significance/chi2 annotations already embedded
    //        in the histogram are repositioned to the bottom-left, clear of
    //        both the new label and the y-axis.
    // ========================================================================
    for (int ic = 0; ic < NCENT; ic++) {
        for (int ip = 0; ip < NPT; ip++) {
            for (int iq = 0; iq < NQ; iq++) {

                cv->Clear();
                cv->SetLeftMargin  (0.16);
                cv->SetRightMargin (0.05);
                cv->SetBottomMargin(0.14);
                cv->SetTopMargin   (0.10);

                TString hname = Form("h_mass_%s_q2bin%d_%s",
                                      centTag[ic], iq, ptTag[ip]);
                TH1D* h = (TH1D*)d->Get(hname);

                if (!h) {
                    Printf("WARNING: %s not found", hname.Data());
                    TLatex miss;
                    miss.SetNDC();
                    miss.SetTextFont(42);
                    miss.SetTextSize(0.07);
                    miss.SetTextAlign(22);
                    miss.DrawLatex(0.5, 0.5, "no data");
                    cv->SaveAs(Form("%s/d_c%d_p%d_q%d.pdf",
                                    dataTmpDir.Data(), ic, ip, iq));
                    continue;
                }

                TH1D* hc = (TH1D*)h->Clone(Form("hc_d_%d_%d_%d", ic, ip, iq));
                hc->SetDirectory(nullptr);

                hc->GetXaxis()->SetTitle("Mass (GeV)");
                hc->GetYaxis()->SetTitle("Entries/5 MeV");
                hc->GetXaxis()->SetTitleSize (0.055);
                hc->GetYaxis()->SetTitleSize (0.055);
                hc->GetXaxis()->SetLabelSize (0.050);
                hc->GetYaxis()->SetLabelSize (0.050);
                hc->GetXaxis()->SetTitleOffset(0.95);
                hc->GetYaxis()->SetTitleOffset(1.15);
                hc->GetXaxis()->SetNdivisions(505);

                hc->SetMarkerStyle(20);
                hc->SetMarkerSize(0.7);
                hc->SetMarkerColor(kBlack);
                hc->SetLineColor(kBlack);

                // Reposition the Yield / significance / chi2 TLatex lines
                // (embedded in the histogram's fFunctions, originally saved
                // at x=0.14, which sits partly behind the y-axis once the
                // frame's left margin is set to 0.16). Move them to the
                // bottom-left of the frame -- the flattest, emptiest part
                // of the mass spectrum -- clear of the axis.
                {
                    const double bx = 0.20;
                    const double by[3] = { 0.30, 0.24, 0.18 };
                    int idx = 0;
                    TIter nextFunc(hc->GetListOfFunctions());
                    TObject* fobj;
                    while ((fobj = nextFunc())) {
                        if (fobj->InheritsFrom("TLatex")) {
                            TLatex* latx = (TLatex*)fobj;
                            latx->SetX(bx);
                            if (idx < 3) latx->SetY(by[idx]);
                            idx++;
                        }
                    }
                }

                hc->Draw("PE"); // + all attached fit curves / repositioned Yield block

                // cent / q-bin / pT label, top-left INSIDE the frame
                TLatex lat1, lat2;
                lat1.SetNDC(); lat2.SetNDC();
                lat1.SetTextFont(42); lat2.SetTextFont(42);
                lat1.SetTextSize(0.050); lat2.SetTextSize(0.050);
                lat1.SetTextAlign(12); lat2.SetTextAlign(12);
                lat1.DrawLatex(0.20, 0.87,
                    Form("Cent : %s (q bin =%d)", centLabel[ic], iq));
                lat2.DrawLatex(0.20, 0.80,
                    Form("%s <p_{T}< %s GeV", ptLo[ip], ptHi[ip]));

                // "CMS Preliminary" (top-left) / "PbPb sqrt(sNN)=5.36 TeV" (top-right)
                DrawCMSHeader(0.16, 0.05, 0.925, /*isMC=*/false);

                cv->Modified();
                cv->Update();
                cv->SaveAs(Form("%s/d_c%d_p%d_q%d.pdf",
                                dataTmpDir.Data(), ic, ip, iq));
                delete hc;
            }
            Printf("Saved data PDFs: cent %s  pT %s", centLabel[ic], ptLabel[ip]);
        }
    }

    // ========================================================================
    // MC: h_mc_sig_v2_<ic>_<ip>_<ip>_cent<CENT>   (45 plots)
    // Style: black points instead of the saved "hist" style, fit curve forced
    //        to blue (saved as red), no stat/fit box. "Cent / pT" legend
    //        sits top-left INSIDE the frame, same style as data.
    // ========================================================================
    for (int ic = 0; ic < NCENT; ic++) {
        for (int ip = 0; ip < NPT; ip++) {

            cv->Clear();
            cv->SetLeftMargin  (0.16);
            cv->SetRightMargin (0.05);
            cv->SetBottomMargin(0.14);
            cv->SetTopMargin   (0.10);

            TString hname = Form("h_mc_sig_v2_%d_%d_%d_%s", ic, ip, ip, centTag[ic]);
            TH1D* h = (TH1D*)d->Get(hname);

            if (!h) {
                Printf("WARNING: %s not found", hname.Data());
                TLatex miss;
                miss.SetNDC();
                miss.SetTextFont(42);
                miss.SetTextSize(0.07);
                miss.SetTextAlign(22);
                miss.DrawLatex(0.5, 0.5, "no data");
                cv->SaveAs(Form("%s/m_c%d_p%d.pdf", mcTmpDir.Data(), ic, ip));
                continue;
            }

            TH1D* hc = (TH1D*)h->Clone(Form("hc_m_%d_%d", ic, ip));
            hc->SetDirectory(nullptr);

            hc->GetXaxis()->SetTitle("Mass (GeV)");
            hc->GetYaxis()->SetTitle("Entries");
            hc->GetXaxis()->SetRangeUser(1.75, 1.98); // match data plot range
            hc->GetXaxis()->SetTitleSize (0.055);
            hc->GetYaxis()->SetTitleSize (0.055);
            hc->GetXaxis()->SetLabelSize (0.050);
            hc->GetYaxis()->SetLabelSize (0.050);
            hc->GetXaxis()->SetTitleOffset(0.95);
            hc->GetYaxis()->SetTitleOffset(1.15);
            hc->GetXaxis()->SetNdivisions(505);

            hc->SetMarkerStyle(20);
            hc->SetMarkerSize(0.7);
            hc->SetMarkerColor(kBlack);
            hc->SetLineColor(kBlack);

            // fit curve in blue (saved fit is red)
            {
                TIter it(hc->GetListOfFunctions());
                TObject* o;
                while ((o = it())) {
                    if (o->InheritsFrom("TF1")) {
                        ((TF1*)o)->SetLineColor(kBlue);
                        ((TF1*)o)->SetLineWidth(2);
                    }
                }
            }

            hc->Draw("PE");

            // cent / pT legend, top-left INSIDE the frame (same style as data)
            TLatex lat1, lat2;
            lat1.SetNDC(); lat2.SetNDC();
            lat1.SetTextFont(42); lat2.SetTextFont(42);
            lat1.SetTextSize(0.050); lat2.SetTextSize(0.050);
            lat1.SetTextAlign(12); lat2.SetTextAlign(12);
            lat1.DrawLatex(0.20, 0.87, Form("Cent : %s", centLabel[ic]));
            lat2.DrawLatex(0.20, 0.80,
                Form("%s <p_{T}< %s GeV", ptLo[ip], ptHi[ip]));

            // "CMS Simulation" (top-left) / "PbPb sqrt(sNN)=5.36 TeV" (top-right)
            DrawCMSHeader(0.16, 0.05, 0.925, /*isMC=*/true);

            cv->Modified();
            cv->Update();
            cv->SaveAs(Form("%s/m_c%d_p%d.pdf", mcTmpDir.Data(), ic, ip));
            delete hc;
        }
        Printf("Saved MC PDFs: cent %s", centLabel[ic]);
    }

    delete cv;
    f->Close();

    // ========================================================================
    // Tile both sets into their final multi-page PDFs
    // ========================================================================
    Printf("Tiling data PDFs (45 pages, 12/page) ...");
    gSystem->Exec("/usr/local/bin/python3.11 tile_mass_pdfs.py data");

    Printf("Tiling MC PDFs (5 pages, 9/page) ...");
    gSystem->Exec("/usr/local/bin/python3.11 tile_mass_pdfs.py mc");

    Printf("Done.");
}
