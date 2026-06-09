// plot_mass_qbins.C
//
// Plot D0 mass histograms (with fits) from the mass_fitted directory.
// Strategy:
//   1) Save each of the 450 individual q-bin plots as a single-pad PDF
//      in mass_pdfs_tmp/.
//   2) Use pdfjam (--nup 5x2 --landscape) to tile 10 q-bin PDFs per page.
//   3) Use gs (Ghostscript) to merge all 45 pages into one final PDF.
//
// This avoids TPaveText / TLegend repositioning issues that arise in
// multi-pad (Divide) layouts — each plot is saved exactly as ROOT
// renders it on a single canvas.
//
// Requirements: pdfjam (part of TeXLive) and gs (Ghostscript) in PATH.
//   macOS:  brew install ghostscript  &&  brew install pdfjam   (or via MacTeX)
//
// Input:  ROOT/Flow_MB0to31_May8_out_combined.root   (mass_fitted dir)
// Output: plot_mass_q2bins_May12.pdf
// Usage:  root -l -q plot_mass_qbins.C

#include "TFile.h"
#include "TDirectoryFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

void setMassStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
}

void plot_mass_qbins()
{
    setMassStyle();
    gROOT->ForceStyle();

    // ── open file ─────────────────────────────────────────────────────────────
    TFile* f = TFile::Open("ROOT/Flow_MB0to31_May8_out_combined.root");
    if (!f || f->IsZombie()) {
        Printf("ERROR: cannot open ROOT/Flow_MB0to31_May8_out_combined.root");
        return;
    }
    TDirectoryFile* d = (TDirectoryFile*)f->Get("mass_fitted");
    if (!d) {
        Printf("ERROR: mass_fitted directory not found");
        f->Close(); return;
    }

    // ── definitions ───────────────────────────────────────────────────────────
    const int   NCENT              = 5;
    const char* centTag  [NCENT]   = { "cent0to10", "cent10to20", "cent20to30",
                                        "cent30to40", "cent40to50" };
    const char* centLabel[NCENT]   = { "0-10%", "10-20%", "20-30%",
                                        "30-40%", "40-50%" };

    const int   NPT                = 9;
    const char* ptTag    [NPT]     = { "pT2to3",  "pT3to4",  "pT4to5",
                                        "pT5to6",  "pT6to8",  "pT8to10",
                                        "pT10to15","pT15to30","pT30to100" };
    const char* ptLabel  [NPT]     = { "2-3",  "3-4",  "4-5",
                                        "5-6",  "6-8",  "8-10",
                                        "10-15","15-30","30-100" };

    const int NQ = 10;

    // ── Step 1: save each q-bin plot as an individual single-pad PDF ──────────
    TString tmpDir = "mass_pdfs_tmp";
    gSystem->mkdir(tmpDir, kTRUE);

    // Single-pad canvas — ROOT places TPaveText/TLegend correctly here
    TCanvas* cv = new TCanvas("cv_single", "", 600, 500);

    for (int ic = 0; ic < NCENT; ic++) {
        for (int ip = 0; ip < NPT; ip++) {
            for (int iq = 0; iq < NQ; iq++) {

                cv->Clear();
                cv->SetLeftMargin  (0.14);
                cv->SetRightMargin (0.05);
                cv->SetBottomMargin(0.13);
                cv->SetTopMargin   (0.16);

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
                    cv->SaveAs(Form("%s/c%d_p%d_q%d.pdf",
                                    tmpDir.Data(), ic, ip, iq));
                    continue;
                }

                TH1D* hc = (TH1D*)h->Clone(Form("hc_%d_%d_%d", ic, ip, iq));
                hc->SetDirectory(nullptr);

                // axis cosmetics
                hc->GetXaxis()->SetTitleSize (0.060);
                hc->GetYaxis()->SetTitleSize (0.060);
                hc->GetXaxis()->SetLabelSize (0.055);
                hc->GetYaxis()->SetLabelSize (0.055);
                hc->GetXaxis()->SetTitleOffset(0.90);
                hc->GetYaxis()->SetTitleOffset(1.05);
                hc->GetXaxis()->SetNdivisions(505);

                // reduce TF1 line widths for PDF legibility
                {
                    TIter it(hc->GetListOfFunctions());
                    TObject* o;
                    while ((o = it())) {
                        if (o->InheritsFrom("TF1"))
                            ((TF1*)o)->SetLineWidth(1);
                    }
                }

                hc->SetMarkerStyle(20);
                hc->SetMarkerSize(0.7);
                hc->Draw("E");
                cv->Update();

                // header: cent + pT + q-bin label (top centre, above frame)
                TLatex lat;
                lat.SetNDC();
                lat.SetTextFont(62);
                lat.SetTextSize(0.048);
                lat.SetTextAlign(22);
                lat.DrawLatex(0.50, 0.96,
                    Form("Cent %s  p_{T} %s GeV/c  q-bin %d",
                         centLabel[ic], ptLabel[ip], iq));

                cv->Modified();
                cv->Update();
                cv->SaveAs(Form("%s/c%d_p%d_q%d.pdf",
                                tmpDir.Data(), ic, ip, iq));
                delete hc;
            }
            Printf("Saved individual PDFs: cent %s  pT %s",
                   centLabel[ic], ptLabel[ip]);
        }
    }

    delete cv;
    f->Close();

    // ── Step 2: tile into 45-page PDF via Python / PyMuPDF ───────────────────
    // tile_mass_pdfs.py arranges 10 individual PDFs per page (5×2 landscape)
    // and writes the final output file.
    Printf("Tiling individual PDFs into 45-page PDF ...");
    // Use python3.11 which has PyMuPDF (fitz) installed
    int pyRet = gSystem->Exec("/usr/local/bin/python3.11 tile_mass_pdfs.py");
    if (pyRet != 0)
        Printf("ERROR: tile_mass_pdfs.py failed (exit %d). "
               "Individual PDFs are in %s/", pyRet, tmpDir.Data());
    else
        Printf("Saved: plot_mass_q2bins_May12.pdf");

    Printf("Done.");
}
