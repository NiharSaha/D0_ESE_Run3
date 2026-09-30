// qn_Compare_2018vs2023.C
// Compares q2 and q3 raw distributions for each 1% centrality bin (0-50%)
// between 2018 PbPb and 2023 PbPb datasets.
//
// Usage (in ROOT):
//   .x qn_Compare_2018vs2023.C
//   .x qn_Compare_2018vs2023.C("ROOT/Quantiles_MB0to31_out_combined_May27.root",
//                               "ROOT/output_PbPb2018_total.root", 10)
//
// plotsPerCanvas: how many 1% bins to show per canvas (default 10)

#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TGaxis.h"

void qn_Compare_2018vs2023(
    TString file2023       = "Quantiles_MB0to31_12NUQbin_out_combined_Sept14.root",
    TString file2018       = "output_PbPb2018_total.root",
    int     plotsPerCanvas = 10
) {
    gStyle->SetOptStat(0);
    TGaxis::SetMaxDigits(4);

    TFile *f23 = TFile::Open(file2023);
    if (!f23 || f23->IsZombie()) {
        std::cerr << "Cannot open 2023 file: " << file2023 << std::endl;
        return;
    }
    TFile *f18 = TFile::Open(file2018);
    if (!f18 || f18->IsZombie()) {
        std::cerr << "Cannot open 2018 file: " << file2018 << std::endl;
        return;
    }

    const int N_CENT = 50;   // 0-50% in 1% steps

    

    // Retrieve raw q-vector histogram from a file.
    // Tries new path "Slices/hRawQ2_c{j}" first, then legacy fallback.
    auto getHist = [&](TFile *f, const char *qUp, int j) -> TH1D* {
        TString qLow(qUp); qLow.ToLower();
        // 1. 2023 new format: Slices/hRawQ2_c{j}
        TH1D *h = (TH1D*)f->Get(Form("Slices/hRaw%s_c%d", qUp, j));
        // 2. 2018 format: q2_hfpm_cen_{j}_{j+1}
        if (!h)
            h = (TH1D*)f->Get(Form("%s_hfpm_cen_%d_%d", qLow.Data(), j, j+1));
        return h;
    };

    TLatex ltx;
    ltx.SetNDC();
    ltx.SetTextFont(42);

    // Canvas grid
    int nCols = (plotsPerCanvas <= 5) ? plotsPerCanvas : 5;
    int nRows = (plotsPerCanvas + nCols - 1) / nCols;

    for (const char *qUp : {"Q2", "Q3"}) {
        TString axTitle = (TString(qUp) == "Q2") ? "q_{2}" : "q_{3}";

        int nCanvases = (N_CENT + plotsPerCanvas - 1) / plotsPerCanvas;

        for (int iC = 0; iC < nCanvases; iC++) {
            int jStart = iC * plotsPerCanvas;
            int jEnd   = std::min(jStart + plotsPerCanvas, N_CENT);
            int nPads  = jEnd - jStart;
            int totalPads = nCols * nRows;

            TCanvas *c = new TCanvas(
                Form("c_%s_%d", qUp, iC),
                Form("2018 vs 2023 - %s  cent %d-%d%%", axTitle.Data(), jStart, jEnd),
                300 * nCols, 300 * nRows);
            c->Divide(nCols, nRows);

            for (int p = 0; p < nPads; p++) {
                int j = jStart + p;
                c->cd(p + 1);
                gPad->SetLeftMargin(0.16);
                gPad->SetBottomMargin(0.16);
                gPad->SetRightMargin(0.04);
                gPad->SetTopMargin(0.06);

                TH1D *h23raw = getHist(f23, qUp, j);
                TH1D *h18raw = getHist(f18, qUp, j);

                bool has23 = h23raw && h23raw->GetEntries() > 0;
                bool has18 = h18raw && h18raw->GetEntries() > 0;

                if (!has23 && !has18) {
                    ltx.SetTextAlign(22); ltx.SetTextSize(0.07);
                    ltx.DrawLatex(0.5, 0.5, Form("No data cent %d%%", j));
                    continue;
                }

                // Clone and rebin 2018 histogram if needed
                TH1D *h23n = nullptr;
                TH1D *h18n = nullptr;
                if (has23) {
                    h23n = (TH1D*)h23raw->Clone(Form("h23n_%s_%d", qUp, j));
                    h23n->SetDirectory(0);
                }
                if (has18) {
                    h18n = (TH1D*)h18raw->Clone(Form("h18n_%s_%d", qUp, j));
                    h18n->SetDirectory(0);
                    // Rebin 2018 to match 2023 if it's finer
                    if (has23 && h18n->GetNbinsX() > h23n->GetNbinsX()) {
                        int factor = h18n->GetNbinsX() / h23n->GetNbinsX();
                        if (factor > 1) {
                            h18n->Rebin(factor);
                        }
                    }
                }

                // Normalize to unit area
                if (h23n && h23n->Integral() > 0) h23n->Scale(1.0 / h23n->Integral());
                if (h18n && h18n->Integral() > 0) h18n->Scale(1.0 / h18n->Integral());

                // Determine combined x and y range so both histograms are visible
                double xMin = 0;
                double xMax = (strcmp(qUp, "Q2") == 0) ? 0.3 : 0.2;
                double yMax = 0;
                auto updateRanges = [&](TH1D *h) {
                    if (!h) return;
                    yMax = std::max(yMax, h->GetMaximum());
                };
                updateRanges(h23n);
                updateRanges(h18n);

                // Style helper
                auto applyStyle = [&](TH1D *h, int col, int lstyle) {
                    if (!h) return;
                    h->SetLineColor(col);
                    h->SetLineWidth(1);
                    h->SetLineStyle(lstyle);
                    h->SetTitle("");
                    h->GetXaxis()->SetTitle(axTitle);
                    h->GetYaxis()->SetTitle("Normalised entries");
                    h->GetXaxis()->SetTitleSize(0.07);
                    h->GetYaxis()->SetTitleSize(0.07);
                    h->GetXaxis()->SetLabelSize(0.063);
                    h->GetYaxis()->SetLabelSize(0.063);
                    h->GetXaxis()->SetTitleOffset(0.95);
                    h->GetYaxis()->SetTitleOffset(1.1);
                    if (strcmp(qUp, "Q3") == 0) {
                        h->GetXaxis()->SetNdivisions(505); // Reduce divisions for smaller range
                    }
                };
                applyStyle(h23n, kRed+1, 1);
                applyStyle(h18n, kBlue+1,  1);

                // Draw: set axis range to cover both before drawing second
                if (h23n) {
                    h23n->GetXaxis()->SetRangeUser(xMin, xMax);
                    h23n->SetMaximum(yMax * 1.25);
                    h23n->Draw("HIST");
                }
                if (h18n) {
                    h18n->GetXaxis()->SetRangeUser(xMin, xMax);
                    if (!h23n) { h18n->SetMaximum(yMax * 1.25); h18n->Draw("HIST"); }
                    else        h18n->Draw("HIST SAME");
                }

                // Per-pad legend (top-right)
                TLegend *leg = new TLegend(0.55, 0.80, 0.98, 0.93);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.065);
                if (h23n) leg->AddEntry(h23n, "2023 PbPb", "l");
                if (h18n) leg->AddEntry(h18n, "2018 PbPb", "l");
                leg->Draw();

                // Centrality label in a box (top-left)
                TPaveText *pt = new TPaveText(0.17, 0.85, 0.5, 0.93, "NDC");
                pt->SetBorderSize(0);
                pt->SetFillColor(0);
                pt->SetFillStyle(0);
                pt->SetTextFont(42);
                pt->SetTextSize(0.065);
                pt->SetTextAlign(12); // Align left
                pt->AddText(Form("Cent %d-%d%%", j, j+1));
                pt->Draw();
            }

            // Any remaining empty pads in the grid — leave blank
            c->SaveAs(Form("Compare_%s_cent%02dto%02d.pdf", qUp, jStart, jEnd - 1));
        }
    }

    std::cout << "Done. PDFs saved for Q2 and Q3." << std::endl;
}
