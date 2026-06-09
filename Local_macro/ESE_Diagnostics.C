#include <iostream>
#include <vector>
#include <fstream>
#include "TFile.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TGaxis.h"

/*
OPTION 1: Overlay 1% raw bins
OPTION 2: Stacked 1% slices
OPTION 3: Global Cut Generation
OPTION 4: Merged Quantiles
*/

void ESE_Diagnostics(
             int option = 4,
             TString inputMergedFile = "ROOT/Quantiles_MB0to31_out_combined_May27.root"
             ) {


  gStyle->SetOptStat(0);
    TGaxis::SetMaxDigits(4);
    TFile *fIn = TFile::Open(inputMergedFile);
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "Error: Cannot open " << inputMergedFile << std::endl;
        return;
    }

    const int N_CENTBIN = 5;
    Int_t min_centbin[N_CENTBIN] = {0, 10, 20, 30, 40};
    Int_t max_centbin[N_CENTBIN] = {10, 20, 30, 40, 50};
    const int N_CENTBIN_1PC = 90;
    const int N_Q_BINS = 10;

    double probs[N_Q_BINS + 1];
    for (int i = 0; i <= N_Q_BINS; ++i) probs[i] = (double)i / N_Q_BINS;

    // The ROOT file was produced with the original 4-bin centrality structure
    // {0-10, 10-30, 30-50, 50-90}.  Histogram names encode that old coarse index
    // (e.g. cent1 covers fine bins 10..29).  Map each 1%-fine-bin j to the
    // correct old cent index so fIn->Get() finds the right histogram.
    auto fileCentIdx = [](int j) -> int {
        if (j <  10) return 0;
        if (j <  30) return 1;
        if (j <  50) return 2;
        return 3;
    };

    TLatex ltx;
    ltx.SetNDC();
    ltx.SetTextFont(42);

    // 10 visually distinct colors shared across all options
    const int qColors[N_Q_BINS] = {
        kRed,       // bin 1
        kBlue,      // bin 2
        kGreen+2,   // bin 3
        kOrange+1,  // bin 4
        kMagenta,   // bin 5
        kCyan+2,    // bin 6
        kViolet,    // bin 7
        kSpring-6,  // bin 8 (olive/dark yellow-green)
        kPink+6,    // bin 9
        kAzure+7    // bin 10
    };

    // =================================================================
    // OPTION 1: Overlay 1% raw bins (Separate q2 and q3 canvases)
    // =================================================================

if (option == 1) {
    // 1. Initialize canvases for Q2 and Q3
    TCanvas *c1_q2 = new TCanvas("c1_q2", "Overlay 1% Bins - Q2 (All Centralities)", 1500, 900);
    c1_q2->Divide(3, 2);
    TCanvas *c1_q3 = new TCanvas("c1_q3", "Overlay 1% Bins - Q3 (All Centralities)", 1500, 900);
    c1_q3->Divide(3, 2);

    // Use the shared color palette
    std::vector<int> baseColors(qColors, qColors + N_Q_BINS);

    for (int i = 0; i < N_CENTBIN; i++) {
        // --- Process Q2 ---
        c1_q2->cd(i + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        
        TLegend *leg2 = new TLegend(0.6, 0.45, 0.93, 0.9);
        leg2->SetHeader(Form("q_{2} for cent: %d-%d%%", min_centbin[i], max_centbin[i]));
        leg2->SetBorderSize(0);
        leg2->SetTextSize(0.035);
        leg2->SetNColumns(2); 

        bool first2 = true;
        for (int j = min_centbin[i]; j < max_centbin[i]; j++) {
            // Try new combined file structure first (Slices/...), then fall back to old layout
            TH1D *h = nullptr;
            h = (TH1D*)fIn->Get(Form("Slices/hRawQ2_c%d", j));
            if (!h) h = (TH1D*)fIn->Get(Form("q2_raw_1pc/hist_q2_tot_cent%d_c%d", fileCentIdx(j), j));
            if (!h || h->GetEntries() == 0) continue;

            // Apply high-contrast colors and axis titles
            int colorIdx = (j - min_centbin[i]) % (int)baseColors.size();
            int color = baseColors[colorIdx];
            if ((j - min_centbin[i]) >= (int)baseColors.size()) color += 2; 

            h->SetLineColor(color);
            h->SetLineWidth(2);
            h->SetTitle(""); // Clear default title
            h->GetXaxis()->SetTitle("q_{2}");
            h->GetYaxis()->SetTitle("Entries");
            h->GetXaxis()->SetTitleSize(0.062);
            h->GetYaxis()->SetTitleSize(0.062);
            h->GetXaxis()->SetLabelSize(0.056);
            h->GetYaxis()->SetLabelSize(0.056);
            h->GetXaxis()->SetTitleOffset(1.0);
            h->GetYaxis()->SetTitleOffset(1.2);
            //h->GetXaxis()->SetRangeUser(0, 0.18);
            
            h->Draw(first2 ? "HIST" : "HIST SAME");
            if (j % 1 == 0) leg2->AddEntry(h, Form("%d-%d%%", j, j+1), "l");
            first2 = false;
        }
        leg2->Draw();
        ltx.SetTextFont(62); ltx.SetTextSize(0.065); ltx.SetTextAlign(11);
        ltx.DrawLatex(0.25, 0.92, "CMS");
        ltx.SetTextFont(42); ltx.SetTextSize(0.055); ltx.SetTextAlign(31);
        ltx.DrawLatex(0.90, 0.92, "PbPb 5.36 TeV");

        // --- Process Q3 ---
        c1_q3->cd(i + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);

        TLegend *leg3 = new TLegend(0.6, 0.45, 0.93, 0.9);
        leg3->SetHeader(Form("q_{3} for cent: %d-%d%%", min_centbin[i], max_centbin[i]));
        leg3->SetBorderSize(0);
        leg3->SetTextSize(0.035);
        leg3->SetNColumns(2);

        bool first3 = true;
        for (int j = min_centbin[i]; j < max_centbin[i]; j++) {
            TH1D *h = nullptr;
            h = (TH1D*)fIn->Get(Form("Slices/hRawQ3_c%d", j));
            if (!h) h = (TH1D*)fIn->Get(Form("q3_raw_1pc/hist_q3_tot_cent%d_c%d", fileCentIdx(j), j));
            if (!h || h->GetEntries() == 0) continue;

            int colorIdx = (j - min_centbin[i]) % (int)baseColors.size();
            int color = baseColors[colorIdx];
            if ((j - min_centbin[i]) >= (int)baseColors.size()) color += 2;

            h->SetLineColor(color);
            h->SetLineWidth(2);
            h->SetTitle("");
            h->GetXaxis()->SetTitle("q_{3}");
            h->GetYaxis()->SetTitle("Entries");
            h->GetXaxis()->SetTitleSize(0.062);
            h->GetYaxis()->SetTitleSize(0.062);
            h->GetXaxis()->SetLabelSize(0.056);
            h->GetYaxis()->SetLabelSize(0.056);
            h->GetXaxis()->SetTitleOffset(1.0);
            h->GetYaxis()->SetTitleOffset(1.2);
            //h->GetXaxis()->SetRangeUser(0, 0.22);
            
            h->Draw(first3 ? "HIST" : "HIST SAME");
            if (j % 1 == 0) leg3->AddEntry(h, Form("%d-%d%%", j, j+1), "l");
            first3 = false;
        }
        leg3->Draw();
        ltx.SetTextFont(62); ltx.SetTextSize(0.065); ltx.SetTextAlign(11);
        ltx.DrawLatex(0.25, 0.92, "CMS");
        ltx.SetTextFont(42); ltx.SetTextSize(0.055); ltx.SetTextAlign(31);
        ltx.DrawLatex(0.90, 0.92, "PbPb 5.36 TeV");
    }
}
  






    

    // =================================================================
    // OPTION 2: Stacked 1% slices (Warning: Generates many windows)
    // =================================================================
    if (option == 2) {
        // To keep q2 and q3 together for the same centrality 1% bin
        for (int j = 0; j < N_CENTBIN_1PC; j++) {
            int i_idx = -1;
            for(int i=0; i<N_CENTBIN; i++) if(j >= min_centbin[i] && j < max_centbin[i]) i_idx = i;
            if (i_idx == -1) continue;

            TH1D *hQ2 = (TH1D*)fIn->Get(Form("Slices/hRawQ2_c%d", j));
            if (!hQ2) hQ2 = (TH1D*)fIn->Get(Form("q2_raw_1pc/hist_q2_tot_cent%d_c%d", fileCentIdx(j), j));
            TH1D *hQ3 = (TH1D*)fIn->Get(Form("Slices/hRawQ3_c%d", j));
            if (!hQ3) hQ3 = (TH1D*)fIn->Get(Form("q3_raw_1pc/hist_q3_tot_cent%d_c%d", fileCentIdx(j), j));
            if (!hQ2 || hQ2->GetEntries() == 0) continue;

            TCanvas *c2 = new TCanvas(Form("c2_cent%d", j), Form("Quantile Slices Cent %d%%", j), 1000, 500);
            c2->Divide(2,1);

            // Left side: Q2
            c2->cd(1);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            hQ2->SetTitle("");
            hQ2->GetXaxis()->SetTitle("q_{2}");
            hQ2->GetYaxis()->SetTitle("Entries");
            hQ2->GetXaxis()->SetTitleSize(0.062);
            hQ2->GetYaxis()->SetTitleSize(0.062);
            hQ2->GetXaxis()->SetLabelSize(0.056);
            hQ2->GetYaxis()->SetLabelSize(0.056);
            hQ2->GetXaxis()->SetTitleOffset(1.0);
            hQ2->GetYaxis()->SetTitleOffset(1.2);
            hQ2->SetLineColor(kBlack);
            hQ2->SetLineWidth(2);
            hQ2->Draw("hist");
            // Prefer precomputed slice histograms if present
            bool hasPrecomputedSlices = true;
            for (int k = 0; k < N_Q_BINS; ++k) {
                if (!fIn->Get(Form("Slices/hSliceQ2_c%d_q%d", j, k+1))) { hasPrecomputedSlices = false; break; }
            }
            if (hasPrecomputedSlices) {
                for (int k = 0; k < N_Q_BINS; k++) {
                    TH1D *s = (TH1D*)fIn->Get(Form("Slices/hSliceQ2_c%d_q%d", j, k+1));
                    if (!s) continue;
                    int color = qColors[k];
                    s->SetLineColor(color);
                    s->SetFillColorAlpha(color, 0.35);
                    s->Draw("hist same");
                }
            } else {
                double cuts2[N_Q_BINS + 1];
                hQ2->GetQuantiles(N_Q_BINS + 1, cuts2, probs);
                for (int k = 0; k < N_Q_BINS; k++) {
                    int color = qColors[k];
                    TH1D *s = (TH1D*)hQ2->Clone(Form("sQ2_j%d_k%d", j, k));
                    for (int b = 1; b <= s->GetNbinsX(); b++) {
                        double center = s->GetBinCenter(b);
                        if (center < cuts2[k] || center >= cuts2[k+1]) {
                            s->SetBinContent(b, 0);
                            s->SetBinError(b, 0);
                        }
                    }
                    s->SetLineColor(color);
                    s->SetFillColorAlpha(color, 0.35);
                    s->Draw("hist same");
                }
            }
            ltx.SetTextFont(62); ltx.SetTextSize(0.065); ltx.SetTextAlign(11);
            ltx.DrawLatex(0.25, 0.92, "CMS");
            ltx.SetTextFont(42); ltx.SetTextSize(0.055); ltx.SetTextAlign(31);
            ltx.DrawLatex(0.90, 0.92, "PbPb 5.36 TeV");

            // Right side: Q3
            c2->cd(2);
            if(hQ3) {
                gPad->SetLeftMargin(0.15);
                gPad->SetBottomMargin(0.15);
                hQ3->SetTitle("");
                hQ3->GetXaxis()->SetTitle("q_{3}");
                hQ3->GetYaxis()->SetTitle("Entries");
                hQ3->GetXaxis()->SetTitleSize(0.062);
                hQ3->GetYaxis()->SetTitleSize(0.062);
                hQ3->GetXaxis()->SetLabelSize(0.056);
                hQ3->GetYaxis()->SetLabelSize(0.056);
                hQ3->GetXaxis()->SetTitleOffset(1.0);
                hQ3->GetYaxis()->SetTitleOffset(1.2);
                hQ3->SetLineColor(kBlack);
                hQ3->SetLineWidth(2);
                hQ3->Draw("hist");
                bool hasPrecomputedSlices3 = true;
                for (int k = 0; k < N_Q_BINS; ++k) {
                    if (!fIn->Get(Form("Slices/hSliceQ3_c%d_q%d", j, k+1))) { hasPrecomputedSlices3 = false; break; }
                }
                if (hasPrecomputedSlices3) {
                    for (int k = 0; k < N_Q_BINS; k++) {
                        TH1D *s = (TH1D*)fIn->Get(Form("Slices/hSliceQ3_c%d_q%d", j, k+1));
                        if (!s) continue;
                        int color = qColors[k];
                        s->SetLineColor(color);
                        s->SetFillColorAlpha(color, 0.35);
                        s->Draw("hist same");
                    }
                } else {
                    double cuts3[N_Q_BINS + 1];
                    hQ3->GetQuantiles(N_Q_BINS + 1, cuts3, probs);
                    for (int k = 0; k < N_Q_BINS; k++) {
                        int color = qColors[k];
                        TH1D *s = (TH1D*)hQ3->Clone(Form("sQ3_j%d_k%d", j, k));
                        for (int b = 1; b <= s->GetNbinsX(); b++) {
                            double center = s->GetBinCenter(b);
                            if (center < cuts3[k] || center >= cuts3[k+1]) {
                                s->SetBinContent(b, 0);
                                s->SetBinError(b, 0);
                            }
                        }
                        s->SetLineColor(color);
                        s->SetFillColorAlpha(color, 0.35);
                        s->Draw("hist same");
                    }
                }
                ltx.SetTextFont(62); ltx.SetTextSize(0.065); ltx.SetTextAlign(11);
                ltx.DrawLatex(0.25, 0.92, "CMS");
                ltx.SetTextFont(42); ltx.SetTextSize(0.055); ltx.SetTextAlign(31);
                ltx.DrawLatex(0.90, 0.92, "PbPb 5.36 TeV");
            }
        }
    }

    // =================================================================
    // OPTION 3: Global Cut Generation (ROOT File + Header File)
    // =================================================================
    if (option == 3) {

    TFile *fOut = new TFile("ESE_GlobalCuts.root", "RECREATE");
    std::ofstream hFile("ESE_Cuts_Jan18.h"); // Create the header file

        // Write Header File Guard and Comments
        hFile << "#ifndef ESE_CUTS_H" << endl;
        hFile << "#define ESE_CUTS_H" << endl << endl;
        hFile << "/* Auto-generated ESE Global Cuts */" << endl;
        hFile << "/* Rows: 90 Centrality Bins (1% each) */" << endl;
        hFile << "/* Columns: 11 Quantile Boundaries (0%, 10%, ..., 100%) */" << endl << endl;

        double q2_cuts_table[N_CENTBIN_1PC][N_Q_BINS + 1];
        double q3_cuts_table[N_CENTBIN_1PC][N_Q_BINS + 1];

        cout << ">>> Extracting Global Quantiles and generating ESE_Cuts.h..." << endl;
        
        for (int i = 0; i < N_CENTBIN; i++) {
            for (int j = min_centbin[i]; j < max_centbin[i]; j++) {
                    TH1D *hQ2 = (TH1D*)fIn->Get(Form("Slices/hRawQ2_c%d", j));
                    if (!hQ2) hQ2 = (TH1D*)fIn->Get(Form("q2_raw_1pc/hist_q2_tot_cent%d_c%d", fileCentIdx(j), j));
                    TH1D *hQ3 = (TH1D*)fIn->Get(Form("Slices/hRawQ3_c%d", j));
                    if (!hQ3) hQ3 = (TH1D*)fIn->Get(Form("q3_raw_1pc/hist_q3_tot_cent%d_c%d", fileCentIdx(j), j));
                
                if (hQ2 && hQ2->GetEntries() > 0) {
                    hQ2->GetQuantiles(N_Q_BINS + 1, q2_cuts_table[j], probs);
                    TAxis *ax = new TAxis(N_Q_BINS, q2_cuts_table[j]);
                    ax->Write(Form("q2_axis_cent%d", j));
                } else {
                    // Try to read precomputed decile slices and derive boundaries from them
                    bool foundSlices = true;
                    for (int k = 0; k < N_Q_BINS; ++k) {
                        TH1D *hs = (TH1D*)fIn->Get(Form("Slices/hSliceQ2_c%d_q%d", j, k+1));
                        if (!hs) { foundSlices = false; break; }
                    }
                    if (foundSlices) {
                        // Reconstruct cut boundaries from precomputed slice histograms
                        for (int k = 0; k < N_Q_BINS; ++k) {
                            TH1D *hs = (TH1D*)fIn->Get(Form("Slices/hSliceQ2_c%d_q%d", j, k+1));
                            double lowEdge = 0.0;
                            if (hs) {
                                int firstBin = 0;
                                for (int b = 1; b <= hs->GetNbinsX(); ++b) {
                                    if (hs->GetBinContent(b) > 0) { firstBin = b; break; }
                                }
                                if (firstBin > 0) lowEdge = hs->GetBinLowEdge(firstBin);
                            }
                            q2_cuts_table[j][k] = lowEdge;
                        }
                        // upper edge of last slice
                        TH1D *hlast = (TH1D*)fIn->Get(Form("Slices/hSliceQ2_c%d_q%d", j, N_Q_BINS));
                        double upEdge = 0.0;
                        if (hlast) {
                            int lastBin = 0;
                            for (int b = hlast->GetNbinsX(); b >= 1; --b) {
                                if (hlast->GetBinContent(b) > 0) { lastBin = b; break; }
                            }
                            if (lastBin > 0) upEdge = hlast->GetBinLowEdge(lastBin) + hlast->GetBinWidth(lastBin);
                        }
                        q2_cuts_table[j][N_Q_BINS] = upEdge;
                    } else {
                        for(int k=0; k<=N_Q_BINS; k++) q2_cuts_table[j][k] = 0.0;
                    }
                }

                if (hQ3 && hQ3->GetEntries() > 0) {
                    hQ3->GetQuantiles(N_Q_BINS + 1, q3_cuts_table[j], probs);
                    TAxis *ax = new TAxis(N_Q_BINS, q3_cuts_table[j]);
                    ax->Write(Form("q3_axis_cent%d", j));
                } else {
                    bool foundSlices = true;
                    for (int k = 0; k < N_Q_BINS; ++k) {
                        TH1D *hs = (TH1D*)fIn->Get(Form("Slices/hSliceQ3_c%d_q%d", j, k+1));
                        if (!hs) { foundSlices = false; break; }
                    }
                    if (foundSlices) {
                        for (int k = 0; k < N_Q_BINS; ++k) {
                            TH1D *hs = (TH1D*)fIn->Get(Form("Slices/hSliceQ3_c%d_q%d", j, k+1));
                            double lowEdge = 0.0;
                            if (hs) {
                                int firstBin = 0;
                                for (int b = 1; b <= hs->GetNbinsX(); ++b) {
                                    if (hs->GetBinContent(b) > 0) { firstBin = b; break; }
                                }
                                if (firstBin > 0) lowEdge = hs->GetBinLowEdge(firstBin);
                            }
                            q3_cuts_table[j][k] = lowEdge;
                        }
                        TH1D *hlast = (TH1D*)fIn->Get(Form("Slices/hSliceQ3_c%d_q%d", j, N_Q_BINS));
                        double upEdge = 0.0;
                        if (hlast) {
                            int lastBin = 0;
                            for (int b = hlast->GetNbinsX(); b >= 1; --b) {
                                if (hlast->GetBinContent(b) > 0) { lastBin = b; break; }
                            }
                            if (lastBin > 0) upEdge = hlast->GetBinLowEdge(lastBin) + hlast->GetBinWidth(lastBin);
                        }
                        q3_cuts_table[j][N_Q_BINS] = upEdge;
                    } else {
                        for(int k=0; k<=N_Q_BINS; k++) q3_cuts_table[j][k] = 0.0;
                    }
                }
            }
        }

        // --- Write Q2 Array to Header ---
        hFile << "double q2_cuts[90][11] = {" << endl;
        for (int j = 0; j < N_CENTBIN_1PC; j++) {
            hFile << "  {";
            for (int k = 0; k <= N_Q_BINS; k++) {
                hFile << Form("%.8f", q2_cuts_table[j][k]) << (k == N_Q_BINS ? "" : ", ");
            }
            hFile << "}" << (j == N_CENTBIN_1PC - 1 ? "" : ",") << " // Cent " << j << "%" << endl;
        }
        hFile << "};" << endl << endl;

        // --- Write Q3 Array to Header ---
        hFile << "double q3_cuts[90][11] = {" << endl;
        for (int j = 0; j < N_CENTBIN_1PC; j++) {
            hFile << "  {";
            for (int k = 0; k <= N_Q_BINS; k++) {
                hFile << Form("%.8f", q3_cuts_table[j][k]) << (k == N_Q_BINS ? "" : ", ");
            }
            hFile << "}" << (j == N_CENTBIN_1PC - 1 ? "" : ",") << " // Cent " << j << "%" << endl;
        }
        hFile << "};" << endl << endl;
        
        hFile << "#endif" << endl;
        hFile.close();
        fOut->Close();
        cout << ">>> Successfully generated ESE_GlobalCuts.root and ESE_Cuts.h" << endl;
    }

    // =================================================================
    // OPTION 4: Merged Quantiles (Sum of 1% slices with Legend)
    // =================================================================


if (option == 4) {
 
    TCanvas *c4_q2 = new TCanvas("c4_q2", "Q2 Merged vs Slices", 1500, 900);
    c4_q2->Divide(3, 2);
    TCanvas *c4_q3 = new TCanvas("c4_q3", "Q3 Merged vs Slices", 1500, 900);
    c4_q3->Divide(3, 2);

    for (int i = 0; i < N_CENTBIN; i++) {
        TH1D *hSliceSumQ2[N_Q_BINS] = {nullptr};
        TH1D *hSliceSumQ3[N_Q_BINS] = {nullptr};
        TH1D *hTotalCoarseQ2 = nullptr;
        TH1D *hTotalCoarseQ3 = nullptr;

        bool initialized = false;

        // --- 1. Accumulate granular 1% slices ---
        for (int j = min_centbin[i]; j < max_centbin[i]; j++) {
            TH1D *h1pcQ2 = (TH1D*)fIn->Get(Form("Slices/hRawQ2_c%d", j));
            if (!h1pcQ2) h1pcQ2 = (TH1D*)fIn->Get(Form("q2_raw_1pc/hist_q2_tot_cent%d_c%d", fileCentIdx(j), j));
            TH1D *h1pcQ3 = (TH1D*)fIn->Get(Form("Slices/hRawQ3_c%d", j));
            if (!h1pcQ3) h1pcQ3 = (TH1D*)fIn->Get(Form("q3_raw_1pc/hist_q3_tot_cent%d_c%d", fileCentIdx(j), j));

            // If no coarse histogram but slice histograms exist, use one slice to obtain binning
            if ((!h1pcQ2 || h1pcQ2->GetEntries() == 0)) {
                for (int kk = 0; kk < N_Q_BINS && !h1pcQ2; ++kk) {
                    TH1D *hs = (TH1D*)fIn->Get(Form("Slices/hSliceQ2_c%d_q%d", j, kk+1));
                    if (hs) {
                        h1pcQ2 = (TH1D*)hs->Clone(Form("tmp_h1pcQ2_c%d", j));
                        h1pcQ2->Reset();
                        h1pcQ2->SetDirectory(0);
                        break;
                    }
                }
            }
            if ((!h1pcQ3 || h1pcQ3->GetEntries() == 0)) {
                for (int kk = 0; kk < N_Q_BINS && !h1pcQ3; ++kk) {
                    TH1D *hs = (TH1D*)fIn->Get(Form("Slices/hSliceQ3_c%d_q%d", j, kk+1));
                    if (hs) {
                        h1pcQ3 = (TH1D*)hs->Clone(Form("tmp_h1pcQ3_c%d", j));
                        h1pcQ3->Reset();
                        h1pcQ3->SetDirectory(0);
                        break;
                    }
                }
            }

            if (!h1pcQ2) continue;

            // Initialize total and slice histograms using the binning of the first valid input
            if (!initialized) {
                hTotalCoarseQ2 = (TH1D*)h1pcQ2->Clone(Form("hTotalQ2_c%d", i));
                hTotalCoarseQ2->Reset();
                hTotalCoarseQ2->SetDirectory(0);

                if (h1pcQ3) {
                    hTotalCoarseQ3 = (TH1D*)h1pcQ3->Clone(Form("hTotalQ3_c%d", i));
                    hTotalCoarseQ3->Reset();
                    hTotalCoarseQ3->SetDirectory(0);
                }

                for(int k=0; k < N_Q_BINS; k++) {
                    hSliceSumQ2[k] = (TH1D*)h1pcQ2->Clone(Form("hSumQ2_c%d_q%d", i, k));
                    hSliceSumQ2[k]->Reset();
                    hSliceSumQ2[k]->SetDirectory(0);

                    if (h1pcQ3) {
                        hSliceSumQ3[k] = (TH1D*)h1pcQ3->Clone(Form("hSumQ3_c%d_q%d", i, k));
                        hSliceSumQ3[k]->Reset();
                        hSliceSumQ3[k]->SetDirectory(0);
                    }
                }
                initialized = true;
            }

	    // Simple bin-by-bin addition (No "Merge" warnings now)
            hTotalCoarseQ2->Add(h1pcQ2);
            if (hTotalCoarseQ3 && h1pcQ3) hTotalCoarseQ3->Add(h1pcQ3);

            // Prefer to use precomputed per-decile slice histograms when available
            bool usedPrecomputedQ2 = false;
            for (int k = 0; k < N_Q_BINS; ++k) {
                TH1D *hs = (TH1D*)fIn->Get(Form("Slices/hSliceQ2_c%d_q%d", j, k+1));
                if (hs) {
                    hSliceSumQ2[k]->Add(hs);
                    usedPrecomputedQ2 = true;
                }
            }
            if (!usedPrecomputedQ2) {
                // Calculate Quantiles for the current 1% slice and fill by bin centers
                double cuts2[N_Q_BINS + 1];
                h1pcQ2->GetQuantiles(N_Q_BINS + 1, cuts2, probs);
                for (int k = 0; k < N_Q_BINS; k++) {
                    for(int b=1; b <= h1pcQ2->GetNbinsX(); b++) {
                        double val2 = h1pcQ2->GetBinCenter(b);
                        if(val2 >= cuts2[k] && val2 < cuts2[k+1]) 
                            hSliceSumQ2[k]->Fill(val2, h1pcQ2->GetBinContent(b));
                    }
                }
            }

            bool usedPrecomputedQ3 = false;
            for (int k = 0; k < N_Q_BINS; ++k) {
                TH1D *hs3 = (TH1D*)fIn->Get(Form("Slices/hSliceQ3_c%d_q%d", j, k+1));
                if (hs3 && hSliceSumQ3[k]) {
                    hSliceSumQ3[k]->Add(hs3);
                    usedPrecomputedQ3 = true;
                }
            }
            if (!usedPrecomputedQ3 && h1pcQ3) {
                double cuts3[N_Q_BINS + 1];
                h1pcQ3->GetQuantiles(N_Q_BINS + 1, cuts3, probs);
                for (int k = 0; k < N_Q_BINS; k++) {
                    for(int b=1; b <= h1pcQ3->GetNbinsX(); b++) {
                        double val3 = h1pcQ3->GetBinCenter(b);
                        if(val3 >= cuts3[k] && val3 < cuts3[k+1]) 
                            hSliceSumQ3[k]->Fill(val3, h1pcQ3->GetBinContent(b));
                    }
                }
            }
        }

        // --- 2. Plotting Q2 ---
        if (hTotalCoarseQ2) {
            c4_q2->cd(i + 1);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            //gPad->SetLogy();
            hTotalCoarseQ2->SetTitle("");
            hTotalCoarseQ2->GetXaxis()->SetTitle("q_{2}");
            hTotalCoarseQ2->GetYaxis()->SetTitle("Counts");
            hTotalCoarseQ2->SetLineColor(kBlack);
            hTotalCoarseQ2->SetLineWidth(3);
            hTotalCoarseQ2->GetXaxis()->SetTitleSize(0.062);
            hTotalCoarseQ2->GetYaxis()->SetTitleSize(0.062);
            hTotalCoarseQ2->GetXaxis()->SetLabelSize(0.056);
            hTotalCoarseQ2->GetYaxis()->SetLabelSize(0.056);
            hTotalCoarseQ2->GetXaxis()->SetTitleOffset(1.0);
            hTotalCoarseQ2->GetYaxis()->SetTitleOffset(1.2);
            hTotalCoarseQ2->Draw("hist");

            TLegend *leg2 = new TLegend(0.55, 0.45, 0.9, 0.9);
            leg2->SetHeader(Form("q2 in cent %d-%d%%", min_centbin[i], max_centbin[i]));
            leg2->SetBorderSize(0);
            leg2->SetNColumns(2);
            leg2->SetTextSize(0.03);
            leg2->AddEntry(hTotalCoarseQ2, "Total q_{2}", "l");

            for(int k=0; k<N_Q_BINS; k++) {
                int color = qColors[k];
                hSliceSumQ2[k]->SetLineColor(color);
                hSliceSumQ2[k]->SetLineWidth(2);
                hSliceSumQ2[k]->SetFillColorAlpha(color, 0.35);
                hSliceSumQ2[k]->Draw("hist same");
                leg2->AddEntry(hSliceSumQ2[k], Form("q_{2} bin %d", k+1), "l");
            }
            leg2->Draw();
            ltx.SetTextFont(62); ltx.SetTextSize(0.065); ltx.SetTextAlign(11);
            ltx.DrawLatex(0.25, 0.92, "CMS");
            ltx.SetTextFont(42); ltx.SetTextSize(0.055); ltx.SetTextAlign(31);
            ltx.DrawLatex(0.90, 0.92, "PbPb 5.36 TeV");
        }

	// --- 3. Plotting Q3 ---
        if (hTotalCoarseQ3) {
            c4_q3->cd(i + 1);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);
            //gPad->SetLogy();
            hTotalCoarseQ3->SetTitle("");
            hTotalCoarseQ3->GetXaxis()->SetTitle("q_{3}");
            hTotalCoarseQ3->GetYaxis()->SetTitle("Counts");
            hTotalCoarseQ3->SetLineColor(kBlack);
            hTotalCoarseQ3->SetLineWidth(3);
            hTotalCoarseQ3->GetXaxis()->SetTitleSize(0.062);
            hTotalCoarseQ3->GetYaxis()->SetTitleSize(0.062);
            hTotalCoarseQ3->GetXaxis()->SetLabelSize(0.056);
            hTotalCoarseQ3->GetYaxis()->SetLabelSize(0.056);
            hTotalCoarseQ3->GetXaxis()->SetTitleOffset(1.0);
            hTotalCoarseQ3->GetYaxis()->SetTitleOffset(1.2);
            hTotalCoarseQ3->Draw("hist");

            TLegend *leg3 = new TLegend(0.55, 0.45, 0.9, 0.9);
            leg3->SetHeader(Form("q3 in cent %d-%d%%", min_centbin[i], max_centbin[i]));
            leg3->SetBorderSize(0);
            leg3->SetNColumns(2);
            leg3->SetTextSize(0.03);
            leg3->AddEntry(hTotalCoarseQ3, "Total q_{3}", "l");

            for(int k=0; k<N_Q_BINS; k++) {
                int color = qColors[k];
                hSliceSumQ3[k]->SetLineColor(color);
                hSliceSumQ3[k]->SetLineWidth(2);
                hSliceSumQ3[k]->SetFillColorAlpha(color, 0.35);
                hSliceSumQ3[k]->Draw("hist same");
                leg3->AddEntry(hSliceSumQ3[k], Form("q_{3} bin %d", k+1), "l");
            }
            leg3->Draw();
            ltx.SetTextFont(62); ltx.SetTextSize(0.065); ltx.SetTextAlign(11);
            ltx.DrawLatex(0.25, 0.92, "CMS");
            ltx.SetTextFont(42); ltx.SetTextSize(0.055); ltx.SetTextAlign(31);
            ltx.DrawLatex(0.90, 0.92, "PbPb 5.36 TeV");
        }
    }
}



    

    
}
