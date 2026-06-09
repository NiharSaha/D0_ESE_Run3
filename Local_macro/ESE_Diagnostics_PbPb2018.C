#include <iostream>
#include <vector>
#include <fstream>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TAxis.h"

void ESE_Diagnostics_PbPb2018(int option = 4, int choice = 1) {
    gStyle->SetOptStat(0);
    
    // File name for 2018 PbPb data
    TString inputFileName = "q2_final_bin_HF_total_PbPb2018.root";
    TFile *fIn = TFile::Open(inputFileName);
    if (!fIn || fIn->IsZombie()) {
        std::cout << "Error: Cannot open " << inputFileName << std::endl;
        return;
    }

    const int N_Q_BINS = 10;
    // Rainbow colors matching Figure 5 style
    int colors[] = {kBlue+2, kBlue, kCyan, kGreen+1, kYellow-3, kOrange-3, kMagenta, kViolet, kGray+2, kBlack};

    // Centrality mapping
    int start_c, end_c;
    TString label;
    if (choice == 1)      { start_c = 0;  end_c = 10; label = "0-10%"; }
    else if (choice == 2) { start_c = 10; end_c = 30; label = "10-30%"; }
    else if (choice == 3) { start_c = 30; end_c = 50; label = "30-50%"; }
    else if (choice == 4) { start_c = 50; end_c = 90; label = "50-90%"; }
    else if (option != 5) { std::cout << "Invalid choice! Use 1, 2, 3, or 4." << std::endl; return; }

    // Probabilities for the Quantiles method
    double probs[11];
    for (int i = 0; i <= 10; ++i) probs[i] = (double)i / 10.0;

    // =================================================================
    // OPTION 1: Overlay individual 1% raw bins
    // =================================================================

if (option == 1) {
    TCanvas *c1 = new TCanvas("c1", "2018 Overlay 1% Bins - Divided Canvas", 1200, 1000);
    // Divide into 2 columns and 2 rows
    c1->Divide(2, 2); 

    // Core color palette
    std::vector<int> baseColors = {kRed, kBlue, kGreen+2, kBlack, kOrange+1, kMagenta, kCyan+1, kAzure+7, kSpring-6, kGray+1};
    std::vector<std::pair<int, int>> ranges = {{0, 10}, {10, 30}, {30, 50}};

    // Map our 3 ranges to specific pad numbers in the 2x2 grid
    // Pad 1 (top-left), Pad 2 (top-right), Pad 3 (bottom-left... we can center this later if needed)
    std::vector<int> padIndices = {1, 2, 3}; 

    for (int i_range = 0; i_range < (int)ranges.size(); i_range++) {
        c1->cd(padIndices[i_range]); // Use standard Divide numbering
        
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetTicks(1,1);

        int start_c = ranges[i_range].first;
        int end_c   = ranges[i_range].second;

        TLegend *leg = new TLegend(0.58, 0.20, 0.93, 0.88);
        leg->SetHeader(Form("q2 for cent: %d-%d%%", start_c, end_c));
        leg->SetBorderSize(0);
        leg->SetTextSize(0.035);
        leg->SetNColumns(1); // Standard column for readability

        bool firstInPad = true;

        for (int j = start_c; j < end_c; j++) {
            TString h_name = Form("h_raw_reconst_%d", j);
            TH1D *h_1pc = new TH1D(h_name, ";q_{2};Entries", 700, 0, 0.35);

            // Reconstruct the 1% bin from the 2018 flat root file
            for (int k = 0; k < 10; k++) {
                TH1D *slice = (TH1D*)fIn->Get(Form("hist_q2_HFtotal_cen%d_%d_q2_%d", j, j+1, k));
                if (slice) h_1pc->Add(slice);
            }

            // Normalizing to probability density
            //if (h_1pc->Integral() > 0) h_1pc->Scale(1.0 / h_1pc->Integral(), "width");

            // Cycle colors for high contrast
            int colorIdx = (j - start_c) % (int)baseColors.size();
            int color = baseColors[colorIdx];
            if ((j - start_c) >= (int)baseColors.size()) color += 2; 

            h_1pc->SetLineColor(color);
            h_1pc->SetLineWidth(2);
            h_1pc->GetXaxis()->SetRangeUser(0, 0.50);
            
            h_1pc->Draw(firstInPad ? "HIST" : "HIST SAME");
            
            // Only add legend entries for every 2nd bin if the range is large (10-30, 30-50)
            if ((end_c - start_c) > 10) {
                if (j % 2 == 0) leg->AddEntry(h_1pc, Form("%d-%d%%", j, j+1), "l");
            } else {
                leg->AddEntry(h_1pc, Form("%d-%d%%", j, j+1), "l");
            }
            
            firstInPad = false;
        }
        leg->Draw();
    }
    
    // Optional: Add a label in the empty 4th pad (bottom-right)
    c1->cd(4);
    TLatex tex;
    tex.SetTextSize(0.05);
    tex.DrawLatex(0.1, 0.5, "PbPb 2018 ESE Diagnostics");
}


    
    // =================================================================
    // OPTION 2: Stacked 1% Quantile Slices (Using Quantiles method)
    // =================================================================
    if (option == 2) {
        for (int j = start_c; j < end_c; j++) {
            TCanvas *c2 = new TCanvas(Form("c2_%d", j), Form("Quantile Method Slices: %d-%d%%", j, j+1), 800, 600);
            
            // Reconstruct total for this 1% bin to calculate quantiles
            TH1D *h_sum = new TH1D("h_temp", "", 700, 0, 0.35);
            for (int k = 0; k < 10; k++) {
                TH1D *s = (TH1D*)fIn->Get(Form("hist_q2_HFtotal_cen%d_%d_q2_%d", j, j+1, k));
                if (s) h_sum->Add(s);
            }

            double q_cuts[11];
            h_sum->GetQuantiles(11, q_cuts, probs); // Quantiles method

            for (int k = 0; k < N_Q_BINS; k++) {
                TH1D *h_slice = (TH1D*)fIn->Get(Form("hist_q2_HFtotal_cen%d_%d_q2_%d", j, j+1, k));
                if (!h_slice) continue;
                if (h_slice->Integral() > 0) h_slice->Scale(1.0 / h_slice->Integral(), "width");
                h_slice->SetLineColor(colors[k]);
                h_slice->Draw(k == 0 ? "HIST" : "HIST SAME");
            }
            delete h_sum;
        }
    }

    // =================================================================
    // OPTION 3: Global Cut Generation (ESE_Cuts_2018.h)
    // =================================================================
    if (option == 3) {
        ofstream hFile("ESE_Cuts_2018.h");
        hFile << "#ifndef ESE_CUTS_2018_H\n#define ESE_CUTS_2018_H\n\n";
        hFile << "double q2_cuts_2018[90][11] = {\n";
        
        for (int j = 0; j < 90; j++) {
            TH1D *h_sum = new TH1D("h_temp", "", 700, 0, 0.35);
            for (int k = 0; k < 10; k++) {
                TH1D *s = (TH1D*)fIn->Get(Form("hist_q2_HFtotal_cen%d_%d_q2_%d", j, j+1, k));
                if (s) h_sum->Add(s);
            }
            double q_cuts[11];
            if (h_sum->Integral() > 0) h_sum->GetQuantiles(11, q_cuts, probs); // Quantiles method
            else for(int p=0; p<11; p++) q_cuts[p] = 0.0;

            hFile << "  {";
            for (int p = 0; p < 11; p++) hFile << Form("%.8f", q_cuts[p]) << (p == 10 ? "" : ", ");
            hFile << "}" << (j == 89 ? "" : ",") << " // Cent " << j << "%\n";
            delete h_sum;
        }
        hFile << "};\n\n#endif";
        hFile.close();
        std::cout << ">>> Generated ESE_Cuts_2018.h using Quantiles method." << std::endl;
    }

    // =================================================================
    // OPTION 4: Merged Overlapping Quantiles (Integrated Reference View)
    // =================================================================

    if (option == 4) {
    TCanvas *c4 = new TCanvas("c4", "2018 Quantile Overlap - Raw Counts", 1200, 1000);
    // 2x2 grid for the three ranges; Pad 3 is bottom-left, Pad 4 is bottom-right
    c4->Divide(2, 2); 

    std::vector<std::pair<int, int>> ranges = {{0, 10}, {10, 30}, {30, 50}};
    // Rainbow colors for the 10 quantiles
    int colors[] = {kBlue+2, kBlue, kCyan, kGreen+1, kYellow-3, kOrange-3, kMagenta, kViolet, kGray+2, kBlack};

    for (int i_range = 0; i_range < (int)ranges.size(); i_range++) {
        c4->cd(i_range + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetTicks(1,1);

        int start_c = ranges[i_range].first;
        int end_c   = ranges[i_range].second;

        TLegend *leg4 = new TLegend(0.55, 0.45, 0.92, 0.88);
        leg4->SetHeader(Form("q_{2} in centrality %d-%d%%", start_c, end_c));
        leg4->SetBorderSize(0);
        leg4->SetTextSize(0.03);
        leg4->SetFillStyle(0);

        // --- 1. Reconstruct Raw Inclusive Background Envelope ---
        TH1D *hTotalRaw = new TH1D(Form("hTotalRaw_%d", i_range), ";q_{2};Counts", 700, 0, 0.35);
        for (int k = 0; k < 10; k++) {
            for (int j = start_c; j < end_c; j++) {
                // Accessing 2018 flat file histograms directly
                TH1D *s = (TH1D*)fIn->Get(Form("hist_q2_HFtotal_cen%d_%d_q2_%d", j, j+1, k));
                if (s) hTotalRaw->Add(s);
            }
        }
        
        // Removed scaling/normalization to keep it as raw distribution
        hTotalRaw->SetFillColorAlpha(kGray, 0.35);
        hTotalRaw->SetLineColor(kGray+1);
        hTotalRaw->GetXaxis()->SetRangeUser(0, 0.16);
        hTotalRaw->GetYaxis()->SetTitleOffset(1.6);
        hTotalRaw->Draw("HIST");
        leg4->AddEntry(hTotalRaw, "Inclusive q2 ", "f");

        // --- 2. Reconstruct and Overlay Raw Quantile Curves ---
        for (int k = 0; k < 10; k++) {
            TH1D *hQRaw = new TH1D(Form("hQRaw_%d_%d", i_range, k), "", 700, 0, 0.35);
            for (int j = start_c; j < end_c; j++) {
                TH1D *s = (TH1D*)fIn->Get(Form("hist_q2_HFtotal_cen%d_%d_q2_%d", j, j+1, k));
                if (s) hQRaw->Add(s);
            }
            
            // Maintaining raw counts to stay within inclusive boundary
            hQRaw->SetLineColor(colors[k]);
            hQRaw->SetLineWidth(2);
            hQRaw->Draw("HIST SAME");
            leg4->AddEntry(hQRaw, Form("q_{2} bin %d", k), "l");
        }
        
        leg4->Draw();
    }
    
    // Label for the fourth empty pad
    c4->cd(4);
    TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.05);
    tex.DrawLatex(0.15, 0.5, "PbPb 2018: Raw q_{2} Quantiles");
}

    // =================================================================
    // OPTION 5: All 4 Inclusive Centralities in One Canvas
    // =================================================================
    if (option == 5) {
        TCanvas *c5 = new TCanvas("c5", "2018 Inclusive q2 - All Centralities", 900, 700);
        TLegend *leg5 = new TLegend(0.6, 0.6, 0.88, 0.88);
        int c_colors[] = {kRed+1, kBlue+1, kGreen+2, kBlack};
        TString c_labels[] = {"0-10%", "10-30%", "30-50%", "50-90%"};
        int c_min[] = {0, 10, 30, 50};
        int c_max[] = {10, 30, 50, 90};

        for (int i = 0; i < 4; i++) {
            TH1D *hC = new TH1D(Form("hC_%d", i), ";q_{2};P(q_{2})", 700, 0, 0.35);
            for (int j = c_min[i]; j < c_max[i]; j++) {
                for (int k = 0; k < 10; k++) {
                    TH1D *s = (TH1D*)fIn->Get(Form("hist_q2_HFtotal_cen%d_%d_q2_%d", j, j+1, k));
                    if (s) hC->Add(s);
                }
            }
            if (hC->Integral() > 0) hC->Scale(1.0 / hC->Integral(), "width");
            hC->SetLineColor(c_colors[i]);
            hC->SetLineWidth(3);
            hC->GetXaxis()->SetRangeUser(0, 0.25);
            hC->Draw(i == 0 ? "HIST" : "HIST SAME");
            leg5->AddEntry(hC, c_labels[i], "l");
        }
        leg5->Draw();
    }
}
