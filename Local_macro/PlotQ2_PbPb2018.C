#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>

void PlotQ2_PbPb2018() {
    // 1. Configuration and Selection
    int choice = 1; // 1: 0-10%, 2: 10-30%, 3: 30-50%
    
    int start_c, end_c;
    TString label;
    
    // Define centrality ranges based on choice
    if (choice == 1) {
        start_c = 0; end_c = 10; label = "0-10%";
    } else if (choice == 2) {
        start_c = 10; end_c = 30; label = "10-30%";
    } else if (choice == 3) {
        start_c = 30; end_c = 50; label = "30-50%";
    } else {
        std::cout << "Invalid selection! Use 1, 2, or 3." << std::endl;
        return;
    }

    gStyle->SetOptStat(0);
    // Updated input file name
    TFile *f = TFile::Open("q2_final_bin_HF_total_PbPb2018.root");
    if (!f || f->IsZombie()) {
        std::cout << "Error: Cannot find q2_final_bin_HF_total_PbPb2018.root" << std::endl;
        return;
    }

    // Rainbow-style colors for the 10 quantiles to match Figure 5
    int colors[] = {kBlue+2, kBlue, kCyan, kGreen+1, kYellow-3, kOrange-3, kMagenta, kViolet, kGray+2, kBlack};

    // 2. Create the overlapping canvas
    TCanvas *c1 = new TCanvas("c_overlap", Form("Overlapping q2 %s", label.Data()), 900, 700);
    c1->SetLeftMargin(0.12);
    
    TLegend *leg = new TLegend(0.58, 0.45, 0.88, 0.88);
    leg->SetHeader(Form("q_{2} in Centrality %s", label.Data()));
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);

    // 3. Build and Draw the INCLUSIVE background (The "Envelope")
    // Binned with 700 bins between 0 to 0.35
    TH1D *h_inclusive = new TH1D("h_inclusive", ";q_{2};P(q_{2})", 700, 0, 0.35);
    for (int i_q = 0; i_q < 10; i_q++) {
        for (int i_c = start_c; i_c < end_c; i_c++) {
            // Summing 1% centrality bins from the TFile keys
            TH1D *h_slice = (TH1D*)f->Get(Form("hist_q2_HFtotal_cen%d_%d_q2_%d", i_c, i_c+1, i_q));
            if (h_slice) h_inclusive->Add(h_slice);
        }
    }
    
    // Normalize to Unit Area
    if (h_inclusive->Integral() > 0) h_inclusive->Scale(1.0 / h_inclusive->Integral(), "width");
    
    h_inclusive->SetFillColorAlpha(kGray, 0.35); // Shaded background envelope
    h_inclusive->SetLineColor(kGray+1);
    h_inclusive->SetLineWidth(1);
    h_inclusive->GetXaxis()->SetRangeUser(0, 0.16); // Focus on peak region per Figure 5
    h_inclusive->Draw("HIST");
    leg->AddEntry(h_inclusive, "Inclusive", "f");

    // 4. Build and Overlay the 10 QUANTILE curves
    for (int i_q = 0; i_q < 10; i_q++) {
        TString h_sum_name = Form("h_quant_q%d", i_q);
        TH1D *h_quant = new TH1D(h_sum_name, "", 700, 0, 0.35);

        for (int i_c = start_c; i_c < end_c; i_c++) {
            TH1D *h_slice = (TH1D*)f->Get(Form("hist_q2_HFtotal_cen%d_%d_q2_%d", i_c, i_c+1, i_q));
            if (h_slice) h_quant->Add(h_slice);
        }

        if (h_quant->Integral() > 0) h_quant->Scale(1.0 / h_quant->Integral(), "width");

        h_quant->SetLineColor(colors[i_q]);
        h_quant->SetLineWidth(2);
        h_quant->Draw("HIST SAME");
        leg->AddEntry(h_quant, Form("q_{2} bin %d", i_q), "l");
    }

    leg->Draw();
    std::cout << ">>> Successfully generated overlapping plot for centrality " << label << std::endl;
}
