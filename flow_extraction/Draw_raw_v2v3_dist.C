#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include <iostream>

#include "Analysis_bin.h"

void Draw_raw_v2v3_dist()
{
    const char* infile = "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/save_vnbinning_outputs_Mar9_v0/build/vnbinning_out_combined.root";

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TFile* f = TFile::Open(infile);
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open file: " << infile << std::endl;
        return;
    }

    const int NCOLS = N_CENTBINS;   // 6
    const int NROWS = N_PTBINS;     // 12
    const int padW  = 280;
    const int padH  = 220;

    TCanvas* c = new TCanvas("c_v2v3", "v2v3 distributions",
                              NCOLS * padW, NROWS * padH);
    c->Divide(NCOLS, NROWS, 0.0005, 0.0005);

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);

    for (int ip = 0; ip < N_PTBINS; ++ip) {
        for (int ic = 0; ic < N_CENTBINS; ++ic) {

            int pad_idx = ip * NCOLS + ic + 1;
            TVirtualPad* pad = c->cd(pad_idx);
            pad->SetLeftMargin(0.18);
            pad->SetRightMargin(0.03);
            pad->SetTopMargin(0.10);
            pad->SetBottomMargin(0.18);

            // ---- retrieve histograms ----
            TString name_v2 = Form("h_v2_dist_%s_%s", cen_name[ic], pt_name[ip]);
            TString name_v3 = Form("h_v3_dist_%s_%s", cen_name[ic], pt_name[ip]);

            TH1D* h_v2 = (TH1D*)f->Get(name_v2);
            TH1D* h_v3 = (TH1D*)f->Get(name_v3);

            if (!h_v2) { std::cerr << "[WARN] Missing: " << name_v2 << "\n"; continue; }
            if (!h_v3) { std::cerr << "[WARN] Missing: " << name_v3 << "\n"; continue; }

            // Clone so we don't modify the originals
            h_v2 = (TH1D*)h_v2->Clone(Form("cv2_%d_%d", ic, ip));
            h_v3 = (TH1D*)h_v3->Clone(Form("cv3_%d_%d", ic, ip));
            h_v2->SetDirectory(nullptr);
            h_v3->SetDirectory(nullptr);

            // ---- normalise to unit area ----
            if (h_v2->Integral() > 0) h_v2->Scale(1.0 / h_v2->Integral());
            if (h_v3->Integral() > 0) h_v3->Scale(1.0 / h_v3->Integral());

            // ---- style: thin lines, NO fill ----
            h_v2->SetLineColor(kRed+1);
            h_v2->SetLineWidth(1.0);
            h_v2->SetFillStyle(0);
            h_v2->SetMarkerSize(0);

            h_v3->SetLineColor(kAzure+2);
            h_v3->SetLineWidth(1.0);
            h_v3->SetFillStyle(0);
            h_v3->SetMarkerSize(0);

            // ---- axis range ----
            double ymax = std::max(h_v2->GetMaximum(), h_v3->GetMaximum());
            h_v2->SetMaximum(1.55 * ymax);
            h_v2->SetMinimum(0.);

            // ---- axis labels ----
            h_v2->GetXaxis()->SetTitle("v_{n}^{SP}");
            h_v2->GetXaxis()->SetTitleSize(0.08);
            h_v2->GetXaxis()->SetLabelSize(0.07);
            h_v2->GetXaxis()->SetTitleOffset(1.0);
            h_v2->GetYaxis()->SetTitle("Norm. entries");
            h_v2->GetYaxis()->SetTitleSize(0.07);
            h_v2->GetYaxis()->SetLabelSize(0.065);
            h_v2->GetYaxis()->SetTitleOffset(1.2);
            h_v2->SetTitle("");

            h_v2->Draw("HIST");
            h_v3->Draw("HIST SAME");

            // ---- legend: shifted more to the right ----
            TLegend* leg = new TLegend(0.70, 0.70, 0.97, 0.89);  // was (0.58, 0.70, 0.96, 0.89)
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.08);
            leg->AddEntry(h_v2, "v_{2}^{SP}", "l");
            leg->AddEntry(h_v3, "v_{3}^{SP}", "l");
            leg->Draw();

            // ---- centrality + pT labels: shifted lower ----
            tex.SetTextAlign(13);
            tex.SetTextSize(0.078);
            tex.DrawLatex(0.24, 0.83,                              
                Form("Cent: %d-%d%%", cen_edges[ic], cen_edges[ic+1]));
            tex.SetTextSize(0.078);
            tex.DrawLatex(0.24, 0.73,                              
                Form("%.0f< p_{T}< %.0f GeV/c", pt_edges[ip], pt_edges[ip+1]));

        } // cent
    } // pT

    c->SaveAs("v2v3_distributions_allbins_new.pdf");
    std::cout << "Saved: v2v3_distributions_allbins.pdf" << std::endl;

    f->Close();
    delete c;
}

