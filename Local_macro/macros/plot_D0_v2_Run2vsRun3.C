// plot_D0_v2_Run2vsRun3.C
//
// Compares D0 v2 vs pT for each q2-bin (Run2 vs Run3).
//
// Run2 file : ROOT/v2_vs_pT_for_q2_bins_5_pT_bins.root
//   Objects : v2_graph_cen{min}_{max}_q2bin_{iq}  (TGraphErrors, v2 vs pT, 5 pts)
//             x = pT midpoint, y = D0 v2,  errors from GetErrorY
//   pT bins : 2-4, 4-6, 6-10, 10-15, 15-30  (midpoints: 3,5,8,12.5,22.5)
//   cent tag: cen_{min}_{max}  e.g. cen_0_10, cen_40_50
//
// Run3 file : ROOT/Flow_MB0to31_May19_out_combined_v3.root
//   Directory: vn_vs_qbin
//   Objects  : v2_vs_q2bin_cent{X}to{Y}_pT{lo}to{hi}  (TGraphErrors)
//   v2 for q2-bin iq → read GetPointY(iq), GetErrorY(iq) from each pT graph
//   pT bins used: 2-3,3-4,4-5,5-6,6-8,8-10,10-15,15-30 (30-100 excluded)
//
// Produces one PDF per (centrality, q2-bin) in comparison_plots/
// Then tiles them: 5 cols × 2 rows per page → 5-page summary PDF.
//
// Usage:
//   root -l -b -q plot_D0_v2_Run2vsRun3.C
//   root -l -b -q 'plot_D0_v2_Run2vsRun3.C("ROOT/run2.root","ROOT/run3.root","v3")'

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"

void plot_D0_v2_Run2vsRun3(
    const char* run2_file = "ROOT/v2_vs_pT_for_q2_bins_5_pT_bins.root",
    const char* run3_file = "ROOT/Flow_MB0to31_May19_out_combined_v3.root",
    const char* tag       = "May21_v3")
{
    // =========================================================================
    // Constants
    // =========================================================================
    const int N_QBINS = 10;

    // Centrality
    const int N_CENT = 5;
    const int minCent[N_CENT] = {  0, 10, 20, 30, 40 };
    const int maxCent[N_CENT] = { 10, 20, 30, 40, 50 };

    // Run3 pT bins (exclude 30-100)
    const int NPT_R3 = 8;
    const char* r3_ptTag[NPT_R3] = {
        "pT2to3","pT3to4","pT4to5","pT5to6",
        "pT6to8","pT8to10","pT10to15","pT15to30"
    };
    // Arithmetic midpoints for x-axis
    const double r3_ptCen[NPT_R3] = { 2.5, 3.5, 4.5, 5.5, 7.0, 9.0, 12.5, 22.5 };
    // Half-widths for x error bars
    const double r3_ptHW[NPT_R3]  = { 0.5, 0.5, 0.5, 0.5, 1.0, 1.0,  2.5,  7.5 };

    // Run2 cent tag format: cen_{min}_{max}  (e.g. cen_0_10, cen_40_50)
    // Run3 cent tag format: cent{min}to{max} (e.g. cent0to10, cent40to50)

    // =========================================================================
    // Style
    // =========================================================================
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // =========================================================================
    // Open files
    // =========================================================================
    TFile* fR2 = TFile::Open(run2_file);
    TFile* fR3 = TFile::Open(run3_file);
    if (!fR2 || fR2->IsZombie()) { std::cerr << "ERROR: Cannot open Run2 file: " << run2_file << "\n"; return; }
    if (!fR3 || fR3->IsZombie()) { std::cerr << "ERROR: Cannot open Run3 file: " << run3_file << "\n"; return; }

    auto* dir_r3 = (TDirectory*)fR3->Get("vn_vs_qbin");
    if (!dir_r3) { std::cerr << "ERROR: vn_vs_qbin directory not found in Run3 file\n"; return; }

    // =========================================================================
    // Output directory
    // =========================================================================
    TString plotDir = "comparison_plots";
    gSystem->mkdir(plotDir, kTRUE);

    TLatex ltx;
    ltx.SetNDC();
    ltx.SetTextFont(42);

    // =========================================================================
    // Pre-load Run3 TGraphErrors for all (cent, pT) combinations
    // gr_r3[ic][ip] = TGraphErrors with N_QBINS points (v2 vs q2bin index)
    // =========================================================================
    TGraphErrors* gr_r3[N_CENT][NPT_R3] = {};
    for (int ic = 0; ic < N_CENT; ++ic) {
        TString cent_r3 = Form("cent%dto%d", minCent[ic], maxCent[ic]);
        for (int ip = 0; ip < NPT_R3; ++ip) {
            TString gname = Form("v2_vs_q2bin_%s_%s", cent_r3.Data(), r3_ptTag[ip]);
            gr_r3[ic][ip] = (TGraphErrors*)dir_r3->Get(gname);
            if (!gr_r3[ic][ip])
                std::cerr << "  WARN: " << gname << " not found\n";
        }
    }

    // =========================================================================
    // Loop: one individual plot per (centrality × q2-bin)
    // =========================================================================
    for (int ic = 0; ic < N_CENT; ++ic)
    {
        TString cent_r2 = Form("cen%d_%d", minCent[ic], maxCent[ic]);

        for (int iq = 0; iq < N_QBINS; ++iq)
        {
            // ── Load Run2 TGraphErrors (v2 vs pT, 5 pts) ─────────────────────
            TString r2name = Form("v2_graph_%s_q2bin_%d", cent_r2.Data(), iq);
            TGraphErrors* gr_r2raw = (TGraphErrors*)fR2->Get(r2name);
            if (!gr_r2raw) {
                std::cerr << "  WARN: " << r2name << " not found\n";
            }

            // ── Build Run3 v2-vs-pT graph for this q2-bin ────────────────────
            // For each pT bin, read the point at index iq
            std::vector<double> r3x, r3ex, r3y, r3ey;
            for (int ip = 0; ip < NPT_R3; ++ip) {
                if (!gr_r3[ic][ip] || iq >= gr_r3[ic][ip]->GetN()) continue;
                double v  = gr_r3[ic][ip]->GetPointY(iq);
                double e  = gr_r3[ic][ip]->GetErrorY(iq);
                if (e <= 0) continue;
                r3x.push_back(r3_ptCen[ip]);
                r3ex.push_back(r3_ptHW[ip]);
                r3y.push_back(v);
                r3ey.push_back(e);
            }

            // ── y-axis range ──────────────────────────────────────────────────
            double ylo =  1e9, yhi = -1e9;
            auto upd = [&](double v, double e) {
                if (e <= 0) return;
                ylo = std::min(ylo, v - e);
                yhi = std::max(yhi, v + e);
            };
            if (gr_r2raw)
                for (int i = 0; i < gr_r2raw->GetN(); ++i)
                    upd(gr_r2raw->GetPointY(i), gr_r2raw->GetErrorY(i));
            for (int i = 0; i < (int)r3y.size(); ++i) upd(r3y[i], r3ey[i]);
            if (ylo > yhi) { ylo = -0.05; yhi = 0.30; }
            double range = std::max(yhi - ylo, 1e-6);
            ylo -= 0.20 * range;
            yhi += 0.35 * range;

            // ── Canvas ────────────────────────────────────────────────────────
            TString baseName = Form("D0v2_comp_cent%dto%d_q2bin%d",
                                    minCent[ic], maxCent[ic], iq);
            TCanvas* cv = new TCanvas(baseName, baseName, 600, 500);
            cv->SetLeftMargin(0.14);
            cv->SetRightMargin(0.05);
            cv->SetBottomMargin(0.14);
            cv->SetTopMargin(0.09);
            cv->SetGridx();
            cv->SetGridy();

            // ── Frame ─────────────────────────────────────────────────────────
            TH2F* hf = new TH2F(Form("hf_%s", baseName.Data()), "",
                                 100, 1.5, 30.0, 100, ylo, yhi);
            hf->GetXaxis()->SetTitle("p_{T}^{D^{0}} (GeV/c)");
            hf->GetYaxis()->SetTitle("D^{0} v_{2}");
            hf->GetXaxis()->SetTitleSize(0.050);
            hf->GetYaxis()->SetTitleSize(0.050);
            hf->GetXaxis()->SetLabelSize(0.040);
            hf->GetYaxis()->SetLabelSize(0.040);
            hf->GetXaxis()->SetTitleOffset(1.10);
            hf->GetYaxis()->SetTitleOffset(1.30);
            hf->GetXaxis()->SetNdivisions(506);
            hf->GetYaxis()->SetNdivisions(506);
            hf->Draw();

            // y = 0 line
            TLine* zero = new TLine(1.5, 0.0, 30.0, 0.0);
            zero->SetLineColor(kGray + 1);
            zero->SetLineStyle(2);
            zero->SetLineWidth(1);
            zero->Draw();

            // ── Run2 graph ────────────────────────────────────────────────────
            TGraphErrors* gr_r2 = nullptr;
            if (gr_r2raw && gr_r2raw->GetN() > 0) {
                gr_r2 = (TGraphErrors*)gr_r2raw->Clone(Form("gr_r2_%s", baseName.Data()));
                gr_r2->SetMarkerStyle(20);
                gr_r2->SetMarkerSize(1.4);
                gr_r2->SetMarkerColor(kBlue + 1);
                gr_r2->SetLineColor(kBlue + 1);
                gr_r2->SetLineWidth(2);
                gr_r2->Draw("P SAME");
            }

            // ── Run3 graph ────────────────────────────────────────────────────
            TGraphErrors* gr_r3plot = nullptr;
            if (!r3y.empty()) {
                gr_r3plot = new TGraphErrors((int)r3y.size(),
                    r3x.data(), r3y.data(), r3ex.data(), r3ey.data());
                gr_r3plot->SetName(Form("gr_r3_%s", baseName.Data()));
                gr_r3plot->SetMarkerStyle(21);
                gr_r3plot->SetMarkerSize(1.4);
                gr_r3plot->SetMarkerColor(kRed + 1);
                gr_r3plot->SetLineColor(kRed + 1);
                gr_r3plot->SetLineWidth(2);
                gr_r3plot->Draw("P SAME");
            }

            // ── Legend (Run2/Run3) – top-right ────────────────────────────────
            TLegend* leg = new TLegend(0.56, 0.78, 0.94, 0.90);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.042);
            if (gr_r2)     leg->AddEntry(gr_r2,     "Run2  5.02 TeV", "pe");
            if (gr_r3plot) leg->AddEntry(gr_r3plot, "Run3  5.36 TeV", "pe");
            leg->Draw();

            // ── Labels ────────────────────────────────────────────────────────
            ltx.SetTextAlign(11);
            ltx.SetTextFont(62);
            ltx.SetTextSize(0.055);
            ltx.DrawLatex(0.14, 0.935, "CMS");
            ltx.SetTextFont(42);
            ltx.SetTextSize(0.040);
            ltx.SetTextAlign(31);
            ltx.DrawLatex(0.95, 0.935, "PbPb");

            // Cent and q-bin – top-left
            ltx.SetTextAlign(11);
            ltx.SetTextSize(0.040);
            ltx.DrawLatex(0.16, 0.86,
                Form("%d < Cent. < %d%%", minCent[ic], maxCent[ic]));
            ltx.DrawLatex(0.16, 0.79,
                Form("q_{2} bin = %d", iq));

            TString outPath = Form("%s/%s.pdf", plotDir.Data(), baseName.Data());
            cv->SaveAs(outPath);
            delete leg; delete hf; delete zero;
            if (gr_r2) delete gr_r2;
            if (gr_r3plot) delete gr_r3plot;
            delete cv;

        } // iq
        std::cout << "  Centrality " << minCent[ic] << "-" << maxCent[ic]
                  << "%: " << N_QBINS << " plots saved\n";
    } // ic

    fR2->Close();
    fR3->Close();

    // =========================================================================
    // Tile individual PDFs into 5-page summary
    // =========================================================================
    std::cout << "\nTiling PDFs...\n";
    gSystem->Exec(Form(
        "/usr/local/bin/python3.11 tile_D0v2_comparison.py %s %s",
        plotDir.Data(), tag));

    std::cout << "\nDone!\n"
              << "  Individual plots: " << plotDir << "/\n"
              << Form("  Summary PDF    : D0_v2_Run2vsRun3_%s.pdf  (5 pages)\n", tag);
}

#ifndef __CINT__
#ifndef __CLING__
int main() { plot_D0_v2_Run2vsRun3(); return 0; }
#endif
#endif
