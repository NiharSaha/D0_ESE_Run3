// plot_vn_vs_qbin_D0Chg_combined.C
//
// Option B: plots D0 v3 AND charged-particle v3 as a function of q-bin index
// on the same canvas with two y-axes (left = D0 v3, right = charged v3).
// A Pearson-r correlation coefficient is also printed on each plot.
//
// Produces one PDF per (centrality x pT) bin, saved to vn_vs_qbin_plots/v3/
// Also tiles them into a summary PDF.
//
// Inputs — same files as make_ESE_scatter.C
//   d0_file  : directory vn_vs_qbin  → TGraphErrors  v3_vs_q3bin_<cent>_<pt>
//   chg_file : directory vsQbin_TProfile → TProfile   hp_v3_vsq3_pt0p5to3_<cent>
//
// Usage:
//   root -l -q plot_vn_vs_qbin_D0Chg_combined.C

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <numeric>

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TPad.h"

// ---------------------------------------------------------------------------
// Pearson correlation coefficient + naive error (Fisher z, n-2 dof)
// ---------------------------------------------------------------------------
static double PearsonR(const std::vector<double>& x, const std::vector<double>& y)
{
    int n = (int)x.size();
    if (n < 2) return 0.0;
    double mx = 0, my = 0;
    for (int i = 0; i < n; ++i) { mx += x[i]; my += y[i]; }
    mx /= n; my /= n;
    double num = 0, dx2 = 0, dy2 = 0;
    for (int i = 0; i < n; ++i) {
        double dx = x[i] - mx, dy = y[i] - my;
        num += dx * dy; dx2 += dx * dx; dy2 += dy * dy;
    }
    double denom = std::sqrt(dx2 * dy2);
    return (denom > 0) ? num / denom : 0.0;
}

void plot_vn_vs_qbin_D0Chg_combined(
    const char* d0_file  = "ROOT/Flow_MB0to31_May19_out_combined_v3.root",
    const char* chg_file = "ROOT/flow_Analysis_chg_out_combined_May11.root",
    const char* tag      = "May21_v0",
    bool        isRun3   = true)   // true -> pt0p5to3 (Run3); false -> pt1to3 (Run2)
{
    // =========================================================================
    // Constants — D0 analysis
    // =========================================================================
    const int N_QBINS = 10;

    // v2: 9 pT bins
    const int N_PTBINS_V2 = 9;
    const double pt_edges_v2[N_PTBINS_V2 + 1] = { 2, 3, 4, 5, 6, 8, 10, 15, 30, 100 };

    // v3: 7 pT bins
    const int N_PTBINS_V3 = 7;
    const double pt_edges_v3[N_PTBINS_V3 + 1] = { 2, 4, 6, 8, 10, 20, 50, 100 };

    // =========================================================================
    // Constants — charged-particle analysis
    // =========================================================================
    const int N_CENTBINS_CHG = 5;
    const int    min_cent_chg[N_CENTBINS_CHG] = {  0, 10, 20, 30, 40 };
    const int    max_cent_chg[N_CENTBINS_CHG] = { 10, 20, 30, 40, 50 };
    const char* cent_lbl_chg[N_CENTBINS_CHG] = {
        "cent0to10","cent10to20","cent20to30","cent30to40","cent40to50"
    };

    // =========================================================================
    // Open files
    // =========================================================================
    TFile* fD0  = TFile::Open(d0_file);
    TFile* fChg = TFile::Open(chg_file);
    if (!fD0  || fD0->IsZombie())  { std::cerr << "ERROR: Cannot open D0 file:  " << d0_file  << "\n"; return; }
    if (!fChg || fChg->IsZombie()) { std::cerr << "ERROR: Cannot open chg file: " << chg_file << "\n"; return; }

    auto* dir_d0_qbin = (TDirectory*)fD0->Get("vn_vs_qbin");
    auto* dir_chg_vsq = (TDirectory*)fChg->Get("vsQbin_TProfile");

    if (!dir_d0_qbin) { std::cerr << "WARNING: D0 vn_vs_qbin dir not found\n"; }
    if (!dir_chg_vsq) { std::cerr << "WARNING: CHG vsQbin_TProfile dir not found\n"; }

    // =========================================================================
    // Output
    // =========================================================================
    gSystem->mkdir("vn_vs_qbin_plots",    kTRUE);
    gSystem->mkdir("vn_vs_qbin_plots/v2", kTRUE);
    gSystem->mkdir("vn_vs_qbin_plots/v3", kTRUE);

    TString out_file = Form("ESE_vn_vs_qbin_%s.root", tag);
    TFile* fOut   = new TFile(out_file.Data(), "RECREATE");
    auto*  dir_v2 = fOut->mkdir("vn_vs_qbin_v2");
    auto*  dir_v3 = fOut->mkdir("vn_vs_qbin_v3");

    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    TLatex ltx;
    ltx.SetNDC();
    ltx.SetTextFont(42);

    // =========================================================================
    // Pre-load charged-particle v2 and v3 vs q-bin TProfiles
    // =========================================================================
    // Run3: 0.5 < pT < 3 GeV/c   Run2: 1 < pT < 3 GeV/c
    const char* chg_pt_tag   = isRun3 ? "pt0p5to3" : "pt1to3";
    const char* chg_pt_range = isRun3 ? "0.5-3"    : "1-3";     // for axis labels

    TProfile* hp_chg_v2[N_CENTBINS_CHG] = {};
    TProfile* hp_chg_v3[N_CENTBINS_CHG] = {};
    for (int ic = 0; ic < N_CENTBINS_CHG; ++ic) {
        hp_chg_v2[ic] = dir_chg_vsq
            ? (TProfile*)dir_chg_vsq->Get(Form("hp_v2_vsq2_%s_%s", chg_pt_tag, cent_lbl_chg[ic]))
            : nullptr;
        hp_chg_v3[ic] = dir_chg_vsq
            ? (TProfile*)dir_chg_vsq->Get(Form("hp_v3_vsq3_%s_%s", chg_pt_tag, cent_lbl_chg[ic]))
            : nullptr;
        if (!hp_chg_v2[ic])
            std::cerr << "WARNING: hp_v2_vsq2_" << chg_pt_tag << "_" << cent_lbl_chg[ic] << " not found\n";
        if (!hp_chg_v3[ic])
            std::cerr << "WARNING: hp_v3_vsq3_" << chg_pt_tag << "_" << cent_lbl_chg[ic] << " not found\n";
    }

    // =========================================================================
    // PRE-SCAN: Compute per-(harmonic, centrality) axis ranges by scanning all
    // pT bins (excluding pT >= 30 which is unreliable).  These are shared
    // across all pT bins within the same (ivn, ic).
    // =========================================================================
    double cent_d0_ylo[2][N_CENTBINS_CHG],  cent_d0_yhi[2][N_CENTBINS_CHG];
    double cent_chg_ylo[2][N_CENTBINS_CHG], cent_chg_yhi[2][N_CENTBINS_CHG];

    for (int ivn = 2; ivn <= 3; ++ivn) {
        int vi             = ivn - 2;
        int n_pt_scan      = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
        const double* pt_e = (ivn == 2) ? pt_edges_v2 : pt_edges_v3;
        TProfile** hp_scan = (ivn == 2) ? hp_chg_v2   : hp_chg_v3;
        for (int ic = 0; ic < N_CENTBINS_CHG; ++ic) {
            double gd0_min =  1e9, gd0_max = -1e9;
            double gch_min =  1e9, gch_max = -1e9;
            // Charged vn range (same for all pT, pT-integrated)
            for (int iq = 0; iq < N_QBINS; ++iq)
                if (hp_scan[ic] && hp_scan[ic]->GetBinEntries(iq+1) > 0) {
                    double cv = hp_scan[ic]->GetBinContent(iq+1);
                    double ce = hp_scan[ic]->GetBinError(iq+1);
                    if (ce <= 0) continue;
                    gch_min = std::min(gch_min, cv - ce);
                    gch_max = std::max(gch_max, cv + ce);
                }
            // D0 vn range: scan pT bins (skip pT >= 30)
            for (int ip = 0; ip < n_pt_scan; ++ip) {
                if (pt_e[ip] >= 30) continue;
                TString ptTag  = Form("pT%dto%d", (int)pt_e[ip], (int)pt_e[ip+1]);
                TString grName = Form("v%d_vs_q%dbin_%s_%s", ivn, ivn, cent_lbl_chg[ic], ptTag.Data());
                TGraphErrors* gr = dir_d0_qbin
                    ? (TGraphErrors*)dir_d0_qbin->Get(grName) : nullptr;
                if (!gr || gr->GetN() == 0) continue;
                for (int iq = 0; iq < N_QBINS && iq < gr->GetN(); ++iq) {
                    double dv = gr->GetPointY(iq), de = gr->GetErrorY(iq);
                    double ce = (hp_scan[ic] && hp_scan[ic]->GetBinEntries(iq+1) > 0)
                                ? hp_scan[ic]->GetBinError(iq+1) : 0.0;
                    if (de <= 0 || ce <= 0) continue;
                    gd0_min = std::min(gd0_min, dv - de);
                    gd0_max = std::max(gd0_max, dv + de);
                }
            }
            if (gd0_min > gd0_max) { gd0_min = 0.0; gd0_max = 0.3; }  // fallback
            if (gch_min > gch_max) { gch_min = 0.0; gch_max = 0.3; }
            double d0_marg  = 0.15 * std::max(gd0_max - gd0_min, 1e-6);
            double chg_marg = 0.15 * std::max(gch_max - gch_min, 1e-6);
            cent_d0_ylo[vi][ic]  = gd0_min  - d0_marg;
            cent_d0_yhi[vi][ic]  = gd0_max  + d0_marg;
            cent_chg_ylo[vi][ic] = std::max(0.0, gch_min - chg_marg);
            cent_chg_yhi[vi][ic] = gch_max  + chg_marg;
            // For v2: unify D0 and charged into a single shared range per centrality
            if (ivn == 2) {
                double unified_lo = std::min(cent_d0_ylo[vi][ic], cent_chg_ylo[vi][ic]);
                double unified_hi = std::max(cent_d0_yhi[vi][ic], cent_chg_yhi[vi][ic]);
                cent_d0_ylo[vi][ic]  = unified_lo;
                cent_d0_yhi[vi][ic]  = unified_hi;
                cent_chg_ylo[vi][ic] = unified_lo;
                cent_chg_yhi[vi][ic] = unified_hi;
            }
            std::cout << Form("Range v%d %s: D0=[%.4f,%.4f] chg=[%.4f,%.4f]\n",
                ivn, cent_lbl_chg[ic],
                cent_d0_ylo[vi][ic], cent_d0_yhi[vi][ic],
                cent_chg_ylo[vi][ic], cent_chg_yhi[vi][ic]);
        }
    }

    // =========================================================================
    // Loop over harmonics (v2 and v3), centrality, and pT bins
    // =========================================================================
    for (int ivn = 2; ivn <= 3; ++ivn)
    {
        int           n_pt     = (ivn == 2) ? N_PTBINS_V2  : N_PTBINS_V3;
        const double* pt_edges = (ivn == 2) ? pt_edges_v2  : pt_edges_v3;
        TProfile**    hp_chg   = (ivn == 2) ? hp_chg_v2    : hp_chg_v3;
        TDirectory*   dir_out  = (ivn == 2) ? dir_v2       : dir_v3;

        for (int ic = 0; ic < N_CENTBINS_CHG; ++ic)
        {
            // Charged-particle vn per q-bin
            double chg_vn_q[N_QBINS] = {}, chg_vn_q_e[N_QBINS] = {};
            for (int iq = 0; iq < N_QBINS; ++iq) {
                if (hp_chg[ic] && hp_chg[ic]->GetBinEntries(iq + 1) > 0) {
                    chg_vn_q[iq]   = hp_chg[ic]->GetBinContent(iq + 1);
                    chg_vn_q_e[iq] = hp_chg[ic]->GetBinError(iq + 1);
                }
            }

            for (int ip = 0; ip < n_pt; ++ip)
            {
                const double pt_lo = pt_edges[ip];
                const double pt_hi = pt_edges[ip + 1];
                TString ptTag  = Form("pT%dto%d", (int)pt_lo, (int)pt_hi);
                TString grName = Form("v%d_vs_q%dbin_%s_%s", ivn, ivn, cent_lbl_chg[ic], ptTag.Data());

                TGraphErrors* gr_d0 = dir_d0_qbin
                    ? (TGraphErrors*)dir_d0_qbin->Get(grName)
                    : nullptr;
                if (!gr_d0 || gr_d0->GetN() == 0) {
                    std::cerr << "SKIP: " << grName << " not found or empty\n";
                    continue;
                }

                // Collect valid q-bin points
                std::vector<double> d0_v, d0_e, chg_v, chg_e, qx_valid;
                double d0_q[N_QBINS]   = {}, d0_q_e[N_QBINS]   = {};
                double chg_q[N_QBINS]  = {}, chg_q_e[N_QBINS]  = {};
                double qx[N_QBINS]     = {}, qex[N_QBINS]       = {};
                int npts = 0;

                for (int iq = 0; iq < N_QBINS && iq < gr_d0->GetN(); ++iq) {
                    double dv  = gr_d0->GetPointY(iq);
                    double de  = gr_d0->GetErrorY(iq);
                    double cv  = chg_vn_q[iq];
                    double ce  = chg_vn_q_e[iq];
                    if (de <= 0 || ce <= 0) continue;
                    d0_q[npts]    = dv;  d0_q_e[npts]   = de;
                    chg_q[npts]   = cv;  chg_q_e[npts]  = ce;
                    qx[npts]      = (double)iq; qex[npts] = 0.0;
                    d0_v.push_back(dv);  d0_e.push_back(de);
                    chg_v.push_back(cv); chg_e.push_back(ce);
                    qx_valid.push_back((double)iq);
                    ++npts;
                }
                if (npts < 2) continue;

                // Pearson r between D0 vn and charged vn
                double rho = PearsonR(d0_v, chg_v);
                (void)rho; // reserved for future use

                // ================================================================
                // Axis ranges: per-centrality, shared across all pT bins
                // ================================================================
                double d0_ylo  = cent_d0_ylo[ivn - 2][ic];
                double d0_yhi  = cent_d0_yhi[ivn - 2][ic];
                double chg_ylo = cent_chg_ylo[ivn - 2][ic];
                double chg_yhi = cent_chg_yhi[ivn - 2][ic];

                // ================================================================
                // Build TGraphErrors
                // ================================================================
                TString base = Form("vn_vs_qbin_v%d_%s_%s", ivn, cent_lbl_chg[ic], ptTag.Data());

                TGraphErrors* gr_d0_plot = new TGraphErrors(npts, qx, d0_q,  qex, d0_q_e);
                gr_d0_plot->SetName(Form("gr_d0_%s",  base.Data()));

                TGraphErrors* gr_chg_plot = new TGraphErrors(npts, qx, chg_q, qex, chg_q_e);
                gr_chg_plot->SetName(Form("gr_chg_%s", base.Data()));

                // D0: filled circles, blue
                gr_d0_plot->SetMarkerStyle(20);
                gr_d0_plot->SetMarkerSize(1.4);
                gr_d0_plot->SetMarkerColor(kBlue + 1);
                gr_d0_plot->SetLineColor(kBlue + 1);
                gr_d0_plot->SetLineWidth(2);

                // charged: solid squares, red
                gr_chg_plot->SetMarkerStyle(21);
                gr_chg_plot->SetMarkerSize(1.4);
                gr_chg_plot->SetMarkerColor(kRed + 1);
                gr_chg_plot->SetLineColor(kRed + 1);
                gr_chg_plot->SetLineWidth(2);

                TString cvName = Form("cv_%s", base.Data());
                TCanvas* cv = new TCanvas(cvName, base, 700, 550);
                cv->SetLeftMargin(0.14);
                cv->SetRightMargin(0.14);
                cv->SetBottomMargin(0.14);
                cv->SetTopMargin(0.08);
                cv->SetGridx();
                cv->SetGridy();

                // Frame drawn on left axis (D0 vn)
                TH2F* hf = new TH2F(Form("hf_%s", base.Data()), "",
                                     10, -0.5, (double)N_QBINS - 0.5,
                                     10, d0_ylo, d0_yhi);
                hf->GetXaxis()->SetTitle(Form("q_{%d} bin index", ivn));
                if (ivn == 2) {
                    hf->GetYaxis()->SetTitle("v_{2}");
                } else {
                    hf->GetYaxis()->SetTitle("v_{3}");
                }
                hf->GetXaxis()->SetTitleSize(0.060); // larger
                hf->GetYaxis()->SetTitleSize(0.060); // larger
                hf->GetXaxis()->SetTitleOffset(0.95); // more centered
                hf->GetYaxis()->SetTitleOffset(0.85); // more centered
                hf->GetXaxis()->SetLabelSize(0.038);
                hf->GetYaxis()->SetLabelSize(0.038);
                hf->GetXaxis()->SetNdivisions(10);
                hf->GetYaxis()->SetTitleColor(kBlack); // black
                hf->GetYaxis()->SetLabelColor(kBlack); // black
                hf->Draw();


                gr_d0_plot->Draw("P SAME");
                gr_chg_plot->Draw("P SAME");

                // Legend — top right
                TLegend* leg = new TLegend(0.52, 0.78, 0.86, 0.90);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.038);
                leg->AddEntry(gr_d0_plot,    Form("D^{0} v_{%d}", ivn),            "pe");
                leg->AddEntry(gr_chg_plot, Form("charged particle v_{%d}", ivn), "pe");
                leg->Draw();

                // Centrality + pT labels — top left, two lines
                ltx.SetTextAlign(11);
                ltx.SetTextFont(42);
                ltx.SetTextSize(0.038);
                ltx.DrawLatex(0.17, 0.86,
                    Form("Centrality: %d-%d%%", min_cent_chg[ic], max_cent_chg[ic]));
                ltx.DrawLatex(0.17, 0.79,
                    Form("%.0f < p_{T}^{D^{0}} < %.0f GeV/c", pt_lo, pt_hi));

                // CMS label
                ltx.SetTextAlign(11);
                ltx.SetTextFont(62);
                ltx.SetTextSize(0.060); // larger CMS
                ltx.DrawLatex(0.14, 0.945, "CMS");
                ltx.SetTextAlign(31);
                ltx.SetTextFont(42);
                ltx.SetTextSize(0.042);
                ltx.DrawLatex(0.86, 0.945, "PbPb 5.36 TeV");
                ltx.SetTextAlign(11);

                cv->SaveAs(Form("vn_vs_qbin_plots/v%d/%s.pdf", ivn, base.Data()));

                dir_out->cd();
                gr_d0_plot->Write();
                gr_chg_plot->Write();
                cv->Write();

                delete leg;
                delete hf;
                delete gr_d0_plot;
                delete gr_chg_plot;
                delete cv;

                std::cout << "Created: " << base << "\n";

            } // ip loop
        } // ic loop
    } // ivn loop

    fOut->Write();
    fOut->Close();
    fD0->Close();
    fChg->Close();

    // =========================================================================
    // Tile all individual PDFs into a summary multi-page PDF via Python
    // =========================================================================
    std::cout << "\nTiling PDFs...\n";
    gSystem->Exec("/usr/local/bin/python3.11 tile_vn_vs_qbin.py 2");
    gSystem->Exec("/usr/local/bin/python3.11 tile_vn_vs_qbin.py 3");
    gSystem->Exec(Form("cp vn_vs_qbin_summary_v2_May13.pdf vn_vs_qbin_summary_v2_%s.pdf", tag));
    gSystem->Exec(Form("cp vn_vs_qbin_summary_v3_May13.pdf vn_vs_qbin_summary_v3_%s.pdf", tag));

    std::cout << "\nDone!\n"
              << "  Individual PDFs -> vn_vs_qbin_plots/v2/  and  vn_vs_qbin_plots/v3/\n"
              << Form("  Summary PDFs    -> vn_vs_qbin_summary_v2_%s.pdf  and  vn_vs_qbin_summary_v3_%s.pdf\n", tag, tag)
              << "  Graphs          -> " << out_file << "  (directories: vn_vs_qbin_v2, vn_vs_qbin_v3)\n";
}
