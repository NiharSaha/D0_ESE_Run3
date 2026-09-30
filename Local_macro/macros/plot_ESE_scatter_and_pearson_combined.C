//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// plot_ESE_scatter_and_pearson_combined.C
//
// Combined macro that:
//   1) Builds scatter plots of D0 vn vs charged-particle vn per q-bin,
//      for every (centrality × D0-pT) bin  [make_ESE_scatter.C]
//   2) Computes the weighted Pearson correlation coefficient r between
//      D0 vn and charged-particle vn per q-bin, as a function of D0 pT,
//      for each centrality class  [plot_pearson_ESE.C]
//
// The scatter TGraphErrors built in Step 1 are passed directly in memory
// to Step 2 — no intermediate file read/write is needed between steps.
//
// Formula (Pearson):
//          sum_i  w_i (x_i - <x_w>) (y_i - <y_w>)
//   r = -----------------------------------------------
//       sqrt[ sum_i w_i(x_i-<x_w>)^2 * sum_i w_i(y_i-<y_w>)^2 ]
//
//   w_i = 1/sigma_i^2,  sigma_i = sqrt((Dx_i)^2 + (Dy_i)^2)
//
// Outputs — all written inside a single tagged output directory:
//   {outdir}/scatter_plots/v{2,3}/*.pdf
//   {outdir}/scatter_plots/v{2,3}/merged_scatter_v{2,3}.pdf
//   {outdir}/pearson_plots/v{2,3}/*.pdf
//   {outdir}/slope_intercept_plots/slope_intercept_v{2,3}_vs_pT.pdf
//   {outdir}/mc_r_plots/v{2,3}/*.pdf
//   {outdir}/{scatter_file}  (ROOT: scatter_v2, scatter_v3)
//   {outdir}/{pearson_file}  (ROOT: pearson_v2, pearson_v3)
//
// outdir is auto-generated as  ESE_YYYYMMDD_{version}  when not provided.
//
// Usage:
//   root -l -b -q 'plot_ESE_scatter_and_pearson_combined.C()'   // -> ESE_20260519_v1/
//   root -l -b -q 'plot_ESE_scatter_and_pearson_combined.C("ROOT/Flow_MB0to31_May19_out_combined_v2.root",
//                                             "ROOT/flow_Analysis_chg_out_combined_May11.root",
//                                             "ESE_scatter_May19.root",
//                                             "pearson_ESE_May19.root",
//                                             "", "v2")'  // -> ESE_20260519_v2/
//   Pass a non-empty 5th arg to use a fully custom directory name.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <random>
#include <string>
#include <cstdio>

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVirtualFitter.h"
#include "TROOT.h"
#include "TMath.h"

// =========================================================================
// Helper: weighted Pearson r
//   w_i = 1 / (ex_i^2 + ey_i^2)
//   Returns NaN if the result is ill-defined.
// =========================================================================
static double WeightedPearson(const std::vector<double>& x,
                              const std::vector<double>& y,
                              const std::vector<double>& ex,
                              const std::vector<double>& ey)
{
    int n = static_cast<int>(x.size());
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();

    std::vector<double> w(n);
    double sum_w = 0.0;
    for (int i = 0; i < n; ++i) {
        double sig2 = ex[i]*ex[i] + ey[i]*ey[i];
        w[i]   = (sig2 > 0.0) ? 1.0/sig2 : 0.0;
        sum_w += w[i];
    }
    if (sum_w <= 0.0) return std::numeric_limits<double>::quiet_NaN();

    double xw = 0.0, yw = 0.0;
    for (int i = 0; i < n; ++i) { xw += w[i]*x[i]; yw += w[i]*y[i]; }
    xw /= sum_w;
    yw /= sum_w;

    double cov = 0.0, varx = 0.0, vary = 0.0;
    for (int i = 0; i < n; ++i) {
        double dx = x[i] - xw;
        double dy = y[i] - yw;
        cov  += w[i]*dx*dy;
        varx += w[i]*dx*dx;
        vary += w[i]*dy*dy;
    }
    double denom = std::sqrt(varx * vary);
    if (denom <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    return cov / denom;
}

// =========================================================================
// Main function
// =========================================================================
void plot_ESE_scatter_and_pearson_combined(
    const char* d0_file      = "ROOT/Flow_MB0to31_May19_out_combined_v3.root",
    const char* chg_file     = "ROOT/flow_Analysis_chg_out_combined_May11.root",
    const char* scatter_file = "ESE_scatter_May20_v3.root",
    const char* pearson_file = "pearson_ESE_May20_v3.root",
    const char* outdir_tag   = "ESE_May20_v2",
    bool        isRun3       = false)   // true -> pt0p5to3 (Run3); false -> pt1to3 (Run2)
{
    // =========================================================================
    // Settings
    // =========================================================================

    // Normalize D0 vn and charged vn by their q-bin mean before scatter
    const bool   NORMALIZE_BY_MEAN     = false;

    // Axis ranges — used only when NORMALIZE_BY_MEAN = true
    const double X_AXIS_MIN            = 0.0;
    const double X_AXIS_MAX            = 2.0;
    const double Y_AXIS_MIN            = -1.0;
    const double Y_AXIS_MAX            = 4.0;

    // MC resampling
    const int    N_MC                  = 50000;

    // Uncertainty interval method:
    //   false — equal-tails percentile [p15.85, p84.15]  (asymmetric)
    //   true  — symmetric expansion around r_meas (Run-2 style)
    const bool   USE_SYMMETRIC_INTERVAL = true;

    // =========================================================================
    // Bin definitions
    // =========================================================================
    const int N_CENTBINS = 5;
    const int N_QBINS    = 10;

    // Centrality labels — first 5 of the D0 bins (overlap with charged-particle file)
    const char* cen_name[N_CENTBINS] = {
        "cent0to10","cent10to20","cent20to30","cent30to40","cent40to50"
    };
    const int min_cent[N_CENTBINS] = { 0, 10, 20, 30, 40};
    const int max_cent[N_CENTBINS] = {10, 20, 30, 40, 50};

    // v2: 9 pT bins
    const int N_PTBINS_V2 = 9;
    const double pt_edges_v2[N_PTBINS_V2 + 1] = {
        2, 3, 4, 5, 6, 8, 10, 15, 30, 100
    };

    // v3: 7 pT bins
    const int N_PTBINS_V3 = 7;
    const double pt_edges_v3[N_PTBINS_V3 + 1] = {
        2, 4, 6, 8, 10, 20, 50, 100
    };

    // Maximum pT bins across both harmonics
    const int N_PTBINS_MAX = 9;   // max(N_PTBINS_V2, N_PTBINS_V3)

    // pT bin centres and half-widths (geometric mean, for log-scale axis)
    double pt_cen_v2[N_PTBINS_V2], pt_elo_v2[N_PTBINS_V2], pt_ehi_v2[N_PTBINS_V2];
    for (int ip = 0; ip < N_PTBINS_V2; ++ip) {
        pt_cen_v2[ip] = std::sqrt(pt_edges_v2[ip] * pt_edges_v2[ip+1]);
        pt_elo_v2[ip] = pt_cen_v2[ip] - pt_edges_v2[ip];
        pt_ehi_v2[ip] = pt_edges_v2[ip+1] - pt_cen_v2[ip];
    }
    double pt_cen_v3[N_PTBINS_V3], pt_elo_v3[N_PTBINS_V3], pt_ehi_v3[N_PTBINS_V3];
    for (int ip = 0; ip < N_PTBINS_V3; ++ip) {
        pt_cen_v3[ip] = std::sqrt(pt_edges_v3[ip] * pt_edges_v3[ip+1]);
        pt_elo_v3[ip] = pt_cen_v3[ip] - pt_edges_v3[ip];
        pt_ehi_v3[ip] = pt_edges_v3[ip+1] - pt_cen_v3[ip];
    }

    // Per-centrality style (Pearson plots)
    const int col   [N_CENTBINS] = {kBlack, kRed+1, kGreen+2, kBlue+1, kMagenta+1};
    const int marker[N_CENTBINS] = {20,     21,     22,       23,      33};

    // =========================================================================
    // In-memory scatter graph storage  [ivn-2][ic][ip]
    // Populated in Section 1, consumed in Section 2.
    // =========================================================================
    TGraphErrors* gr_sc_mem[2][N_CENTBINS][N_PTBINS_MAX] = {};

    // =========================================================================
    // Open input files
    // =========================================================================
    TFile* fD0  = TFile::Open(d0_file);
    TFile* fChg = TFile::Open(chg_file);
    if (!fD0  || fD0->IsZombie())  { std::cerr << "ERROR: Cannot open D0 file:  " << d0_file  << "\n"; return; }
    if (!fChg || fChg->IsZombie()) { std::cerr << "ERROR: Cannot open chg file: " << chg_file << "\n"; return; }

    auto* dir_d0_qbin = (TDirectory*)fD0->Get("vn_vs_qbin");
    auto* dir_chg_vsq = (TDirectory*)fChg->Get("vsQbin_TProfile");
    if (!dir_d0_qbin) std::cerr << "WARNING: D0 vn_vs_qbin dir not found\n";
    if (!dir_chg_vsq) std::cerr << "WARNING: CHG vsQbin_TProfile dir not found\n";

    // =========================================================================
    // Output directory — use the tag directly as provided
    // =========================================================================
    TString outbase = outdir_tag;
    std::cout << "\n=== Output directory: " << outbase << " ===\n";

    // Resolved paths for the ROOT output files
    TString scatter_path = Form("%s/%s", outbase.Data(), gSystem->BaseName(scatter_file));
    TString pearson_path = Form("%s/%s", outbase.Data(), gSystem->BaseName(pearson_file));

    // =========================================================================
    // Create output directories
    // =========================================================================
    gSystem->mkdir(outbase,                                             kTRUE);
    gSystem->mkdir(Form("%s/scatter_plots",            outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/scatter_plots/v2",         outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/scatter_plots/v3",         outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/pearson_plots",            outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/pearson_plots/v2",         outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/pearson_plots/v3",         outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/slope_intercept_plots",    outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/mc_r_plots",               outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/mc_r_plots/v2",            outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/mc_r_plots/v3",            outbase.Data()), kTRUE);

    // =========================================================================
    // Create output ROOT files
    // =========================================================================
    TFile* fSc  = new TFile(scatter_path, "RECREATE");
    auto*  sc_v2 = fSc->mkdir("scatter_v2");
    auto*  sc_v3 = fSc->mkdir("scatter_v3");

    TFile* fPr   = new TFile(pearson_path, "RECREATE");
    TDirectory* pr_v2 = fPr->mkdir("pearson_v2");
    TDirectory* pr_v3 = fPr->mkdir("pearson_v3");

    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    std::mt19937_64 rng(2025);

    TLatex ltx;
    ltx.SetNDC();
    ltx.SetTextFont(42);

    // =========================================================================
    // Pre-load charged-particle vn vs q-bin TProfiles
    // =========================================================================
    // Run3: 0.5 < pT < 3 GeV/c   Run2: 1 < pT < 3 GeV/c
    const char* chg_pt_tag   = isRun3 ? "pt0p5to3" : "pt1to3";
    const char* chg_pt_range = isRun3 ? "0.5-3"    : "1-3";     // for axis labels

    TProfile* hp_chg_v2[N_CENTBINS] = {};
    TProfile* hp_chg_v3[N_CENTBINS] = {};
    for (int ic = 0; ic < N_CENTBINS; ++ic) {
        hp_chg_v2[ic] = dir_chg_vsq
            ? (TProfile*)dir_chg_vsq->Get(Form("hp_v2_vsq2_%s_%s", chg_pt_tag, cen_name[ic]))
            : nullptr;
        hp_chg_v3[ic] = dir_chg_vsq
            ? (TProfile*)dir_chg_vsq->Get(Form("hp_v3_vsq3_%s_%s", chg_pt_tag, cen_name[ic]))
            : nullptr;
        if (!hp_chg_v2[ic]) std::cerr << "WARNING: hp_v2_vsq2_" << chg_pt_tag << "_" << cen_name[ic] << " not found\n";
        if (!hp_chg_v3[ic]) std::cerr << "WARNING: hp_v3_vsq3_" << chg_pt_tag << "_" << cen_name[ic] << " not found\n";
    }
    (void)chg_pt_range; // used in axis labels later if needed

    // =========================================================================
    // PRE-SCAN: Compute per-(harmonic, centrality) axis ranges by scanning
    // all pT bins and taking the global min/max (including ±error bars).
    // These are used uniformly for every pT bin within the same (ivn, ic).
    // =========================================================================
    double cent_xlo[2][N_CENTBINS], cent_xhi[2][N_CENTBINS];
    double cent_ylo[2][N_CENTBINS], cent_yhi[2][N_CENTBINS];

    for (int ivn = 2; ivn <= 3; ++ivn) {
        int vi = ivn - 2;
        int n_pt_scan       = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
        const double* pt_e  = (ivn == 2) ? pt_edges_v2 : pt_edges_v3;
        for (int ic = 0; ic < N_CENTBINS; ++ic) {
            double gxmin =  1e9, gxmax = -1e9;
            double gymin =  1e9, gymax = -1e9;
            TProfile* hp_scan = (ivn == 2) ? hp_chg_v2[ic] : hp_chg_v3[ic];
            double vchg[N_QBINS] = {}, echg[N_QBINS] = {};
            for (int iq = 0; iq < N_QBINS; ++iq)
                if (hp_scan && hp_scan->GetBinEntries(iq+1) > 0) {
                    vchg[iq] = hp_scan->GetBinContent(iq+1);
                    echg[iq] = hp_scan->GetBinError(iq+1);
                }
            for (int ip = 0; ip < n_pt_scan; ++ip) {
                if (pt_e[ip] >= 30) continue;  // skip unreliable high-pT bin (30-100) from range scan
                TString ptTag = Form("pT%dto%d", (int)pt_e[ip], (int)pt_e[ip+1]);
                TGraphErrors* gr = dir_d0_qbin
                    ? (TGraphErrors*)dir_d0_qbin->Get(
                          Form("v%d_vs_q%dbin_%s_%s", ivn, ivn, cen_name[ic], ptTag.Data()))
                    : nullptr;
                if (!gr || gr->GetN() == 0) continue;
                for (int iq = 0; iq < N_QBINS && iq < gr->GetN(); ++iq) {
                    double cv = vchg[iq], ce = echg[iq];
                    double dv = gr->GetPointY(iq), de = gr->GetErrorY(iq);
                    if (ce <= 0 || de <= 0) continue;
                    gxmin = std::min(gxmin, cv - ce);  gxmax = std::max(gxmax, cv + ce);
                    gymin = std::min(gymin, dv - de);  gymax = std::max(gymax, dv + de);
                }
            }
            if (gxmin > gxmax) { gxmin = 0.0; gxmax = 0.3; }  // fallback
            if (gymin > gymax) { gymin = 0.0; gymax = 0.3; }
            double xmarg = 0.15 * std::max(gxmax - gxmin, 1e-6);
            double ymarg = 0.15 * std::max(gymax - gymin, 1e-6);
            cent_xlo[vi][ic] = std::max(0.0, gxmin - xmarg);
            cent_xhi[vi][ic] = gxmax + xmarg;
            cent_ylo[vi][ic] = gymin - ymarg;
            cent_yhi[vi][ic] = gymax + ymarg;
            std::cout << Form("Axis range v%d %s: x=[%.4f,%.4f] y=[%.4f,%.4f]\n",
                ivn, cen_name[ic],
                cent_xlo[vi][ic], cent_xhi[vi][ic],
                cent_ylo[vi][ic], cent_yhi[vi][ic]);
        }
    }

    // =========================================================================
    // SECTION 1: Build scatter graphs
    // Loop: centrality × harmonic × pT
    // =========================================================================
    std::cout << "\n=== SECTION 1: Building scatter graphs ===\n";

    for (int ic = 0; ic < N_CENTBINS; ++ic)
    {
        // Charged-particle vn per q-bin (pT-integrated 0.5-3.0 GeV/c)
        double chg_v2_q[N_QBINS] = {}, chg_v2_q_e[N_QBINS] = {};
        double chg_v3_q[N_QBINS] = {}, chg_v3_q_e[N_QBINS] = {};
        for (int iq = 0; iq < N_QBINS; ++iq) {
            if (hp_chg_v2[ic] && hp_chg_v2[ic]->GetBinEntries(iq + 1) > 0) {
                chg_v2_q[iq]   = hp_chg_v2[ic]->GetBinContent(iq + 1);
                chg_v2_q_e[iq] = hp_chg_v2[ic]->GetBinError(iq + 1);
            }
            if (hp_chg_v3[ic] && hp_chg_v3[ic]->GetBinEntries(iq + 1) > 0) {
                chg_v3_q[iq]   = hp_chg_v3[ic]->GetBinContent(iq + 1);
                chg_v3_q_e[iq] = hp_chg_v3[ic]->GetBinError(iq + 1);
            }
        }

        for (int ivn = 2; ivn <= 3; ++ivn)
        {
            int           n_pt       = (ivn == 2) ? N_PTBINS_V2  : N_PTBINS_V3;
            const double* pt_edges   = (ivn == 2) ? pt_edges_v2  : pt_edges_v3;
            double*       vn_chg_q   = (ivn == 2) ? chg_v2_q     : chg_v3_q;
            double*       vn_chg_q_e = (ivn == 2) ? chg_v2_q_e   : chg_v3_q_e;

            for (int ip = 0; ip < n_pt; ++ip)
            {
                const double pt_lo = pt_edges[ip];
                const double pt_hi = pt_edges[ip + 1];
                TString ptTag = Form("pT%dto%d", (int)pt_lo, (int)pt_hi);

                TGraphErrors* gr_d0_vn_q = dir_d0_qbin
                    ? (TGraphErrors*)dir_d0_qbin->Get(
                          Form("v%d_vs_q%dbin_%s_%s", ivn, ivn, cen_name[ic], ptTag.Data()))
                    : nullptr;

                if (!gr_d0_vn_q || gr_d0_vn_q->GetN() == 0) continue;

                // Build scatter points (one per q-bin)
                std::vector<double> xv, yv, exv, eyv;
                for (int iq = 0; iq < N_QBINS && iq < gr_d0_vn_q->GetN(); ++iq) {
                    double chg_v = vn_chg_q[iq];
                    double chg_e = vn_chg_q_e[iq];
                    double d0_v  = gr_d0_vn_q->GetPointY(iq);
                    double d0_e  = gr_d0_vn_q->GetErrorY(iq);
                    if (chg_e <= 0 || d0_e <= 0) continue;
                    xv.push_back(chg_v);  exv.push_back(chg_e);
                    yv.push_back(d0_v);   eyv.push_back(d0_e);
                }
                if ((int)xv.size() < 2) continue;

                // Normalize by mean vn if requested
                if (NORMALIZE_BY_MEAN) {
                    double mean_x = 0.0, mean_y = 0.0;
                    for (int i = 0; i < (int)xv.size(); ++i) { mean_x += xv[i]; mean_y += yv[i]; }
                    mean_x /= (double)xv.size();
                    mean_y /= (double)xv.size();
                    if (std::abs(mean_x) > 1e-9 && std::abs(mean_y) > 1e-9) {
                        for (int i = 0; i < (int)xv.size(); ++i) {
                            exv[i] /= std::abs(mean_x);  xv[i] /= mean_x;
                            eyv[i] /= std::abs(mean_y);  yv[i] /= mean_y;
                        }
                    }
                }

                TString sc_name = Form("scatter_v%d_%s_%s", ivn, cen_name[ic], ptTag.Data());

                TGraphErrors* gr_sc = new TGraphErrors(
                    (int)xv.size(), xv.data(), yv.data(), exv.data(), eyv.data());
                gr_sc->SetName(sc_name);
                gr_sc->SetTitle(sc_name);

                // Axis ranges: fixed when normalizing, per-centrality data-driven otherwise
                double xlo, xhi, ylo, yhi;
                if (NORMALIZE_BY_MEAN) {
                    xlo = X_AXIS_MIN;
                    xhi = X_AXIS_MAX;
                    ylo = Y_AXIS_MIN;
                    yhi = Y_AXIS_MAX;
                } else {
                    xlo = cent_xlo[ivn - 2][ic];
                    xhi = cent_xhi[ivn - 2][ic];
                    ylo = cent_ylo[ivn - 2][ic];
                    yhi = cent_yhi[ivn - 2][ic];
                }

                // Linear fit  y = p0 + p1*x
                TF1* flin = new TF1(Form("flin_%s", sc_name.Data()), "[0]+[1]*x", xlo, xhi);
                flin->SetParameters(0.0, 1.0);
                flin->SetLineColor(kRed);
                flin->SetLineWidth(2);
                gr_sc->Fit(flin, "QR");

                // 1-sigma confidence interval band
                TH1D* h_ci = new TH1D(Form("hci_%s", sc_name.Data()), "", 400, xlo, xhi);
                h_ci->SetDirectory(nullptr);
                TVirtualFitter* vfitter = TVirtualFitter::GetFitter();
                if (vfitter) vfitter->GetConfidenceIntervals(h_ci, 0.68);
                h_ci->SetFillColorAlpha(kRed, 0.35);
                h_ci->SetFillStyle(1001);
                h_ci->SetMarkerSize(0);

                // Canvas
                TCanvas* cv = new TCanvas(Form("cv_%s", sc_name.Data()), sc_name, 600, 600);
                cv->SetLeftMargin(0.14);
                cv->SetBottomMargin(0.14);
                cv->SetRightMargin(0.05);
                cv->SetTopMargin(0.07);
                cv->SetGridx();
                cv->SetGridy();

                TH2F* hf = new TH2F(Form("hf_%s", sc_name.Data()), "",
                                    10, xlo, xhi, 10, ylo, yhi);
                const char* x_title = NORMALIZE_BY_MEAN
                    ? Form("charged particle v_{%d} / <v_{%d}>", ivn, ivn)
                    : Form("charged particle v_{%d}", ivn);
                const char* y_title = NORMALIZE_BY_MEAN
                    ? Form("D^{0} v_{%d} / <v_{%d}>", ivn, ivn)
                    : Form("D^{0} v_{%d}", ivn);
                hf->GetXaxis()->SetTitle(x_title);
                hf->GetYaxis()->SetTitle(y_title);
                hf->GetXaxis()->SetTitleSize(0.050);
                hf->GetYaxis()->SetTitleSize(0.050);
                hf->GetXaxis()->SetTitleOffset(1.10);
                hf->GetYaxis()->SetTitleOffset(1.25);
                hf->GetXaxis()->SetLabelSize(0.040);
                hf->GetYaxis()->SetLabelSize(0.040);
                hf->Draw();

                h_ci->Draw("E3 SAME");
                flin->Draw("SAME");
                gr_sc->SetMarkerStyle(20);
                gr_sc->SetMarkerSize(0.9);
                gr_sc->SetMarkerColor(kBlack);
                gr_sc->SetLineColor(kBlack);
                gr_sc->Draw("P SAME");

                TLegend* leg = new TLegend(0.10, 0.81, 0.45, 0.92);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.040);
                leg->AddEntry((TObject*)nullptr,
                    Form("Centrality: %d-%d%%", min_cent[ic], max_cent[ic]), "");
                leg->AddEntry((TObject*)nullptr,
                    Form("%.0f < p_{T} < %.0f GeV/c", pt_lo, pt_hi), "");
                leg->Draw();

                ltx.SetTextAlign(31);
                ltx.SetTextSize(0.038);
                ltx.DrawLatex(0.93, 0.89,
                    Form("slope = %.3f #pm %.3f",
                         flin->GetParameter(1), flin->GetParError(1)));
                ltx.DrawLatex(0.93, 0.83,
                    Form("intercept = %.3f #pm %.3f",
                         flin->GetParameter(0), flin->GetParError(0)));
                ltx.DrawLatex(0.93, 0.77,
                    Form("#chi^{2}/NDF = %.2f",
                         flin->GetChisquare() / std::max(1.0, (double)flin->GetNDF())));
                ltx.SetTextAlign(11);
                ltx.SetTextFont(62);
                ltx.SetTextSize(0.048);
                ltx.DrawLatex(0.14, 0.945, "CMS");
                ltx.SetTextAlign(31);
                ltx.SetTextFont(42);
                ltx.SetTextSize(0.042);
                ltx.DrawLatex(0.95, 0.945, "PbPb 5.36 TeV");
                ltx.SetTextAlign(11);

                const char* subdir = (ivn == 2) ? "v2" : "v3";
                cv->SaveAs(Form("%s/scatter_plots/%s/%s.pdf", outbase.Data(), subdir, sc_name.Data()));

                // Write scatter graph and canvas to scatter ROOT file
                if (ivn == 2) sc_v2->cd(); else sc_v3->cd();
                gr_sc->Write();
                cv->Write();

                // Store graph pointer for Section 2 (NOT deleted here)
                gr_sc_mem[ivn - 2][ic][ip] = gr_sc;

                delete leg;
                delete hf;
                delete h_ci;
                delete flin;
                delete cv;

                std::cout << "Created scatter: " << sc_name << "\n";
            } // ip
        } // ivn
    } // ic

    // Write and close scatter file; input files no longer needed
    fSc->Write();
    fSc->Close();
    fD0->Close();
    fChg->Close();

    // =========================================================================
    // Merged scatter PDFs — tile the already-saved individual PDFs using
    // Python/PyMuPDF.  One page per centrality, 4-column grid.
    // Output is pixel-identical to the individual plots (no ROOT re-draw).
    // =========================================================================
    std::cout << "=== Generating merged scatter PDFs (Python/PyMuPDF) ===\n";
    {
        TString py_script = Form("%s/_merge_scatter_tmp.py", outbase.Data());
        FILE* fpy = fopen(py_script.Data(), "w");
        if (!fpy) {
            std::cerr << "WARNING: Cannot write temp Python script; merged PDFs skipped.\n";
        } else {
            fprintf(fpy,
                "import sys\n"
                "try:\n"
                "    import fitz\n"
                "except ModuleNotFoundError:\n"
                "    sys.exit('ERROR: PyMuPDF not found.  Install: pip install pymupdf')\n"
                "import os\n"
                "from collections import defaultdict\n"
                "\n"
                "CENT_ORDER=['0to10','10to20','20to30','30to40','40to50']\n"
                "\n"
                "def pt_key(f):\n"
                "    return int(f.split('_pT')[1].split('to')[0])\n"
                "\n"
                "def merge(inp, out, cols=4):\n"
                "    pdfs=sorted([f for f in os.listdir(inp)\n"
                "                 if f.endswith('.pdf') and not f.startswith('merged')],\n"
                "                key=pt_key)\n"
                "    if not pdfs: print('No PDFs in',inp); return\n"
                "    groups=defaultdict(list)\n"
                "    for f in pdfs:\n"
                "        try: cent=f.split('_cent')[1].split('_')[0]\n"
                "        except IndexError: continue\n"
                "        groups[cent].append(f)\n"
                "    with fitz.open(os.path.join(inp,pdfs[0])) as s:\n"
                "        w,h=s[0].rect.width,s[0].rect.height\n"
                "    doc=fitz.open()\n"
                "    for cent in [c for c in CENT_ORDER if c in groups]:\n"
                "        files=sorted(groups[cent],key=pt_key)\n"
                "        rows=-(-len(files)//cols)\n"
                "        page=doc.new_page(width=w*cols,height=h*rows)\n"
                "        for i,fn in enumerate(files[:cols*rows]):\n"
                "            ci,ri=i%%cols,i//cols\n"
                "            with fitz.open(os.path.join(inp,fn)) as src:\n"
                "                pix=src[0].get_pixmap(dpi=150)\n"
                "            page.insert_image(\n"
                "                fitz.Rect(ci*w,ri*h,(ci+1)*w,(ri+1)*h),pixmap=pix)\n"
                "        print(f'  cent{cent}: {len(files)} plots')\n"
                "    doc.save(out,garbage=4,deflate=True)\n"
                "    doc.close()\n"
                "    print('Saved ->',out)\n"
                "\n");
            fprintf(fpy,
                "merge('%s/scatter_plots/v2','%s/scatter_plots/v2/merged_scatter_v2.pdf')\n",
                outbase.Data(), outbase.Data());
            fprintf(fpy,
                "merge('%s/scatter_plots/v3','%s/scatter_plots/v3/merged_scatter_v3.pdf')\n",
                outbase.Data(), outbase.Data());
            fclose(fpy);

            // Auto-detect the Python interpreter that has PyMuPDF (fitz).
            // Tries each candidate in order; picks the first one where
            //   python -c 'import fitz'  exits 0.
            TString py_exe;
            {
                const char* cands[] = {
                    "/usr/local/bin/python3.11",
                    "/opt/homebrew/bin/python3",
                    "/usr/bin/python3",
                    "python3",
                    "python",
                    nullptr
                };
                for (int k = 0; cands[k] && py_exe.IsNull(); ++k) {
                    TString probe = Form("%s -c 'import fitz' 2>/dev/null", cands[k]);
                    if (gSystem->Exec(probe.Data()) == 0)
                        py_exe = cands[k];
                }
            }

            if (py_exe.IsNull()) {
                std::cerr << "WARNING: No Python interpreter with PyMuPDF found.\n"
                          << "  Install with:  pip install pymupdf\n"
                          << "  Merged PDFs skipped; individual PDFs are still in "
                          << outbase << "/scatter_plots/\n";
            } else {
                std::cout << "  Using Python: " << py_exe << "\n";
                int ret = gSystem->Exec(Form("%s '%s'", py_exe.Data(), py_script.Data()));
                if (ret != 0)
                    std::cerr << "WARNING: Python merge script failed (exit " << ret
                              << "). Individual PDFs are still in "
                              << outbase << "/scatter_plots/\n";
            }
            gSystem->Unlink(py_script.Data());
        }
    }

    std::cout << "  Individual scatter PDFs -> " << outbase << "/scatter_plots/v{2,3}/\n"
              << "  Scatter ROOT -> " << scatter_path << "\n";

    // =========================================================================
    // SECTION 2: Compute Pearson r from the in-memory scatter graphs
    // =========================================================================
    std::cout << "\n=== SECTION 2: Computing Pearson r ===\n";

    for (int ivn = 2; ivn <= 3; ++ivn)
    {
        TDirectory* outDir    = (ivn == 2) ? pr_v2 : pr_v3;
        int           n_pt        = (ivn == 2) ? N_PTBINS_V2  : N_PTBINS_V3;
        const double* pt_edges_vn = (ivn == 2) ? pt_edges_v2  : pt_edges_v3;
        double*       pt_cen_vn   = (ivn == 2) ? pt_cen_v2    : pt_cen_v3;
        double*       pt_elo_vn   = (ivn == 2) ? pt_elo_v2    : pt_elo_v3;
        double*       pt_ehi_vn   = (ivn == 2) ? pt_ehi_v2    : pt_ehi_v3;

        double r_val[N_CENTBINS][N_PTBINS_MAX] = {};
        double r_elo[N_CENTBINS][N_PTBINS_MAX] = {};
        double r_ehi[N_CENTBINS][N_PTBINS_MAX] = {};
        bool   r_ok [N_CENTBINS][N_PTBINS_MAX] = {};
        double s_val[N_CENTBINS][N_PTBINS_MAX] = {};
        double s_err[N_CENTBINS][N_PTBINS_MAX] = {};
        bool   s_ok [N_CENTBINS][N_PTBINS_MAX] = {};
        double b_val[N_CENTBINS][N_PTBINS_MAX] = {};
        double b_err[N_CENTBINS][N_PTBINS_MAX] = {};
        bool   b_ok [N_CENTBINS][N_PTBINS_MAX] = {};
        TH1D*  h_mc_arr[N_CENTBINS][N_PTBINS_MAX] = {};

        // =====================================================================
        // Compute r and MC uncertainty for every (centrality, pT) cell
        // =====================================================================
        for (int ic = 0; ic < N_CENTBINS; ++ic)
        {
            for (int ip = 0; ip < n_pt; ++ip)
            {
                // Use the scatter graph built in Section 1 directly from memory
                TGraphErrors* gr = gr_sc_mem[ivn - 2][ic][ip];
                if (!gr || gr->GetN() < 2) continue;

                TString ptTag_ip = Form("pT%dto%d",
                    (int)pt_edges_vn[ip], (int)pt_edges_vn[ip+1]);

                int np = gr->GetN();
                std::vector<double> xv(np), yv(np), exv(np), eyv(np);
                for (int i = 0; i < np; ++i) {
                    xv[i]  = gr->GetPointX(i);
                    yv[i]  = gr->GetPointY(i);
                    exv[i] = gr->GetErrorX(i);
                    eyv[i] = gr->GetErrorY(i);
                }

                // Measured Pearson r
                double r_meas = WeightedPearson(xv, yv, exv, eyv);
                if (std::isnan(r_meas)) continue;

                // MC resampling: smear D0 y-values only
                std::vector<double> r_mc;
                r_mc.reserve(N_MC);
                for (int imc = 0; imc < N_MC; ++imc) {
                    std::vector<double> yv_mc(np);
                    for (int i = 0; i < np; ++i) {
                        std::normal_distribution<double> gaus(yv[i], eyv[i]);
                        yv_mc[i] = gaus(rng);
                    }
                    double r_tmp = WeightedPearson(xv, yv_mc, exv, eyv);
                    if (!std::isnan(r_tmp)) r_mc.push_back(r_tmp);
                }
                if (static_cast<int>(r_mc.size()) < 100) continue;

                std::sort(r_mc.begin(), r_mc.end());
                int    n_mc = static_cast<int>(r_mc.size());
                double r_elo_v, r_ehi_v;

                if (!USE_SYMMETRIC_INTERVAL) {
                    // Equal-tails percentile interval
                    double p16 = r_mc[static_cast<int>(0.1585 * n_mc)];
                    double p84 = r_mc[static_cast<int>(0.8415 * n_mc)];
                    r_elo_v = r_meas - p16;
                    r_ehi_v = p84 - r_meas;
                } else {
                    // Symmetric expansion around r_meas (Run-2 style)
                    double lo_d = 0.0, hi_d = 2.0;
                    for (int iter = 0; iter < 60; ++iter) {
                        double mid = 0.5 * (lo_d + hi_d);
                        auto it_lo = std::lower_bound(r_mc.begin(), r_mc.end(), r_meas - mid);
                        auto it_hi = std::upper_bound(r_mc.begin(), r_mc.end(), r_meas + mid);
                        double frac = static_cast<double>(
                            std::distance(it_lo, it_hi)) / n_mc;
                        if (frac >= 0.683) hi_d = mid;
                        else               lo_d = mid;
                    }
                    r_elo_v = hi_d;
                    r_ehi_v = hi_d;
                }

                r_val[ic][ip] = r_meas;
                r_elo[ic][ip] = r_elo_v;
                r_ehi[ic][ip] = r_ehi_v;
                r_ok [ic][ip] = true;

                // Linear fit  y = p0 + p1*x  (for slope and intercept)
                {
                    double xmin_f = *std::min_element(xv.begin(), xv.end());
                    double xmax_f = *std::max_element(xv.begin(), xv.end());
                    double xmarg  = 0.3 * std::max(xmax_f - xmin_f, 1e-6);
                    TGraphErrors* gr_fit = new TGraphErrors(
                        np, xv.data(), yv.data(), exv.data(), eyv.data());
                    TF1* flin_si = new TF1(
                        Form("flin_si_%d_%d_%d", ivn, ic, ip),
                        "[0]+[1]*x", xmin_f - xmarg, xmax_f + xmarg);
                    flin_si->SetParameters(0.0, 1.0);
                    gr_fit->Fit(flin_si, "QR0");
                    double perr0 = flin_si->GetParError(0);
                    double perr1 = flin_si->GetParError(1);
                    if (perr0 > 0 && perr1 > 0) {
                        s_val[ic][ip] = flin_si->GetParameter(1);
                        s_err[ic][ip] = perr1;
                        s_ok [ic][ip] = true;
                        b_val[ic][ip] = flin_si->GetParameter(0);
                        b_err[ic][ip] = perr0;
                        b_ok [ic][ip] = true;
                    }
                    delete flin_si;
                    delete gr_fit;
                }

                std::cout << Form("v%d %s %s : r = %.3f  +%.3f -%.3f  (N_q=%d)\n",
                    ivn, cen_name[ic], ptTag_ip.Data(),
                    r_meas, r_ehi[ic][ip], r_elo[ic][ip], np);

                // MC r-distribution histogram (stored for combined canvas)
                TH1D* h_mc = new TH1D(
                    Form("h_mc_v%d_%s_%s", ivn, cen_name[ic], ptTag_ip.Data()),
                    Form(";r (MC);counts"), 200, -1.0, 1.0);
                h_mc->SetDirectory(nullptr);
                for (double rv : r_mc) h_mc->Fill(rv);
                h_mc->SetFillColor(kBlue-9);
                h_mc->SetFillStyle(1001);
                h_mc->SetLineColor(kBlue+1);
                h_mc_arr[ic][ip] = h_mc;
                outDir->cd();
                h_mc->Write();

            } // ip
        } // ic

        // =====================================================================
        // Build TGraphAsymmErrors per centrality (Pearson r vs pT)
        // =====================================================================
        TGraphAsymmErrors* gr_r[N_CENTBINS] = {};
        for (int ic = 0; ic < N_CENTBINS; ++ic)
        {
            std::vector<double> xv, yv, exlo, exhi, eylo, eyhi;
            for (int ip = 0; ip < n_pt; ++ip) {
                if (!r_ok[ic][ip]) continue;
                xv  .push_back(pt_cen_vn[ip]);
                yv  .push_back(r_val[ic][ip]);
                exlo.push_back(pt_elo_vn[ip]);
                exhi.push_back(pt_ehi_vn[ip]);
                eylo.push_back(r_elo[ic][ip]);
                eyhi.push_back(r_ehi[ic][ip]);
            }
            if (xv.empty()) continue;

            gr_r[ic] = new TGraphAsymmErrors(
                static_cast<int>(xv.size()),
                xv.data(), yv.data(),
                exlo.data(), exhi.data(),
                eylo.data(), eyhi.data());
            gr_r[ic]->SetName(Form("pearson_v%d_%s", ivn, cen_name[ic]));
            gr_r[ic]->SetTitle(Form("Pearson v%d %s", ivn, cen_name[ic]));
            gr_r[ic]->SetMarkerStyle(marker[ic]);
            gr_r[ic]->SetMarkerSize(1.1);
            gr_r[ic]->SetMarkerColor(col[ic]);
            gr_r[ic]->SetLineColor(col[ic]);
            gr_r[ic]->SetLineWidth(1);
            outDir->cd();
            gr_r[ic]->Write();
        }

        // =====================================================================
        // Build TGraphAsymmErrors for slope and intercept vs pT
        // =====================================================================
        TGraphAsymmErrors* gr_slope[N_CENTBINS] = {};
        TGraphAsymmErrors* gr_inter[N_CENTBINS] = {};
        for (int ic = 0; ic < N_CENTBINS; ++ic) {
            std::vector<double> xs, ys, exlos, exhis, eylos, eyhis;
            std::vector<double> xb, yb, exlob, exhib, eylob, eyhib;
            for (int ip = 0; ip < n_pt; ++ip) {
                if (s_ok[ic][ip]) {
                    xs   .push_back(pt_cen_vn[ip]);
                    ys   .push_back(s_val[ic][ip]);
                    exlos.push_back(pt_elo_vn[ip]);
                    exhis.push_back(pt_ehi_vn[ip]);
                    eylos.push_back(s_err[ic][ip]);
                    eyhis.push_back(s_err[ic][ip]);
                }
                if (b_ok[ic][ip]) {
                    xb   .push_back(pt_cen_vn[ip]);
                    yb   .push_back(b_val[ic][ip]);
                    exlob.push_back(pt_elo_vn[ip]);
                    exhib.push_back(pt_ehi_vn[ip]);
                    eylob.push_back(b_err[ic][ip]);
                    eyhib.push_back(b_err[ic][ip]);
                }
            }
            if (!xs.empty()) {
                gr_slope[ic] = new TGraphAsymmErrors(
                    static_cast<int>(xs.size()),
                    xs.data(), ys.data(),
                    exlos.data(), exhis.data(), eylos.data(), eyhis.data());
                gr_slope[ic]->SetName(Form("slope_v%d_%s", ivn, cen_name[ic]));
                gr_slope[ic]->SetTitle(Form("Slope v%d %s", ivn, cen_name[ic]));
                gr_slope[ic]->SetMarkerStyle(marker[ic]);
                gr_slope[ic]->SetMarkerSize(1.1);
                gr_slope[ic]->SetMarkerColor(col[ic]);
                gr_slope[ic]->SetLineColor(col[ic]);
                gr_slope[ic]->SetLineWidth(1);
                outDir->cd();
                gr_slope[ic]->Write();
            }
            if (!xb.empty()) {
                gr_inter[ic] = new TGraphAsymmErrors(
                    static_cast<int>(xb.size()),
                    xb.data(), yb.data(),
                    exlob.data(), exhib.data(), eylob.data(), eyhib.data());
                gr_inter[ic]->SetName(Form("inter_v%d_%s", ivn, cen_name[ic]));
                gr_inter[ic]->SetTitle(Form("Intercept v%d %s", ivn, cen_name[ic]));
                gr_inter[ic]->SetMarkerStyle(marker[ic]);
                gr_inter[ic]->SetMarkerSize(1.1);
                gr_inter[ic]->SetMarkerColor(col[ic]);
                gr_inter[ic]->SetLineColor(col[ic]);
                gr_inter[ic]->SetLineWidth(1);
                outDir->cd();
                gr_inter[ic]->Write();
            }
        }

        // =====================================================================
        // Combined two-panel Pearson canvas  (left: 0-20%, right: 20-50%)
        // =====================================================================
        const int PANEL_IC[2][3] = { {0, 1, -1}, {2, 3, 4} };
        const int PANEL_NC[2]    = { 2, 3 };

        const double Y_LO =  -1.0;
        const double Y_HI =   2.5;
        const double X_LO =   1.0;
        const double X_HI = 110.0;

        TCanvas* cv_comb = new TCanvas(
            Form("cv_pearson_v%d", ivn),
            Form("Pearson r vs pT (v%d)", ivn),
            1200, 600);
        cv_comb->SetFillColor(0);
        cv_comb->SetBorderSize(0);

        const double bm    = 0.15;
        const double tm    = 0.09;
        const double split = 0.50;

        TPad* pad[2];
        pad[0] = new TPad("pad0","", 0.00, 0.00, split, 1.00);
        pad[1] = new TPad("pad1","", split, 0.00, 1.00,  1.00);

        pad[0]->SetLeftMargin(0.14);  pad[0]->SetRightMargin(0.000);
        pad[0]->SetBottomMargin(bm);  pad[0]->SetTopMargin(tm);
        pad[0]->SetLogx();            pad[0]->SetFillColor(0);

        pad[1]->SetLeftMargin(0.000); pad[1]->SetRightMargin(0.05);
        pad[1]->SetBottomMargin(bm);  pad[1]->SetTopMargin(tm);
        pad[1]->SetLogx();            pad[1]->SetFillColor(0);

        cv_comb->cd(); pad[0]->Draw();
        cv_comb->cd(); pad[1]->Draw();

        for (int ipanel = 0; ipanel < 2; ++ipanel)
        {
            pad[ipanel]->cd();

            TH2F* hf = new TH2F(Form("hf_v%d_p%d", ivn, ipanel), "",
                                 100, X_LO, X_HI, 100, Y_LO, Y_HI);
            hf->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            hf->GetXaxis()->SetTitleSize(0.060);
            hf->GetXaxis()->SetTitleOffset(0.95);
            hf->GetXaxis()->SetLabelSize(0.050);
            hf->GetXaxis()->SetMoreLogLabels();
            hf->GetXaxis()->SetNoExponent();
            hf->GetXaxis()->CenterTitle(kTRUE);
            if (ipanel == 0) {
                hf->GetYaxis()->SetTitle("r");
                hf->GetYaxis()->SetTitleSize(0.072);
                hf->GetYaxis()->SetTitleOffset(0.82);
                hf->GetYaxis()->SetLabelSize(0.050);
                hf->GetYaxis()->CenterTitle(kTRUE);
            } else {
                hf->GetYaxis()->SetTitle("");
                hf->GetYaxis()->SetLabelSize(0.0);
                hf->GetYaxis()->SetTickLength(0.03);
            }
            hf->Draw("AXIS");

            TLine* lref = new TLine(X_LO, 1.0, X_HI, 1.0);
            lref->SetLineStyle(2);
            lref->SetLineWidth(1);
            lref->SetLineColor(kGray+2);
            lref->Draw();

            for (int k = 0; k < PANEL_NC[ipanel]; ++k) {
                int ic = PANEL_IC[ipanel][k];
                if (ic < 0 || !gr_r[ic]) continue;
                gr_r[ic]->Draw("P Z SAME");
            }

            // Centrality legend — top right of each panel
            if (ipanel == 0) {
                const double leg_x2 = 0.97;
                const double leg_x1 = leg_x2 - 0.38;
                const double leg_y2 = 0.88;
                const double leg_y1 = leg_y2 - 0.082 * (PANEL_NC[ipanel] + 1);
                TLegend* leg = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.050);
                leg->AddEntry((TObject*)nullptr, "Centrality", "");
                for (int k = 0; k < PANEL_NC[ipanel]; ++k) {
                    int ic = PANEL_IC[ipanel][k];
                    if (ic < 0 || !gr_r[ic]) continue;
                    leg->AddEntry(gr_r[ic],
                        Form("%d-%d%%", min_cent[ic], max_cent[ic]), "p");
                }
                leg->Draw();
            } else {
                const double leg_x2 = 0.94;
                const double leg_x1 = leg_x2 - 0.40;
                const double leg_y2 = 0.90;
                const double leg_y1 = leg_y2 - 0.082 * (PANEL_NC[ipanel] + 1);
                TLegend* leg = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.050);
                leg->AddEntry((TObject*)nullptr, "Centrality", "");
                for (int k = 0; k < PANEL_NC[ipanel]; ++k) {
                    int ic = PANEL_IC[ipanel][k];
                    if (ic < 0 || !gr_r[ic]) continue;
                    leg->AddEntry(gr_r[ic],
                        Form("%d-%d%%", min_cent[ic], max_cent[ic]), "p");
                }
                leg->Draw();
            }

            TLatex ltx2;
            ltx2.SetNDC();
            if (ipanel == 0) {
                ltx2.SetTextFont(62);
                ltx2.SetTextSize(0.065);
                ltx2.SetTextAlign(11);
                ltx2.DrawLatex(0.17, 0.930, "CMS");
            }
            ltx2.SetTextFont(42);
            ltx2.SetTextSize(0.050);
            ltx2.SetTextAlign(11);
            ltx2.DrawLatex((ipanel == 0) ? 0.17 : 0.05, 0.860, "|y| < 1");

        } // ipanel

        cv_comb->cd();
        TLatex ltx_comb;
        ltx_comb.SetNDC();
        ltx_comb.SetTextFont(42);
        ltx_comb.SetTextSize(0.050);
        ltx_comb.SetTextAlign(31);
        ltx_comb.DrawLatex(0.975, 0.940, "PbPb 5.36 TeV");

        cv_comb->SaveAs(Form("%s/pearson_plots/v%d/pearson_v%d_vs_pT_combined.pdf", outbase.Data(), ivn, ivn));
        outDir->cd();
        cv_comb->Write();
        delete cv_comb;

        // =====================================================================
        // Slope + Intercept combined 4-pad canvas
        // Top row: slope vs pT;  Bottom row: intercept vs pT
        // Left column: 0-20%;   Right column: 20-50%
        // =====================================================================
        {
        const double SL_LO = -1.0, SL_HI = 4.0;
        const double IN_LO = -2.5, IN_HI = 2.5;
        const double y_lo_si[2] = {SL_LO, IN_LO};
        const double y_hi_si[2] = {SL_HI, IN_HI};
        const double y_ref_si[2] = {1.0, 0.0};
        const char*  y_tit_si[2] = {"Slope", "Intercept"};

        TCanvas* cv_si = new TCanvas(
            Form("cv_si_v%d", ivn),
            Form("Slope & Intercept vs pT (v%d)", ivn),
            1200, 800);
        cv_si->SetFillColor(0);
        cv_si->SetBorderSize(0);

        TPad* padsi[4];
        padsi[0] = new TPad("padsi0","", 0.00, 0.50, 0.50, 1.00);
        padsi[1] = new TPad("padsi1","", 0.50, 0.50, 1.00, 1.00);
        padsi[2] = new TPad("padsi2","", 0.00, 0.00, 0.50, 0.50);
        padsi[3] = new TPad("padsi3","", 0.50, 0.00, 1.00, 0.50);

        padsi[0]->SetLeftMargin(0.14);  padsi[0]->SetRightMargin(0.000);
        padsi[0]->SetBottomMargin(0.00); padsi[0]->SetTopMargin(0.12);
        padsi[1]->SetLeftMargin(0.000); padsi[1]->SetRightMargin(0.05);
        padsi[1]->SetBottomMargin(0.00); padsi[1]->SetTopMargin(0.12);
        padsi[2]->SetLeftMargin(0.14);  padsi[2]->SetRightMargin(0.000);
        padsi[2]->SetBottomMargin(0.28); padsi[2]->SetTopMargin(0.00);
        padsi[3]->SetLeftMargin(0.000); padsi[3]->SetRightMargin(0.05);
        padsi[3]->SetBottomMargin(0.28); padsi[3]->SetTopMargin(0.00);

        for (int i = 0; i < 4; ++i) {
            padsi[i]->SetLogx();
            padsi[i]->SetFillColor(0);
            cv_si->cd();
            padsi[i]->Draw();
        }

        for (int irow = 0; irow < 2; ++irow) {
            TGraphAsymmErrors** gr_cur = (irow == 0) ? gr_slope : gr_inter;
            for (int ipanel = 0; ipanel < 2; ++ipanel) {
                int  ipad    = irow * 2 + ipanel;
                bool is_top  = (irow   == 0);
                bool is_left = (ipanel == 0);
                padsi[ipad]->cd();

                TH2F* hfsi = new TH2F(
                    Form("hfsi_v%d_%d_%d", ivn, irow, ipanel), "",
                    100, X_LO, X_HI, 100, y_lo_si[irow], y_hi_si[irow]);
                if (!is_top) {
                    hfsi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                    hfsi->GetXaxis()->SetTitleSize(0.090);
                    hfsi->GetXaxis()->SetTitleOffset(0.80);
                    hfsi->GetXaxis()->SetLabelSize(0.075);
                    hfsi->GetXaxis()->CenterTitle(kTRUE);
                } else {
                    hfsi->GetXaxis()->SetTitle("");
                    hfsi->GetXaxis()->SetLabelSize(0.0);
                    hfsi->GetXaxis()->SetTickLength(0.06);
                }
                hfsi->GetXaxis()->SetMoreLogLabels();
                hfsi->GetXaxis()->SetNoExponent();
                if (is_left) {
                    hfsi->GetYaxis()->SetTitle(y_tit_si[irow]);
                    hfsi->GetYaxis()->SetTitleSize(0.090);
                    hfsi->GetYaxis()->SetTitleOffset(0.72);
                    hfsi->GetYaxis()->SetLabelSize(0.075);
                    hfsi->GetYaxis()->CenterTitle(kTRUE);
                } else {
                    hfsi->GetYaxis()->SetTitle("");
                    hfsi->GetYaxis()->SetLabelSize(0.0);
                    hfsi->GetYaxis()->SetTickLength(0.03);
                }
                hfsi->Draw("AXIS");

                TLine* lref_si = new TLine(X_LO, y_ref_si[irow], X_HI, y_ref_si[irow]);
                lref_si->SetLineStyle(2);
                lref_si->SetLineWidth(1);
                lref_si->SetLineColor(kGray+2);
                lref_si->Draw();

                for (int k = 0; k < PANEL_NC[ipanel]; ++k) {
                    int ic = PANEL_IC[ipanel][k];
                    if (ic < 0 || !gr_cur[ic]) continue;
                    gr_cur[ic]->Draw("P Z SAME");
                }

                if (is_top) {
                    TLatex ltxsi;
                    ltxsi.SetNDC();
                    if (is_left) {
                        ltxsi.SetTextFont(62);
                        ltxsi.SetTextSize(0.095);
                        ltxsi.SetTextAlign(11);
                        ltxsi.DrawLatex(0.17, 0.915, "CMS");
                    }
                    const double sleg_x2 = is_left ? 0.97 : 0.94;
                    const double sleg_x1 = sleg_x2 - 0.44;
                    const double sleg_y2 = 0.87;
                    const double sleg_y1 = sleg_y2 - 0.095 * (PANEL_NC[ipanel] + 1);
                    TLegend* sleg = new TLegend(sleg_x1, sleg_y1, sleg_x2, sleg_y2);
                    sleg->SetBorderSize(0);
                    sleg->SetFillStyle(0);
                    sleg->SetTextFont(42);
                    sleg->SetTextSize(0.078);
                    sleg->AddEntry((TObject*)nullptr, "Centrality", "");
                    for (int k = 0; k < PANEL_NC[ipanel]; ++k) {
                        int ic = PANEL_IC[ipanel][k];
                        if (ic < 0 || !gr_cur[ic]) continue;
                        sleg->AddEntry(gr_cur[ic],
                            Form("%d-%d%%", min_cent[ic], max_cent[ic]), "p");
                    }
                    sleg->Draw();
                    ltxsi.SetTextFont(42);
                    ltxsi.SetTextSize(0.078);
                    ltxsi.SetTextAlign(11);
                    ltxsi.DrawLatex(is_left ? 0.17 : 0.05, 0.80, "|y| < 1");
                }
            } // ipanel
        } // irow

        cv_si->cd();
        TLatex ltx_si_top;
        ltx_si_top.SetNDC();
        ltx_si_top.SetTextFont(42);
        ltx_si_top.SetTextSize(0.035);
        ltx_si_top.SetTextAlign(31);
        ltx_si_top.DrawLatex(0.975, 0.955, "PbPb 5.36 TeV");

        cv_si->SaveAs(Form("%s/slope_intercept_plots/slope_intercept_v%d_vs_pT.pdf", outbase.Data(), ivn));
        outDir->cd();
        cv_si->Write();
        delete cv_si;
        } // slope+intercept canvas scope

        // =====================================================================
        // MC r-distribution combined canvases (4 x N_ROWS grid per centrality)
        // =====================================================================
        for (int ic = 0; ic < N_CENTBINS; ++ic)
        {
            bool any_mc = false;
            for (int ip = 0; ip < n_pt; ++ip)
                if (h_mc_arr[ic][ip]) { any_mc = true; break; }
            if (!any_mc) continue;

            const int N_COLS_MC = 4;
            const int N_ROWS_MC = (n_pt + N_COLS_MC - 1) / N_COLS_MC;
            TCanvas* cv_mc = new TCanvas(
                Form("cv_mc_v%d_%s", ivn, cen_name[ic]),
                Form("MC r v%d %s", ivn, cen_name[ic]),
                1600, N_ROWS_MC * 300);
            cv_mc->Divide(N_COLS_MC, N_ROWS_MC, 0.001, 0.001);

            for (int ip = 0; ip < n_pt; ++ip)
            {
                cv_mc->cd(ip + 1);
                gPad->SetLeftMargin(0.13);
                gPad->SetRightMargin(0.03);
                gPad->SetBottomMargin(0.16);
                gPad->SetTopMargin(0.10);

                TH1D* hh = h_mc_arr[ic][ip];
                if (!hh) {
                    TH2F* hblank = new TH2F(
                        Form("hblank_%d_%d_%d", ivn, ic, ip), "",
                        10, -1.0, 1.0, 10, 0, 10);
                    hblank->GetXaxis()->SetTitle("r");
                    hblank->GetXaxis()->SetRangeUser(-1.0, 1.0);
                    hblank->Draw();
                    continue;
                }

                TH1D* hdraw = (TH1D*)hh->Clone(Form("hdraw_%d_%d_%d", ivn, ic, ip));
                hdraw->SetDirectory(nullptr);
                hdraw->GetXaxis()->SetTitle("r");
                hdraw->GetXaxis()->SetTitleSize(0.070);
                hdraw->GetXaxis()->SetTitleOffset(0.90);
                hdraw->GetXaxis()->SetLabelSize(0.060);
                bool wide_range = (ip == 0 || pt_edges_vn[ip] >= 10.0);
                hdraw->GetXaxis()->SetRangeUser(wide_range ? -1.0 : 0.0, 1.0);
                hdraw->GetYaxis()->SetTitle("Counts");
                hdraw->GetYaxis()->SetTitleSize(0.070);
                hdraw->GetYaxis()->SetTitleOffset(0.85);
                hdraw->GetYaxis()->SetLabelSize(0.060);
                hdraw->GetYaxis()->SetMaxDigits(3);
                gStyle->SetStripDecimals(kFALSE);
                hdraw->Draw("HIST");

                double ymax = hdraw->GetMaximum();

                if (r_ok[ic][ip]) {
                    TLine* lr = new TLine(r_val[ic][ip], 0, r_val[ic][ip], ymax);
                    lr->SetLineColor(kRed);
                    lr->SetLineWidth(2);
                    lr->Draw();
                    double p16v = r_val[ic][ip] - r_elo[ic][ip];
                    double p84v = r_val[ic][ip] + r_ehi[ic][ip];
                    TLine* ll16 = new TLine(p16v, 0, p16v, ymax * 0.65);
                    ll16->SetLineColor(kRed); ll16->SetLineStyle(2); ll16->SetLineWidth(1);
                    ll16->Draw();
                    TLine* ll84 = new TLine(p84v, 0, p84v, ymax * 0.65);
                    ll84->SetLineColor(kRed); ll84->SetLineStyle(2); ll84->SetLineWidth(1);
                    ll84->Draw();
                }

                TLegend* leg_mc = new TLegend(0.11, 0.73, 0.48, 0.9);
                leg_mc->SetBorderSize(0);
                leg_mc->SetFillStyle(0);
                leg_mc->SetTextFont(42);
                leg_mc->SetTextSize(0.072);
                leg_mc->AddEntry((TObject*)nullptr,
                    Form("Cent: %d-%d%%", min_cent[ic], max_cent[ic]), "");
                leg_mc->AddEntry((TObject*)nullptr,
                    Form("%.0f< p_{T} < %.0f GeV/c", pt_edges_vn[ip], pt_edges_vn[ip+1]), "");
                leg_mc->Draw();

            } // ip

            cv_mc->cd();
            TLatex ltx_cen_mc;
            ltx_cen_mc.SetNDC();
            ltx_cen_mc.SetTextFont(42);
            ltx_cen_mc.SetTextSize(0.022);
            ltx_cen_mc.SetTextAlign(11);
            ltx_cen_mc.DrawLatex(0.01, 0.003,
                Form("v_{%d}  |  Centrality %d-%d%%  |  PbPb 5.36 TeV",
                     ivn, min_cent[ic], max_cent[ic]));

            cv_mc->SaveAs(Form("%s/mc_r_plots/v%d/mc_r_v%d_%s.pdf", outbase.Data(), ivn, ivn, cen_name[ic]));
            outDir->cd();
            cv_mc->Write();
            delete cv_mc;

        } // ic (MC r canvases)

        // Delete stored MC histograms
        for (int ic = 0; ic < N_CENTBINS; ++ic)
            for (int ip = 0; ip < n_pt; ++ip)
                { delete h_mc_arr[ic][ip]; h_mc_arr[ic][ip] = nullptr; }

        // =====================================================================
        // Per-centrality individual Pearson r canvas
        // =====================================================================
        for (int ic = 0; ic < N_CENTBINS; ++ic)
        {
            if (!gr_r[ic]) continue;

            TCanvas* cv2 = new TCanvas(
                Form("cv_r_v%d_%s", ivn, cen_name[ic]), "", 640, 580);
            cv2->SetLeftMargin(0.14);
            cv2->SetRightMargin(0.05);
            cv2->SetBottomMargin(0.14);
            cv2->SetTopMargin(0.09);
            cv2->SetLogx();
            cv2->SetGridx();
            cv2->SetGridy();

            TH2F* hf2 = new TH2F(
                Form("hf2_v%d_%s", ivn, cen_name[ic]), "",
                100, X_LO, X_HI, 100, Y_LO, Y_HI);
            hf2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            hf2->GetXaxis()->SetTitleSize(0.055);
            hf2->GetXaxis()->SetTitleOffset(1.05);
            hf2->GetXaxis()->SetLabelSize(0.045);
            hf2->GetXaxis()->SetMoreLogLabels();
            hf2->GetXaxis()->SetNoExponent();
            hf2->GetYaxis()->SetTitle("r");
            hf2->GetYaxis()->SetTitleSize(0.055);
            hf2->GetYaxis()->SetTitleOffset(1.20);
            hf2->GetYaxis()->SetLabelSize(0.045);
            hf2->Draw("AXIS");

            TLine* lref2 = new TLine(X_LO, 1.0, X_HI, 1.0);
            lref2->SetLineStyle(2);
            lref2->SetLineWidth(1);
            lref2->SetLineColor(kGray+2);
            lref2->Draw();

            gr_r[ic]->Draw("P Z SAME");

            TLegend* leg2 = new TLegend(0.52, 0.82, 0.93, 0.90);
            leg2->SetBorderSize(0);
            leg2->SetFillStyle(0);
            leg2->SetTextFont(42);
            leg2->SetTextSize(0.045);
            leg2->AddEntry(gr_r[ic],
                Form("Centrality: %d-%d%%", min_cent[ic], max_cent[ic]), "p");
            leg2->Draw();

            TLatex ltx2;
            ltx2.SetNDC();
            ltx2.SetTextAlign(11);
            ltx2.SetTextFont(62);
            ltx2.SetTextSize(0.058);
            ltx2.DrawLatex(0.15, 0.935, "CMS");
            ltx2.SetTextFont(42);
            ltx2.SetTextSize(0.045);
            ltx2.DrawLatex(0.15, 0.872, "PbPb 5.36 TeV");
            ltx2.SetTextAlign(31);
            ltx2.DrawLatex(0.94, 0.935, "|y| < 1");

            cv2->SaveAs(Form("%s/pearson_plots/v%d/pearson_v%d_%s.pdf",
                             outbase.Data(), ivn, ivn, cen_name[ic]));
            outDir->cd();
            cv2->Write();

            delete hf2;
            delete lref2;
            delete leg2;
            delete cv2;
        } // ic

        std::cout << "v" << ivn << " Pearson plots complete.\n";

    } // ivn

    // =========================================================================
    // Clean up in-memory scatter graphs
    // =========================================================================
    for (int iv = 0; iv < 2; ++iv)
        for (int ic = 0; ic < N_CENTBINS; ++ic)
            for (int ip = 0; ip < N_PTBINS_MAX; ++ip)
                { delete gr_sc_mem[iv][ic][ip]; gr_sc_mem[iv][ic][ip] = nullptr; }

    fPr->Write();
    fPr->Close();

    std::cout << "\n=== Done ===\n"
              << "  Output directory        -> " << outbase << "/\n"
              << "  Individual scatter PDFs -> " << outbase << "/scatter_plots/v{2,3}/\n"
              << "  Merged scatter PDFs     -> " << outbase << "/scatter_plots/v2/merged_scatter_v2.pdf\n"
              << "                         -> " << outbase << "/scatter_plots/v3/merged_scatter_v3.pdf\n"
              << "  Pearson PDFs            -> " << outbase << "/pearson_plots/v{2,3}/\n"
              << "  Slope/Int.              -> " << outbase << "/slope_intercept_plots/\n"
              << "  MC r plots              -> " << outbase << "/mc_r_plots/v{2,3}/\n"
              << "  Scatter ROOT            -> " << scatter_path << "\n"
              << "  Pearson ROOT            -> " << pearson_path << "\n";
}


