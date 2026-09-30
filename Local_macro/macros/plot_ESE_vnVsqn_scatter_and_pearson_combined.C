//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// plot_ESE_scatter_and_pearson_combined.C
//
// Combined macro that generates ALL ESE D0-vs-charged-particle vn plots in
// one run, from the same pair of input ROOT files:
//   0) vn vs q-bin-index and vn vs mean-q, D0 & charged particle overlaid,
//      per (centrality x D0-pT) bin  [ported from plot_vn_vs_qbin_D0Chg_combined.C]
//   1) Scatter plots of D0 vn vs charged-particle vn per q-bin,
//      for every (centrality × D0-pT) bin  [make_ESE_scatter.C]
//   2) Weighted Pearson correlation coefficient r between D0 vn and
//      charged-particle vn per q-bin, as a function of D0 pT,
//      for each centrality class  [plot_pearson_ESE.C]
//
// Each numbered step (SECTION 0, 1, 2) is a clearly delimited block in the
// code below, so the individual steps stay easy to follow even though they
// all now run from a single macro call. The scatter TGraphErrors built in
// Section 1 are passed directly in memory to Section 2 — no intermediate
// file read/write is needed between those two steps.
//
// Formula (Pearson):
//          sum_i  w_i (x_i - <x_w>) (y_i - <y_w>)
//   r = -----------------------------------------------
//       sqrt[ sum_i w_i(x_i-<x_w>)^2 * sum_i w_i(y_i-<y_w>)^2 ]
//
//   w_i = 1/sigma_i^2,  sigma_i = sqrt((Dx_i)^2 + (Dy_i)^2)
//
// Outputs — all written inside a single tagged output directory. Subdirectories
// are numbered 1_.. through 6_.. so they list in pipeline order in any
// alphabetically-sorting file browser (Finder/VSCode/ls), not just in this comment:
//   {outdir}/1_vn_vs_qbin_plots/v{2,3}/*.pdf
//   {outdir}/1_vn_vs_qbin_plots/summary_v{2,3}_qbin_{vn_qbin_tag}.pdf   (merged/tiled)
//   {outdir}/2_vn_vs_meanq_plots/v{2,3}/*.pdf
//   {outdir}/2_vn_vs_meanq_plots/summary_v{2,3}_meanq_{vn_qbin_tag}.pdf (merged/tiled)
//   {outdir}/ESE_vn_vs_qbin_{vn_qbin_tag}.root
//   {outdir}/ESE_vn_vs_meanq_{vn_qbin_tag}.root
//   {outdir}/3_scatter_plots/v{2,3}/*.pdf
//   {outdir}/3_scatter_plots/v{2,3}/merged_scatter_v{2,3}.pdf
//   {outdir}/4_pearson_plots/v{2,3}/*.pdf
//   {outdir}/5_slope_intercept_plots/slope_intercept_v{2,3}_vs_pT.pdf
//   {outdir}/6_mc_r_plots/v{2,3}/*.pdf
//   {outdir}/6_mc_r_plots/v{2,3}/merged_mc_r_v{2,3}.pdf
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

// #########################################################################
// ####  INPUT FILES / ARGUMENTS  —  edit these, no need to scroll down  ####
// ####  (these are just the defaults for plot_ESE_scatter_and_pearson_  ####
// ####  combined()'s arguments, defined below; passing explicit args    ####
// ####  to the macro call still overrides them)                        ####
// #########################################################################

// D0 and charged-particle flow input files (must share the same q-bin scheme)
//static const char* INPUT_D0_FILE  = "D0_Flow_12NUQbin_diffq2q3_out_combined_Aug24.root"; //w/o Eff correction
//static const char* INPUT_CHG_FILE = "Charge_flow_Analysis_12NUQbin_diffq2q3_out_combined_Aug24.root"; //w/o Eff correction

static const char* INPUT_CHG_FILE = "Charge_Flow_wANDwo_Eff_combined_Sept25.root"; //contains both w/ & w/o Eff correction

//static const char* INPUT_D0_FILE  = "D0_Flow_woEff_out_combined_Sept29.root"; // w/o Eff correction
//static const char* INPUT_D0_FILE  = "D0_Flow_wEff_combined_Sept25.root"; //with Eff correction
static const char* INPUT_D0_FILE  = "D0_Flow_wEFF_combined_Sept29.root";

static const bool  INPUT_USE_EFFW = true;

static const std::string INPUT_TAG = "Sept29_wEff_v1";

// Output ROOT files for the scatter (Section 1) and Pearson (Section 2) results — derived from INPUT_TAG
static const std::string INPUT_SCATTER_FILE_OUT_STR = "ESE_scatter_" + INPUT_TAG + ".root";
static const std::string INPUT_PEARSON_FILE_OUT_STR = "pearson_ESE_" + INPUT_TAG + ".root";
static const char* INPUT_SCATTER_FILE_OUT = INPUT_SCATTER_FILE_OUT_STR.c_str();
static const char* INPUT_PEARSON_FILE_OUT = INPUT_PEARSON_FILE_OUT_STR.c_str();

// Output directory tag — derived from INPUT_TAG; everything from all sections is written under this directory
static const std::string INPUT_OUTDIR_TAG_STR = "ESE_12NUQbin_diffq2q3_" + INPUT_TAG;
static const char* INPUT_OUTDIR_TAG = INPUT_OUTDIR_TAG_STR.c_str();

// --- SECTION 0 (vn vs q-bin / vn vs mean-q) settings ---
// true -> charged particle 0.5 < pT < 3 GeV/c (Run3); false -> 1 < pT < 3 GeV/c (Run2)
static const bool  INPUT_IS_RUN3 = true;
static const bool  INPUT_SHOW_MEAN_ERR = true;   // x-axis error = mean error (mean-q pass only)
static const bool  INPUT_SHOW_STDDEV   = false;  // x-axis error = std-dev (overrides showMeanErr if true)


// For Inclusive centrality (0-50%)
static const char* INCL_CENT_TAG   = "cent0to50";
static const char* INCL_CENT_LABEL = "0-50%";
static const int   INCL_CENT_MIN   = 0;
static const int   INCL_CENT_MAX   = 50;
static const char* INCL_SUFFIX     = "_inclusiveCent";

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
// Helper: locate a Python interpreter that has PyMuPDF (fitz) installed.
// Tries common candidates in order; returns "" if none is found.
// =========================================================================
static TString FindPythonWithFitz()
{
    const char* cands[] = {
        "/usr/local/bin/python3.11",
        "/opt/homebrew/bin/python3",
        "/usr/bin/python3",
        "python3",
        "python",
        nullptr
    };
    for (int k = 0; cands[k]; ++k) {
        TString probe = Form("%s -c 'import fitz' 2>/dev/null", cands[k]);
        if (gSystem->Exec(probe.Data()) == 0) return TString(cands[k]);
    }
    return TString();
}

// =========================================================================
// SECTION 0 helper: vn vs q-bin-index and vn vs mean-q, D0 & charged-particle
// overlay  [ported from plot_vn_vs_qbin_D0Chg_combined.C]
//
// Fully self-contained: opens/closes its own handles on d0_file & chg_file,
// and writes everything under its own sub-directories of `outbase`
// (numbered 1_/2_ so they sort into pipeline order in any file browser):
//   {outbase}/1_vn_vs_qbin_plots/v{2,3}/*.pdf
//   {outbase}/1_vn_vs_qbin_plots/summary_v{2,3}_qbin_{tag}.pdf   (merged/tiled)
//   {outbase}/2_vn_vs_meanq_plots/v{2,3}/*.pdf
//   {outbase}/2_vn_vs_meanq_plots/summary_v{2,3}_meanq_{tag}.pdf (merged/tiled)
//   {outbase}/ESE_vn_vs_qbin_{tag}.root
//   {outbase}/ESE_vn_vs_meanq_{tag}.root
// =========================================================================
static void RunVnVsQbinD0Chg(const TString& outbase,
                              const char* d0_file,
                              const char* chg_file,
                              const char* tag,
                              bool        isRun3,
                              bool        showMeanErr,
                              bool        showStdDev,
                              bool        useEffW)
{
    const int N_QBINS = 12;
    const int N_PTBINS_V2 = 9;
    const double pt_edges_v2[N_PTBINS_V2 + 1] = { 2, 3, 4, 5, 6, 8, 10, 15, 30, 100 };
    const int N_PTBINS_V3 = 6;
    const double pt_edges_v3[N_PTBINS_V3 + 1] = { 2, 4, 6, 8, 10, 20, 50 };

    // Index 5 (INCL_CENT_TAG = "cent0to50") is the inclusive-centrality bin,
    // appended after the 5 differential ones — see the INCL_CENT_* comment
    // near the top of the file.
    const int N_CENTBINS_CHG = 6;
    const int    min_cent_chg[N_CENTBINS_CHG] = {  0, 10, 20, 30, 40, INCL_CENT_MIN };
    const int    max_cent_chg[N_CENTBINS_CHG] = { 10, 20, 30, 40, 50, INCL_CENT_MAX };
    const char* cent_lbl_chg[N_CENTBINS_CHG] = {
        "cent0to10","cent10to20","cent20to30","cent30to40","cent40to50", INCL_CENT_TAG
    };

    TString py_exe = FindPythonWithFitz();
    if (py_exe.IsNull())
        std::cerr << "WARNING: No Python interpreter with PyMuPDF found; "
                     "merged vn-vs-qbin/meanq PDFs will be skipped.\n";

    for (int pass = 0; pass < 2; ++pass)
    {
        bool doMeanQ = (pass == 1);
        TString typeTag = doMeanQ ? "meanq" : "qbin";

        std::cout << "\n--- vn-vs-" << typeTag << " pass ---\n";

        TFile* fD0  = TFile::Open(d0_file);
        TFile* fChg = TFile::Open(chg_file);
        if (!fD0  || fD0->IsZombie())  { std::cerr << "ERROR: Cannot open D0 file.\n"; return; }
        if (!fChg || fChg->IsZombie()) { std::cerr << "ERROR: Cannot open chg file.\n"; return; }

        TString d0_dir_name   = "vn_vs_qbin";
        TString chg_dir_name  = useEffW ? "vsQbin_TProfile_effW" : "vsQbin_TProfile";
        TString d0_qdist_dir  = "q_distributions";
        TString chg_qdist_dir = "Qdist_perBin";  // no _effW variant

        auto* dir_d0_qbin   = (TDirectory*)fD0->Get(d0_dir_name);
        auto* dir_chg_vsq   = (TDirectory*)fChg->Get(chg_dir_name);
        auto* dir_d0_qdist  = (TDirectory*)fD0->Get(d0_qdist_dir);
        auto* dir_chg_qdist = (TDirectory*)fChg->Get(chg_qdist_dir);

        if (!dir_d0_qbin) std::cerr << "WARNING: D0 " << d0_dir_name << " dir not found\n";
        if (!dir_chg_vsq) std::cerr << "WARNING: CHG " << chg_dir_name << " dir not found\n";
        if (doMeanQ && !dir_d0_qdist)  std::cerr << "WARNING: D0 " << d0_qdist_dir << " dir not found\n";
        if (doMeanQ && !dir_chg_qdist) std::cerr << "WARNING: CHG " << chg_qdist_dir << " dir not found\n";

        // Numeric prefix (1_, 2_) forces pipeline order in alphabetically-sorting
        // file browsers (Finder/VSCode/ls), which otherwise show folders A-Z.
        TString plot_base_dir = Form("%s/%d_vn_vs_%s_plots", outbase.Data(), pass + 1, typeTag.Data());
        gSystem->mkdir(plot_base_dir, kTRUE);
        gSystem->mkdir(plot_base_dir + "/v2", kTRUE);
        gSystem->mkdir(plot_base_dir + "/v3", kTRUE);

        TString out_file = Form("%s/ESE_vn_vs_%s_%s.root", outbase.Data(), typeTag.Data(), tag);
        TFile* fOut  = new TFile(out_file.Data(), "RECREATE");
        auto*  dir_v2 = fOut->mkdir("vn_vs_q_v2");
        auto*  dir_v3 = fOut->mkdir("vn_vs_q_v3");

        TLatex ltx;
        ltx.SetNDC();
        ltx.SetTextFont(42);

        const char* chg_pt_tag = isRun3 ? "pt0p5to3" : "pt1to3";
        TString chg_prof_base2 = "hp_v2_vsq2";
        TString chg_prof_base3 = "hp_v3_vsq3";

        TProfile* hp_chg_v2[N_CENTBINS_CHG] = {};
        TProfile* hp_chg_v3[N_CENTBINS_CHG] = {};

        TH1D* hq2_dist_d0[N_CENTBINS_CHG][N_QBINS]  = {};
        TH1D* hq3_dist_d0[N_CENTBINS_CHG][N_QBINS]  = {};
        TH1D* hq2_dist_chg[N_CENTBINS_CHG][N_QBINS] = {};
        TH1D* hq3_dist_chg[N_CENTBINS_CHG][N_QBINS] = {};

        for (int ic = 0; ic < N_CENTBINS_CHG; ++ic) {
            hp_chg_v2[ic] = dir_chg_vsq ? (TProfile*)dir_chg_vsq->Get(Form("%s_%s_%s", chg_prof_base2.Data(), chg_pt_tag, cent_lbl_chg[ic])) : nullptr;
            hp_chg_v3[ic] = dir_chg_vsq ? (TProfile*)dir_chg_vsq->Get(Form("%s_%s_%s", chg_prof_base3.Data(), chg_pt_tag, cent_lbl_chg[ic])) : nullptr;

            if (doMeanQ) {
                for (int iq = 0; iq < N_QBINS; ++iq) {
                    if (dir_d0_qdist) {
                        hq2_dist_d0[ic][iq] = (TH1D*)dir_d0_qdist->Get(Form("hq2_dist_%s_q2bin%d", cent_lbl_chg[ic], iq));
                        hq3_dist_d0[ic][iq] = (TH1D*)dir_d0_qdist->Get(Form("hq3_dist_%s_q3bin%d", cent_lbl_chg[ic], iq));
                    }
                    if (dir_chg_qdist) {
                        hq2_dist_chg[ic][iq] = (TH1D*)dir_chg_qdist->Get(Form("hq2_dist_cen%d_q2bin%d", ic, iq));
                        hq3_dist_chg[ic][iq] = (TH1D*)dir_chg_qdist->Get(Form("hq3_dist_cen%d_q3bin%d", ic, iq));
                    }
                }
            }
        }

        double cent_d0_ylo[2][N_CENTBINS_CHG],  cent_d0_yhi[2][N_CENTBINS_CHG];
        double cent_chg_ylo[2][N_CENTBINS_CHG], cent_chg_yhi[2][N_CENTBINS_CHG];

        for (int ivn = 2; ivn <= 3; ++ivn) {
            int vi             = ivn - 2;
            int n_pt_scan      = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
            const double* pt_e = (ivn == 2) ? pt_edges_v2 : pt_edges_v3;
            TProfile** hp_scan = (ivn == 2) ? hp_chg_v2   : hp_chg_v3;
            TString d0_gr_base = Form("v%d_vs_q%dbin", ivn, ivn);

            for (int ic = 0; ic < N_CENTBINS_CHG; ++ic) {
                double gd0_min =  1e9, gd0_max = -1e9;
                double gch_min =  1e9, gch_max = -1e9;

                for (int iq = 0; iq < N_QBINS; ++iq) {
                    if (hp_scan[ic] && hp_scan[ic]->GetBinEntries(iq+1) > 0) {
                        double cv = hp_scan[ic]->GetBinContent(iq+1);
                        double ce = hp_scan[ic]->GetBinError(iq+1);
                        if (ce <= 0) continue;
                        gch_min = std::min(gch_min, cv - ce);
                        gch_max = std::max(gch_max, cv + ce);
                    }
                }

                for (int ip = 0; ip < n_pt_scan; ++ip) {
                    if (ivn==2) {if (pt_e[ip] >= 30) continue;}
                    if (ivn==3) {if ( pt_e[ip] >= 20) continue;}
                    TString ptTag  = Form("pT%dto%d", (int)pt_e[ip], (int)pt_e[ip+1]);
                    TString grName = Form("%s_%s_%s", d0_gr_base.Data(), cent_lbl_chg[ic], ptTag.Data());
                    TGraphErrors* gr = dir_d0_qbin ? (TGraphErrors*)dir_d0_qbin->Get(grName) : nullptr;

                    if (!gr || gr->GetN() == 0) continue;
                    for (int iq = 0; iq < N_QBINS && iq < gr->GetN(); ++iq) {
                        double dv = gr->GetPointY(iq), de = gr->GetErrorY(iq);
                        double ce = (hp_scan[ic] && hp_scan[ic]->GetBinEntries(iq+1) > 0) ? hp_scan[ic]->GetBinError(iq+1) : 0.0;
                        if (de <= 0 || ce <= 0) continue;
                        gd0_min = std::min(gd0_min, dv - de);
                        gd0_max = std::max(gd0_max, dv + de);
                    }
                }
                if (gd0_min > gd0_max) { gd0_min = 0.0; gd0_max = 0.3; }
                if (gch_min > gch_max) { gch_min = 0.0; gch_max = 0.3; }
                double d0_marg  = 0.25 * std::max(gd0_max - gd0_min, 1e-6);
                double chg_marg = 0.25 * std::max(gch_max - gch_min, 1e-6);
                cent_d0_ylo[vi][ic]  = gd0_min  - d0_marg;
                cent_d0_yhi[vi][ic]  = gd0_max  + d0_marg;
                cent_chg_ylo[vi][ic] = std::max(0.0, gch_min - chg_marg);
                cent_chg_yhi[vi][ic] = gch_max  + chg_marg;

                if (ivn == 2) {
                    double unified_lo = std::min(cent_d0_ylo[vi][ic], cent_chg_ylo[vi][ic]);
                    double unified_hi = std::max(cent_d0_yhi[vi][ic], cent_chg_yhi[vi][ic]);
                    cent_d0_ylo[vi][ic]  = unified_lo;
                    cent_d0_yhi[vi][ic]  = unified_hi;
                    cent_chg_ylo[vi][ic] = unified_lo;
                    cent_chg_yhi[vi][ic] = unified_hi;
                }
            }
        }

        // --- Plot Generation ---
        for (int ivn = 2; ivn <= 3; ++ivn)
        {
            int           n_pt       = (ivn == 2) ? N_PTBINS_V2  : N_PTBINS_V3;
            const double* pt_edges   = (ivn == 2) ? pt_edges_v2  : pt_edges_v3;
            TProfile** hp_chg        = (ivn == 2) ? hp_chg_v2    : hp_chg_v3;
            TDirectory* dir_out      = (ivn == 2) ? dir_v2       : dir_v3;
            TString d0_gr_base       = Form("v%d_vs_q%dbin", ivn, ivn);

            for (int ic = 0; ic < N_CENTBINS_CHG; ++ic)
            {
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
                    TString grName = Form("%s_%s_%s", d0_gr_base.Data(), cent_lbl_chg[ic], ptTag.Data());

                    TGraphErrors* gr_d0 = dir_d0_qbin ? (TGraphErrors*)dir_d0_qbin->Get(grName) : nullptr;
                    if (!gr_d0 || gr_d0->GetN() == 0) continue;

                    double d0_q[N_QBINS]  = {}, d0_q_e[N_QBINS]  = {};
                    double chg_q[N_QBINS] = {}, chg_q_e[N_QBINS] = {};

                    double qx_d0[N_QBINS] = {}, qex_d0[N_QBINS] = {};
                    double qx_chg[N_QBINS] = {}, qex_chg[N_QBINS] = {};
                    int npts = 0;

                    for (int iq = 0; iq < N_QBINS && iq < gr_d0->GetN(); ++iq) {
                        double dv  = gr_d0->GetPointY(iq);
                        double de  = gr_d0->GetErrorY(iq);
                        double cv  = chg_vn_q[iq];
                        double ce  = chg_vn_q_e[iq];
                        if (de <= 0 || ce <= 0) continue;

                        double current_qx_d0, current_qex_d0;
                        double current_qx_chg, current_qex_chg;

                        if (doMeanQ) {
                            TH1D* hdist_d0 = (ivn == 2) ? hq2_dist_d0[ic][iq] : hq3_dist_d0[ic][iq];
                            if (hdist_d0 && hdist_d0->GetEntries() > 0) {
                                current_qx_d0  = hdist_d0->GetMean();
                                if (showStdDev)       current_qex_d0 = hdist_d0->GetStdDev();
                                else if (showMeanErr) current_qex_d0 = hdist_d0->GetMeanError();
                                else                  current_qex_d0 = 0.0;
                            } else {
                                current_qx_d0  = gr_d0->GetPointX(iq);
                                current_qex_d0 = 0.0;
                            }

                            TH1D* hdist_chg = (ivn == 2) ? hq2_dist_chg[ic][iq] : hq3_dist_chg[ic][iq];
                            if (hdist_chg && hdist_chg->GetEntries() > 0) {
                                current_qx_chg  = hdist_chg->GetMean();
                                if (showStdDev)       current_qex_chg = hdist_chg->GetStdDev();
                                else if (showMeanErr) current_qex_chg = hdist_chg->GetMeanError();
                                else                  current_qex_chg = 0.0;
                            } else {
                                current_qx_chg  = current_qx_d0;
                                current_qex_chg = 0.0;
                            }
                        } else {
                            current_qx_d0  = (double)iq;  current_qex_d0  = 0.0;
                            current_qx_chg = (double)iq;  current_qex_chg = 0.0;
                        }

                        d0_q[npts]    = dv;  d0_q_e[npts]   = de;
                        chg_q[npts]   = cv;  chg_q_e[npts]  = ce;

                        qx_d0[npts]   = current_qx_d0;   qex_d0[npts]  = current_qex_d0;
                        qx_chg[npts]  = current_qx_chg;  qex_chg[npts] = current_qex_chg;

                        ++npts;
                    }
                    if (npts < 2) continue;

                    double d0_ylo  = cent_d0_ylo[ivn - 2][ic];
                    double d0_yhi  = cent_d0_yhi[ivn - 2][ic];

                    double actual_xmin = std::min(*std::min_element(qx_d0, qx_d0 + npts), *std::min_element(qx_chg, qx_chg + npts));
                    double actual_xmax = std::max(*std::max_element(qx_d0, qx_d0 + npts), *std::max_element(qx_chg, qx_chg + npts));
                    double x_range = actual_xmax - actual_xmin;

                    double frame_xmin = doMeanQ ? std::max(0.0, actual_xmin - 0.15 * x_range) : -0.5;
                    double frame_xmax = doMeanQ ? (actual_xmax + 0.15 * x_range) : ((double)N_QBINS - 0.5);

                    TString base = Form("vn_vs_%s_v%d_%s_%s", typeTag.Data(), ivn, cent_lbl_chg[ic], ptTag.Data());

                    TGraphErrors* gr_d0_plot = new TGraphErrors(npts, qx_d0, d0_q, qex_d0, d0_q_e);
                    gr_d0_plot->SetName(Form("gr_d0_%s",  base.Data()));
                    gr_d0_plot->SetMarkerStyle(20);
                    gr_d0_plot->SetMarkerSize(1.4);
                    gr_d0_plot->SetMarkerColor(kBlue + 1);
                    gr_d0_plot->SetLineColor(kBlue + 1);
                    gr_d0_plot->SetLineWidth(2);

                    TGraphErrors* gr_chg_plot = new TGraphErrors(npts, qx_chg, chg_q, qex_chg, chg_q_e);
                    gr_chg_plot->SetName(Form("gr_chg_%s", base.Data()));
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

                    TH2F* hf = new TH2F(Form("hf_%s", base.Data()), "",
                                         10, frame_xmin, frame_xmax,
                                         10, d0_ylo, d0_yhi);
                    hf->GetXaxis()->SetTitle(doMeanQ ? Form("<q_{%d}>", ivn) : Form("q_{%d} bin index", ivn));
                    hf->GetYaxis()->SetTitle(ivn == 2 ? "v_{2}" : "v_{3}");
                    hf->GetXaxis()->SetTitleSize(0.060);
                    hf->GetYaxis()->SetTitleSize(0.060);
                    hf->GetXaxis()->SetTitleOffset(0.95);
                    hf->GetYaxis()->SetTitleOffset(0.85);
                    hf->GetXaxis()->SetLabelSize(0.038);
                    hf->GetYaxis()->SetLabelSize(0.038);
                    if (!doMeanQ) hf->GetXaxis()->SetNdivisions(10);
                    hf->GetYaxis()->SetTitleColor(kBlack);
                    hf->GetYaxis()->SetLabelColor(kBlack);
                    hf->Draw();

                    gr_d0_plot->Draw("P SAME");
                    gr_chg_plot->Draw("P SAME");

                    TLegend* leg = new TLegend(0.52, 0.78, 0.86, 0.90);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->SetTextFont(42);
                    leg->SetTextSize(0.038);
                    leg->AddEntry(gr_d0_plot,  Form("D^{0} v_{%d}", ivn), "pe");
                    leg->AddEntry(gr_chg_plot, Form("charged particle v_{%d}", ivn), "pe");
                    leg->Draw();

                    ltx.SetTextAlign(11);
                    ltx.SetTextFont(42);
                    ltx.SetTextSize(0.038);
                    ltx.DrawLatex(0.17, 0.86, Form("Centrality: %d-%d%%", min_cent_chg[ic], max_cent_chg[ic]));
                    ltx.DrawLatex(0.17, 0.79, Form("%.0f < p_{T}^{D^{0}} < %.0f GeV/c", pt_lo, pt_hi));

                    ltx.SetTextAlign(11);
                    ltx.SetTextFont(42);
                    ltx.SetTextSize(0.060);
                    ltx.DrawLatex(0.14, 0.945, "#bf{CMS} #it{Preliminary}");
                    ltx.SetTextAlign(31);
                    ltx.SetTextFont(42);
                    ltx.SetTextSize(0.052);
                    ltx.DrawLatex(0.86, 0.945, "PbPb 5.36 TeV");
                    ltx.SetTextAlign(11);

                    cv->SaveAs(Form("%s/v%d/%s.pdf", plot_base_dir.Data(), ivn, base.Data()));

                    dir_out->cd();
                    gr_d0_plot->Write();
                    gr_chg_plot->Write();
                    cv->Write();

                    delete leg;
                    delete hf;
                    delete gr_d0_plot;
                    delete gr_chg_plot;
                    delete cv;
                }
            }
        }

        fOut->Write();
        fOut->Close();
        fD0->Close();
        fChg->Close();

        if (!py_exe.IsNull()) {
            TString sum_v2 = Form("%s/summary_v2_%s_%s", plot_base_dir.Data(), typeTag.Data(), tag);
            TString sum_v3 = Form("%s/summary_v3_%s_%s", plot_base_dir.Data(), typeTag.Data(), tag);
            std::cout << "Tiling " << typeTag << " PDFs...\n";
            gSystem->Exec(Form("%s tile_vn_vs_qbin.py 2 %s %s", py_exe.Data(), plot_base_dir.Data(), sum_v2.Data()));
            gSystem->Exec(Form("%s tile_vn_vs_qbin.py 3 %s %s", py_exe.Data(), plot_base_dir.Data(), sum_v3.Data()));

            // Inclusive centrality (0-50%) — single-page summary, kept under a
            // clearly distinct name so it never mixes with the 5-page
            // differential-centrality summary above.
            TString sum_v2_incl = sum_v2 + INCL_SUFFIX;
            TString sum_v3_incl = sum_v3 + INCL_SUFFIX;
            std::cout << "Tiling " << typeTag << " PDFs (inclusive centrality)...\n";
            gSystem->Exec(Form("%s tile_vn_vs_qbin.py 2 %s %s %s", py_exe.Data(), plot_base_dir.Data(), sum_v2_incl.Data(), INCL_CENT_TAG));
            gSystem->Exec(Form("%s tile_vn_vs_qbin.py 3 %s %s %s", py_exe.Data(), plot_base_dir.Data(), sum_v3_incl.Data(), INCL_CENT_TAG));
        }
    } // pass

    std::cout << "vn-vs-qbin / vn-vs-meanq plots complete -> "
              << outbase << "/1_vn_vs_qbin_plots/, " << outbase << "/2_vn_vs_meanq_plots/\n";
}

// =========================================================================
// Main function
// =========================================================================
void plot_ESE_vnVsqn_scatter_and_pearson_combined(
    const char* d0_file          = INPUT_D0_FILE,
    const char* chg_file         = INPUT_CHG_FILE,
    const char* scatter_file_out = INPUT_SCATTER_FILE_OUT,
    const char* pearson_file_out = INPUT_PEARSON_FILE_OUT,
    const char* outdir_tag       = INPUT_OUTDIR_TAG,        // Output directory tag
    bool        isRun3           = INPUT_IS_RUN3,           // true -> pt0p5to3 (Run3); false -> pt1to3 (Run2)
    // --- SECTION 0 (vn vs q-bin / vn vs mean-q) settings ---
    const char* vn_qbin_tag      = INPUT_TAG.c_str(),       // tag used in ESE_vn_vs_{qbin,meanq}_<tag>.root / summary PDFs
    bool        showMeanErr      = INPUT_SHOW_MEAN_ERR,     // x-axis error = mean error (mean-q pass only)
    bool        showStdDev       = INPUT_SHOW_STDDEV,       // x-axis error = std-dev (overrides showMeanErr if true)
    bool        useEffW           = INPUT_USE_EFFW)         // true -> "*_effW" charged-particle dirs; false -> uncorrected dirs
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
    // Index 5 (INCL_CENT_TAG = "cent0to50") is the inclusive-centrality bin,
    // appended after the 5 differential ones — see the INCL_CENT_* comment
    // near the top of the file. It flows through every per-centrality loop
    // below automatically; only the differential-only combined/panel plots
    // (2-panel Pearson, 4-pad slope/intercept) deliberately ignore it, each
    // getting its own inclusive-only counterpart further down.
    const int N_CENTBINS = 6;
    const int N_QBINS    = 12;

    // Centrality labels — first 5 differential + 1 inclusive (0-50%)
    const char* cen_name[N_CENTBINS] = {
        "cent0to10","cent10to20","cent20to30","cent30to40","cent40to50", INCL_CENT_TAG
    };
    const int min_cent[N_CENTBINS] = { 0, 10, 20, 30, 40, INCL_CENT_MIN};
    const int max_cent[N_CENTBINS] = {10, 20, 30, 40, 50, INCL_CENT_MAX};

    // v2: 9 pT bins
    const int N_PTBINS_V2 = 9;
    const double pt_edges_v2[N_PTBINS_V2 + 1] = {
        2, 3, 4, 5, 6, 8, 10, 15, 30, 100
    };

    // v3: 6 pT bins
    const int N_PTBINS_V3 = 6;
    const double pt_edges_v3[N_PTBINS_V3 + 1] = {
        2, 4, 6, 8, 10, 20, 50
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

    // Per-centrality style (Pearson/slope/intercept plots) — solid circle for every
    // centrality, distinguished by color only. Last entry (kOrange+7) is the
    // inclusive-centrality bin.
    const int col        [N_CENTBINS] = {kBlack, kRed+1, kGreen+2, kBlue+1, kMagenta+1, kOrange+7};
    const int MARKER_STYLE = 20;
    const double MARKER_SIZE = 1.4;

    // =========================================================================
    // In-memory scatter graph storage  [ivn-2][ic][ip]
    // Populated in Section 1, consumed in Section 2.
    // =========================================================================
    TGraphErrors* gr_sc_mem[2][N_CENTBINS][N_PTBINS_MAX] = {};

    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    // =========================================================================
    // Output directory — use the tag directly as provided
    // =========================================================================
    TString outbase = outdir_tag;
    std::cout << "\n=== Output directory: " << outbase << " ===\n";

    // Sub-directory names, numbered 1_.. through 6_.. so any alphabetically-sorting
    // file browser (Finder/VSCode/ls) lists them in pipeline order. SECTION 0's own
    // 1_vn_vs_qbin_plots / 2_vn_vs_meanq_plots are built inside RunVnVsQbinD0Chg().
    const char* DIR_SCATTER         = "3_scatter_plots";
    const char* DIR_PEARSON         = "4_pearson_plots";
    const char* DIR_SLOPE_INTERCEPT = "5_slope_intercept_plots";
    const char* DIR_MC_R            = "6_mc_r_plots";

    // =========================================================================
    // SECTION 0: vn vs q-bin-index / vn vs mean-q  (D0 & charged particle overlay)
    // [ported from plot_vn_vs_qbin_D0Chg_combined.C]
    // Fully self-contained — runs to completion before SECTION 1 starts.
    // =========================================================================
    std::cout << "\n=== SECTION 0: vn vs q-bin / vn vs mean-q plots ===\n";
    RunVnVsQbinD0Chg(outbase, d0_file, chg_file, vn_qbin_tag, isRun3, showMeanErr, showStdDev, useEffW);

    // =========================================================================
    // Open input files
    // =========================================================================
    TFile* fD0  = TFile::Open(d0_file);
    TFile* fChg = TFile::Open(chg_file);
    if (!fD0  || fD0->IsZombie())  { std::cerr << "ERROR: Cannot open D0 file:  " << d0_file  << "\n"; return; }
    if (!fChg || fChg->IsZombie()) { std::cerr << "ERROR: Cannot open chg file: " << chg_file << "\n"; return; }

    TString chg_dir_name = useEffW ? "vsQbin_TProfile_effW" : "vsQbin_TProfile";
    auto* dir_d0_qbin = (TDirectory*)fD0->Get("vn_vs_qbin");
    auto* dir_chg_vsq = (TDirectory*)fChg->Get(chg_dir_name);
    if (!dir_d0_qbin) std::cerr << "WARNING: D0 vn_vs_qbin dir not found\n";
    if (!dir_chg_vsq) std::cerr << "WARNING: CHG " << chg_dir_name << " dir not found\n";

    // Resolved paths for the ROOT output files
    TString scatter_path = Form("%s/%s", outbase.Data(), gSystem->BaseName(scatter_file_out));
    TString pearson_path = Form("%s/%s", outbase.Data(), gSystem->BaseName(pearson_file_out));

    // =========================================================================
    // Create output directories
    // =========================================================================
    gSystem->mkdir(outbase,                                                   kTRUE);
    gSystem->mkdir(Form("%s/%s",    outbase.Data(), DIR_SCATTER),             kTRUE);
    gSystem->mkdir(Form("%s/%s/v2", outbase.Data(), DIR_SCATTER),             kTRUE);
    gSystem->mkdir(Form("%s/%s/v3", outbase.Data(), DIR_SCATTER),             kTRUE);
    gSystem->mkdir(Form("%s/%s",    outbase.Data(), DIR_PEARSON),             kTRUE);
    gSystem->mkdir(Form("%s/%s/v2", outbase.Data(), DIR_PEARSON),             kTRUE);
    gSystem->mkdir(Form("%s/%s/v3", outbase.Data(), DIR_PEARSON),             kTRUE);
    gSystem->mkdir(Form("%s/%s",    outbase.Data(), DIR_SLOPE_INTERCEPT),     kTRUE);
    gSystem->mkdir(Form("%s/%s",    outbase.Data(), DIR_MC_R),                kTRUE);
    gSystem->mkdir(Form("%s/%s/v2", outbase.Data(), DIR_MC_R),                kTRUE);
    gSystem->mkdir(Form("%s/%s/v3", outbase.Data(), DIR_MC_R),                kTRUE);

    // =========================================================================
    // Create output ROOT files
    // =========================================================================
    TFile* fSc  = new TFile(scatter_path, "RECREATE");
    auto*  sc_v2 = fSc->mkdir("scatter_v2");
    auto*  sc_v3 = fSc->mkdir("scatter_v3");

    TFile* fPr   = new TFile(pearson_path, "RECREATE");
    TDirectory* pr_v2 = fPr->mkdir("pearson_v2");
    TDirectory* pr_v3 = fPr->mkdir("pearson_v3");

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
                if (ivn==2)if (pt_e[ip] >= 30) continue;  // skip unreliable high-pT bin (30-100) from range scan
                if (ivn==3)if (pt_e[ip] >= 20) continue;  
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
                ltx.DrawLatex(0.92, 0.32,
                    Form("slope = %.3f #pm %.3f",
                         flin->GetParameter(1), flin->GetParError(1)));
                ltx.DrawLatex(0.92, 0.27,
                    Form("intercept = %.3f #pm %.3f",
                         flin->GetParameter(0), flin->GetParError(0)));
                ltx.DrawLatex(0.92, 0.22,
                    Form("#chi^{2}/NDF = %.2f",
                         flin->GetChisquare() / std::max(1.0, (double)flin->GetNDF())));
                ltx.SetTextAlign(11);
                ltx.SetTextFont(42);
                ltx.SetTextSize(0.048);
                ltx.DrawLatex(0.14, 0.945, "#bf{CMS} #it{Preliminary}");
                ltx.SetTextAlign(31);
                ltx.SetTextFont(42);
                ltx.SetTextSize(0.042);
                ltx.DrawLatex(0.95, 0.945, "PbPb 5.36 TeV");
                ltx.SetTextAlign(11);

                const char* subdir = (ivn == 2) ? "v2" : "v3";
                cv->SaveAs(Form("%s/%s/%s/%s.pdf", outbase.Data(), DIR_SCATTER, subdir, sc_name.Data()));

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
                "INCL_CENT_ORDER=['0to50']\n"
                "\n"
                "def pt_key(f):\n"
                "    return int(f.split('_pT')[1].split('to')[0])\n"
                "\n"
                "def merge(inp, out, cent_order, cols=3):\n"
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
                "    for cent in [c for c in cent_order if c in groups]:\n"
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
                "    if len(doc)==0:\n"
                "        print('No matching centralities for',out); doc.close(); return\n"
                "    doc.save(out,garbage=4,deflate=True)\n"
                "    doc.close()\n"
                "    print('Saved ->',out)\n"
                "\n");
            fprintf(fpy,
                "merge('%s/%s/v2','%s/%s/v2/merged_scatter_v2.pdf', CENT_ORDER)\n",
                outbase.Data(), DIR_SCATTER, outbase.Data(), DIR_SCATTER);
            fprintf(fpy,
                "merge('%s/%s/v3','%s/%s/v3/merged_scatter_v3.pdf', CENT_ORDER)\n",
                outbase.Data(), DIR_SCATTER, outbase.Data(), DIR_SCATTER);
            // Inclusive centrality (0-50%) — separate 1-page PDF, kept under a
            // clearly distinct name so it never mixes with the differential
            // merged_scatter_v{2,3}.pdf above.
            fprintf(fpy,
                "merge('%s/%s/v2','%s/%s/v2/merged_scatter_v2%s.pdf', INCL_CENT_ORDER)\n",
                outbase.Data(), DIR_SCATTER, outbase.Data(), DIR_SCATTER, INCL_SUFFIX);
            fprintf(fpy,
                "merge('%s/%s/v3','%s/%s/v3/merged_scatter_v3%s.pdf', INCL_CENT_ORDER)\n",
                outbase.Data(), DIR_SCATTER, outbase.Data(), DIR_SCATTER, INCL_SUFFIX);
            fclose(fpy);

            TString py_exe = FindPythonWithFitz();

            if (py_exe.IsNull()) {
                std::cerr << "WARNING: No Python interpreter with PyMuPDF found.\n"
                          << "  Install with:  pip install pymupdf\n"
                          << "  Merged PDFs skipped; individual PDFs are still in "
                          << outbase << "/" << DIR_SCATTER << "/\n";
            } else {
                std::cout << "  Using Python: " << py_exe << "\n";
                int ret = gSystem->Exec(Form("%s '%s'", py_exe.Data(), py_script.Data()));
                if (ret != 0)
                    std::cerr << "WARNING: Python merge script failed (exit " << ret
                              << "). Individual PDFs are still in "
                              << outbase << "/" << DIR_SCATTER << "/\n";
            }
            gSystem->Unlink(py_script.Data());
        }
    }

    std::cout << "  Individual scatter PDFs -> " << outbase << "/" << DIR_SCATTER << "/v{2,3}/\n"
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
                    // Clamp to the physically valid range of a Pearson r ([-1,1]) —
                    // the bisection search can converge on a half-width that would
                    // otherwise push r_meas ± hi_d outside that range (most visible
                    // for noisy/weak correlations at peripheral centrality / high pT).
                    r_elo_v = std::min(hi_d, r_meas + 1.0);   // ensures r_meas - r_elo_v >= -1
                    r_ehi_v = std::min(hi_d, 1.0 - r_meas);   // ensures r_meas + r_ehi_v <= 1
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
            gr_r[ic]->SetMarkerStyle(MARKER_STYLE);
            gr_r[ic]->SetMarkerSize(MARKER_SIZE);
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
                gr_slope[ic]->SetMarkerStyle(MARKER_STYLE);
                gr_slope[ic]->SetMarkerSize(MARKER_SIZE);
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
                gr_inter[ic]->SetMarkerStyle(MARKER_STYLE);
                gr_inter[ic]->SetMarkerSize(MARKER_SIZE);
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
            hf->GetXaxis()->SetTitleOffset(1.15);
            hf->GetXaxis()->SetLabelSize(0.050);
            hf->GetXaxis()->SetMoreLogLabels();
            hf->GetXaxis()->SetNoExponent();
            hf->GetXaxis()->CenterTitle(kTRUE);
            if (ipanel == 0) {
                hf->GetYaxis()->SetTitle("r");
                hf->GetYaxis()->SetTitleSize(0.072);
                hf->GetYaxis()->SetTitleOffset(1.05);
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
                ltx2.SetTextFont(42);
                ltx2.SetTextSize(0.065);
                ltx2.SetTextAlign(11);
                ltx2.DrawLatex(0.17, 0.930, "#bf{CMS} #it{Preliminary}");
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
        ltx_comb.SetTextSize(0.060);
        ltx_comb.SetTextAlign(31);
        ltx_comb.DrawLatex(0.975, 0.940, "PbPb 5.36 TeV");

        cv_comb->SaveAs(Form("%s/%s/v%d/pearson_v%d_vs_pT_combined.pdf", outbase.Data(), DIR_PEARSON, ivn, ivn));
        outDir->cd();
        cv_comb->Write();
        delete cv_comb;

        // =====================================================================
        // Slope + Intercept combined 4-pad canvas
        // Top row: slope vs pT;  Bottom row: intercept vs pT
        // Left column: 0-20%;   Right column: 20-50%
        // =====================================================================
        {
        const double SL_LO = (ivn == 2) ? -1.0 : -3.0;
        const double SL_HI = (ivn == 2) ?  4.0 :  6.0;
        double IN_LO = (ivn == 2) ? -0.5 : -1.0;
        double IN_HI = (ivn == 2) ?  0.5 :  1.0;
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

        padsi[0]->SetLeftMargin(0.16);  padsi[0]->SetRightMargin(0.000);
        padsi[0]->SetBottomMargin(0.00); padsi[0]->SetTopMargin(0.12);
        padsi[1]->SetLeftMargin(0.000); padsi[1]->SetRightMargin(0.05);
        padsi[1]->SetBottomMargin(0.00); padsi[1]->SetTopMargin(0.12);
        padsi[2]->SetLeftMargin(0.16);  padsi[2]->SetRightMargin(0.000);
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
                    hfsi->GetYaxis()->SetNdivisions(505);
                if (!is_top) {
                    hfsi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                    hfsi->GetXaxis()->SetTitleSize(0.090);
                    hfsi->GetXaxis()->SetTitleOffset(1.00);
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
                    hfsi->GetYaxis()->SetTitleOffset(0.82);
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
                        ltxsi.SetTextFont(42);
                        ltxsi.SetTextSize(0.095);
                        ltxsi.SetTextAlign(11);
                        ltxsi.DrawLatex(0.17, 0.915, "#bf{CMS} #it{Preliminary}");
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
        ltx_si_top.SetTextSize(0.045);
        ltx_si_top.SetTextAlign(31);
        ltx_si_top.DrawLatex(0.975, 0.955, "PbPb 5.36 TeV");

        cv_si->SaveAs(Form("%s/%s/slope_intercept_v%d_vs_pT.pdf", outbase.Data(), DIR_SLOPE_INTERCEPT, ivn));
        outDir->cd();
        cv_si->Write();
        delete cv_si;
        } // slope+intercept canvas scope

        // =====================================================================
        // Inclusive-centrality (0-50%) Slope + Intercept vs pT — single-panel
        // counterpart to the 4-pad differential-only canvas above (that one
        // is hardcoded to the 5-differential-centrality left/right panel
        // split via PANEL_IC/PANEL_NC and has no slot for a 6th curve).
        // Saved under a distinct filename so it never overwrites/duplicates
        // slope_intercept_v{ivn}_vs_pT.pdf above.
        // =====================================================================
        {
        const int inclIc = N_CENTBINS - 1; // INCL_CENT_TAG entry
        if (gr_slope[inclIc] || gr_inter[inclIc]) {
            const double SL_LO = -1.0, SL_HI = 4.0;
            const double IN_LO = (ivn == 2) ? -0.5 : -1.0;
            const double IN_HI = (ivn == 2) ?  0.5 :  1.0;
            const double y_lo_si[2]  = {SL_LO, IN_LO};
            const double y_hi_si[2]  = {SL_HI, IN_HI};
            const double y_ref_si[2] = {1.0, 0.0};
            const char*  y_tit_si[2] = {"Slope", "Intercept"};
            TGraphAsymmErrors* gr_cur_si[2] = { gr_slope[inclIc], gr_inter[inclIc] };

            TCanvas* cv_si_incl = new TCanvas(
                Form("cv_si_incl_v%d", ivn),
                Form("Slope & Intercept vs pT (v%d, inclusive cent)", ivn),
                1200, 500);
            cv_si_incl->Divide(2, 1, 0.015, 0.001);

            for (int icol = 0; icol < 2; ++icol) {
                cv_si_incl->cd(icol + 1);
                gPad->SetLeftMargin(0.16);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.09);
                gPad->SetLogx();

                TH2F* hfi = new TH2F(
                    Form("hfi_v%d_%d", ivn, icol), "",
                    100, X_LO, X_HI, 100, y_lo_si[icol], y_hi_si[icol]);
                hfi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                hfi->GetXaxis()->SetTitleSize(0.050);
                hfi->GetXaxis()->SetTitleOffset(1.15);
                hfi->GetXaxis()->SetLabelSize(0.042);
                hfi->GetXaxis()->SetMoreLogLabels();
                hfi->GetXaxis()->SetNoExponent();
                hfi->GetYaxis()->SetTitle(y_tit_si[icol]);
                hfi->GetYaxis()->SetTitleSize(0.050);
                hfi->GetYaxis()->SetTitleOffset(1.35);
                hfi->GetYaxis()->SetLabelSize(0.042);
                hfi->GetYaxis()->SetNdivisions(505);
                hfi->Draw("AXIS");

                TLine* lref_i = new TLine(X_LO, y_ref_si[icol], X_HI, y_ref_si[icol]);
                lref_i->SetLineStyle(2);
                lref_i->SetLineWidth(1);
                lref_i->SetLineColor(kGray+2);
                lref_i->Draw();

                if (gr_cur_si[icol]) gr_cur_si[icol]->Draw("P Z SAME");

                if (icol == 0) {
                    TLatex ltxi;
                    ltxi.SetNDC();
                    ltxi.SetTextFont(42);
                    ltxi.SetTextSize(0.050);
                    ltxi.SetTextAlign(11);
                    ltxi.DrawLatex(0.18, 0.935, "#bf{CMS} #it{Preliminary}");
                }
                TLatex ltxlab;
                ltxlab.SetNDC();
                ltxlab.SetTextFont(42);
                ltxlab.SetTextSize(0.045);
                ltxlab.SetTextAlign(11);
                ltxlab.DrawLatex(0.20, 0.86, Form("Centrality: %s", INCL_CENT_LABEL));
                ltxlab.DrawLatex(0.20, 0.80, "|y| < 1");
            }

            cv_si_incl->cd();
            TLatex ltx_si_incl_top;
            ltx_si_incl_top.SetNDC();
            ltx_si_incl_top.SetTextFont(42);
            ltx_si_incl_top.SetTextSize(0.035);
            ltx_si_incl_top.SetTextAlign(31);
            ltx_si_incl_top.DrawLatex(0.99, 0.965, "PbPb 5.36 TeV");

            cv_si_incl->SaveAs(Form("%s/%s/slope_intercept_v%d_vs_pT_%s.pdf",
                                     outbase.Data(), DIR_SLOPE_INTERCEPT, ivn, INCL_CENT_TAG));
            outDir->cd();
            cv_si_incl->Write();
            delete cv_si_incl;
        }
        } // inclusive slope+intercept canvas scope

        // =====================================================================
        // MC r-distribution combined canvases (4 x N_ROWS grid per centrality)
        // =====================================================================
        for (int ic = 0; ic < N_CENTBINS; ++ic)
        {
            bool any_mc = false;
            for (int ip = 0; ip < n_pt; ++ip)
                if (h_mc_arr[ic][ip]) { any_mc = true; break; }
            if (!any_mc) continue;

            // 3 columns: gives an exact 3x3 grid for v2 (9 pT bins) and 3x2 for v3 (6 pT bins)
            const int N_COLS_MC = 3;
            const int N_ROWS_MC = (n_pt + N_COLS_MC - 1) / N_COLS_MC;
            TCanvas* cv_mc = new TCanvas(
                Form("cv_mc_v%d_%s", ivn, cen_name[ic]),
                Form("MC r v%d %s", ivn, cen_name[ic]),
                1600, N_ROWS_MC * 400);
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
                hdraw->GetXaxis()->SetRangeUser(-1.0, 1.0);
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

            cv_mc->SaveAs(Form("%s/%s/v%d/mc_r_v%d_%s.pdf", outbase.Data(), DIR_MC_R, ivn, ivn, cen_name[ic]));
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
            hf2->GetXaxis()->SetTitleOffset(1.15);
            hf2->GetXaxis()->SetLabelSize(0.045);
            hf2->GetXaxis()->SetMoreLogLabels();
            hf2->GetXaxis()->SetNoExponent();
            hf2->GetYaxis()->SetTitle("r");
            hf2->GetYaxis()->SetTitleSize(0.055);
            hf2->GetYaxis()->SetTitleOffset(1.30);
            hf2->GetYaxis()->SetLabelSize(0.045);
            hf2->Draw("AXIS");

            TLine* lref2 = new TLine(X_LO, 1.0, X_HI, 1.0);
            lref2->SetLineStyle(2);
            lref2->SetLineWidth(1);
            lref2->SetLineColor(kGray+2);
            lref2->Draw();

            gr_r[ic]->Draw("P Z SAME");

            TLegend* leg2 = new TLegend(0.52, 0.78, 0.93, 0.90);
            leg2->SetBorderSize(0);
            leg2->SetFillStyle(0);
            leg2->SetTextFont(42);
            leg2->SetTextSize(0.045);
            leg2->AddEntry(gr_r[ic],
                Form("Centrality: %d-%d%%", min_cent[ic], max_cent[ic]), "p");
            leg2->AddEntry((TObject*)nullptr, "|y| < 1", "");
            leg2->Draw();

            TLatex ltx2;
            ltx2.SetNDC();
            ltx2.SetTextAlign(11);
            ltx2.SetTextFont(42);
            ltx2.SetTextSize(0.058);
            ltx2.DrawLatex(0.15, 0.935, "#bf{CMS} #it{Preliminary}");
            ltx2.SetTextAlign(31);
            ltx2.SetTextFont(42);
            ltx2.SetTextSize(0.052);
            ltx2.DrawLatex(0.94, 0.935, "PbPb 5.36 TeV");

            cv2->SaveAs(Form("%s/%s/v%d/pearson_v%d_%s.pdf",
                             outbase.Data(), DIR_PEARSON, ivn, ivn, cen_name[ic]));
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
    // Merged MC r PDFs — concatenate the already-saved per-centrality PDFs
    // (each already a full pT-bin grid for one centrality) into a single
    // multi-page PDF, one page per centrality, using Python/PyMuPDF.
    // Kept in the same directory as the individual plots.
    // =========================================================================
    std::cout << "=== Generating merged MC r PDFs (Python/PyMuPDF) ===\n";
    {
        TString py_script = Form("%s/_merge_mcr_tmp.py", outbase.Data());
        FILE* fpy = fopen(py_script.Data(), "w");
        if (!fpy) {
            std::cerr << "WARNING: Cannot write temp Python script; merged MC r PDFs skipped.\n";
        } else {
            fprintf(fpy,
                "import sys\n"
                "try:\n"
                "    import fitz\n"
                "except ModuleNotFoundError:\n"
                "    sys.exit('ERROR: PyMuPDF not found.  Install: pip install pymupdf')\n"
                "import os\n"
                "\n"
                "CENT_ORDER=['cent0to10','cent10to20','cent20to30','cent30to40','cent40to50']\n"
                "\n"
                "def concat(inp, out, prefix):\n"
                "    doc=fitz.open()\n"
                "    found=0\n"
                "    for cent in CENT_ORDER:\n"
                "        fname=os.path.join(inp, f'{prefix}_{cent}.pdf')\n"
                "        if not os.path.isfile(fname):\n"
                "            print(f'  MISSING: {fname}')\n"
                "            continue\n"
                "        with fitz.open(fname) as src:\n"
                "            doc.insert_pdf(src)\n"
                "        found+=1\n"
                "    if found==0:\n"
                "        print('No PDFs found in',inp); return\n"
                "    doc.save(out,garbage=4,deflate=True)\n"
                "    doc.close()\n"
                "    print(f'Saved {found}-page summary to',out)\n"
                "\n");
            fprintf(fpy,
                "concat('%s/%s/v2','%s/%s/v2/merged_mc_r_v2.pdf','mc_r_v2')\n",
                outbase.Data(), DIR_MC_R, outbase.Data(), DIR_MC_R);
            fprintf(fpy,
                "concat('%s/%s/v3','%s/%s/v3/merged_mc_r_v3.pdf','mc_r_v3')\n",
                outbase.Data(), DIR_MC_R, outbase.Data(), DIR_MC_R);
            fclose(fpy);

            TString py_exe = FindPythonWithFitz();

            if (py_exe.IsNull()) {
                std::cerr << "WARNING: No Python interpreter with PyMuPDF found.\n"
                          << "  Install with:  pip install pymupdf\n"
                          << "  Merged MC r PDFs skipped; individual PDFs are still in "
                          << outbase << "/" << DIR_MC_R << "/\n";
            } else {
                std::cout << "  Using Python: " << py_exe << "\n";
                int ret = gSystem->Exec(Form("%s '%s'", py_exe.Data(), py_script.Data()));
                if (ret != 0)
                    std::cerr << "WARNING: Python merge script failed (exit " << ret
                              << "). Individual PDFs are still in "
                              << outbase << "/" << DIR_MC_R << "/\n";
            }
            gSystem->Unlink(py_script.Data());
        }
    }

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
              << "  vn-vs-qbin PDFs         -> " << outbase << "/1_vn_vs_qbin_plots/v{2,3}/  (merged: summary_v{2,3}_qbin_" << vn_qbin_tag << ".pdf, inclusive cent: summary_v{2,3}_qbin_" << vn_qbin_tag << INCL_SUFFIX << ".pdf)\n"
              << "  vn-vs-meanq PDFs        -> " << outbase << "/2_vn_vs_meanq_plots/v{2,3}/  (merged: summary_v{2,3}_meanq_" << vn_qbin_tag << ".pdf, inclusive cent: summary_v{2,3}_meanq_" << vn_qbin_tag << INCL_SUFFIX << ".pdf)\n"
              << "  Individual scatter PDFs -> " << outbase << "/" << DIR_SCATTER << "/v{2,3}/  (inclusive cent tag: " << INCL_CENT_TAG << ")\n"
              << "  Merged scatter PDFs     -> " << outbase << "/" << DIR_SCATTER << "/v2/merged_scatter_v2.pdf\n"
              << "                         -> " << outbase << "/" << DIR_SCATTER << "/v3/merged_scatter_v3.pdf\n"
              << "                         -> " << outbase << "/" << DIR_SCATTER << "/v2/merged_scatter_v2" << INCL_SUFFIX << ".pdf  (inclusive cent)\n"
              << "                         -> " << outbase << "/" << DIR_SCATTER << "/v3/merged_scatter_v3" << INCL_SUFFIX << ".pdf  (inclusive cent)\n"
              << "  Pearson PDFs            -> " << outbase << "/" << DIR_PEARSON << "/v{2,3}/  (inclusive cent: pearson_v{2,3}_" << INCL_CENT_TAG << ".pdf)\n"
              << "  Slope/Int.              -> " << outbase << "/" << DIR_SLOPE_INTERCEPT << "/  (inclusive cent: slope_intercept_v{2,3}_vs_pT_" << INCL_CENT_TAG << ".pdf)\n"
              << "  MC r plots              -> " << outbase << "/" << DIR_MC_R << "/v{2,3}/  (inclusive cent: mc_r_v{2,3}_" << INCL_CENT_TAG << ".pdf)\n"
              << "  Merged MC r PDFs        -> " << outbase << "/" << DIR_MC_R << "/v2/merged_mc_r_v2.pdf\n"
              << "                         -> " << outbase << "/" << DIR_MC_R << "/v3/merged_mc_r_v3.pdf\n"
              << "  Scatter ROOT            -> " << scatter_path << "\n"
              << "  Pearson ROOT            -> " << pearson_path << "\n"
              << "  vn-vs-qbin ROOT         -> " << outbase << "/ESE_vn_vs_qbin_" << vn_qbin_tag << ".root\n"
              << "  vn-vs-meanq ROOT        -> " << outbase << "/ESE_vn_vs_meanq_" << vn_qbin_tag << ".root\n";
}


