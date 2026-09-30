//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// plot_ESE_systematics.C
//
// Systematic-uncertainty study for the Pearson r / slope / intercept results computed by
// plot_ESE_scatter_and_pearson_combined.C.
//
// ALL SETTINGS YOU MIGHT WANT TO CHANGE — input files, source labels, the 6 on/off switches
// below, thresholds, MC settings — live in ONE block ("USER CONFIGURATION") at the very top
// of this file, right after the #include lines. Edit that block, then just run:
//     root -l -b -q 'plot_ESE_systematics.C()'
// (no arguments needed — nothing below the user-configuration block should need editing).
//
// INPUT: for the nominal result AND for each of 7 systematic sources, the SAME PAIR of raw
// analysis files that plot_ESE_scatter_and_pearson_combined.C itself takes (d0_file, chg_file).
// Each source has EITHER one variation (e.g. a single alternative cut) OR two variations
// (e.g. "Loose"/"Tight" or "Down"/"Up") to compare against nominal — see source_d0_file[][2]
// in the user-configuration block, where the second slot is left nullptr for a one-variation
// source. For a two-variation source, at every (harmonic, centrality, pT) bin this macro
// compares BOTH variations to nominal and adopts whichever one deviates more
// (max(|variantA-nominal|, |variantB-nominal|), sign kept) as that source's systematic at
// that bin — the dominant side is allowed to differ from bin to bin. This macro re-derives
// the D0-vn-vs-charged-vn scatter per q-bin, the weighted Pearson r (with the same
// MC-resampling uncertainty as the nominal macro) and the linear-fit slope/intercept,
// independently for each raw-file variation — it does not need any pre-computed
// pearson_ESE_*.root file. It does NOT reproduce the per-bin scatter/Pearson diagnostic PDFs
// that plot_ESE_scatter_and_pearson_combined.C makes for a single dataset; it only needs
// (and keeps in memory) the final r / slope / intercept vs pT graphs.
//
// WHERE TO PUT YOUR FILES:
//   Put the nominal (D0_Flow_*.root, Charge_flow_*.root) pair and each source's 1-or-2
//   variation pairs into the "syst_inputs/" subdirectory (created next to this macro), then
//   edit the user-configuration block at the top of this file so nominal_d0_file /
//   nominal_chg_file / source_d0_file[][2] / source_chg_file[][2] / source_var_label[][2] /
//   source_label[] point at your actual files and carry the actual systematic-source names
//   (source_d0_file[isrc][1] left nullptr means that source has only one variation).
//
// TEN INDEPENDENT CALCULATIONS (each gated by its own on/off switch in the user-configuration
// block so you can run/review them one at a time):
//
// IMPORTANT: Pearson r's systematic uncertainty is NOT combined from the 7 sources by a
// per-source quadratic sum, unlike Slope/Intercept. Since r is bounded to [-1,1] by
// definition, an envelope built from adding per-source r-differences in quadrature is not a
// valid way to propagate its uncertainty. Instead (STEP 8, following the same procedure as
// the Run 2 analysis) each D0 v_n value feeding the Pearson-r calculation is Gaussian-smeared
// with sigma = the total per-q-bin systematic uncertainty on that D0 v_n (quadratic sum across
// sources, from STEP 7's ComputeScatterTotalSyst), while the weights used in the r calculation
// stay the STATISTICAL uncertainties — exactly the same MC-resampling machinery already used
// for r's statistical uncertainty (see WeightedPearson / ComputeDatasetGraphs, ported from
// plot_ESE_vnVsqn_scatter_and_pearson_combined.C), just with the smearing sigma swapped from
// stat to total-systematic. Because of this, Pearson r is EXCLUDED from STEPs 2/3/4/5/6 below
// (which all still use the per-source quadratic-sum method, valid for Slope/Intercept) — its
// own systematic result lives only in STEP 8's separate pearson_syst_plots/ directory and text
// table, kept apart from Slope/Intercept's totalunc_plots/ so the two methods are never mixed
// in the same output.
//
//   STEP 1 (DO_STEP1_OVERLAY) — nominal vs. each systematic source, overlaid (still includes r
//       as a raw diagnostic comparison; this step computes no uncertainty, so the r/slope/
//       intercept methodology split above doesn't apply here).
//       Per source: 4 combined canvases (pearson_v2, pearson_v3, slope+intercept_v2,
//       slope+intercept_v3), each already covering all 5 centrality classes in sub-panels,
//       exactly like the combined canvases in the nominal macro.  4 x 7 = 28 PDFs.
//   STEP 2 (DO_STEP2_ABSDIFF) — Slope/Intercept only: absolute difference (source - nominal)
//       vs pT.  2-canvas-per-source layout (slope+intercept_v2, slope+intercept_v3).  2 x 7 = 14 PDFs.
//   STEP 3 (DO_STEP3_RELDIFF) — Slope/Intercept only: relative difference
//       (source-nominal)/nominal x 100% vs pT.  Same 2-canvas-per-source layout; points with
//       |nominal| below a threshold are omitted (kept in STEP 2's absolute-difference plots).
//       2 x 7 = 14 PDFs.
//   STEP 4 (DO_STEP4_SUMMARY) — Slope/Intercept only: one canvas per (harmonic x quantity) =
//       4 canvases, each with one sub-pad per centrality, overlaying all 7 sources'
//       relative-difference curves vs pT plus a shaded min/max envelope band per pT bin.
//   STEP 5 (DO_STEP5_TOTALUNC) — Slope/Intercept only: vs pT with statistical error bars
//       (from the nominal file) and a total-systematic box (quadratic sum of the 7 source
//       differences, per pT and centrality bin).  2 PDFs (slope+intercept_v2, slope+intercept_v3).
//   STEP 6 (DO_STEP6_TEXTDUMP) — Slope/Intercept only: writes the full numeric breakdown
//       (nominal, each source's absolute/relative shift, and the total quadrature-sum
//       systematic) to plain-text tables in text_output/, for pasting into an analysis note.
//       Each quantity's table is followed by a min/max-systematic-source-per-bin summary.
//   STEP 7 (DO_STEP7_SCATTER_SYST) — the nominal D0-vn-vs-charged-vn scatter plots (one per
//       harmonic/centrality/pT bin, same layout as plot_ESE_scatter_and_pearson_combined.C's
//       own scatter plots: linear fit + 1-sigma confidence band, legend, CMS headers), with
//       an added per-q-bin total-systematic box drawn behind each statistical point: the box's
//       y-half-height is the real propagated D0-v_n systematic (quadratic sum across sources);
//       its x-half-width is a FIXED fraction of the point's charged-vn value (SYST_BOX_XFRAC),
//       not a separately computed charged-vn systematic. Also computes (regardless of this
//       switch) the per-q-bin total-systematic values STEP 8/9/10 all need.
//   STEP 8 (DO_STEP8_PEARSON_SYST) — Pearson r's own systematic uncertainty via Gaussian MC
//       resampling (see the IMPORTANT note above), in its own pearson_syst_plots/ directory
//       and text_output/pearson_systematic_v{2,3}.txt table (Nominal, StatLo/Hi, SystLo/Hi;
//       no per-source breakdown, since the per-source method doesn't apply here).  2 PDFs.
//       Also saves the underlying r_mc distribution histograms (one sub-pad per pT bin, one
//       PDF per centrality, same layout as the nominal macro's own mc_r_plots) to
//       mc_r_syst_plots/v{2,3}/, plus a merged_mc_r_syst_v{2,3}.pdf per harmonic (concatenated
//       via Python/PyMuPDF, same technique as the nominal macro's merged mc_r_plots).
//   STEP 9 (DO_STEP9_VN_VS_QBIN_SYST) — D0 v_n vs q-bin index (0-11), one combined canvas per
//       (harmonic, centrality) with a sub-pad per pT bin (3-column grid, same layout as STEP 8's
//       mc_r_syst canvases), nominal points + stat error bars, and a total-systematic box per
//       point: y-half-height = the real propagated systematic (reused from STEP 7's
//       ComputeScatterTotalSyst, .sysy); x-half-width = SYST_BOX_FIXED_HALFWIDTH (a single
//       hand-picked absolute number, since a q-bin index has no physical systematic of its
//       own). Saved to vn_vs_qbin_syst_plots/v{2,3}/, plus a merged PDF per harmonic
//       (Python/PyMuPDF concat).
//   STEP 10 (DO_STEP10_VN_VS_MEANQ_SYST) — plots the exact same points at the exact same
//       q-bin-index x-positions as STEP 9 (same box, same SYST_BOX_FIXED_HALFWIDTH); the only
//       difference is that each tick is relabeled with that bin's actual mean-q value (read
//       once from the nominal D0 file's "vn_vs_qmean" directory -- no source-variant files
//       need to be re-read for this) instead of its plain index. Saved to
//       vn_vs_meanq_syst_plots/v{2,3}/, plus a merged PDF per harmonic.
//
// Usage (after editing the user-configuration block at the top of this file):
//   root -l -b -q 'plot_ESE_systematics.C()'
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
#include <cstdio>
#include <random>

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TBox.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TString.h"

// #########################################################################
// ####  USER CONFIGURATION — EDIT THIS SECTION, THEN JUST RUN:          ####
// ####    root -l -b -q 'plot_ESE_systematics.C()'                      ####
// ####  Everything below this section is internal and should not need  ####
// ####  to change for a normal run.                                    ####
// #########################################################################

// Number of systematic sources. Only change this if you are actually adding/removing a
// source slot — if you do, also resize g_src_col/g_src_marker further below.
static const int N_SOURCES = 7;

// --- Input files ---------------------------------------------------------
static const char* nominal_d0_file  = "syst_inputs/D0_Flow_12NUQbin_diffq2q3_out_combined_Aug24.root";
static const char* nominal_chg_file = "syst_inputs/Charge_flow_Analysis_12NUQbin_diffq2q3_out_combined_Aug24.root";

// Put each source's (D0 file, charged-particle file) pair(s) into syst_inputs/ (created
// next to this macro). Every file must have the exact same layout as the nominal pair
// (same directory names / q-bin,pT,centrality tags) that
// plot_ESE_scatter_and_pearson_combined.C reads.
//
// Each source has a SECOND slot for sources evaluated with two variations (e.g.
// "Loose"/"Tight" or "Down"/"Up"): fill BOTH [0] and [1] for such a source, and leave [1]
// as nullptr for a source with only one variation. When both are filled, this macro
// compares each to nominal per (harmonic, centrality, pT) bin and adopts whichever
// deviates more as that source's systematic at that bin.
static const char* source_d0_file[N_SOURCES][2] = {
    { "syst_inputs/D0_Flow_SigPDF_combined_Sept11.root",       nullptr },
    { "syst_inputs/D0_Flow_out_BkgPDF_combined_Sept11.root",       nullptr },
    { "syst_inputs/D0_Flow_BDTloose_combined_Sept11.root",  "syst_inputs/D0_Flow_BDTtight_combined_Sept11.root" },
    { "syst_inputs/D0_Flow_SysCentDown_combined_Sept11.root",  "syst_inputs/D0_Flow_SysCentUp_combined_Sept11.root" },
    { "syst_inputs/D0_Flow_qnDOWN_combined.root",  "syst_inputs/D0_Flow_qnUP_combined.root"},
    { nullptr,       nullptr },
    { nullptr,       nullptr },
};
static const char* source_chg_file[N_SOURCES][2] = {
    { "syst_inputs/Charge_flow_Analysis_12NUQbin_diffq2q3_out_combined_Aug24.root",      nullptr },
    { "syst_inputs/Charge_flow_Analysis_12NUQbin_diffq2q3_out_combined_Aug24.root",      nullptr },
    { "syst_inputs/Charge_flow_Analysis_12NUQbin_diffq2q3_out_combined_Aug24.root", "syst_inputs/Charge_flow_Analysis_12NUQbin_diffq2q3_out_combined_Aug24.root" }, 
    { "syst_inputs/flow_Analysis_chg_centrailityDOWN_combined.root", "syst_inputs/flow_Analysis_chg_centrailityUP_combined.root" },
    { "syst_inputs/flow_Analysis_chg_qnDOWN_combined.root", "syst_inputs/flow_Analysis_chg_qnUP_combined.root" },
    { nullptr,      nullptr },
    { nullptr,      nullptr },
};
// Sub-labels for the plot annotations of two-variation sources (ignored — and may stay
// {nullptr, nullptr} — for a one-variation source). E.g. {"Loose","Tight"} or {"Down","Up"}.
static const char* source_var_label[N_SOURCES][2] = {
    { nullptr, nullptr },
    { nullptr, nullptr },
    { "Loose", "Tight" },
    { "Down", "Up" },
    { "Down", "Up" },
    { "Down", "Up" },
    { "Down", "Up" },
};
// Each source's own display name — used in filenames, legends and the text-dump columns.
static const char* source_label[N_SOURCES] = {
    "SigPDF", "BkgPDF", "BDT", "Centrality", "qn quantile", "PromptFraction", "Syst7"
};

// A source marked true here gets its systematic from a precomputed per-(centrality, pT) table
// instead of rerunning the analysis on varied files — its source_d0_file[]/source_chg_file[]
// entries above are ignored and should stay {nullptr, nullptr}. Currently only the prompt-
// fraction source (see prompt_frac_shift_file below) uses this. Treated as an intercept-only
// shift (Delta r = 0, Delta Slope = 0, Delta Intercept = the table's vn_shift), with a
// separately measured shift value per harmonic (see LoadPromptFractionShiftTable).
static const bool source_is_table_shift[N_SOURCES] = {
    false, false, false, false, false, true, false
};
// Prompt-fraction shift table produced by compute_promptfrac_vn_shift.C: vn_shift = factor *
// (v_n,1 - v_n,2), with (v_n,1 - v_n,2) measured directly from data (nominal vs. loose-DCA
// "SysPromFrac" file) rather than assumed. v2 and v3 rows each use that harmonic's own native
// pT binning (g_pt_edges_v2 / g_pt_edges_v3) — no cross-harmonic pT-bin borrowing needed.
static const char* prompt_frac_shift_file = "PromptFraction_vn_factor_measured.txt";

// --- Output ----------------------------------------------------------------
// Output directory tag. Leave "" to auto-generate "ESE_systematics_YYYYMMDD".
static const char* outdir_tag = "";

// --- Which of the 8 calculations to run (toggle independently) -------------
// NOTE: STEPs 2-6 only ever produce Slope/Intercept output now -- Pearson r's systematic
// uncertainty uses a different method (Gaussian MC resampling, not a per-source quadratic
// sum) and lives entirely in STEP 8; see the IMPORTANT note in the header comment above.
static const bool DO_STEP1_OVERLAY   = true;  // nominal vs each source, overlaid (r + slope/intercept)
static const bool DO_STEP2_ABSDIFF   = true;  // absolute difference vs pT (slope/intercept only)
static const bool DO_STEP3_RELDIFF   = true;  // relative difference (%) vs pT (slope/intercept only)
static const bool DO_STEP4_SUMMARY   = true;  // all-sources summary + envelope (slope/intercept only)
static const bool DO_STEP5_TOTALUNC  = true;  // total systematic (quadratic sum) box plot (slope/intercept only)
static const bool DO_STEP6_TEXTDUMP  = true;  // numeric table for the analysis note (slope/intercept only)
static const bool DO_STEP7_SCATTER_SYST = true;  // D0-vn-vs-charged-vn scatter plots with a total-syst box per point
static const bool DO_STEP8_PEARSON_SYST = true;  // Pearson r systematic uncertainty via Gaussian MC resampling
static const bool DO_STEP9_VN_VS_QBIN_SYST  = true;  // D0 v_n vs q-bin index, with a total-syst box per point
static const bool DO_STEP10_VN_VS_MEANQ_SYST = true;  // D0 v_n vs mean-q value, with a total-syst box per point

// Fixed (not data-derived) systematic-box x-half-widths, used wherever the x-axis has no
// independent "systematic uncertainty" of its own to propagate (STEP 7's charged-vn axis;
// STEP 9's q-bin-index axis; STEP 10's mean-q axis). The box's y-half-height, by contrast, is
// always the real propagated D0-v_n systematic (quadratic sum across sources).
static const double SYST_BOX_XFRAC = 0.02; // STEP 7: fixed 2% of the point's own charged-vn value
// STEP 9 and STEP 10 both plot v_n at the same q-bin-index x-positions (0..N_QBINS-1) -- STEP
// 10 only relabels those same tick positions with each bin's mean-q value, it does not use a
// differently-scaled x-axis. So one plain, hand-picked absolute half-width (in q-bin-index
// units) renders identically in both.
static const double SYST_BOX_FIXED_HALFWIDTH = 0.20;

// --- Misc analysis settings --------------------------------------------------
static const bool isRun3 = true;  // true -> pt0p5to3 (Run3); false -> pt1to3 (Run2)

// Relative-difference points are omitted when |nominal| falls below these thresholds
// (kept in the absolute-difference plots / text dump regardless).
static const double REL_ABS_THRESHOLD_R     = 0.02;
static const double REL_ABS_THRESHOLD_SLOPE = 0.05;
static const double REL_ABS_THRESHOLD_INTER = 0.02;

// Pearson-r MC-resampling uncertainty settings — same meaning as in
// plot_ESE_scatter_and_pearson_combined.C. N_MC is applied independently to every raw file
// loaded, i.e. (1 nominal + sum of each source's 1-or-2 variations) times, so with 4
// two-variation sources that's 1+3x1+4x2 = 12 (vs 8 with none), ~50% more run time. Lower
// this while iterating on the macro itself, then set back to 50000 for the final result.
static const int    N_MC                   = 50000;
static const bool   USE_SYMMETRIC_INTERVAL = true;

// #########################################################################
// ####  END OF USER CONFIGURATION                                       ####
// #########################################################################

// =========================================================================
// Bin definitions — MUST match plot_ESE_scatter_and_pearson_combined.C
// =========================================================================
static const int N_CENTBINS = 5;
static const char* g_cen_name[N_CENTBINS] = {
    "cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"
};
static const int g_min_cent[N_CENTBINS] = { 0, 10, 20, 30, 40};
static const int g_max_cent[N_CENTBINS] = {10, 20, 30, 40, 50};

static const int N_PTBINS_V2 = 9;
static const double g_pt_edges_v2[N_PTBINS_V2 + 1] = { 2, 3, 4, 5, 6, 8, 10, 15, 30, 100 };
static const int N_PTBINS_V3 = 6;
static const double g_pt_edges_v3[N_PTBINS_V3 + 1] = { 2, 4, 6, 8, 10, 20, 50 };
static const int N_PTBINS_MAX = 9;

// Number of q-bins the ESE event-shape classification splits each event into per q2/q3 —
// MUST match plot_ESE_scatter_and_pearson_combined.C (used for the STEP 7 scatter plots).
static const int N_QBINS = 12;

// Per-centrality style (filled = nominal, open = systematic-source overlay)
static const int g_col[N_CENTBINS]        = {kBlack, kRed+1, kGreen+2, kBlue+1, kMagenta+1};
static const int g_marker[N_CENTBINS]     = {20,     21,     22,       23,      33};
static const int g_marker_open[N_CENTBINS]= {24,     25,     26,       32,      27};

// Per-source plotting style (summary/diff plots). N_SOURCES is defined in the user
// configuration section above.
static const int g_src_col[N_SOURCES]    = {kOrange+7, kAzure+2, kSpring-1, kViolet+1, kCyan+2, kPink+6, kGray+2};
static const int g_src_marker[N_SOURCES] = {20, 21, 22, 23, 33, 34, 29};

// Panel grouping used throughout (matches the nominal macro's combined canvases)
static const int PANEL_IC[2][3] = { {0, 1, -1}, {2, 3, 4} };
static const int PANEL_NC[2]    = { 2, 3 };
static const double X_LO = 1.0;
static const double X_HI = 110.0;

// =========================================================================
// Small helpers
// =========================================================================
static void ComputePtCenters(const double* edges, int n, double* cen, double* elo, double* ehi)
{
    for (int i = 0; i < n; ++i) {
        cen[i] = std::sqrt(edges[i] * edges[i + 1]);
        elo[i] = cen[i] - edges[i];
        ehi[i] = edges[i + 1] - cen[i];
    }
}

static int FindPointIndex(TGraphAsymmErrors* g, double xtarget, double tol = 1e-4)
{
    if (!g) return -1;
    for (int i = 0; i < g->GetN(); ++i)
        if (std::fabs(g->GetPointX(i) - xtarget) < tol) return i;
    return -1;
}

struct PtVal { double y = 0.0, elo = 0.0, ehi = 0.0; bool ok = false; };

// One q-bin point of the D0-vn-vs-charged-vn scatter (STEP 7): x = charged vn, y = D0 vn,
// with their statistical errors, exactly as plot_ESE_scatter_and_pearson_combined.C builds
// the scatter that feeds the Pearson r / slope / intercept fit.
struct ScatterPoint { double x = 0.0, y = 0.0, ex = 0.0, ey = 0.0; bool ok = false; };

// Total systematic (quadratic sum across sources) on one scatter point's x and y, separately.
struct ScatterSyst { double sysx = 0.0, sysy = 0.0; bool ok = false; };

// Returns lbl, or fallback if lbl is null — guards against a source's source_var_label[][]
// entry being left {nullptr, nullptr} in the user configuration while its source_d0_file[][]
// entry has a real second file (making it a two-variation source): printing/using a null
// const char* directly (e.g. via std::cout <<) crashes, so every read of source_var_label
// goes through this instead of being used raw.
static const char* SafeLabel(const char* lbl, const char* fallback)
{
    return lbl ? lbl : fallback;
}

static PtVal GetPtVal(TGraphAsymmErrors* g, double xtarget, double tol = 1e-4)
{
    PtVal p;
    int idx = FindPointIndex(g, xtarget, tol);
    if (idx < 0) return p;
    p.y   = g->GetPointY(idx);
    p.elo = g->GetErrorYlow(idx);
    p.ehi = g->GetErrorYhigh(idx);
    p.ok  = true;
    return p;
}

// Build a (source - nominal) difference graph, optionally as a percentage relative
// difference. Points where |nominal| < rel_abs_threshold are dropped when relative=true.
static TGraphAsymmErrors* MakeDiffGraph(
    TGraphAsymmErrors* g_src, TGraphAsymmErrors* g_nom,
    const double* pt_cen, const double* pt_elo, const double* pt_ehi, int n_pt,
    bool relative, double rel_abs_threshold,
    int color, int marker, int linestyle = 1)
{
    std::vector<double> xv, yv, exlo, exhi, eylo, eyhi;
    for (int ip = 0; ip < n_pt; ++ip) {
        PtVal pn = GetPtVal(g_nom, pt_cen[ip]);
        PtVal ps = GetPtVal(g_src, pt_cen[ip]);
        if (!pn.ok || !ps.ok) continue;
        if (relative && std::fabs(pn.y) < rel_abs_threshold) continue;
        double d = ps.y - pn.y;
        if (relative) d = d / pn.y * 100.0;
        xv.push_back(pt_cen[ip]);   yv.push_back(d);
        exlo.push_back(pt_elo[ip]); exhi.push_back(pt_ehi[ip]);
        eylo.push_back(0.0);        eyhi.push_back(0.0);
    }
    if (xv.empty()) return nullptr;
    auto* g = new TGraphAsymmErrors(
        (int)xv.size(), xv.data(), yv.data(), exlo.data(), exhi.data(), eylo.data(), eyhi.data());
    g->SetMarkerColor(color); g->SetLineColor(color);
    g->SetMarkerStyle(marker); g->SetLineStyle(linestyle);
    g->SetMarkerSize(1.0); g->SetLineWidth(2);
    return g;
}

// For a source with one or two raw-file variations, builds the "adopted" value-space graph:
// at each pT bin, if both variants have a point, keep whichever deviates more from nominal
// (max(|value-nominal|), sign/value kept, not just the delta); if only one variant has a
// point there, use it as-is. g_varB may be nullptr (or all-null-per-bin) for a
// one-variation source, in which case this reduces exactly to "use g_varA". Y-errors on the
// result are set to 0 — this synthesized graph is only ever read via GetPtVal(...).y
// downstream (MakeDiffGraph / MakeEnvelopeBand / ComputeTotalSyst), never drawn directly.
static TGraphAsymmErrors* BuildEffectiveGraph(
    TGraphAsymmErrors* g_nom, TGraphAsymmErrors* g_varA, TGraphAsymmErrors* g_varB,
    const double* pt_cen, const double* pt_elo, const double* pt_ehi, int n_pt)
{
    std::vector<double> xv, yv, exlo, exhi, eylo, eyhi;
    for (int ip = 0; ip < n_pt; ++ip) {
        PtVal pn = GetPtVal(g_nom, pt_cen[ip]);
        if (!pn.ok) continue;
        PtVal pa = GetPtVal(g_varA, pt_cen[ip]);
        PtVal pb = GetPtVal(g_varB, pt_cen[ip]);

        double eff_y; bool any = false;
        if (pa.ok && pb.ok) {
            eff_y = (std::fabs(pa.y - pn.y) >= std::fabs(pb.y - pn.y)) ? pa.y : pb.y;
            any = true;
        } else if (pa.ok) {
            eff_y = pa.y; any = true;
        } else if (pb.ok) {
            eff_y = pb.y; any = true;
        }
        if (!any) continue;

        xv.push_back(pt_cen[ip]);   yv.push_back(eff_y);
        exlo.push_back(pt_elo[ip]); exhi.push_back(pt_ehi[ip]);
        eylo.push_back(0.0);        eyhi.push_back(0.0);
    }
    if (xv.empty()) return nullptr;
    return new TGraphAsymmErrors(
        (int)xv.size(), xv.data(), yv.data(), exlo.data(), exhi.data(), eylo.data(), eyhi.data());
}

// Min/max envelope (relative difference, %) across all N_SOURCES sources, per pT bin,
// for a single centrality. Bins with no valid source point, or |nominal| below
// rel_abs_threshold, are dropped.
static TGraphAsymmErrors* MakeEnvelopeBand(
    TGraphAsymmErrors* g_nom, TGraphAsymmErrors* g_src[N_SOURCES],
    const double* pt_cen, const double* pt_elo, const double* pt_ehi, int n_pt,
    double rel_abs_threshold)
{
    std::vector<double> xv, yv, exlo, exhi, eylo, eyhi;
    for (int ip = 0; ip < n_pt; ++ip) {
        PtVal pn = GetPtVal(g_nom, pt_cen[ip]);
        if (!pn.ok || std::fabs(pn.y) < rel_abs_threshold) continue;
        double mn = 1e18, mx = -1e18; bool any = false;
        for (int isrc = 0; isrc < N_SOURCES; ++isrc) {
            PtVal ps = GetPtVal(g_src[isrc], pt_cen[ip]);
            if (!ps.ok) continue;
            double d = (ps.y - pn.y) / pn.y * 100.0;
            mn = std::min(mn, d); mx = std::max(mx, d); any = true;
        }
        if (!any) continue;
        double mid = 0.5 * (mn + mx);
        xv.push_back(pt_cen[ip]);   yv.push_back(mid);
        exlo.push_back(pt_elo[ip]); exhi.push_back(pt_ehi[ip]);
        eylo.push_back(mid - mn);   eyhi.push_back(mx - mid);
    }
    if (xv.empty()) return nullptr;
    auto* g = new TGraphAsymmErrors(
        (int)xv.size(), xv.data(), yv.data(), exlo.data(), exhi.data(), eylo.data(), eyhi.data());
    g->SetFillColorAlpha(kGray + 1, 0.45);
    g->SetLineColor(kGray + 2);
    g->SetLineWidth(1);
    return g;
}

// =========================================================================
// Weighted Pearson r  (identical to plot_ESE_scatter_and_pearson_combined.C)
//   w_i = 1 / (ex_i^2 + ey_i^2);  returns NaN if ill-defined.
// =========================================================================
static double WeightedPearson(const std::vector<double>& x, const std::vector<double>& y,
                               const std::vector<double>& ex, const std::vector<double>& ey)
{
    int n = static_cast<int>(x.size());
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();

    std::vector<double> w(n);
    double sum_w = 0.0;
    for (int i = 0; i < n; ++i) {
        double sig2 = ex[i] * ex[i] + ey[i] * ey[i];
        w[i] = (sig2 > 0.0) ? 1.0 / sig2 : 0.0;
        sum_w += w[i];
    }
    if (sum_w <= 0.0) return std::numeric_limits<double>::quiet_NaN();

    double xw = 0.0, yw = 0.0;
    for (int i = 0; i < n; ++i) { xw += w[i] * x[i]; yw += w[i] * y[i]; }
    xw /= sum_w; yw /= sum_w;

    double cov = 0.0, varx = 0.0, vary = 0.0;
    for (int i = 0; i < n; ++i) {
        double dx = x[i] - xw, dy = y[i] - yw;
        cov += w[i] * dx * dy; varx += w[i] * dx * dx; vary += w[i] * dy * dy;
    }
    double denom = std::sqrt(varx * vary);
    if (denom <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    return cov / denom;
}

// =========================================================================
// Re-derives, for ONE dataset (one D0 file + one charged-particle file), the same
// weighted-Pearson-r / slope / intercept vs pT graphs that
// plot_ESE_scatter_and_pearson_combined.C computes — reusing its exact scatter-building,
// MC-resampling-uncertainty and linear-fit logic — but keeping only the final graphs in
// memory (no scatter/Pearson/MC diagnostic PDFs or ROOT files are written here).
// out_r/out_slope/out_inter are each TGraphAsymmErrors*[2][N_CENTBINS] (index 0 -> v2, 1 -> v3).
// =========================================================================
static void ComputeDatasetGraphs(
    const char* d0_file, const char* chg_file, bool isRun3,
    int N_MC, bool USE_SYMMETRIC_INTERVAL, std::mt19937_64& rng,
    TGraphAsymmErrors* out_r[2][N_CENTBINS],
    TGraphAsymmErrors* out_slope[2][N_CENTBINS],
    TGraphAsymmErrors* out_inter[2][N_CENTBINS],
    ScatterPoint out_scatter[2][N_CENTBINS][N_PTBINS_MAX][N_QBINS])
{

    TFile* fD0  = TFile::Open(d0_file);
    TFile* fChg = TFile::Open(chg_file);
    if (!fD0 || fD0->IsZombie() || !fChg || fChg->IsZombie()) {
        std::cerr << "ERROR: cannot open D0/chg file pair: " << d0_file << " , " << chg_file << "\n";
        if (fD0)  fD0->Close();
        if (fChg) fChg->Close();
        return;
    }
    auto* dir_d0_qbin = (TDirectory*)fD0->Get("vn_vs_qbin");
    auto* dir_chg_vsq = (TDirectory*)fChg->Get("vsQbin_TProfile");
    if (!dir_d0_qbin || !dir_chg_vsq) {
        std::cerr << "ERROR: vn_vs_qbin / vsQbin_TProfile directory missing in "
                   << d0_file << " / " << chg_file << "\n";
        fD0->Close(); fChg->Close();
        return;
    }

    const char* chg_pt_tag = isRun3 ? "pt0p5to3" : "pt1to3";
    TProfile* hp_chg_v2[N_CENTBINS] = {};
    TProfile* hp_chg_v3[N_CENTBINS] = {};
    for (int ic = 0; ic < N_CENTBINS; ++ic) {
        hp_chg_v2[ic] = (TProfile*)dir_chg_vsq->Get(Form("hp_v2_vsq2_%s_%s", chg_pt_tag, g_cen_name[ic]));
        hp_chg_v3[ic] = (TProfile*)dir_chg_vsq->Get(Form("hp_v3_vsq3_%s_%s", chg_pt_tag, g_cen_name[ic]));
    }

    for (int ivn = 2; ivn <= 3; ++ivn) {
        int vi = ivn - 2;
        int n_pt = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
        const double* pt_edges = (ivn == 2) ? g_pt_edges_v2 : g_pt_edges_v3;

        double pt_cen[N_PTBINS_MAX], pt_elo[N_PTBINS_MAX], pt_ehi[N_PTBINS_MAX];
        ComputePtCenters(pt_edges, n_pt, pt_cen, pt_elo, pt_ehi);

        double r_val[N_CENTBINS][N_PTBINS_MAX] = {}, r_elo[N_CENTBINS][N_PTBINS_MAX] = {}, r_ehi[N_CENTBINS][N_PTBINS_MAX] = {};
        bool   r_ok [N_CENTBINS][N_PTBINS_MAX] = {};
        double s_val[N_CENTBINS][N_PTBINS_MAX] = {}, s_err[N_CENTBINS][N_PTBINS_MAX] = {};
        bool   s_ok [N_CENTBINS][N_PTBINS_MAX] = {};
        double b_val[N_CENTBINS][N_PTBINS_MAX] = {}, b_err[N_CENTBINS][N_PTBINS_MAX] = {};
        bool   b_ok [N_CENTBINS][N_PTBINS_MAX] = {};

        for (int ic = 0; ic < N_CENTBINS; ++ic) {
            double chg_v2_q[N_QBINS] = {}, chg_v2_q_e[N_QBINS] = {};
            double chg_v3_q[N_QBINS] = {}, chg_v3_q_e[N_QBINS] = {};
            for (int iq = 0; iq < N_QBINS; ++iq) {
                if (hp_chg_v2[ic] && hp_chg_v2[ic]->GetBinEntries(iq + 1) > 0) {
                    chg_v2_q[iq] = hp_chg_v2[ic]->GetBinContent(iq + 1);
                    chg_v2_q_e[iq] = hp_chg_v2[ic]->GetBinError(iq + 1);
                }
                if (hp_chg_v3[ic] && hp_chg_v3[ic]->GetBinEntries(iq + 1) > 0) {
                    chg_v3_q[iq] = hp_chg_v3[ic]->GetBinContent(iq + 1);
                    chg_v3_q_e[iq] = hp_chg_v3[ic]->GetBinError(iq + 1);
                }
            }
            double* vn_chg_q   = (ivn == 2) ? chg_v2_q   : chg_v3_q;
            double* vn_chg_q_e = (ivn == 2) ? chg_v2_q_e : chg_v3_q_e;

            for (int ip = 0; ip < n_pt; ++ip) {
                TString ptTag = Form("pT%dto%d", (int)pt_edges[ip], (int)pt_edges[ip + 1]);
                auto* gr_d0 = (TGraphErrors*)dir_d0_qbin->Get(
                    Form("v%d_vs_q%dbin_%s_%s", ivn, ivn, g_cen_name[ic], ptTag.Data()));
                if (!gr_d0 || gr_d0->GetN() == 0) continue;

                std::vector<double> xv, yv, exv, eyv;
                for (int iq = 0; iq < N_QBINS && iq < gr_d0->GetN(); ++iq) {
                    double chg_v = vn_chg_q[iq], chg_e = vn_chg_q_e[iq];
                    double d0_v = gr_d0->GetPointY(iq), d0_e = gr_d0->GetErrorY(iq);
                    if (chg_e <= 0 || d0_e <= 0) continue;
                    xv.push_back(chg_v); exv.push_back(chg_e);
                    yv.push_back(d0_v);  eyv.push_back(d0_e);
                    out_scatter[vi][ic][ip][iq] = {chg_v, d0_v, chg_e, d0_e, true};
                }
                if ((int)xv.size() < 2) continue;

                double r_meas = WeightedPearson(xv, yv, exv, eyv);
                if (std::isnan(r_meas)) continue;

                int np = (int)xv.size();
                std::vector<double> r_mc; r_mc.reserve(N_MC);
                for (int imc = 0; imc < N_MC; ++imc) {
                    std::vector<double> yv_mc(np);
                    for (int i = 0; i < np; ++i) {
                        std::normal_distribution<double> gaus(yv[i], eyv[i]);
                        yv_mc[i] = gaus(rng);
                    }
                    double r_tmp = WeightedPearson(xv, yv_mc, exv, eyv);
                    if (!std::isnan(r_tmp)) r_mc.push_back(r_tmp);
                }
                if ((int)r_mc.size() < 100) continue;
                std::sort(r_mc.begin(), r_mc.end());
                int n_mc = (int)r_mc.size();
                double r_elo_v, r_ehi_v;
                if (!USE_SYMMETRIC_INTERVAL) {
                    double p16 = r_mc[(int)(0.1585 * n_mc)];
                    double p84 = r_mc[(int)(0.8415 * n_mc)];
                    r_elo_v = r_meas - p16; r_ehi_v = p84 - r_meas;
                } else {
                    double lo_d = 0.0, hi_d = 2.0;
                    for (int iter = 0; iter < 60; ++iter) {
                        double mid = 0.5 * (lo_d + hi_d);
                        auto it_lo = std::lower_bound(r_mc.begin(), r_mc.end(), r_meas - mid);
                        auto it_hi = std::upper_bound(r_mc.begin(), r_mc.end(), r_meas + mid);
                        double frac = (double)std::distance(it_lo, it_hi) / n_mc;
                        if (frac >= 0.683) hi_d = mid; else lo_d = mid;
                    }
                    r_elo_v = hi_d; r_ehi_v = hi_d;
                }
                r_val[ic][ip] = r_meas; r_elo[ic][ip] = r_elo_v; r_ehi[ic][ip] = r_ehi_v; r_ok[ic][ip] = true;

                {
                    double xmin_f = *std::min_element(xv.begin(), xv.end());
                    double xmax_f = *std::max_element(xv.begin(), xv.end());
                    double xmarg  = 0.3 * std::max(xmax_f - xmin_f, 1e-6);
                    auto* gr_fit = new TGraphErrors(np, xv.data(), yv.data(), exv.data(), eyv.data());
                    auto* flin = new TF1(Form("flin_sys_%d_%d_%d_%p", ivn, ic, ip, (void*)d0_file),
                                          "[0]+[1]*x", xmin_f - xmarg, xmax_f + xmarg);
                    flin->SetParameters(0.0, 1.0);
                    gr_fit->Fit(flin, "QR0");
                    double perr0 = flin->GetParError(0), perr1 = flin->GetParError(1);
                    if (perr0 > 0 && perr1 > 0) {
                        s_val[ic][ip] = flin->GetParameter(1); s_err[ic][ip] = perr1; s_ok[ic][ip] = true;
                        b_val[ic][ip] = flin->GetParameter(0); b_err[ic][ip] = perr0; b_ok[ic][ip] = true;
                    }
                    delete flin; delete gr_fit;
                }
            } // ip
        } // ic

        for (int ic = 0; ic < N_CENTBINS; ++ic) {
            std::vector<double> xr, yr, exlo, exhi, eylo, eyhi;
            for (int ip = 0; ip < n_pt; ++ip) {
                if (!r_ok[ic][ip]) continue;
                xr.push_back(pt_cen[ip]);   yr.push_back(r_val[ic][ip]);
                exlo.push_back(pt_elo[ip]); exhi.push_back(pt_ehi[ip]);
                eylo.push_back(r_elo[ic][ip]); eyhi.push_back(r_ehi[ic][ip]);
            }
            if (!xr.empty())
                out_r[vi][ic] = new TGraphAsymmErrors((int)xr.size(), xr.data(), yr.data(),
                    exlo.data(), exhi.data(), eylo.data(), eyhi.data());

            std::vector<double> xs, ys, exlos, exhis, eylos, eyhis;
            std::vector<double> xb, yb, exlob, exhib, eylob, eyhib;
            for (int ip = 0; ip < n_pt; ++ip) {
                if (s_ok[ic][ip]) {
                    xs.push_back(pt_cen[ip]);    ys.push_back(s_val[ic][ip]);
                    exlos.push_back(pt_elo[ip]); exhis.push_back(pt_ehi[ip]);
                    eylos.push_back(s_err[ic][ip]); eyhis.push_back(s_err[ic][ip]);
                }
                if (b_ok[ic][ip]) {
                    xb.push_back(pt_cen[ip]);    yb.push_back(b_val[ic][ip]);
                    exlob.push_back(pt_elo[ip]); exhib.push_back(pt_ehi[ip]);
                    eylob.push_back(b_err[ic][ip]); eyhib.push_back(b_err[ic][ip]);
                }
            }
            if (!xs.empty())
                out_slope[vi][ic] = new TGraphAsymmErrors((int)xs.size(), xs.data(), ys.data(),
                    exlos.data(), exhis.data(), eylos.data(), eyhis.data());
            if (!xb.empty())
                out_inter[vi][ic] = new TGraphAsymmErrors((int)xb.size(), xb.data(), yb.data(),
                    exlob.data(), exhib.data(), eylob.data(), eyhib.data());
        }
    } // ivn

    fD0->Close();
    fChg->Close();
}

// =========================================================================
// STEP 10: reads the nominal D0 file's "vn_vs_qmean" directory to get each q-bin's actual
// mean-q value (as opposed to STEP 9's plain q-bin index) -- these graphs carry the exact same
// v_n values as "vn_vs_qbin" (only the x-axis representation differs), so this is purely a
// per-q-bin x-position lookup; the y-values/systematics already computed for scatter_nom /
// scatter_syst_all apply unchanged. Only the nominal file is read (not every source variant --
// some don't even have this directory, e.g. the current SigPDF file), since the mean-q
// x-position of a point is a fixed reference, not something with its own propagated
// uncertainty here.
// =========================================================================
static void LoadMeanQNominal(const char* d0_file,
    double out_meanq[2][N_CENTBINS][N_PTBINS_MAX][N_QBINS],
    bool out_ok[2][N_CENTBINS][N_PTBINS_MAX][N_QBINS])
{
    for (int vi = 0; vi < 2; ++vi)
        for (int ic = 0; ic < N_CENTBINS; ++ic)
            for (int ip = 0; ip < N_PTBINS_MAX; ++ip)
                for (int iq = 0; iq < N_QBINS; ++iq)
                    out_ok[vi][ic][ip][iq] = false;

    TFile* fD0 = TFile::Open(d0_file);
    if (!fD0 || fD0->IsZombie()) {
        std::cerr << "WARNING: cannot open " << d0_file << " for mean-q lookup (STEP 10 skipped).\n";
        return;
    }
    auto* dir = (TDirectory*)fD0->Get("vn_vs_qmean");
    if (!dir) {
        std::cerr << "WARNING: no vn_vs_qmean directory in " << d0_file << " (STEP 10 skipped).\n";
        fD0->Close();
        return;
    }

    for (int ivn = 2; ivn <= 3; ++ivn) {
        int vi = ivn - 2;
        int n_pt = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
        const double* pt_edges = (ivn == 2) ? g_pt_edges_v2 : g_pt_edges_v3;
        for (int ic = 0; ic < N_CENTBINS; ++ic) {
            for (int ip = 0; ip < n_pt; ++ip) {
                TString name = Form("v%d_vs_meanq%d_%s_pT%dto%d", ivn, ivn, g_cen_name[ic],
                                     (int)pt_edges[ip], (int)pt_edges[ip + 1]);
                auto* g = (TGraphErrors*)dir->Get(name.Data());
                if (!g) continue;
                for (int iq = 0; iq < N_QBINS && iq < g->GetN(); ++iq) {
                    out_meanq[vi][ic][ip][iq] = g->GetPointX(iq);
                    out_ok[vi][ic][ip][iq] = true;
                }
            }
        }
    }
    fD0->Close();
}

// =========================================================================
// Prompt-fraction systematic: reads a precomputed per-(harmonic, centrality, pT) shift table
// (produced by compute_promptfrac_vn_shift.C; see that file / the table's own header for the
// formula) instead of rerunning the analysis on varied files. Each row is tagged with its
// harmonic (leading "2"/"3" column) and is matched against that harmonic's own native pT
// binning (g_pt_edges_v2 / g_pt_edges_v3); any row that doesn't match a known harmonic/
// centrality/pT tag is skipped.
// =========================================================================
static bool LoadPromptFractionShiftTable(const char* path,
    double shift_v2[N_CENTBINS][N_PTBINS_V2], double shift_v3[N_CENTBINS][N_PTBINS_V3])
{
    for (int ic = 0; ic < N_CENTBINS; ++ic) {
        for (int ip = 0; ip < N_PTBINS_V2; ++ip) shift_v2[ic][ip] = 0.0;
        for (int ip = 0; ip < N_PTBINS_V3; ++ip) shift_v3[ic][ip] = 0.0;
    }

    std::ifstream fin(path);
    if (!fin.is_open()) {
        std::cerr << "WARNING: cannot open prompt-fraction shift table: " << path << "\n";
        return false;
    }

    TString ptTag_v2[N_PTBINS_V2];
    for (int ip = 0; ip < N_PTBINS_V2; ++ip)
        ptTag_v2[ip] = Form("pT%dto%d", (int)g_pt_edges_v2[ip], (int)g_pt_edges_v2[ip + 1]);
    TString ptTag_v3[N_PTBINS_V3];
    for (int ip = 0; ip < N_PTBINS_V3; ++ip)
        ptTag_v3[ip] = Form("pT%dto%d", (int)g_pt_edges_v3[ip], (int)g_pt_edges_v3[ip + 1]);

    std::string line;
    int nfilled = 0;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        int vn;
        std::string centTok, ptTok;
        double f1, f1e, f2, f2e, factor, vn1_avg, vn2_avg, dvn_meas, shift;
        if (!(iss >> vn >> centTok >> ptTok >> f1 >> f1e >> f2 >> f2e >> factor
                   >> vn1_avg >> vn2_avg >> dvn_meas >> shift)) continue; // header/malformed

        int ic = -1;
        for (int k = 0; k < N_CENTBINS; ++k) if (centTok == g_cen_name[k]) { ic = k; break; }
        if (ic < 0) {
            std::cerr << "WARNING: prompt-fraction table row not recognized (cent='" << centTok
                      << "') -- skipped.\n";
            continue;
        }
        if (vn == 2) {
            int ip = -1;
            for (int k = 0; k < N_PTBINS_V2; ++k) if (ptTok == ptTag_v2[k].Data()) { ip = k; break; }
            if (ip < 0) {
                std::cerr << "WARNING: prompt-fraction table row not recognized (vn=2, pT='"
                          << ptTok << "') -- skipped.\n";
                continue;
            }
            shift_v2[ic][ip] = shift;
            ++nfilled;
        } else if (vn == 3) {
            int ip = -1;
            for (int k = 0; k < N_PTBINS_V3; ++k) if (ptTok == ptTag_v3[k].Data()) { ip = k; break; }
            if (ip < 0) {
                std::cerr << "WARNING: prompt-fraction table row not recognized (vn=3, pT='"
                          << ptTok << "') -- skipped.\n";
                continue;
            }
            shift_v3[ic][ip] = shift;
            ++nfilled;
        } else {
            std::cerr << "WARNING: prompt-fraction table row has unrecognized vn='" << vn << "' -- skipped.\n";
        }
    }
    std::cout << "Loaded " << nfilled << " prompt-fraction shift entries from " << path << "\n";
    return nfilled > 0;
}

// Builds the single "variant 0" for a table-shift source: r and slope are untouched clones of
// nominal (Delta r = Delta Slope = 0 by construction), intercept is nominal's intercept plus
// the table's shift at that (harmonic, centrality, pT) bin — v2 and v3 each use their own
// native-binned shift table, no cross-harmonic borrowing.
static void BuildPromptFractionVariant(
    double shift_v2[N_CENTBINS][N_PTBINS_V2], double shift_v3[N_CENTBINS][N_PTBINS_V3],
    TGraphAsymmErrors* gr_r_nom[2][N_CENTBINS], TGraphAsymmErrors* gr_slope_nom[2][N_CENTBINS],
    TGraphAsymmErrors* gr_inter_nom[2][N_CENTBINS],
    TGraphAsymmErrors* out_r[2][N_CENTBINS], TGraphAsymmErrors* out_slope[2][N_CENTBINS],
    TGraphAsymmErrors* out_inter[2][N_CENTBINS])
{
    double pt_cen_v2[N_PTBINS_V2], pt_elo_v2[N_PTBINS_V2], pt_ehi_v2[N_PTBINS_V2];
    ComputePtCenters(g_pt_edges_v2, N_PTBINS_V2, pt_cen_v2, pt_elo_v2, pt_ehi_v2);
    double pt_cen_v3[N_PTBINS_V3], pt_elo_v3[N_PTBINS_V3], pt_ehi_v3[N_PTBINS_V3];
    ComputePtCenters(g_pt_edges_v3, N_PTBINS_V3, pt_cen_v3, pt_elo_v3, pt_ehi_v3);

    for (int vi = 0; vi < 2; ++vi) {
        int n_pt = (vi == 0) ? N_PTBINS_V2 : N_PTBINS_V3;
        const double* pt_cen = (vi == 0) ? pt_cen_v2 : pt_cen_v3;
        for (int ic = 0; ic < N_CENTBINS; ++ic) {
            if (gr_r_nom[vi][ic])     out_r[vi][ic]     = (TGraphAsymmErrors*)gr_r_nom[vi][ic]->Clone();
            if (gr_slope_nom[vi][ic]) out_slope[vi][ic] = (TGraphAsymmErrors*)gr_slope_nom[vi][ic]->Clone();

            TGraphAsymmErrors* gnom_b = gr_inter_nom[vi][ic];
            if (!gnom_b) continue;
            auto* g = (TGraphAsymmErrors*)gnom_b->Clone();
            for (int i = 0; i < g->GetN(); ++i) {
                double x, y;
                g->GetPoint(i, x, y);
                int ip_bin = -1;
                for (int ip = 0; ip < n_pt; ++ip)
                    if (std::fabs(pt_cen[ip] - x) < 1e-4) { ip_bin = ip; break; }
                double shift = 0.0;
                if (ip_bin >= 0) shift = (vi == 0) ? shift_v2[ic][ip_bin] : shift_v3[ic][ip_bin];
                g->SetPoint(i, x, y + shift);
            }
            out_inter[vi][ic] = g;
        }
    }
}

// Same idea as BuildPromptFractionVariant but at the raw q-bin scatter-point level (STEP 7):
// every q-bin in a (centrality, pT) bin gets the same table shift added to its D0 vn (y),
// with charged vn (x) and both statistical errors left untouched.
static void BuildPromptFractionScatterVariant(
    double shift_v2[N_CENTBINS][N_PTBINS_V2], double shift_v3[N_CENTBINS][N_PTBINS_V3],
    ScatterPoint scatter_nom[2][N_CENTBINS][N_PTBINS_MAX][N_QBINS],
    ScatterPoint out_scatter[2][N_CENTBINS][N_PTBINS_MAX][N_QBINS])
{
    for (int vi = 0; vi < 2; ++vi) {
        int n_pt = (vi == 0) ? N_PTBINS_V2 : N_PTBINS_V3;
        for (int ic = 0; ic < N_CENTBINS; ++ic) {
            for (int ip = 0; ip < n_pt; ++ip) {
                double shift = (vi == 0) ? shift_v2[ic][ip] : shift_v3[ic][ip];
                for (int iq = 0; iq < N_QBINS; ++iq) {
                    ScatterPoint p = scatter_nom[vi][ic][ip][iq];
                    if (p.ok) p.y += shift;
                    out_scatter[vi][ic][ip][iq] = p;
                }
            }
        }
    }
}

// Per-pT-bin bookkeeping for the total-systematic calculation (STEP 5 / STEP 6).
struct BinSyst {
    bool   nom_ok = false;
    double nom_y = 0.0, nom_elo = 0.0, nom_ehi = 0.0;
    double syst_abs = 0.0;
    bool   src_ok[N_SOURCES] = {};
    double src_diff[N_SOURCES] = {};
};

// Total systematic = quadratic sum of (source_i - nominal) over all sources that have a
// valid point in this bin. Sources missing a point contribute 0 (not double counted).
static void ComputeTotalSyst(
    TGraphAsymmErrors* g_nom, TGraphAsymmErrors* g_src[N_SOURCES],
    const double* pt_cen, int n_pt, BinSyst* out)
{
    for (int ip = 0; ip < n_pt; ++ip) {
        BinSyst& b = out[ip];
        b = BinSyst();
        PtVal pn = GetPtVal(g_nom, pt_cen[ip]);
        if (!pn.ok) continue;
        b.nom_ok = true; b.nom_y = pn.y; b.nom_elo = pn.elo; b.nom_ehi = pn.ehi;
        double sumsq = 0.0;
        for (int isrc = 0; isrc < N_SOURCES; ++isrc) {
            PtVal ps = GetPtVal(g_src[isrc], pt_cen[ip]);
            b.src_ok[isrc] = ps.ok;
            if (ps.ok) { double d = ps.y - pn.y; b.src_diff[isrc] = d; sumsq += d * d; }
        }
        b.syst_abs = std::sqrt(sumsq);
    }
}

// STEP 7's per-q-bin analogue of ComputeTotalSyst: for one (harmonic, centrality, pT) bin,
// the quadratic sum across sources of each source's adopted (max-|diff|) deviation from
// nominal, computed independently for x (charged vn) and y (D0 vn) — same "treat each axis
// as its own 1D quantity, pick whichever of a source's up-to-two variants deviates more"
// rule already used for r/slope/intercept, just applied to the raw scatter points that feed
// that fit instead of to the fit's own outputs.
static void ComputeScatterTotalSyst(
    int vi, int ic, int ip,
    ScatterPoint scatter_nom[2][N_CENTBINS][N_PTBINS_MAX][N_QBINS],
    ScatterPoint scatter_var[N_SOURCES][2][2][N_CENTBINS][N_PTBINS_MAX][N_QBINS],
    const bool src_present[N_SOURCES],
    ScatterSyst out[N_QBINS])
{
    ScatterPoint* nom = scatter_nom[vi][ic][ip];
    for (int iq = 0; iq < N_QBINS; ++iq) {
        out[iq] = ScatterSyst();
        if (!nom[iq].ok) continue;
        out[iq].ok = true;
        double sumsqx = 0.0, sumsqy = 0.0;
        for (int isrc = 0; isrc < N_SOURCES; ++isrc) {
            if (!src_present[isrc]) continue;
            ScatterPoint& a = scatter_var[isrc][0][vi][ic][ip][iq];
            ScatterPoint& b = scatter_var[isrc][1][vi][ic][ip][iq];
            double dx = 0.0, dy = 0.0; bool any = false;
            if (a.ok && b.ok) {
                double dxa = a.x - nom[iq].x, dxb = b.x - nom[iq].x;
                dx = (std::fabs(dxa) >= std::fabs(dxb)) ? dxa : dxb;
                double dya = a.y - nom[iq].y, dyb = b.y - nom[iq].y;
                dy = (std::fabs(dya) >= std::fabs(dyb)) ? dya : dyb;
                any = true;
            } else if (a.ok) {
                dx = a.x - nom[iq].x; dy = a.y - nom[iq].y; any = true;
            } else if (b.ok) {
                dx = b.x - nom[iq].x; dy = b.y - nom[iq].y; any = true;
            }
            if (!any) continue;
            sumsqx += dx * dx; sumsqy += dy * dy;
        }
        out[iq].sysx = std::sqrt(sumsqx);
        out[iq].sysy = std::sqrt(sumsqy);
    }
}

// =========================================================================
// STEP 8: Pearson r's own systematic uncertainty, via Gaussian MC resampling rather than a
// per-source quadratic sum (see the IMPORTANT note in the header comment for the physics
// reasoning -- r is bounded to [-1,1], so an envelope built from per-source r-differences
// isn't valid). Procedure, identical to the STATISTICAL MC resampling already used for r
// (WeightedPearson + the N_MC loop inside ComputeDatasetGraphs, itself ported from
// plot_ESE_vnVsqn_scatter_and_pearson_combined.C), except the per-q-bin Gaussian smearing
// applied to each D0 v_n value uses sigma = that point's TOTAL SYSTEMATIC uncertainty
// (ComputeScatterTotalSyst's sysy) instead of its statistical uncertainty -- while the
// weights used inside WeightedPearson for every MC draw remain the point's STATISTICAL
// uncertainties, unchanged. The spread of the resulting r distribution (same symmetric-
// interval or 16/84-percentile extraction as the statistical case) IS the systematic
// uncertainty on r.
// =========================================================================
static void ComputePearsonSystUncertainty(
    ScatterPoint nom[N_QBINS], ScatterSyst syst[N_QBINS],
    int N_MC, bool USE_SYMMETRIC_INTERVAL, std::mt19937_64& rng,
    double& r_meas, double& r_syst_elo, double& r_syst_ehi, bool& ok,
    std::vector<double>* out_r_mc = nullptr)
{
    ok = false;
    std::vector<double> xv, yv, exv, eyv, sysyv;
    for (int iq = 0; iq < N_QBINS; ++iq) {
        if (!nom[iq].ok || nom[iq].ex <= 0 || nom[iq].ey <= 0) continue;
        xv.push_back(nom[iq].x);  yv.push_back(nom[iq].y);
        exv.push_back(nom[iq].ex); eyv.push_back(nom[iq].ey);
        sysyv.push_back(syst[iq].ok ? syst[iq].sysy : 0.0);
    }
    if ((int)xv.size() < 2) return;

    r_meas = WeightedPearson(xv, yv, exv, eyv);
    if (std::isnan(r_meas)) return;

    int np = (int)xv.size();
    std::vector<double> r_mc; r_mc.reserve(N_MC);
    for (int imc = 0; imc < N_MC; ++imc) {
        std::vector<double> yv_mc(np);
        for (int i = 0; i < np; ++i) {
            if (sysyv[i] > 0) {
                std::normal_distribution<double> gaus(yv[i], sysyv[i]);
                yv_mc[i] = gaus(rng);
            } else {
                yv_mc[i] = yv[i]; // no systematic at this point -> no smear
            }
        }
        double r_tmp = WeightedPearson(xv, yv_mc, exv, eyv); // weights stay statistical
        if (!std::isnan(r_tmp)) r_mc.push_back(r_tmp);
    }
    if ((int)r_mc.size() < 100) return;
    std::sort(r_mc.begin(), r_mc.end());
    int n_mc = (int)r_mc.size();

    if (!USE_SYMMETRIC_INTERVAL) {
        double p16 = r_mc[(int)(0.1585 * n_mc)];
        double p84 = r_mc[(int)(0.8415 * n_mc)];
        r_syst_elo = r_meas - p16; r_syst_ehi = p84 - r_meas;
    } else {
        double lo_d = 0.0, hi_d = 2.0;
        for (int iter = 0; iter < 60; ++iter) {
            double mid = 0.5 * (lo_d + hi_d);
            auto it_lo = std::lower_bound(r_mc.begin(), r_mc.end(), r_meas - mid);
            auto it_hi = std::upper_bound(r_mc.begin(), r_mc.end(), r_meas + mid);
            double frac = (double)std::distance(it_lo, it_hi) / n_mc;
            if (frac >= 0.683) hi_d = mid; else lo_d = mid;
        }
        r_syst_elo = hi_d; r_syst_ehi = hi_d;
    }
    if (out_r_mc) *out_r_mc = r_mc;
    ok = true;
}

static void DrawSystBoxesForGraph(TGraphAsymmErrors* gsyst, int color, double halfwidth_frac = 0.05)
{
    if (!gsyst) return;
    for (int i = 0; i < gsyst->GetN(); ++i) {
        double x = gsyst->GetPointX(i), y = gsyst->GetPointY(i);
        double elo = gsyst->GetErrorYlow(i), ehi = gsyst->GetErrorYhigh(i);
        double xlo = x * (1.0 - halfwidth_frac), xhi = x * (1.0 + halfwidth_frac);
        TBox* b = new TBox(xlo, y - elo, xhi, y + ehi);
        b->SetFillColorAlpha(color, 0.30);
        b->SetLineColor(color);
        b->SetLineWidth(1);
        b->Draw();
    }
}

// =========================================================================
// STEP 7: draws the nominal D0-vn-vs-charged-vn scatter for one (harmonic, centrality, pT)
// bin — same TH2F frame, linear fit + 1-sigma confidence band, legend and CMS headers as
// plot_ESE_scatter_and_pearson_combined.C's own scatter plots — with an added per-q-bin
// total-systematic box (gray, from ComputeScatterTotalSyst) drawn behind the statistical
// point and error bar.
// =========================================================================
static void DrawScatterWithSystBox(
    const TString& outpath, int ivn, int ic,
    double pt_lo, double pt_hi,
    ScatterPoint* nom, ScatterSyst* syst,
    double xlo, double xhi, double ylo, double yhi)
{
    std::vector<double> xv, yv, exv, eyv, sysx, sysy;
    for (int iq = 0; iq < N_QBINS; ++iq) {
        if (!nom[iq].ok) continue;
        xv.push_back(nom[iq].x);   yv.push_back(nom[iq].y);
        exv.push_back(nom[iq].ex); eyv.push_back(nom[iq].ey);
        // Box x-half-width is a fixed fraction of the point's own x-value, not the computed
        // charged-vn systematic (sysx) -- the charged-vn axis has no independent systematic
        // to propagate here; only the box's y-half-height (D0-vn systematic) is data-derived.
        sysx.push_back(nom[iq].x * SYST_BOX_XFRAC);
        sysy.push_back(syst[iq].ok ? syst[iq].sysy : 0.0);
    }
    if ((int)xv.size() < 2) return;

    TString sc_name = Form("scatter_syst_v%d_%s_pT%dto%d", ivn, g_cen_name[ic], (int)pt_lo, (int)pt_hi);

    auto* gr_sc = new TGraphErrors((int)xv.size(), xv.data(), yv.data(), exv.data(), eyv.data());
    gr_sc->SetName(sc_name); gr_sc->SetTitle(sc_name);

    // Linear fit y = p0 + p1*x (same as the nominal macro's scatter plots)
    auto* flin = new TF1(Form("flin_%s", sc_name.Data()), "[0]+[1]*x", xlo, xhi);
    flin->SetParameters(0.0, 1.0);
    flin->SetLineColor(kRed);
    flin->SetLineWidth(2);
    gr_sc->Fit(flin, "QR");

    // 1-sigma confidence interval band
    auto* h_ci = new TH1D(Form("hci_%s", sc_name.Data()), "", 400, xlo, xhi);
    h_ci->SetDirectory(nullptr);
    TVirtualFitter* vfitter = TVirtualFitter::GetFitter();
    if (vfitter) vfitter->GetConfidenceIntervals(h_ci, 0.68);
    h_ci->SetFillColorAlpha(kRed, 0.35);
    h_ci->SetFillStyle(1001);
    h_ci->SetMarkerSize(0);

    TCanvas* cv = new TCanvas(Form("cv_%s", sc_name.Data()), sc_name, 600, 600);
    cv->SetLeftMargin(0.14);
    cv->SetBottomMargin(0.14);
    cv->SetRightMargin(0.05);
    cv->SetTopMargin(0.07);
    cv->SetGridx();
    cv->SetGridy();

    auto* hf = new TH2F(Form("hf_%s", sc_name.Data()), "", 10, xlo, xhi, 10, ylo, yhi);
    hf->GetXaxis()->SetTitle(Form("charged particle v_{%d}", ivn));
    hf->GetYaxis()->SetTitle(Form("D^{0} v_{%d}", ivn));
    hf->GetXaxis()->SetTitleSize(0.050);
    hf->GetYaxis()->SetTitleSize(0.050);
    hf->GetXaxis()->SetTitleOffset(1.10);
    hf->GetYaxis()->SetTitleOffset(1.25);
    hf->GetXaxis()->SetLabelSize(0.040);
    hf->GetYaxis()->SetLabelSize(0.040);
    hf->Draw();

    h_ci->Draw("E3 SAME");

    // Total-systematic box per point, drawn behind the fit line and the statistical point.
    for (size_t i = 0; i < xv.size(); ++i) {
        auto* b = new TBox(xv[i] - sysx[i], yv[i] - sysy[i], xv[i] + sysx[i], yv[i] + sysy[i]);
        b->SetFillColorAlpha(kGray + 2, 0.35);
        b->SetLineColor(kGray + 2);
        b->SetLineWidth(1);
        b->Draw();
    }

    flin->Draw("SAME");
    gr_sc->SetMarkerStyle(20);
    gr_sc->SetMarkerSize(0.9);
    gr_sc->SetMarkerColor(kBlack);
    gr_sc->SetLineColor(kBlack);
    gr_sc->Draw("P SAME");

    TLegend* leg = new TLegend(0.10, 0.76, 0.45, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.038);
    leg->AddEntry((TObject*)nullptr, Form("Centrality: %d-%d%%", g_min_cent[ic], g_max_cent[ic]), "");
    leg->AddEntry((TObject*)nullptr, Form("%.0f < p_{T} < %.0f GeV/c", pt_lo, pt_hi), "");
    auto* boxProxy = new TBox();
    boxProxy->SetFillColorAlpha(kGray + 2, 0.35);
    boxProxy->SetLineColor(kGray + 2);
    leg->AddEntry(boxProxy, "Total syst. unc.", "f");
    leg->Draw();

    TLatex ltx;
    ltx.SetNDC();
    ltx.SetTextAlign(31);
    ltx.SetTextSize(0.038);
    ltx.DrawLatex(0.92, 0.32, Form("slope = %.3f #pm %.3f", flin->GetParameter(1), flin->GetParError(1)));
    ltx.DrawLatex(0.92, 0.27, Form("intercept = %.3f #pm %.3f", flin->GetParameter(0), flin->GetParError(0)));
    ltx.DrawLatex(0.92, 0.22, Form("#chi^{2}/NDF = %.2f",
                                    flin->GetChisquare() / std::max(1.0, (double)flin->GetNDF())));
    ltx.SetTextAlign(11);
    ltx.SetTextFont(42);
    ltx.SetTextSize(0.048);
    ltx.DrawLatex(0.14, 0.945, "#bf{CMS} #it{Preliminary}");
    ltx.SetTextAlign(31);
    ltx.SetTextFont(42);
    ltx.SetTextSize(0.042);
    ltx.DrawLatex(0.95, 0.945, "PbPb 5.36 TeV");

    cv->SaveAs(outpath);
    delete leg;
    delete hf;
    delete h_ci;
    delete flin;
    delete gr_sc;
    delete cv;
}

// =========================================================================
// STEP 9 / STEP 10: draws D0 v_n vs the q-bin index (0..N_QBINS-1) -- one combined canvas per
// (harmonic, centrality), 3-column grid of pT-bin sub-pads, same layout as STEP 8's mc_r_syst
// canvases. STEP 9 and STEP 10 plot the exact same points at the exact same x-positions (the
// q-bin index); STEP 10 only relabels those same tick positions with each bin's mean-q value
// instead of its plain index (pass meanq_label / meanq_label_ok for that; leave both nullptr
// for STEP 9's plain numeric axis). Because the underlying x-axis geometry is identical
// between the two, a single fixed box_halfwidth (in q-bin-index units) renders identically
// sized in both -- there is no separate "mean-q axis scale" for it to look wrong against.
// The box's y-half-height is always the real propagated D0-v_n systematic (from
// ComputeScatterTotalSyst, passed in via scatter_syst_ic); the box's x-half-width is that
// single fixed, hand-picked box_halfwidth -- the q-bin-index axis has no systematic of its own.
// =========================================================================
static void DrawVnVsQAxisSystGrid(
    const TString& outpath, int ivn, int ic, int n_pt, const double* pt_edges,
    ScatterPoint scatter_nom_ic[N_PTBINS_MAX][N_QBINS],
    ScatterSyst  scatter_syst_ic[N_PTBINS_MAX][N_QBINS],
    double box_halfwidth, const char* xaxis_title,
    double meanq_label[N_PTBINS_MAX][N_QBINS] = nullptr,
    bool   meanq_label_ok[N_PTBINS_MAX][N_QBINS] = nullptr)
{
    const int N_COLS = 3;
    const int N_ROWS = (n_pt + N_COLS - 1) / N_COLS;
    TCanvas* cv = new TCanvas(Form("cv_vnq_%d_%s_%p", ivn, g_cen_name[ic], (void*)&outpath), "",
                               1600, N_ROWS * 400);
    cv->Divide(N_COLS, N_ROWS, 0.001, 0.001);

    for (int ip = 0; ip < n_pt; ++ip) {
        cv->cd(ip + 1);
        gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.03);
        gPad->SetBottomMargin(0.24); gPad->SetTopMargin(0.10);

        std::vector<double> xv, yv, exv, eyv, bxlo, bxhi, bylo, byhi;
        for (int iq = 0; iq < N_QBINS; ++iq) {
            if (!scatter_nom_ic[ip][iq].ok) continue;
            double x = iq; // plain q-bin index, identical for STEP 9 and STEP 10
            double y = scatter_nom_ic[ip][iq].y, ey = scatter_nom_ic[ip][iq].ey;
            double sy = scatter_syst_ic[ip][iq].ok ? scatter_syst_ic[ip][iq].sysy : 0.0;
            xv.push_back(x); yv.push_back(y); exv.push_back(0.0); eyv.push_back(ey);
            bxlo.push_back(x - box_halfwidth); bxhi.push_back(x + box_halfwidth);
            bylo.push_back(y - sy); byhi.push_back(y + sy);
        }
        if (xv.empty()) {
            TH2F* hblank = new TH2F(Form("hblank_vnq_%d_%d_%d", ivn, ic, ip), "",
                                     N_QBINS, -0.5, N_QBINS - 0.5, 10, -0.2, 0.2);
            hblank->GetXaxis()->SetTitle(xaxis_title);
            hblank->Draw();
            continue;
        }

        double ymin = 1e9, ymax = -1e9;
        for (size_t i = 0; i < yv.size(); ++i) {
            ymin = std::min(ymin, std::min(yv[i] - eyv[i], bylo[i]));
            ymax = std::max(ymax, std::max(yv[i] + eyv[i], byhi[i]));
        }
        double ymarg = 0.15 * std::max(ymax - ymin, 1e-6);

        // Fixed x-range [-0.5, N_QBINS-0.5] always -- exactly N_QBINS bins, one per q-bin
        // index, so bin centers land exactly on 0, 1, ..., N_QBINS-1 (matching the plotted
        // points) whether this panel is STEP 9's or STEP 10's. The half-unit margin on each
        // side comfortably fits box_halfwidth without clipping.
        TH2F* hf = new TH2F(Form("hf_vnq_%d_%d_%d", ivn, ic, ip), "",
                             N_QBINS, -0.5, N_QBINS - 0.5, 10, ymin - ymarg, ymax + ymarg);
        hf->GetXaxis()->SetTitle(xaxis_title);
        hf->GetXaxis()->SetTitleSize(0.065); hf->GetXaxis()->SetTitleOffset(1.25);
        hf->GetYaxis()->SetTitle(Form("D^{0} v_{%d}", ivn));
        hf->GetYaxis()->SetTitleSize(0.065); hf->GetYaxis()->SetTitleOffset(1.05);
        hf->GetYaxis()->SetLabelSize(0.055);

        if (meanq_label) {
            // Relabel only ~5 evenly-spaced tick positions (not all N_QBINS) with that bin's
            // mean-q value -- the points/boxes above are unaffected, only the axis text; bins
            // left unset show no label at all (ROOT leaves them blank once any bin label is
            // set, rather than falling back to a numeric tick), avoiding a crowded axis.
            const int N_SHOWN = 5;
            for (int k = 0; k < N_SHOWN; ++k) {
                int iq = (int)std::lround(k * (double)(N_QBINS - 1) / (N_SHOWN - 1));
                if (meanq_label_ok && meanq_label_ok[ip][iq])
                    hf->GetXaxis()->SetBinLabel(iq + 1, Form("%.3f", meanq_label[ip][iq]));
            }
            hf->GetXaxis()->SetLabelSize(0.055);
        } else {
            hf->GetXaxis()->SetLabelSize(0.055);
        }
        hf->Draw("AXIS");

        for (size_t i = 0; i < xv.size(); ++i) {
            TBox* b = new TBox(bxlo[i], bylo[i], bxhi[i], byhi[i]);
            b->SetFillColorAlpha(kGray + 2, 0.35);
            b->SetLineColor(kGray + 2);
            b->SetLineWidth(1);
            b->Draw();
        }

        auto* gr = new TGraphErrors((int)xv.size(), xv.data(), yv.data(), exv.data(), eyv.data());
        gr->SetMarkerStyle(20); gr->SetMarkerSize(0.9);
        gr->SetMarkerColor(kBlack); gr->SetLineColor(kBlack);
        gr->Draw("P SAME");

        TLegend* leg = new TLegend(0.17, 0.75, 0.60, 0.90);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.060);
        leg->AddEntry((TObject*)nullptr, Form("Cent: %d-%d%%", g_min_cent[ic], g_max_cent[ic]), "");
        leg->AddEntry((TObject*)nullptr, Form("%.0f< p_{T} < %.0f GeV/c", pt_edges[ip], pt_edges[ip + 1]), "");
        leg->Draw();
    }

    cv->SaveAs(outpath);
    delete cv;
}

// =========================================================================
// Generic 2-panel (low-centrality | high-centrality) canvas for a single quantity
// (r, or an abs/rel-difference of r, slope or intercept) vs pT.
//   g_main      : required, per-centrality graph (filled marker / solid line)
//   g_over      : optional per-centrality overlay graph (open marker / dashed line)
//   g_over2     : optional second overlay graph (open marker / dotted line) — used for a
//                 second systematic-source variation (e.g. "Up" alongside g_over's "Down")
//   g_syst_box  : optional per-centrality graph whose Y-errors are drawn as shaded boxes
//                 (used for the total-systematic plot)
// =========================================================================
static void DrawQuantityPanelCanvas(
    const TString& outpath, int ivn, const char* ytitle,
    double Y_LO, double Y_HI, bool draw_ref, double yref,
    TGraphAsymmErrors** g_main,
    TGraphAsymmErrors** g_over = nullptr,
    const char* over_label = nullptr,
    TGraphAsymmErrors** g_syst_box = nullptr,
    const char* corner_label = nullptr,
    TGraphAsymmErrors** g_over2 = nullptr,
    const char* over_label2 = nullptr)
{
    TCanvas* cv = new TCanvas(Form("cv_qp_%d_%p", ivn, (void*)&outpath), "", 1200, 600);
    cv->SetFillColor(0); cv->SetBorderSize(0);

    const double bm = 0.20, tm = 0.09, split = 0.50;
    TPad* pad[2];
    pad[0] = new TPad("p0", "", 0.00, 0.00, split, 1.00);
    pad[1] = new TPad("p1", "", split, 0.00, 1.00, 1.00);
    pad[0]->SetLeftMargin(0.14);  pad[0]->SetRightMargin(0.000);
    pad[0]->SetBottomMargin(bm);  pad[0]->SetTopMargin(tm);
    pad[0]->SetLogx();            pad[0]->SetFillColor(0);
    pad[1]->SetLeftMargin(0.000); pad[1]->SetRightMargin(0.05);
    pad[1]->SetBottomMargin(bm);  pad[1]->SetTopMargin(tm);
    pad[1]->SetLogx();            pad[1]->SetFillColor(0);
    cv->cd(); pad[0]->Draw();
    cv->cd(); pad[1]->Draw();

    for (int ipanel = 0; ipanel < 2; ++ipanel) {
        pad[ipanel]->cd();

        TH2F* hf = new TH2F(Form("hf_qp_%d_%d_%p", ivn, ipanel, (void*)&outpath), "",
                             100, X_LO, X_HI, 100, Y_LO, Y_HI);
        hf->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hf->GetXaxis()->SetTitleSize(0.060); hf->GetXaxis()->SetTitleOffset(0.95);
        hf->GetXaxis()->SetLabelSize(0.050); hf->GetXaxis()->SetMoreLogLabels();
        hf->GetXaxis()->SetNoExponent();     hf->GetXaxis()->CenterTitle(kTRUE);
        if (ipanel == 0) {
            hf->GetYaxis()->SetTitle(ytitle);
            hf->GetYaxis()->SetTitleSize(0.065); hf->GetYaxis()->SetTitleOffset(0.95);
            hf->GetYaxis()->SetLabelSize(0.050); hf->GetYaxis()->CenterTitle(kTRUE);
        } else {
            hf->GetYaxis()->SetTitle(""); hf->GetYaxis()->SetLabelSize(0.0);
            hf->GetYaxis()->SetTickLength(0.03);
        }
        hf->Draw("AXIS");

        if (draw_ref) {
            TLine* lref = new TLine(X_LO, yref, X_HI, yref);
            lref->SetLineStyle(2); lref->SetLineWidth(1); lref->SetLineColor(kGray + 2);
            lref->Draw();
        }

        for (int k = 0; k < PANEL_NC[ipanel]; ++k) {
            int ic = PANEL_IC[ipanel][k];
            if (ic < 0) continue;
            if (g_syst_box && g_syst_box[ic]) DrawSystBoxesForGraph(g_syst_box[ic], g_col[ic]);
            if (g_over2 && g_over2[ic]) g_over2[ic]->Draw("P Z SAME");
            if (g_over && g_over[ic]) g_over[ic]->Draw("P Z SAME");
            if (g_main[ic]) g_main[ic]->Draw("P Z SAME");
        }

        {
            double leg_x2 = (ipanel == 0) ? 0.97 : 0.94;
            double leg_x1 = leg_x2 - 0.38;
            double leg_y2 = 0.88;
            double leg_y1 = leg_y2 - 0.082 * (PANEL_NC[ipanel] + 1);
            TLegend* leg = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.050);
            leg->AddEntry((TObject*)nullptr, "Centrality", "");
            for (int k = 0; k < PANEL_NC[ipanel]; ++k) {
                int ic = PANEL_IC[ipanel][k];
                if (ic < 0 || !g_main[ic]) continue;
                leg->AddEntry(g_main[ic], Form("%d-%d%%", g_min_cent[ic], g_max_cent[ic]), "p");
            }
            leg->Draw();
        }

        TLatex ltx; ltx.SetNDC();
        if (ipanel == 0) {
            ltx.SetTextFont(42); ltx.SetTextSize(0.055); ltx.SetTextAlign(11);
            ltx.DrawLatex(0.17, 0.930, "#bf{CMS} #it{Preliminary}");
        }
        ltx.SetTextFont(42); ltx.SetTextSize(0.045); ltx.SetTextAlign(11);
        ltx.DrawLatex((ipanel == 0) ? 0.17 : 0.05, 0.860, "|y| < 1");
        if (over_label && ipanel == 1) {
            ltx.SetTextSize(0.040); ltx.SetTextAlign(31);
            if (over_label2)
                ltx.DrawLatex(0.94, 0.860, Form("open (dashed) = %s, open (dotted) = %s", over_label, over_label2));
            else
                ltx.DrawLatex(0.94, 0.860, Form("open = %s", over_label));
        }
    }

    cv->cd();
    TLatex ltxc; ltxc.SetNDC(); ltxc.SetTextFont(42); ltxc.SetTextSize(0.045); ltxc.SetTextAlign(31);
    ltxc.DrawLatex(0.975, 0.940, "PbPb 5.36 TeV");
    if (corner_label) {
        TLatex ltxl; ltxl.SetNDC(); ltxl.SetTextFont(42); ltxl.SetTextSize(0.032); ltxl.SetTextAlign(11);
        ltxl.DrawLatex(0.02, 0.015, corner_label);
    }

    cv->SaveAs(outpath);
    delete cv;
}

// =========================================================================
// Generic 4-pad (slope top / intercept bottom, low-cent | high-cent) canvas.
// Mirrors DrawQuantityPanelCanvas but for the slope+intercept pair, with independently
// settable reference lines (1.0/0.0 for absolute values, 0.0/0.0 for differences).
// =========================================================================
static void DrawSlopeInterceptPanelCanvas(
    const TString& outpath, int ivn,
    double SL_LO, double SL_HI, double IN_LO, double IN_HI,
    double SL_REF, double IN_REF,
    TGraphAsymmErrors** g_slope_main, TGraphAsymmErrors** g_inter_main,
    TGraphAsymmErrors** g_slope_over = nullptr, TGraphAsymmErrors** g_inter_over = nullptr,
    const char* over_label = nullptr,
    TGraphAsymmErrors** g_slope_syst = nullptr, TGraphAsymmErrors** g_inter_syst = nullptr,
    const char* corner_label = nullptr,
    TGraphAsymmErrors** g_slope_over2 = nullptr, TGraphAsymmErrors** g_inter_over2 = nullptr,
    const char* over_label2 = nullptr,
    const char* slope_title = "Slope", const char* inter_title = "Intercept")
{
    const double y_lo_si[2]  = {SL_LO, IN_LO};
    const double y_hi_si[2]  = {SL_HI, IN_HI};
    const double y_ref_si[2] = {SL_REF, IN_REF};
    const char*  y_tit_si[2] = {slope_title, inter_title};

    TCanvas* cv = new TCanvas(Form("cv_si_%d_%p", ivn, (void*)&outpath), "", 1200, 800);
    cv->SetFillColor(0); cv->SetBorderSize(0);

    TPad* pad[4];
    pad[0] = new TPad("p0", "", 0.00, 0.50, 0.50, 1.00);
    pad[1] = new TPad("p1", "", 0.50, 0.50, 1.00, 1.00);
    pad[2] = new TPad("p2", "", 0.00, 0.00, 0.50, 0.50);
    pad[3] = new TPad("p3", "", 0.50, 0.00, 1.00, 0.50);
    // A strictly-zero bottom margin here clips the bottom-most y-axis tick label (e.g. a "-3")
    // right at the pad edge -- a small non-zero margin fixes that while still reading as
    // visually flush against the Intercept row below.
    pad[0]->SetLeftMargin(0.14);  pad[0]->SetRightMargin(0.000); pad[0]->SetBottomMargin(0.02); pad[0]->SetTopMargin(0.12);
    pad[1]->SetLeftMargin(0.000); pad[1]->SetRightMargin(0.05);  pad[1]->SetBottomMargin(0.02); pad[1]->SetTopMargin(0.12);
    pad[2]->SetLeftMargin(0.14);  pad[2]->SetRightMargin(0.000); pad[2]->SetBottomMargin(0.28); pad[2]->SetTopMargin(0.00);
    pad[3]->SetLeftMargin(0.000); pad[3]->SetRightMargin(0.05);  pad[3]->SetBottomMargin(0.28); pad[3]->SetTopMargin(0.00);
    for (int i = 0; i < 4; ++i) { pad[i]->SetLogx(); pad[i]->SetFillColor(0); cv->cd(); pad[i]->Draw(); }

    for (int irow = 0; irow < 2; ++irow) {
        TGraphAsymmErrors** gmain  = (irow == 0) ? g_slope_main  : g_inter_main;
        TGraphAsymmErrors** gover  = (irow == 0) ? g_slope_over  : g_inter_over;
        TGraphAsymmErrors** gover2 = (irow == 0) ? g_slope_over2 : g_inter_over2;
        TGraphAsymmErrors** gsyst  = (irow == 0) ? g_slope_syst  : g_inter_syst;
        for (int ipanel = 0; ipanel < 2; ++ipanel) {
            int  ipad    = irow * 2 + ipanel;
            bool is_top  = (irow == 0);
            bool is_left = (ipanel == 0);
            pad[ipad]->cd();

            TH2F* hfsi = new TH2F(Form("hfsi_%d_%d_%d_%p", ivn, irow, ipanel, (void*)&outpath), "",
                                   100, X_LO, X_HI, 100, y_lo_si[irow], y_hi_si[irow]);
            hfsi->GetYaxis()->SetNdivisions(505);
            if (!is_top) {
                hfsi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                hfsi->GetXaxis()->SetTitleSize(0.090); hfsi->GetXaxis()->SetTitleOffset(0.80);
                hfsi->GetXaxis()->SetLabelSize(0.075); hfsi->GetXaxis()->CenterTitle(kTRUE);
            } else {
                hfsi->GetXaxis()->SetTitle(""); hfsi->GetXaxis()->SetLabelSize(0.0);
                hfsi->GetXaxis()->SetTickLength(0.06);
            }
            hfsi->GetXaxis()->SetMoreLogLabels(); hfsi->GetXaxis()->SetNoExponent();
            if (is_left) {
                hfsi->GetYaxis()->SetTitle(y_tit_si[irow]);
                hfsi->GetYaxis()->SetTitleSize(0.090); hfsi->GetYaxis()->SetTitleOffset(0.72);
                hfsi->GetYaxis()->SetLabelSize(0.075); hfsi->GetYaxis()->CenterTitle(kTRUE);
            } else {
                hfsi->GetYaxis()->SetTitle(""); hfsi->GetYaxis()->SetLabelSize(0.0);
                hfsi->GetYaxis()->SetTickLength(0.03);
            }
            hfsi->Draw("AXIS");

            TLine* lref = new TLine(X_LO, y_ref_si[irow], X_HI, y_ref_si[irow]);
            lref->SetLineStyle(2); lref->SetLineWidth(1); lref->SetLineColor(kGray + 2); lref->Draw();

            for (int k = 0; k < PANEL_NC[ipanel]; ++k) {
                int ic = PANEL_IC[ipanel][k];
                if (ic < 0) continue;
                if (gsyst && gsyst[ic]) DrawSystBoxesForGraph(gsyst[ic], g_col[ic]);
                if (gover2 && gover2[ic]) gover2[ic]->Draw("P Z SAME");
                if (gover && gover[ic]) gover[ic]->Draw("P Z SAME");
                if (gmain[ic]) gmain[ic]->Draw("P Z SAME");
            }

            if (is_top) {
                TLatex ltxsi; ltxsi.SetNDC();
                if (is_left) {
                    ltxsi.SetTextFont(42); ltxsi.SetTextSize(0.095); ltxsi.SetTextAlign(11);
                    ltxsi.DrawLatex(0.17, 0.915, "#bf{CMS} #it{Preliminary}");
                }
                double sleg_x2 = is_left ? 0.97 : 0.94;
                double sleg_x1 = sleg_x2 - 0.44;
                double sleg_y2 = 0.87;
                double sleg_y1 = sleg_y2 - 0.095 * (PANEL_NC[ipanel] + 1);
                TLegend* sleg = new TLegend(sleg_x1, sleg_y1, sleg_x2, sleg_y2);
                sleg->SetBorderSize(0); sleg->SetFillStyle(0); sleg->SetTextFont(42); sleg->SetTextSize(0.078);
                sleg->AddEntry((TObject*)nullptr, "Centrality", "");
                for (int k = 0; k < PANEL_NC[ipanel]; ++k) {
                    int ic = PANEL_IC[ipanel][k];
                    if (ic < 0 || !gmain[ic]) continue;
                    sleg->AddEntry(gmain[ic], Form("%d-%d%%", g_min_cent[ic], g_max_cent[ic]), "p");
                }
                sleg->Draw();
                ltxsi.SetTextFont(42); ltxsi.SetTextSize(0.078); ltxsi.SetTextAlign(11);
                ltxsi.DrawLatex(is_left ? 0.17 : 0.05, 0.80, "|y| < 1");
                if (over_label && !is_left) {
                    ltxsi.SetTextSize(0.070); ltxsi.SetTextAlign(31);
                    if (over_label2)
                        ltxsi.DrawLatex(0.94, 0.80, Form("dashed = %s, dotted = %s", over_label, over_label2));
                    else
                        ltxsi.DrawLatex(0.94, 0.80, Form("open = %s", over_label));
                }
            }
        }
    }

    cv->cd();
    TLatex ltxtop; ltxtop.SetNDC(); ltxtop.SetTextFont(42); ltxtop.SetTextSize(0.035); ltxtop.SetTextAlign(31);
    ltxtop.DrawLatex(0.975, 0.955, "PbPb 5.36 TeV");
    if (corner_label) {
        TLatex ltxl; ltxl.SetNDC(); ltxl.SetTextFont(42); ltxl.SetTextSize(0.022); ltxl.SetTextAlign(11);
        ltxl.DrawLatex(0.01, 0.008, corner_label);
    }

    cv->SaveAs(outpath);
    delete cv;
}

// =========================================================================
// Main function
// =========================================================================
void plot_ESE_systematics()
{
    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    // =====================================================================
    // Output directory
    // =====================================================================
    TString outbase = outdir_tag;
    if (outbase.IsNull()) {
        TDatime now;
        outbase = Form("ESE_systematics_%d", now.GetDate());
    }
    std::cout << "\n=== Output directory: " << outbase << " ===\n";

    gSystem->mkdir(outbase, kTRUE);
    gSystem->mkdir(Form("%s/overlay_plots/v2", outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/overlay_plots/v3", outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/absdiff_plots/v2", outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/absdiff_plots/v3", outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/reldiff_plots/v2", outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/reldiff_plots/v3", outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/summary_plots",    outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/totalunc_plots",   outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/pearson_syst_plots", outbase.Data()), kTRUE);
    gSystem->mkdir(Form("%s/text_output",      outbase.Data()), kTRUE);

    // =====================================================================
    // pT bin centers
    // =====================================================================
    double pt_cen_v2[N_PTBINS_V2], pt_elo_v2[N_PTBINS_V2], pt_ehi_v2[N_PTBINS_V2];
    double pt_cen_v3[N_PTBINS_V3], pt_elo_v3[N_PTBINS_V3], pt_ehi_v3[N_PTBINS_V3];
    ComputePtCenters(g_pt_edges_v2, N_PTBINS_V2, pt_cen_v2, pt_elo_v2, pt_ehi_v2);
    ComputePtCenters(g_pt_edges_v3, N_PTBINS_V3, pt_cen_v3, pt_elo_v3, pt_ehi_v3);

    // =====================================================================
    // Compute Pearson r / slope / intercept, directly from raw D0/charged-particle files:
    //  - nominal (ids=0 in gr_r/gr_slope/gr_inter)
    //  - for each of the 7 sources (ids=1..7): every raw-file variation it has is computed
    //    into gr_*_var[isrc][variant]; the "adopted" value per (harmonic,centrality,pT) bin
    //    — the one variation deviating more from nominal, or the only variation if just one
    //    exists — is then synthesized into gr_*[isrc+1], exactly where a one-variation
    //    source's single result lived before. Everything below this loop (STEP4/5/6's
    //    "Adopted" figures) reads only gr_*[isrc+1] and is unaffected by variation count.
    // =====================================================================
    const int N_DS = 1 + N_SOURCES;
    TGraphAsymmErrors* gr_r    [N_DS][2][N_CENTBINS] = {};
    TGraphAsymmErrors* gr_slope[N_DS][2][N_CENTBINS] = {};
    TGraphAsymmErrors* gr_inter[N_DS][2][N_CENTBINS] = {};
    bool ds_loaded[N_DS] = {};

    // Raw per-variation graphs, kept around for STEP1/STEP2-3/STEP6's variant-level detail.
    // [source][variant 0/1][harmonic][centrality]; variant 1 stays all-nullptr for a
    // one-variation source.
    TGraphAsymmErrors* gr_r_var    [N_SOURCES][2][2][N_CENTBINS] = {};
    TGraphAsymmErrors* gr_slope_var[N_SOURCES][2][2][N_CENTBINS] = {};
    TGraphAsymmErrors* gr_inter_var[N_SOURCES][2][2][N_CENTBINS] = {};
    bool var_loaded[N_SOURCES][2] = {};

    // Raw q-bin scatter points (STEP 7) — nominal, and each source's variant(s). static:
    // large enough (~650 KB combined) that keeping them off the call stack is simple insurance.
    static ScatterPoint scatter_nom[2][N_CENTBINS][N_PTBINS_MAX][N_QBINS];
    static ScatterPoint scatter_var[N_SOURCES][2][2][N_CENTBINS][N_PTBINS_MAX][N_QBINS];

    // Nominal mean-q x-position per q-bin (STEP 10 only) -- see LoadMeanQNominal's comment.
    static double meanq_nom[2][N_CENTBINS][N_PTBINS_MAX][N_QBINS];
    static bool   meanq_nom_ok[2][N_CENTBINS][N_PTBINS_MAX][N_QBINS];

    std::mt19937_64 rng(2025);

    // --- Nominal ---
    std::cout << "\n--- Computing Pearson r/slope/intercept for dataset 'nominal' <- "
              << nominal_d0_file << " , " << nominal_chg_file << " ---\n";
    ComputeDatasetGraphs(nominal_d0_file, nominal_chg_file, isRun3, N_MC, USE_SYMMETRIC_INTERVAL, rng,
                          gr_r[0], gr_slope[0], gr_inter[0], scatter_nom);
    if (DO_STEP10_VN_VS_MEANQ_SYST) LoadMeanQNominal(nominal_d0_file, meanq_nom, meanq_nom_ok);
    {
        bool any = false;
        for (int vi = 0; vi < 2; ++vi)
            for (int ic = 0; ic < N_CENTBINS; ++ic)
                if (gr_r[0][vi][ic]) any = true;
        ds_loaded[0] = any;
        if (!any) std::cerr << "WARNING: no valid graphs computed for the nominal dataset.\n";
        for (int vi = 0; vi < 2; ++vi)
            for (int ic = 0; ic < N_CENTBINS; ++ic)
                for (TGraphAsymmErrors* g : {gr_r[0][vi][ic], gr_slope[0][vi][ic], gr_inter[0][vi][ic]}) {
                    if (!g) continue;
                    g->SetMarkerColor(g_col[ic]); g->SetLineColor(g_col[ic]);
                    g->SetMarkerStyle(g_marker[ic]); g->SetLineStyle(1);
                    g->SetMarkerSize(1.1); g->SetLineWidth(1);
                }
        std::cout << "  -> dataset 'nominal' ready.\n";
    }

    // --- Prompt-fraction shift table (used only by source slots with source_is_table_shift=true) ---
    double promptfrac_shift_v2[N_CENTBINS][N_PTBINS_V2] = {};
    double promptfrac_shift_v3[N_CENTBINS][N_PTBINS_V3] = {};
    LoadPromptFractionShiftTable(prompt_frac_shift_file, promptfrac_shift_v2, promptfrac_shift_v3);

    // --- Each of the 7 sources ---
    for (int isrc = 0; isrc < N_SOURCES; ++isrc) {
        bool is_dual = (source_d0_file[isrc][1] != nullptr);

        if (source_is_table_shift[isrc]) {
            std::cout << "\n--- Building source '" << source_label[isrc]
                      << "' from the prompt-fraction shift table (intercept-only) ---\n";
            BuildPromptFractionVariant(promptfrac_shift_v2, promptfrac_shift_v3, gr_r[0], gr_slope[0], gr_inter[0],
                                        gr_r_var[isrc][0], gr_slope_var[isrc][0], gr_inter_var[isrc][0]);
            BuildPromptFractionScatterVariant(promptfrac_shift_v2, promptfrac_shift_v3, scatter_nom, scatter_var[isrc][0]);
            bool any = false;
            for (int vi = 0; vi < 2; ++vi)
                for (int ic = 0; ic < N_CENTBINS; ++ic)
                    if (gr_inter_var[isrc][0][vi][ic]) any = true;
            var_loaded[isrc][0] = any;
            var_loaded[isrc][1] = false;
            if (!any) std::cerr << "WARNING: no valid graphs built for source '" << source_label[isrc] << "'.\n";
            for (int vi = 0; vi < 2; ++vi)
                for (int ic = 0; ic < N_CENTBINS; ++ic)
                    for (TGraphAsymmErrors* g : {gr_r_var[isrc][0][vi][ic], gr_slope_var[isrc][0][vi][ic], gr_inter_var[isrc][0][vi][ic]}) {
                        if (!g) continue;
                        g->SetMarkerColor(g_col[ic]); g->SetLineColor(g_col[ic]);
                        g->SetMarkerStyle(g_marker_open[ic]); g->SetLineStyle(2);
                        g->SetMarkerSize(1.1); g->SetLineWidth(1);
                    }
            std::cout << "  -> source '" << source_label[isrc] << "' ready.\n";
        } else {
        for (int iv = 0; iv < 2; ++iv) {
            if (!source_d0_file[isrc][iv]) continue; // no such variation for this source
            const char* d0f  = source_d0_file[isrc][iv];
            const char* chgf = source_chg_file[isrc][iv];
            const char* vlabel = is_dual ? SafeLabel(source_var_label[isrc][iv], (iv == 0) ? "Var1" : "Var2") : source_label[isrc];

            std::cout << "\n--- Computing Pearson r/slope/intercept for source '" << source_label[isrc]
                      << "' variation '" << vlabel << "' <- " << d0f << " , " << chgf << " ---\n";
            ComputeDatasetGraphs(d0f, chgf, isRun3, N_MC, USE_SYMMETRIC_INTERVAL, rng,
                                  gr_r_var[isrc][iv], gr_slope_var[isrc][iv], gr_inter_var[isrc][iv],
                                  scatter_var[isrc][iv]);

            bool any = false;
            for (int vi = 0; vi < 2; ++vi)
                for (int ic = 0; ic < N_CENTBINS; ++ic)
                    if (gr_r_var[isrc][iv][vi][ic]) any = true;
            var_loaded[isrc][iv] = any;
            if (!any) std::cerr << "WARNING: no valid graphs computed for source '" << source_label[isrc]
                                 << "' variation '" << vlabel << "'.\n";

            int lstyle = (iv == 0) ? 2 : 3; // dashed for variation 0, dotted for variation 1
            for (int vi = 0; vi < 2; ++vi)
                for (int ic = 0; ic < N_CENTBINS; ++ic)
                    for (TGraphAsymmErrors* g : {gr_r_var[isrc][iv][vi][ic], gr_slope_var[isrc][iv][vi][ic], gr_inter_var[isrc][iv][vi][ic]}) {
                        if (!g) continue;
                        g->SetMarkerColor(g_col[ic]); g->SetLineColor(g_col[ic]);
                        g->SetMarkerStyle(g_marker_open[ic]); g->SetLineStyle(lstyle);
                        g->SetMarkerSize(1.1); g->SetLineWidth(1);
                    }
            std::cout << "  -> source '" << source_label[isrc] << "' variation '" << vlabel << "' ready.\n";
        }
        } // end of the else branch (non-table-shift sources)

        if (is_dual && var_loaded[isrc][0] != var_loaded[isrc][1] && (var_loaded[isrc][0] || var_loaded[isrc][1]))
            std::cerr << "WARNING: source '" << source_label[isrc]
                      << "' has only one of its two variations loaded -- its adopted diff will fall back "
                         "to that single surviving variation everywhere.\n";

        // Synthesize the "adopted" (max-|diff|) graph per (harmonic, centrality) bin.
        for (int vi = 0; vi < 2; ++vi) {
            int n_pt = (vi == 0) ? N_PTBINS_V2 : N_PTBINS_V3;
            const double* pt_cen = (vi == 0) ? pt_cen_v2 : pt_cen_v3;
            const double* pt_elo = (vi == 0) ? pt_elo_v2 : pt_elo_v3;
            const double* pt_ehi = (vi == 0) ? pt_ehi_v2 : pt_ehi_v3;
            for (int ic = 0; ic < N_CENTBINS; ++ic) {
                gr_r[isrc+1][vi][ic] = BuildEffectiveGraph(
                    gr_r[0][vi][ic], gr_r_var[isrc][0][vi][ic], gr_r_var[isrc][1][vi][ic],
                    pt_cen, pt_elo, pt_ehi, n_pt);
                gr_slope[isrc+1][vi][ic] = BuildEffectiveGraph(
                    gr_slope[0][vi][ic], gr_slope_var[isrc][0][vi][ic], gr_slope_var[isrc][1][vi][ic],
                    pt_cen, pt_elo, pt_ehi, n_pt);
                gr_inter[isrc+1][vi][ic] = BuildEffectiveGraph(
                    gr_inter[0][vi][ic], gr_inter_var[isrc][0][vi][ic], gr_inter_var[isrc][1][vi][ic],
                    pt_cen, pt_elo, pt_ehi, n_pt);
            }
        }

        bool any_eff = false;
        for (int vi = 0; vi < 2; ++vi)
            for (int ic = 0; ic < N_CENTBINS; ++ic)
                if (gr_r[isrc+1][vi][ic]) any_eff = true;
        ds_loaded[isrc+1] = any_eff;
        if (!any_eff) {
            std::cerr << "WARNING: no valid adopted graphs for source '" << source_label[isrc]
                      << "' -- it will be skipped everywhere below.\n";
            continue;
        }
        for (int vi = 0; vi < 2; ++vi)
            for (int ic = 0; ic < N_CENTBINS; ++ic)
                for (TGraphAsymmErrors* g : {gr_r[isrc+1][vi][ic], gr_slope[isrc+1][vi][ic], gr_inter[isrc+1][vi][ic]}) {
                    if (!g) continue;
                    g->SetMarkerColor(g_col[ic]); g->SetLineColor(g_col[ic]);
                    g->SetMarkerStyle(g_marker_open[ic]); g->SetLineStyle(2);
                    g->SetMarkerSize(1.1); g->SetLineWidth(1);
                }
    }

    // =====================================================================
    // Per-q-bin total systematic uncertainty (quadratic sum across sources, independently in
    // x = charged vn and y = D0 vn) — computed once here since both STEP 7 (scatter+box plots)
    // and STEP 8 (Pearson r's own MC-based systematic) need it.
    // =====================================================================
    bool src_present[N_SOURCES];
    for (int isrc = 0; isrc < N_SOURCES; ++isrc) src_present[isrc] = ds_loaded[isrc + 1];

    static ScatterSyst scatter_syst_all[2][N_CENTBINS][N_PTBINS_MAX][N_QBINS];
    for (int vi = 0; vi < 2; ++vi) {
        int n_pt = (vi == 0) ? N_PTBINS_V2 : N_PTBINS_V3;
        for (int ic = 0; ic < N_CENTBINS; ++ic)
            for (int ip = 0; ip < n_pt; ++ip)
                ComputeScatterTotalSyst(vi, ic, ip, scatter_nom, scatter_var, src_present,
                                         scatter_syst_all[vi][ic][ip]);
    }

    // =====================================================================
    // STEP 1 — Overlay: nominal vs. each systematic source
    // =====================================================================
    if (DO_STEP1_OVERLAY) {
        std::cout << "\n=== STEP 1: overlay plots (nominal vs each source's variation(s)) ===\n";
        const double Y_LO_R = -1.0, Y_HI_R = 2.5;
        for (int isrc = 0; isrc < N_SOURCES; ++isrc) {
            bool is_dual = (source_d0_file[isrc][1] != nullptr);
            const char* lbl0 = is_dual ? SafeLabel(source_var_label[isrc][0], "Var1") : source_label[isrc];
            const char* lbl1 = is_dual ? SafeLabel(source_var_label[isrc][1], "Var2") : nullptr;
            if (!var_loaded[isrc][0] && !var_loaded[isrc][1]) continue;

            for (int ivn = 2; ivn <= 3; ++ivn) {
                int vi = ivn - 2;

                TString path_r = Form("%s/overlay_plots/v%d/pearson_v%d_overlay_%s.pdf",
                                       outbase.Data(), ivn, ivn, source_label[isrc]);
                DrawQuantityPanelCanvas(path_r, ivn, "r", Y_LO_R, Y_HI_R, true, 1.0,
                                         gr_r[0][vi], gr_r_var[isrc][0][vi], lbl0,
                                         nullptr, nullptr, gr_r_var[isrc][1][vi], lbl1);

                double SL_LO, SL_HI, IN_LO, IN_HI;
                if (ivn == 2) { SL_LO = -1.0; SL_HI = 4.0; IN_LO = -0.5; IN_HI = 0.5; }
                else          { SL_LO = -6.0; SL_HI = 6.0; IN_LO = -1.0; IN_HI = 2.0; }
                TString path_si = Form("%s/overlay_plots/v%d/slope_intercept_v%d_overlay_%s.pdf",
                                        outbase.Data(), ivn, ivn, source_label[isrc]);
                DrawSlopeInterceptPanelCanvas(path_si, ivn, SL_LO, SL_HI, IN_LO, IN_HI, 1.0, 0.0,
                                              gr_slope[0][vi], gr_inter[0][vi],
                                              gr_slope_var[isrc][0][vi], gr_inter_var[isrc][0][vi], lbl0,
                                              nullptr, nullptr, nullptr,
                                              gr_slope_var[isrc][1][vi], gr_inter_var[isrc][1][vi], lbl1);
            }
        }
        std::cout << "  -> " << outbase << "/overlay_plots/v{2,3}/  (up to 28 PDFs)\n";
    } else {
        std::cout << "\n=== STEP 1 skipped (DO_STEP1_OVERLAY = false) ===\n";
    }

    // =====================================================================
    // STEP 2 / STEP 3 — Absolute / relative difference vs pT
    // (share the same code path; `relative` selects which one)
    // =====================================================================
    auto RunDiffStep = [&](bool relative, const char* subdir, const char* calc_word)
    {
        std::cout << "\n=== " << calc_word << " difference plots ===\n";
        for (int ivn = 2; ivn <= 3; ++ivn) {
            int vi = ivn - 2;
            int n_pt = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
            const double* pt_cen = (ivn == 2) ? pt_cen_v2 : pt_cen_v3;
            const double* pt_elo = (ivn == 2) ? pt_elo_v2 : pt_elo_v3;
            const double* pt_ehi = (ivn == 2) ? pt_ehi_v2 : pt_ehi_v3;

            // Fixed, common y-range for all sources/centralities/harmonics — the previous
            // auto-scaled range could be blown up by a handful of near-threshold-nominal
            // outlier bins (where dividing by a tiny nominal value gives a huge relative %),
            // squashing every normal-sized point invisibly close to zero.
            // Slope/Intercept only -- Pearson r's systematic no longer comes from a per-source
            // diff (see STEP 8), so it's excluded here; see the header comment's IMPORTANT note.
            double Y_S, Y_B;
            const char *ytit_s, *ytit_b;
            if (relative) {
                Y_S = Y_B = 100.0;
                ytit_s = "#Delta Slope / Slope_{nom} (%)";
                ytit_b = "#Delta Intercept / Intercept_{nom} (%)";
            } else {
                Y_S = Y_B = 0.15;
                ytit_s = "#Delta Slope";
                ytit_b = "#Delta Intercept";
            }

            for (int isrc = 0; isrc < N_SOURCES; ++isrc) {
                if (!var_loaded[isrc][0] && !var_loaded[isrc][1]) continue;
                bool is_dual = (source_d0_file[isrc][1] != nullptr);

                TGraphAsymmErrors* gdiff_s[N_CENTBINS] = {}, *gdiff_s2[N_CENTBINS] = {};
                TGraphAsymmErrors* gdiff_b[N_CENTBINS] = {}, *gdiff_b2[N_CENTBINS] = {};
                for (int ic = 0; ic < N_CENTBINS; ++ic) {
                    gdiff_s[ic] = MakeDiffGraph(gr_slope_var[isrc][0][vi][ic], gr_slope[0][vi][ic], pt_cen, pt_elo, pt_ehi, n_pt,
                                                 relative, REL_ABS_THRESHOLD_SLOPE, g_col[ic], g_marker[ic]);
                    gdiff_b[ic] = MakeDiffGraph(gr_inter_var[isrc][0][vi][ic], gr_inter[0][vi][ic], pt_cen, pt_elo, pt_ehi, n_pt,
                                                 relative, REL_ABS_THRESHOLD_INTER, g_col[ic], g_marker[ic]);
                    if (is_dual) {
                        gdiff_s2[ic] = MakeDiffGraph(gr_slope_var[isrc][1][vi][ic], gr_slope[0][vi][ic], pt_cen, pt_elo, pt_ehi, n_pt,
                                                      relative, REL_ABS_THRESHOLD_SLOPE, g_col[ic], g_marker_open[ic], 2);
                        gdiff_b2[ic] = MakeDiffGraph(gr_inter_var[isrc][1][vi][ic], gr_inter[0][vi][ic], pt_cen, pt_elo, pt_ehi, n_pt,
                                                      relative, REL_ABS_THRESHOLD_INTER, g_col[ic], g_marker_open[ic], 2);
                    }
                }

                const char* corner = is_dual
                    ? Form("solid = %s, open = %s", SafeLabel(source_var_label[isrc][0], "Var1"),
                                                     SafeLabel(source_var_label[isrc][1], "Var2"))
                    : nullptr;

                TString path_si = Form("%s/%s/v%d/slope_intercept_v%d_%s_%s.pdf",
                    outbase.Data(), subdir, ivn, ivn, relative ? "reldiff" : "absdiff", source_label[isrc]);
                DrawSlopeInterceptPanelCanvas(path_si, ivn, -Y_S, Y_S, -Y_B, Y_B, 0.0, 0.0,
                                              gdiff_s, gdiff_b, gdiff_s2, gdiff_b2, nullptr,
                                              nullptr, nullptr, corner,
                                              nullptr, nullptr, nullptr, ytit_s, ytit_b);

                for (int ic = 0; ic < N_CENTBINS; ++ic) {
                    delete gdiff_s[ic]; delete gdiff_b[ic];
                    delete gdiff_s2[ic]; delete gdiff_b2[ic];
                }
            }
        }
        std::cout << "  -> " << outbase << "/" << subdir << "/v{2,3}/  (up to 14 PDFs; slope/intercept only)\n";
    };

    if (DO_STEP2_ABSDIFF) RunDiffStep(false, "absdiff_plots", "STEP 2: absolute");
    else std::cout << "\n=== STEP 2 skipped (DO_STEP2_ABSDIFF = false) ===\n";

    if (DO_STEP3_RELDIFF) RunDiffStep(true, "reldiff_plots", "STEP 3: relative");
    else std::cout << "\n=== STEP 3 skipped (DO_STEP3_RELDIFF = false) ===\n";

    // =====================================================================
    // STEP 4 — Summary: all 7 sources' relative differences + min/max envelope,
    // one canvas per (harmonic x quantity), one sub-pad per centrality.
    // =====================================================================
    if (DO_STEP4_SUMMARY) {
        std::cout << "\n=== STEP 4: summary (all-sources + envelope) plots ===\n";

        for (int ivn = 2; ivn <= 3; ++ivn) {
            int vi = ivn - 2;
            int n_pt = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
            const double* pt_cen = (ivn == 2) ? pt_cen_v2 : pt_cen_v3;
            const double* pt_elo = (ivn == 2) ? pt_elo_v2 : pt_elo_v3;
            const double* pt_ehi = (ivn == 2) ? pt_ehi_v2 : pt_ehi_v3;

            // Slope/Intercept only (iq starts at 1) -- Pearson r's systematic no longer comes
            // from a per-source envelope; see STEP 8 and the header comment's IMPORTANT note.
            for (int iq = 1; iq < 3; ++iq) {
                const char* qname = (iq == 1) ? "slope" : "intercept";
                const char* ytit  = Form("#Delta%s / %s_{nom} (%%)", qname, qname);
                double thr = (iq == 1) ? REL_ABS_THRESHOLD_SLOPE : REL_ABS_THRESHOLD_INTER;

                auto pick = [&](int ids, int ic) -> TGraphAsymmErrors* {
                    return (iq == 1) ? gr_slope[ids][vi][ic] : gr_inter[ids][vi][ic];
                };

                // Fixed y-range (relative difference, %) — a handful of near-threshold-
                // nominal outlier bins could otherwise blow up an auto-scaled range and
                // squash every normal-sized point invisibly close to zero.
                double Y = 100.0;

                TCanvas* cv = new TCanvas(Form("cv_summary_v%d_%s", ivn, qname), "", 1500, 900);
                cv->Divide(3, 2, 0.001, 0.001);

                // Drawn graphs must outlive the canvas's SaveAs() below — ROOT's automatic
                // primitive cleanup silently strips a pad's drawn object from the display the
                // moment that object is deleted, even before the canvas has ever been painted.
                // So these are only deleted after SaveAs()/delete cv, not inside the ic loop.
                TGraphAsymmErrors* gband_keep[N_CENTBINS] = {};
                TGraphAsymmErrors* gsrc_keep[N_CENTBINS][N_SOURCES] = {};

                for (int ic = 0; ic < N_CENTBINS; ++ic) {
                    cv->cd(ic + 1);
                    gPad->SetLogx();
                    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.03);
                    gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.10);

                    TH2F* hf = new TH2F(Form("hf_sum_%d_%s_%d", ivn, qname, ic), "",
                                         100, X_LO, X_HI, 100, -Y, Y);
                    hf->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                    hf->GetXaxis()->SetTitleSize(0.050); hf->GetXaxis()->SetLabelSize(0.045);
                    hf->GetXaxis()->SetMoreLogLabels(); hf->GetXaxis()->SetNoExponent();
                    hf->GetYaxis()->SetTitle(ytit);
                    hf->GetYaxis()->SetTitleSize(0.050); hf->GetYaxis()->SetTitleOffset(1.20);
                    hf->GetYaxis()->SetLabelSize(0.045);
                    hf->Draw("AXIS");

                    TLine* lref = new TLine(X_LO, 0.0, X_HI, 0.0);
                    lref->SetLineStyle(2); lref->SetLineColor(kGray + 2); lref->Draw();

                    TGraphAsymmErrors* srcArr[N_SOURCES] = {};
                    for (int isrc = 0; isrc < N_SOURCES; ++isrc)
                        if (ds_loaded[isrc + 1]) srcArr[isrc] = pick(isrc + 1, ic);

                    TGraphAsymmErrors* gband = MakeEnvelopeBand(pick(0, ic), srcArr, pt_cen, pt_elo, pt_ehi, n_pt, thr);
                    if (gband) gband->Draw("3 SAME");

                    TGraphAsymmErrors* gsrc[N_SOURCES] = {};
                    for (int isrc = 0; isrc < N_SOURCES; ++isrc) {
                        if (!ds_loaded[isrc + 1]) continue;
                        gsrc[isrc] = MakeDiffGraph(srcArr[isrc], pick(0, ic), pt_cen, pt_elo, pt_ehi, n_pt,
                                                    true, thr, g_src_col[isrc], g_src_marker[isrc]);
                        if (gsrc[isrc]) gsrc[isrc]->Draw("PL SAME");
                    }

                    TLatex ltx; ltx.SetNDC(); ltx.SetTextFont(42); ltx.SetTextSize(0.050); ltx.SetTextAlign(11);
                    ltx.DrawLatex(0.17, 0.93, Form("Cent: %d-%d%%", g_min_cent[ic], g_max_cent[ic]));

                    gband_keep[ic] = gband;
                    for (int isrc = 0; isrc < N_SOURCES; ++isrc) gsrc_keep[ic][isrc] = gsrc[isrc];
                }

                // 6th pad: shared legend
                cv->cd(6);
                gPad->Clear();
                TLegend* leg = new TLegend(0.05, 0.05, 0.95, 0.95);
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.060);
                leg->SetHeader("Systematic source", "C");
                for (int isrc = 0; isrc < N_SOURCES; ++isrc) {
                    if (!ds_loaded[isrc + 1]) continue;
                    TGraphAsymmErrors* dummy = new TGraphAsymmErrors(1);
                    dummy->SetPoint(0, -1, -1);
                    dummy->SetMarkerColor(g_src_col[isrc]); dummy->SetLineColor(g_src_col[isrc]);
                    dummy->SetMarkerStyle(g_src_marker[isrc]); dummy->SetMarkerSize(1.2);
                    leg->AddEntry(dummy, source_label[isrc], "pl");
                }
                TGraphAsymmErrors* dummyband = new TGraphAsymmErrors(1);
                dummyband->SetPoint(0, -1, -1);
                dummyband->SetFillColorAlpha(kGray + 1, 0.45); dummyband->SetLineColor(kGray + 2);
                leg->AddEntry(dummyband, "min-max envelope", "f");
                leg->Draw();
                TLatex ltxh; ltxh.SetNDC(); ltxh.SetTextFont(42); ltxh.SetTextSize(0.06); ltxh.SetTextAlign(11);
                ltxh.DrawLatex(0.05, 0.03, Form("v_{%d}, %s : PbPb 5.36 TeV, |y|<1", ivn, qname));

                TString path = Form("%s/summary_plots/reldiff_summary_%s_v%d.pdf", outbase.Data(), qname, ivn);
                cv->SaveAs(path);
                delete cv;

                for (int ic = 0; ic < N_CENTBINS; ++ic) {
                    delete gband_keep[ic];
                    for (int isrc = 0; isrc < N_SOURCES; ++isrc) delete gsrc_keep[ic][isrc];
                }
            }
        }
        std::cout << "  -> " << outbase << "/summary_plots/  (4 PDFs; slope/intercept only)\n";
    } else {
        std::cout << "\n=== STEP 4 skipped (DO_STEP4_SUMMARY = false) ===\n";
    }

    // =====================================================================
    // STEP 5 — Slope/Intercept total systematic uncertainty (quadratic sum of the 7 sources)
    // shown as a shaded box behind the nominal point + its statistical error. Pearson r is
    // excluded here -- its own systematic uncertainty is computed differently, in STEP 8.
    // =====================================================================
    if (DO_STEP5_TOTALUNC) {
        std::cout << "\n=== STEP 5: total-systematic plots (slope/intercept only) ===\n";
        for (int ivn = 2; ivn <= 3; ++ivn) {
            int vi = ivn - 2;
            int n_pt = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
            const double* pt_cen = (ivn == 2) ? pt_cen_v2 : pt_cen_v3;
            const double* pt_elo = (ivn == 2) ? pt_elo_v2 : pt_elo_v3;
            const double* pt_ehi = (ivn == 2) ? pt_ehi_v2 : pt_ehi_v3;

            TGraphAsymmErrors* gbox_s[N_CENTBINS] = {};
            TGraphAsymmErrors* gbox_b[N_CENTBINS] = {};
            for (int ic = 0; ic < N_CENTBINS; ++ic) {
                TGraphAsymmErrors* src_s[N_SOURCES] = {};
                TGraphAsymmErrors* src_b[N_SOURCES] = {};
                for (int isrc = 0; isrc < N_SOURCES; ++isrc) {
                    if (!ds_loaded[isrc + 1]) continue;
                    src_s[isrc] = gr_slope[isrc+1][vi][ic];
                    src_b[isrc] = gr_inter[isrc+1][vi][ic];
                }
                BinSyst bs_s[N_PTBINS_MAX] = {}, bs_b[N_PTBINS_MAX] = {};
                ComputeTotalSyst(gr_slope[0][vi][ic], src_s, pt_cen, n_pt, bs_s);
                ComputeTotalSyst(gr_inter[0][vi][ic], src_b, pt_cen, n_pt, bs_b);

                auto buildBox = [&](BinSyst* bs) -> TGraphAsymmErrors* {
                    std::vector<double> xv, yv, exlo, exhi, eylo, eyhi;
                    for (int ip = 0; ip < n_pt; ++ip) {
                        if (!bs[ip].nom_ok) continue;
                        xv.push_back(pt_cen[ip]);       yv.push_back(bs[ip].nom_y);
                        exlo.push_back(pt_elo[ip]);     exhi.push_back(pt_ehi[ip]);
                        eylo.push_back(bs[ip].syst_abs); eyhi.push_back(bs[ip].syst_abs);
                    }
                    if (xv.empty()) return nullptr;
                    return new TGraphAsymmErrors((int)xv.size(), xv.data(), yv.data(),
                                                  exlo.data(), exhi.data(), eylo.data(), eyhi.data());
                };
                gbox_s[ic] = buildBox(bs_s);
                gbox_b[ic] = buildBox(bs_b);
            }

            double SL_LO, SL_HI, IN_LO, IN_HI;
            if (ivn == 2) { SL_LO = -1.0; SL_HI = 4.0; IN_LO = -0.5; IN_HI = 0.5; }
            else          { SL_LO = -6.0; SL_HI = 6.0; IN_LO = -1.0; IN_HI = 2.0; }
            TString path_si = Form("%s/totalunc_plots/slope_intercept_v%d_with_totalsyst.pdf", outbase.Data(), ivn);
            DrawSlopeInterceptPanelCanvas(path_si, ivn, SL_LO, SL_HI, IN_LO, IN_HI, 1.0, 0.0,
                                          gr_slope[0][vi], gr_inter[0][vi],
                                          nullptr, nullptr, nullptr, gbox_s, gbox_b,
                                          nullptr);

            for (int ic = 0; ic < N_CENTBINS; ++ic) { delete gbox_s[ic]; delete gbox_b[ic]; }
        }
        std::cout << "  -> " << outbase << "/totalunc_plots/  (2 PDFs; slope/intercept only)\n";
    } else {
        std::cout << "\n=== STEP 5 skipped (DO_STEP5_TOTALUNC = false) ===\n";
    }

    // =====================================================================
    // STEP 6 — Text dump of the full numeric breakdown, for the analysis note.
    // =====================================================================
    if (DO_STEP6_TEXTDUMP) {
        std::cout << "\n=== STEP 6: text dump (slope/intercept only) ===\n";
        for (int ivn = 2; ivn <= 3; ++ivn) {
            int vi = ivn - 2;
            int n_pt = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
            const double* pt_edges = (ivn == 2) ? g_pt_edges_v2 : g_pt_edges_v3;
            const double* pt_cen   = (ivn == 2) ? pt_cen_v2 : pt_cen_v3;

            TString path = Form("%s/text_output/systematic_uncertainty_v%d.txt", outbase.Data(), ivn);
            FILE* fp = fopen(path.Data(), "w");
            if (!fp) { std::cerr << "ERROR: cannot open " << path << " for writing\n"; continue; }

            fprintf(fp, "Systematic uncertainty summary : v%d (Slope, Intercept)\n", ivn);
            fprintf(fp, "Total systematic = quadratic sum of each source's Adopted value (source_i - nominal) over the %d sources present in each bin.\n", N_SOURCES);
            fprintf(fp, "For a source with two variations, Adopted = whichever variation deviates more from nominal at that bin (VarB is 'n/a' for a one-variation source).\n");
            fprintf(fp, "'*' next to a relative value marks bins where |nominal| is below the omission threshold used in the reldiff/summary plots.\n");
            fprintf(fp, "Pearson r is NOT in this table -- its systematic uncertainty uses a different, MC-based method (see pearson_syst_plots/ and text_output/pearson_systematic_v%d.txt).\n\n", ivn);

            // Slope/Intercept only (iq starts at 1) -- see the note above.
            for (int iq = 1; iq < 3; ++iq) {
                const char* qname = (iq == 1) ? "Slope" : "Intercept";
                double thr = (iq == 1) ? REL_ABS_THRESHOLD_SLOPE : REL_ABS_THRESHOLD_INTER;
                auto pick = [&](int ids, int ic) -> TGraphAsymmErrors* {
                    return (iq == 1) ? gr_slope[ids][vi][ic] : gr_inter[ids][vi][ic];
                };
                auto pick_var = [&](int isrc, int iv, int ic) -> TGraphAsymmErrors* {
                    return (iq == 1) ? gr_slope_var[isrc][iv][vi][ic] : gr_inter_var[isrc][iv][vi][ic];
                };

                fprintf(fp, "==================================================================================\n");
                fprintf(fp, "Quantity: %s\n", qname);
                fprintf(fp, "==================================================================================\n");
                fprintf(fp, "%-10s %-12s %10s %10s %10s", "Cent", "pT(GeV/c)", "Nominal", "StatLo", "StatHi");
                for (int isrc = 0; isrc < N_SOURCES; ++isrc)
                    fprintf(fp, " %14s_VarA_abs %11s_VarA_rel%% %14s_VarB_abs %11s_VarB_rel%% %13s_Adopted_abs %10s_Adopted_rel%%",
                            source_label[isrc], source_label[isrc], source_label[isrc],
                            source_label[isrc], source_label[isrc], source_label[isrc]);
                fprintf(fp, " %14s %12s\n", "TotalSyst_abs", "TotalSyst_rel%");

                BinSyst bs_all[N_CENTBINS][N_PTBINS_MAX] = {};
                for (int ic = 0; ic < N_CENTBINS; ++ic) {
                    TGraphAsymmErrors* src[N_SOURCES] = {};
                    for (int isrc = 0; isrc < N_SOURCES; ++isrc)
                        if (ds_loaded[isrc + 1]) src[isrc] = pick(isrc + 1, ic);

                    BinSyst* bs = bs_all[ic];
                    ComputeTotalSyst(pick(0, ic), src, pt_cen, n_pt, bs);

                    for (int ip = 0; ip < n_pt; ++ip) {
                        if (!bs[ip].nom_ok) continue;
                        fprintf(fp, "%-10s %4.0f-%-7.0f %10.4f %10.4f %10.4f",
                                g_cen_name[ic], pt_edges[ip], pt_edges[ip + 1],
                                bs[ip].nom_y, bs[ip].nom_elo, bs[ip].nom_ehi);
                        for (int isrc = 0; isrc < N_SOURCES; ++isrc) {
                            bool below = std::fabs(bs[ip].nom_y) < thr;

                            for (int iv = 0; iv < 2; ++iv) {
                                PtVal pv = GetPtVal(pick_var(isrc, iv, ic), pt_cen[ip]);
                                if (!pv.ok) { fprintf(fp, " %20s %17s", "n/a", "n/a"); continue; }
                                double dv = pv.y - bs[ip].nom_y;
                                if (below) fprintf(fp, " %20.5f %16s", dv, "n/a*");
                                else       fprintf(fp, " %20.5f %17.2f", dv, dv / bs[ip].nom_y * 100.0);
                            }

                            if (!bs[ip].src_ok[isrc]) { fprintf(fp, " %19s %18s", "n/a", "n/a"); continue; }
                            double d = bs[ip].src_diff[isrc];
                            if (below) fprintf(fp, " %19.5f %18s", d, "n/a*");
                            else       fprintf(fp, " %19.5f %18.2f", d, d / bs[ip].nom_y * 100.0);
                        }
                        bool below_tot = std::fabs(bs[ip].nom_y) < thr;
                        fprintf(fp, " %14.5f", bs[ip].syst_abs);
                        if (below_tot) fprintf(fp, " %11s*\n", "n/a");
                        else           fprintf(fp, " %12.2f\n", bs[ip].syst_abs / std::fabs(bs[ip].nom_y) * 100.0);
                    }
                }
                fprintf(fp, "\n");

                // --- Min/max systematic-source contribution per (centrality, pT) bin ---
                fprintf(fp, "--- Min/Max systematic-source contribution per bin (%s) ---\n", qname);
                fprintf(fp, "'*' on a rel%% marks bins where |nominal| is below the omission threshold (same rule as above).\n");
                fprintf(fp, "%-10s %-12s %-16s %12s %9s   %-16s %12s %9s\n",
                        "Cent", "pT(GeV/c)", "Min_source", "Min_abs", "Min_rel%", "Max_source", "Max_abs", "Max_rel%");
                for (int ic = 0; ic < N_CENTBINS; ++ic) {
                    for (int ip = 0; ip < n_pt; ++ip) {
                        BinSyst& b = bs_all[ic][ip];
                        if (!b.nom_ok) continue;

                        int imin = -1, imax = -1;
                        double vmin = std::numeric_limits<double>::max(), vmax = -1.0;
                        for (int isrc = 0; isrc < N_SOURCES; ++isrc) {
                            if (!b.src_ok[isrc]) continue;
                            double ad = std::fabs(b.src_diff[isrc]);
                            if (ad < vmin) { vmin = ad; imin = isrc; }
                            if (ad > vmax) { vmax = ad; imax = isrc; }
                        }
                        if (imin < 0) {
                            fprintf(fp, "%-10s %4.0f-%-7.0f %-16s %12s %9s   %-16s %12s %9s\n",
                                    g_cen_name[ic], pt_edges[ip], pt_edges[ip + 1],
                                    "n/a", "n/a", "n/a", "n/a", "n/a", "n/a");
                            continue;
                        }

                        bool below = std::fabs(b.nom_y) < thr;
                        double dmin = b.src_diff[imin], dmax = b.src_diff[imax];
                        char minrel[16], maxrel[16];
                        if (below) {
                            snprintf(minrel, sizeof(minrel), "n/a*");
                            snprintf(maxrel, sizeof(maxrel), "n/a*");
                        } else {
                            snprintf(minrel, sizeof(minrel), "%.2f", dmin / b.nom_y * 100.0);
                            snprintf(maxrel, sizeof(maxrel), "%.2f", dmax / b.nom_y * 100.0);
                        }
                        fprintf(fp, "%-10s %4.0f-%-7.0f %-16s %12.5f %9s   %-16s %12.5f %9s\n",
                                g_cen_name[ic], pt_edges[ip], pt_edges[ip + 1],
                                source_label[imin], dmin, minrel,
                                source_label[imax], dmax, maxrel);
                    }
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
            std::cout << "  -> " << path << "\n";
        }
    } else {
        std::cout << "\n=== STEP 6 skipped (DO_STEP6_TEXTDUMP = false) ===\n";
    }

    // =====================================================================
    // STEP 7 — D0-vn-vs-charged-vn scatter plots with a total-systematic box per q-bin point
    // (same layout as plot_ESE_scatter_and_pearson_combined.C's own scatter plots).
    // =====================================================================
    if (DO_STEP7_SCATTER_SYST) {
        std::cout << "\n=== STEP 7: scatter plots with total-systematic box ===\n";
        gSystem->mkdir(Form("%s/scatter_syst_plots/v2", outbase.Data()), kTRUE);
        gSystem->mkdir(Form("%s/scatter_syst_plots/v3", outbase.Data()), kTRUE);

        // src_present[] and scatter_syst_all[] are computed once, above, right after the
        // per-source loop (STEP 8 needs them too).

        // Per-(harmonic, centrality) axis-range prescan, matching the nominal macro's own
        // technique (global min/max over all pT bins ± errors, with the same unreliable-
        // high-pT-bin exclusion), but also widened to include the systematic-box extents so
        // a box is never clipped by the frame.
        double cent_xlo[2][N_CENTBINS], cent_xhi[2][N_CENTBINS];
        double cent_ylo[2][N_CENTBINS], cent_yhi[2][N_CENTBINS];
        for (int vi = 0; vi < 2; ++vi) {
            int n_pt = (vi == 0) ? N_PTBINS_V2 : N_PTBINS_V3;
            const double* pt_e = (vi == 0) ? g_pt_edges_v2 : g_pt_edges_v3;
            for (int ic = 0; ic < N_CENTBINS; ++ic) {
                double gxmin = 1e9, gxmax = -1e9, gymin = 1e9, gymax = -1e9;
                for (int ip = 0; ip < n_pt; ++ip) {
                    if (vi == 0 && pt_e[ip] >= 30) continue; // skip unreliable high-pT bin from range scan
                    if (vi == 1 && pt_e[ip] >= 20) continue;
                    for (int iq = 0; iq < N_QBINS; ++iq) {
                        ScatterPoint& p = scatter_nom[vi][ic][ip][iq];
                        if (!p.ok) continue;
                        ScatterSyst& s = scatter_syst_all[vi][ic][ip][iq];
                        double hx = std::max(p.ex, p.x * SYST_BOX_XFRAC);
                        double hy = std::max(p.ey, s.ok ? s.sysy : 0.0);
                        gxmin = std::min(gxmin, p.x - hx); gxmax = std::max(gxmax, p.x + hx);
                        gymin = std::min(gymin, p.y - hy); gymax = std::max(gymax, p.y + hy);
                    }
                }
                if (gxmin > gxmax) { gxmin = 0.0; gxmax = 0.3; }
                if (gymin > gymax) { gymin = 0.0; gymax = 0.3; }
                double xmarg = 0.15 * std::max(gxmax - gxmin, 1e-6);
                double ymarg = 0.15 * std::max(gymax - gymin, 1e-6);
                cent_xlo[vi][ic] = std::max(0.0, gxmin - xmarg);
                cent_xhi[vi][ic] = gxmax + xmarg;
                cent_ylo[vi][ic] = gymin - ymarg;
                cent_yhi[vi][ic] = gymax + ymarg;
            }
        }

        int nplots = 0;
        for (int vi = 0; vi < 2; ++vi) {
            int ivn = vi + 2;
            int n_pt = (vi == 0) ? N_PTBINS_V2 : N_PTBINS_V3;
            const double* pt_e = (vi == 0) ? g_pt_edges_v2 : g_pt_edges_v3;
            for (int ic = 0; ic < N_CENTBINS; ++ic) {
                for (int ip = 0; ip < n_pt; ++ip) {
                    bool any = false;
                    for (int iq = 0; iq < N_QBINS; ++iq)
                        if (scatter_nom[vi][ic][ip][iq].ok) any = true;
                    if (!any) continue;

                    TString path = Form("%s/scatter_syst_plots/v%d/scatter_v%d_%s_pT%dto%d.pdf",
                        outbase.Data(), ivn, ivn, g_cen_name[ic], (int)pt_e[ip], (int)pt_e[ip + 1]);
                    DrawScatterWithSystBox(path, ivn, ic, pt_e[ip], pt_e[ip + 1],
                                            scatter_nom[vi][ic][ip], scatter_syst_all[vi][ic][ip],
                                            cent_xlo[vi][ic], cent_xhi[vi][ic],
                                            cent_ylo[vi][ic], cent_yhi[vi][ic]);
                    ++nplots;
                }
            }
        }

        // Merged scatter+syst PDFs — tile the already-saved individual PDFs with Python/
        // PyMuPDF, exactly like plot_ESE_scatter_and_pearson_combined.C's own merged scatter
        // PDFs (one page per centrality, 3-column grid, pixel-identical to the individual
        // plots — no ROOT re-draw). Reuses that same script verbatim, just pointed at
        // scatter_syst_plots/ instead of scatter_plots/, since the filenames follow the
        // identical "scatter_v<ivn>_cent<lo>to<hi>_pT<lo>to<hi>.pdf" convention.
        std::cout << "  Generating merged scatter PDFs (Python/PyMuPDF)...\n";
        {
            TString py_script = Form("%s/_merge_scatter_syst_tmp.py", outbase.Data());
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
                    "def merge(inp, out, cols=3):\n"
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
                    "merge('%s/scatter_syst_plots/v2','%s/scatter_syst_plots/v2/merged_scatter_syst_v2.pdf')\n",
                    outbase.Data(), outbase.Data());
                fprintf(fpy,
                    "merge('%s/scatter_syst_plots/v3','%s/scatter_syst_plots/v3/merged_scatter_syst_v3.pdf')\n",
                    outbase.Data(), outbase.Data());
                fclose(fpy);

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
                              << outbase << "/scatter_syst_plots/\n";
                } else {
                    std::cout << "  Using Python: " << py_exe << "\n";
                    int ret = gSystem->Exec(Form("%s '%s'", py_exe.Data(), py_script.Data()));
                    if (ret != 0)
                        std::cerr << "WARNING: Python merge script failed (exit " << ret
                                  << "). Individual PDFs are still in "
                                  << outbase << "/scatter_syst_plots/\n";
                }
                gSystem->Unlink(py_script.Data());
            }
        }
        std::cout << "  -> " << outbase << "/scatter_syst_plots/v{2,3}/  (" << nplots << " individual PDFs"
                     " + merged_scatter_syst_v{2,3}.pdf)\n";
    } else {
        std::cout << "\n=== STEP 7 skipped (DO_STEP7_SCATTER_SYST = false) ===\n";
    }

    // =====================================================================
    // STEP 8 — Pearson r's own systematic uncertainty via Gaussian MC resampling (NOT a
    // per-source quadratic sum -- see the header comment's IMPORTANT note and
    // ComputePearsonSystUncertainty's own comment for the physics reasoning). Each D0 v_n
    // value is Gaussian-smeared with sigma = its total per-q-bin systematic uncertainty
    // (scatter_syst_all, from STEP 7's ComputeScatterTotalSyst -- computed above regardless
    // of DO_STEP7_SCATTER_SYST), while the weights used in the Pearson-r calculation stay
    // statistical. Kept in its own directory/table, separate from totalunc_plots/ and
    // text_output/systematic_uncertainty_v{2,3}.txt (Slope/Intercept's quadrature-sum method).
    // =====================================================================
    if (DO_STEP8_PEARSON_SYST) {
        std::cout << "\n=== STEP 8: Pearson r systematic uncertainty (Gaussian MC resampling) ===\n";
        const double Y_LO_R = -1.0, Y_HI_R = 2.5;
        gSystem->mkdir(Form("%s/mc_r_syst_plots/v2", outbase.Data()), kTRUE);
        gSystem->mkdir(Form("%s/mc_r_syst_plots/v3", outbase.Data()), kTRUE);

        for (int ivn = 2; ivn <= 3; ++ivn) {
            int vi = ivn - 2;
            int n_pt = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
            const double* pt_edges = (ivn == 2) ? g_pt_edges_v2 : g_pt_edges_v3;
            const double* pt_cen = (ivn == 2) ? pt_cen_v2 : pt_cen_v3;
            const double* pt_elo = (ivn == 2) ? pt_elo_v2 : pt_elo_v3;
            const double* pt_ehi = (ivn == 2) ? pt_ehi_v2 : pt_ehi_v3;

            TString path_txt = Form("%s/text_output/pearson_systematic_v%d.txt", outbase.Data(), ivn);
            FILE* fp = fopen(path_txt.Data(), "w");
            if (!fp) std::cerr << "ERROR: cannot open " << path_txt << " for writing\n";
            if (fp) {
                fprintf(fp, "Pearson r systematic uncertainty summary : v%d\n", ivn);
                fprintf(fp, "Method: Gaussian MC resampling (N_MC=%d) of each D0 v_n point with sigma = its total\n", N_MC);
                fprintf(fp, "per-q-bin systematic uncertainty (quadratic sum across the %d sources); weights used in the\n", N_SOURCES);
                fprintf(fp, "Pearson-r calculation stay statistical. NOT a per-source quadratic sum (r is bounded to\n");
                fprintf(fp, "[-1,1] by definition, so that method -- used for Slope/Intercept -- is not valid for r).\n\n");
                fprintf(fp, "%-10s %-12s %10s %10s %10s %10s %10s\n",
                        "Cent", "pT(GeV/c)", "Nominal", "StatLo", "StatHi", "SystLo", "SystHi");
            }

            TGraphAsymmErrors* gbox_r[N_CENTBINS] = {};
            for (int ic = 0; ic < N_CENTBINS; ++ic) {
                std::vector<double> xv, yv, exlo, exhi, eylo, eyhi;

                // r_mc distribution histograms (one per pT bin), same layout/style as the
                // nominal macro's own mc_r_plots -- built here so STEP 8 doesn't need to
                // recompute the MC resampling a second time.
                TH1D*  h_mc_arr[N_PTBINS_MAX] = {};
                bool   r_ok_ip [N_PTBINS_MAX] = {};
                double r_val_ip[N_PTBINS_MAX] = {}, r_elo_ip[N_PTBINS_MAX] = {}, r_ehi_ip[N_PTBINS_MAX] = {};

                for (int ip = 0; ip < n_pt; ++ip) {
                    double r_meas, r_elo, r_ehi; bool ok;
                    std::vector<double> r_mc;
                    ComputePearsonSystUncertainty(scatter_nom[vi][ic][ip], scatter_syst_all[vi][ic][ip],
                                                   N_MC, USE_SYMMETRIC_INTERVAL, rng, r_meas, r_elo, r_ehi, ok, &r_mc);
                    if (!ok) continue;
                    xv.push_back(pt_cen[ip]);   yv.push_back(r_meas);
                    exlo.push_back(pt_elo[ip]); exhi.push_back(pt_ehi[ip]);
                    eylo.push_back(r_elo);      eyhi.push_back(r_ehi);

                    if (fp) {
                        PtVal pn = GetPtVal(gr_r[0][vi][ic], pt_cen[ip]);
                        fprintf(fp, "%-10s %4.0f-%-7.0f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                                g_cen_name[ic], pt_edges[ip], pt_edges[ip + 1],
                                r_meas, pn.ok ? pn.elo : 0.0, pn.ok ? pn.ehi : 0.0, r_elo, r_ehi);
                    }

                    TH1D* h_mc = new TH1D(
                        Form("h_mc_syst_v%d_%s_pT%dto%d", ivn, g_cen_name[ic], (int)pt_edges[ip], (int)pt_edges[ip + 1]),
                        ";r (MC);counts", 200, -1.0, 1.0);
                    h_mc->SetDirectory(nullptr);
                    for (double rv : r_mc) h_mc->Fill(rv);
                    h_mc->SetFillColor(kBlue - 9);
                    h_mc->SetFillStyle(1001);
                    h_mc->SetLineColor(kBlue + 1);
                    h_mc_arr[ip] = h_mc;
                    r_ok_ip[ip] = true; r_val_ip[ip] = r_meas; r_elo_ip[ip] = r_elo; r_ehi_ip[ip] = r_ehi;
                }
                if (!xv.empty())
                    gbox_r[ic] = new TGraphAsymmErrors((int)xv.size(), xv.data(), yv.data(),
                                                         exlo.data(), exhi.data(), eylo.data(), eyhi.data());

                // --- Combined mc_r canvas for this centrality: one PDF, 3-column grid of
                // pT-bin sub-pads, exactly like the nominal macro's mc_r_plots. ---
                bool any_mc = false;
                for (int ip = 0; ip < n_pt; ++ip) if (h_mc_arr[ip]) any_mc = true;
                if (any_mc) {
                    const int N_COLS_MC = 3;
                    const int N_ROWS_MC = (n_pt + N_COLS_MC - 1) / N_COLS_MC;
                    TCanvas* cv_mc = new TCanvas(Form("cv_mc_syst_v%d_%s", ivn, g_cen_name[ic]), "",
                                                  1600, N_ROWS_MC * 400);
                    cv_mc->Divide(N_COLS_MC, N_ROWS_MC, 0.001, 0.001);

                    for (int ip = 0; ip < n_pt; ++ip) {
                        cv_mc->cd(ip + 1);
                        gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03);
                        gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.10);

                        TH1D* hh = h_mc_arr[ip];
                        if (!hh) {
                            TH2F* hblank = new TH2F(Form("hblank_syst_%d_%d_%d", ivn, ic, ip), "",
                                                     10, -1.0, 1.0, 10, 0, 10);
                            hblank->GetXaxis()->SetTitle("r");
                            hblank->Draw();
                            continue;
                        }

                        hh->GetXaxis()->SetTitle("r");
                        hh->GetXaxis()->SetTitleSize(0.070); hh->GetXaxis()->SetTitleOffset(0.90);
                        hh->GetXaxis()->SetLabelSize(0.060); hh->GetXaxis()->SetRangeUser(-1.0, 1.0);
                        hh->GetYaxis()->SetTitle("Counts");
                        hh->GetYaxis()->SetTitleSize(0.070); hh->GetYaxis()->SetTitleOffset(0.85);
                        hh->GetYaxis()->SetLabelSize(0.060); hh->GetYaxis()->SetMaxDigits(3);
                        gStyle->SetStripDecimals(kFALSE);
                        hh->Draw("HIST");

                        double ymax = hh->GetMaximum();
                        if (r_ok_ip[ip]) {
                            TLine* lr = new TLine(r_val_ip[ip], 0, r_val_ip[ip], ymax);
                            lr->SetLineColor(kRed); lr->SetLineWidth(2); lr->Draw();
                            double lo = r_val_ip[ip] - r_elo_ip[ip], hi = r_val_ip[ip] + r_ehi_ip[ip];
                            TLine* ll_lo = new TLine(lo, 0, lo, ymax * 0.65);
                            ll_lo->SetLineColor(kRed); ll_lo->SetLineStyle(2); ll_lo->SetLineWidth(1); ll_lo->Draw();
                            TLine* ll_hi = new TLine(hi, 0, hi, ymax * 0.65);
                            ll_hi->SetLineColor(kRed); ll_hi->SetLineStyle(2); ll_hi->SetLineWidth(1); ll_hi->Draw();
                        }

                        TLegend* leg_mc = new TLegend(0.11, 0.73, 0.48, 0.9);
                        leg_mc->SetBorderSize(0); leg_mc->SetFillStyle(0);
                        leg_mc->SetTextFont(42); leg_mc->SetTextSize(0.072);
                        leg_mc->AddEntry((TObject*)nullptr, Form("Cent: %d-%d%%", g_min_cent[ic], g_max_cent[ic]), "");
                        leg_mc->AddEntry((TObject*)nullptr, Form("%.0f< p_{T} < %.0f GeV/c", pt_edges[ip], pt_edges[ip + 1]), "");
                        leg_mc->Draw();
                    }

                    cv_mc->cd();
                    TLatex ltx_cen_mc; ltx_cen_mc.SetNDC(); ltx_cen_mc.SetTextFont(42);
                    ltx_cen_mc.SetTextSize(0.022); ltx_cen_mc.SetTextAlign(11);
                    ltx_cen_mc.DrawLatex(0.01, 0.003, Form("v_{%d}  |  Centrality %d-%d%%  |  PbPb 5.36 TeV (systematic MC)",
                                                            ivn, g_min_cent[ic], g_max_cent[ic]));

                    cv_mc->SaveAs(Form("%s/mc_r_syst_plots/v%d/mc_r_syst_v%d_%s.pdf",
                                        outbase.Data(), ivn, ivn, g_cen_name[ic]));
                    delete cv_mc;
                }
                for (int ip = 0; ip < n_pt; ++ip) delete h_mc_arr[ip];
            }
            if (fp) { fclose(fp); std::cout << "  -> " << path_txt << "\n"; }

            TString path_r = Form("%s/pearson_syst_plots/pearson_v%d_with_totalsyst.pdf", outbase.Data(), ivn);
            DrawQuantityPanelCanvas(path_r, ivn, "r", Y_LO_R, Y_HI_R, true, 1.0,
                                     gr_r[0][vi], nullptr, nullptr, gbox_r,
                                     nullptr);

            for (int ic = 0; ic < N_CENTBINS; ++ic) delete gbox_r[ic];
        }
        std::cout << "  -> " << outbase << "/pearson_syst_plots/  (2 PDFs)\n";
        std::cout << "  -> " << outbase << "/mc_r_syst_plots/v{2,3}/  (r-distribution diagnostic plots, one PDF per centrality)\n";

        // Merged mc_r_syst PDFs — concatenate the already-saved per-centrality PDFs (each
        // already a full pT-bin grid) into a single multi-page PDF, one page per centrality,
        // via Python/PyMuPDF. Reuses the exact same "concat" approach (not the 3-column tiling
        // used for STEP 7's scatter PDFs) as plot_ESE_vnVsqn_scatter_and_pearson_combined.C's
        // own merged mc_r_plots, since here each individual PDF is already a full grid.
        std::cout << "  Generating merged MC r PDFs (Python/PyMuPDF)...\n";
        {
            TString py_script = Form("%s/_merge_mcr_syst_tmp.py", outbase.Data());
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
                    "concat('%s/mc_r_syst_plots/v2','%s/mc_r_syst_plots/v2/merged_mc_r_syst_v2.pdf','mc_r_syst_v2')\n",
                    outbase.Data(), outbase.Data());
                fprintf(fpy,
                    "concat('%s/mc_r_syst_plots/v3','%s/mc_r_syst_plots/v3/merged_mc_r_syst_v3.pdf','mc_r_syst_v3')\n",
                    outbase.Data(), outbase.Data());
                fclose(fpy);

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
                              << "  Merged MC r PDFs skipped; individual PDFs are still in "
                              << outbase << "/mc_r_syst_plots/\n";
                } else {
                    std::cout << "  Using Python: " << py_exe << "\n";
                    int ret = gSystem->Exec(Form("%s '%s'", py_exe.Data(), py_script.Data()));
                    if (ret != 0)
                        std::cerr << "WARNING: Python merge script failed (exit " << ret
                                  << "). Individual PDFs are still in "
                                  << outbase << "/mc_r_syst_plots/\n";
                }
                gSystem->Unlink(py_script.Data());
            }
        }
        std::cout << "  -> " << outbase << "/mc_r_syst_plots/v{2,3}/merged_mc_r_syst_v{2,3}.pdf\n";
    } else {
        std::cout << "\n=== STEP 8 skipped (DO_STEP8_PEARSON_SYST = false) ===\n";
    }

    // =====================================================================
    // Shared helper for STEP 9 / STEP 10: concatenates the already-saved per-centrality PDFs
    // (each already a full pT-bin grid) into one multi-page PDF per harmonic, one page per
    // centrality -- identical technique to STEP 8's merged mc_r_syst PDFs.
    // =====================================================================
    auto RunConcatMerge = [&](const char* subdir, const char* base_name) {
        std::cout << "  Generating merged " << base_name << " PDFs (Python/PyMuPDF)...\n";
        TString py_script = Form("%s/_merge_%s_tmp.py", outbase.Data(), base_name);
        FILE* fpy = fopen(py_script.Data(), "w");
        if (!fpy) {
            std::cerr << "WARNING: Cannot write temp Python script; merged " << base_name << " PDFs skipped.\n";
            return;
        }
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
        for (int ivn = 2; ivn <= 3; ++ivn) {
            fprintf(fpy,
                "concat('%s/%s/v%d','%s/%s/v%d/merged_%s_v%d.pdf','%s_v%d')\n",
                outbase.Data(), subdir, ivn, outbase.Data(), subdir, ivn, base_name, ivn, base_name, ivn);
        }
        fclose(fpy);

        TString py_exe;
        {
            const char* cands[] = {
                "/usr/local/bin/python3.11", "/opt/homebrew/bin/python3",
                "/usr/bin/python3", "python3", "python", nullptr
            };
            for (int k = 0; cands[k] && py_exe.IsNull(); ++k) {
                TString probe = Form("%s -c 'import fitz' 2>/dev/null", cands[k]);
                if (gSystem->Exec(probe.Data()) == 0) py_exe = cands[k];
            }
        }
        if (py_exe.IsNull()) {
            std::cerr << "WARNING: No Python interpreter with PyMuPDF found.\n"
                      << "  Install with:  pip install pymupdf\n"
                      << "  Merged " << base_name << " PDFs skipped; individual PDFs are still in "
                      << outbase << "/" << subdir << "/\n";
        } else {
            std::cout << "  Using Python: " << py_exe << "\n";
            int ret = gSystem->Exec(Form("%s '%s'", py_exe.Data(), py_script.Data()));
            if (ret != 0)
                std::cerr << "WARNING: Python merge script failed (exit " << ret
                          << "). Individual PDFs are still in " << outbase << "/" << subdir << "/\n";
        }
        gSystem->Unlink(py_script.Data());
    };

    // =====================================================================
    // STEP 9 — D0 v_n vs q-bin index, with a total-systematic box per point (fixed absolute
    // x-half-width -- the q-bin index itself has no physical systematic uncertainty; the
    // box's y-half-height is the real propagated D0-v_n systematic, reused from STEP 7).
    // =====================================================================
    if (DO_STEP9_VN_VS_QBIN_SYST) {
        std::cout << "\n=== STEP 9: D0 v_n vs q-bin index (with total-syst box) ===\n";
        gSystem->mkdir(Form("%s/vn_vs_qbin_syst_plots/v2", outbase.Data()), kTRUE);
        gSystem->mkdir(Form("%s/vn_vs_qbin_syst_plots/v3", outbase.Data()), kTRUE);

        for (int vi = 0; vi < 2; ++vi) {
            int ivn = vi + 2;
            int n_pt = (vi == 0) ? N_PTBINS_V2 : N_PTBINS_V3;
            const double* pt_e = (vi == 0) ? g_pt_edges_v2 : g_pt_edges_v3;
            for (int ic = 0; ic < N_CENTBINS; ++ic) {
                bool any = false;
                for (int ip = 0; ip < n_pt && !any; ++ip)
                    for (int iq = 0; iq < N_QBINS; ++iq)
                        if (scatter_nom[vi][ic][ip][iq].ok) { any = true; break; }
                if (!any) continue;

                TString path = Form("%s/vn_vs_qbin_syst_plots/v%d/vn_vs_qbin_syst_v%d_%s.pdf",
                                     outbase.Data(), ivn, ivn, g_cen_name[ic]);
                DrawVnVsQAxisSystGrid(path, ivn, ic, n_pt, pt_e,
                                       scatter_nom[vi][ic], scatter_syst_all[vi][ic],
                                       SYST_BOX_FIXED_HALFWIDTH,
                                       "q-bin index");
            }
        }

        RunConcatMerge("vn_vs_qbin_syst_plots", "vn_vs_qbin_syst");
        std::cout << "  -> " << outbase << "/vn_vs_qbin_syst_plots/v{2,3}/merged_vn_vs_qbin_syst_v{2,3}.pdf\n";
    } else {
        std::cout << "\n=== STEP 9 skipped (DO_STEP9_VN_VS_QBIN_SYST = false) ===\n";
    }

    // =====================================================================
    // STEP 10 — D0 v_n vs q-bin index (exact same points, positions and box as STEP 9), with
    // the tick labels replaced by each bin's actual mean-q value instead of its plain index.
    // =====================================================================
    if (DO_STEP10_VN_VS_MEANQ_SYST) {
        std::cout << "\n=== STEP 10: D0 v_n vs mean-q (with total-syst box) ===\n";
        gSystem->mkdir(Form("%s/vn_vs_meanq_syst_plots/v2", outbase.Data()), kTRUE);
        gSystem->mkdir(Form("%s/vn_vs_meanq_syst_plots/v3", outbase.Data()), kTRUE);

        for (int vi = 0; vi < 2; ++vi) {
            int ivn = vi + 2;
            int n_pt = (vi == 0) ? N_PTBINS_V2 : N_PTBINS_V3;
            const double* pt_e = (vi == 0) ? g_pt_edges_v2 : g_pt_edges_v3;
            for (int ic = 0; ic < N_CENTBINS; ++ic) {
                bool any = false;
                for (int ip = 0; ip < n_pt && !any; ++ip)
                    for (int iq = 0; iq < N_QBINS; ++iq)
                        if (scatter_nom[vi][ic][ip][iq].ok) { any = true; break; }
                if (!any) continue;

                TString path = Form("%s/vn_vs_meanq_syst_plots/v%d/vn_vs_meanq_syst_v%d_%s.pdf",
                                     outbase.Data(), ivn, ivn, g_cen_name[ic]);
                DrawVnVsQAxisSystGrid(path, ivn, ic, n_pt, pt_e,
                                       scatter_nom[vi][ic], scatter_syst_all[vi][ic],
                                       SYST_BOX_FIXED_HALFWIDTH,
                                       Form("<q_{%d}>", ivn),
                                       meanq_nom[vi][ic], meanq_nom_ok[vi][ic]);
            }
        }

        RunConcatMerge("vn_vs_meanq_syst_plots", "vn_vs_meanq_syst");
        std::cout << "  -> " << outbase << "/vn_vs_meanq_syst_plots/v{2,3}/merged_vn_vs_meanq_syst_v{2,3}.pdf\n";
    } else {
        std::cout << "\n=== STEP 10 skipped (DO_STEP10_VN_VS_MEANQ_SYST = false) ===\n";
    }

    std::cout << "\n=== Done ===\n  Output directory -> " << outbase << "/\n";
}
