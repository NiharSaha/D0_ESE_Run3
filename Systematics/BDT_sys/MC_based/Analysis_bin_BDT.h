#pragma once
// =============================================================================
//  Analysis_bin_BDT.h
//  Shared configuration for the BDT systematic study on prompt-D0 MC.
//
//  Self-contained: no dependency on the auto-generated vnbinning / quantile
//  headers used by the data analysis.  Centrality and pT bin edges are kept
//  identical to Analysis_bin.h so the MC result can be compared 1:1 with data.
//
//  Design (mirrors the charged-particle analysis charge_vn/for_20Qbin/ and the
//  MC treatment in Systematics/Prompt_frac/getDCA_fromMC.C):
//    * resolution is computed in a SEPARATE step (Calculate_Resolution_MC.C),
//      q-inclusive and centrality-inclusive over the WHOLE MC sample
//      (there is no meaningful centrality in MC -> pool everything for stats);
//    * flow_MC_BDT_sys.C then reads the resolution file and does the whole
//      vn calculation in one pass, filling TProfiles (no intermediate ntuple);
//    * the BDT cut is applied exactly as in getDCA_fromMC.C / getEfficiency.C:
//      the candidate's true event centrality is IGNORED; instead the result is
//      produced once per analysis centrality class, applying that class's BDT
//      cut (evaluated at 2 x the class bin-center) to the full inclusive sample.
// =============================================================================
#include <cmath>
#include <cstdlib>
#include "TString.h"

// ---------------------------------------------------------------------------
//  Centrality classes
//  Only used to pick which BDT cut to apply -- NOT to bin events.  The vn
//  numerator and the resolution are fully centrality-inclusive.
// ---------------------------------------------------------------------------
static const int N_CENTBINS = 5;
[[maybe_unused]] static const char *cen_name[N_CENTBINS] = {
    "cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"};

// centrality argument passed to BDTHandler::getBDTCut for each class:
// 2 x bin-center, identical to getDCA_fromMC.C (centApprox) and
// getEfficiency.C (2*cent_case).
[[maybe_unused]] static const int cent_bdt_arg[N_CENTBINS] = {10, 30, 50, 70, 90};

// ---------------------------------------------------------------------------
//  pT bins  (same edges as Analysis_bin.h)
// ---------------------------------------------------------------------------
static const int N_PTBINS_V2 = 9;   // 2,3,4,5,6,8,10,15,30,100
static const double pt_edges_v2[N_PTBINS_V2 + 1] = {2, 3, 4, 5, 6, 8, 10, 15, 30, 100};

static const int N_PTBINS_V3 = 7;   // 2,4,6,8,10,20,50,100
static const double pt_edges_v3[N_PTBINS_V3 + 1] = {2, 4, 6, 8, 10, 20, 50, 100};

// ---------------------------------------------------------------------------
//  Analysis acceptance
// ---------------------------------------------------------------------------
static const double MIN_PT_ANA = 2.0;
static const double MAX_PT_ANA = 100.0;
static const double MAX_Y_ANA  = 1.0;   // final measurement fiducial region (|y|<1)

// ---------------------------------------------------------------------------
//  MC event weighting
//  Cross-section stitching weights for the nested "Dpt > X GeV" generator-level
//  filtered prompt-D0 samples that make up inputFiles_Prompt_MC_Sept5_noDpt8.txt.
//  Numbers copied verbatim from Systematics/AccxEff/Eff_factor/getEfficiency.C
//  (updated there 2026-09-24) so the MC event mixture is treated identically to
//  the Acc x Eff analysis:
//    * thr=10..80 xsec/filterEff from GenXsecAnalyzer; thr=0,1 still McM values;
//    * thr=8 dropped -- its GEN fragment is identical to thr=1 (MomMinPt=1.0),
//      so its files are also removed from the input list.  With the old thr=8
//      numbers its L_eff was ~9x (L_0+L_1), so events with p_gen in [8,10) were
//      effectively represented by that one sample only -> stats hole at 8-10 GeV.
//
//  The SAME weight is applied to the resolution ingredients and to the flow
//  numerator, so it cancels in the with-BDT / without-BDT systematic ratio and
//  only fixes the relative normalisation of the nested pT samples.
//  Set USE_MC_WEIGHT = false to fall back to a plain (unweighted) pooling,
//  like the data / charged-particle resolution step.
//  (getDCA_fromMC.C fills raw counts with no stitch weight; toggle this off to
//   match that convention.)
// ---------------------------------------------------------------------------
static const bool USE_MC_WEIGHT = true;

// previous (McM-only, 2026-09-06, incl. thr=8):
// N=9; thr = {0, 1, 8, 10, 20, 30, 40, 60, 80}
// xsec_mb   = {4.63e3, 3.67e3, 3.08e1, 1.15e1, 7.1e-1, 1.15e-1, 3.5e-2, 4.57e-3, 1.11e-3}
// filterEff = {7.43e-3, 5.9e-3, 5.8e-3, 3.9e-4, 3.9e-4, 6.3e-5, 2.3e-4, 1.42e-4, 1.11e-4}
// Nevts     = {1.0e6, 2.0e6, 1.5e5, 1.0e5, 5.0e4, 2.0e4, 2.0e4, 2.0e4, 2.0e4}
static const int N_DPT_SAMPLES = 8;
static const double Dpt_threshold[N_DPT_SAMPLES] = {0, 1, 10, 20, 30, 40, 60, 80};
static const double Dpt_xsec_mb[N_DPT_SAMPLES]   = {4.63e3, 3.67e3, 29.36, 1.809, 1.810, 0.1539, 0.03214, 0.01002};
static const double Dpt_filterEff[N_DPT_SAMPLES] = {7.43e-3, 5.9e-3, 3.400e-4, 3.800e-4, 5.500e-5, 1.200e-4, 1.700e-4, 7.000e-5};
static const double Dpt_Nevts[N_DPT_SAMPLES]     = {1.0e6, 2.0e6, 1.0e5, 5.0e4, 2.0e4, 2.0e4, 2.0e4, 2.0e4};

inline double getStitchWeight(double p_gen)
{
  static double Leff[N_DPT_SAMPLES];
  static bool init = false;
  if (!init)
  {
    for (int i = 0; i < N_DPT_SAMPLES; i++)
      Leff[i] = Dpt_Nevts[i] / (Dpt_xsec_mb[i] * Dpt_filterEff[i]);
    init = true;
  }
  double sumL = 0.0;
  for (int i = 0; i < N_DPT_SAMPLES; i++)
    if (Dpt_threshold[i] <= p_gen)
      sumL += Leff[i];
  return (sumL > 0.0) ? 1.0 / sumL : 0.0;
}

// Index into Dpt_threshold[] of the sample a file belongs to, from the
// "..._PT-<thr>_..." dataset name in its path; -1 if not one of the samples.
inline int dptSampleIndex(const TString &fname)
{
  TString tag = "_PT-";
  Ssiz_t i = fname.Index(tag);
  if (i < 0)
    return -1;
  Ssiz_t j = fname.Index("_", i + tag.Length());
  if (j < 0)
    return -1;
  TString thr = fname(i + tag.Length(), j - i - tag.Length());
  if (!thr.IsDigit())
    return -1;
  for (int k = 0; k < N_DPT_SAMPLES; k++)
    if (thr.Atoi() == (int)Dpt_threshold[k])
      return k;
  return -1;
}

// ---------------------------------------------------------------------------
//  Nominal BDT cut table (same file BDTHandler's constructor loads).  The macro
//  validates it and reloads it explicitly, since loadCuts() does no checking.
// ---------------------------------------------------------------------------
static const TString BDT_CUTS_FILE =
    "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/BDT_sys/MC_based/bdt_cuts.csv";

// ---------------------------------------------------------------------------
//  Random-cut control: N_RAND replicas, each keeping a candidate with
//  probability = BDT (cent0to10 cut) pass fraction of its v2 pT bin, taken
//  from the per-sample profiles of a previous hadd-ed flow output.  A random
//  subset has no real effect, so its dvn pulls test the Barlow errors.
// ---------------------------------------------------------------------------
static const int N_RAND = 20;
static const TString RAND_FRAC_REF =
    "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/"
    "BDT_sys_Flow_MC_prompt_Sept29_v2/ROOT/flow_MC_BDT_out_combined.root";

// ---------------------------------------------------------------------------
//  Default location of the hadd-ed resolution file produced by
//  Calculate_Resolution_MC.C  (edit to match your DATE_tag / scratch path).
// ---------------------------------------------------------------------------
static const TString RES_FILE_DEFAULT =
    "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/"
    "BDT_sys_Resolution_MC_prompt_Sept10_v2/ROOT/Resolution_MC_combined.root";
