#pragma once
// =============================================================================
//  SampleCombination.h
//
//  Per-pT-bin combination of the UNWEIGHTED per-Dpt-sample profiles written by
//  flow_MC_BDT_sys.C  (pv{n}_woBDT_Dpt<thr>, pv{n}_wBDT_<cen>_Dpt<thr>).
//
//  In every sample, separately:
//      dvn_s  = vn_s(BDT) - vn_s(noBDT)
//      sig_s  = sqrt| sig_s(BDT)^2 - sig_s(noBDT)^2 |   (Barlow; exact for an
//                                                      unweighted subsample)
//  then inverse-variance average over samples:
//      dvn    = sum(dvn_s / sig_s^2) / sum(1 / sig_s^2),  sig = 1/sqrt(sum 1/sig_s^2)
//  vn(noBDT) is combined the same way with its own errors (closure check).
//
//  Valid because the true vn of the embedded D0 is 0 in every sample, so each
//  sample is an independent closure test of the same quantity.  Replaces the
//  cross-section stitch weights, whose ~80x jump at the Dpt>10 threshold lets a
//  handful of gen-pT<10 candidates migrating to reco pT>10 dominate sum(w^2).
//
//  A sample enters a bin only with >= MIN_ENTRIES_WO noBDT and
//  >= MIN_ENTRIES_W BDT-selected candidates (profile errors are unreliable
//  for a handful of entries) and a non-zero Barlow error.
// =============================================================================
#include <cmath>
#include "TFile.h"
#include "TProfile.h"
#include "Analysis_bin_BDT.h"

static const double MIN_ENTRIES_WO = 50;
static const double MIN_ENTRIES_W  = 20;

struct CombinedBin
{
  double vwo = 0, ewo = 0; // vn(noBDT), inverse-variance combined
  double d = 0, de = 0;    // dvn and its error, inverse-variance combined
  double nwo = 0;          // noBDT candidates in the samples used
  int nsamples = 0;        // samples used
};

// true if the per-sample profiles exist in f
inline bool hasSampleProfiles(TFile *f)
{
  return f->Get(Form("pv2_woBDT_Dpt%d", (int)Dpt_threshold[0])) != nullptr;
}

// n = 2 or 3, sel = selected-profile tag ("wBDT_<cen>" or "wRAND<r>"),
// ib = profile bin (1-based)
inline CombinedBin combineSamplesSel(TFile *f, int n, const TString &sel, int ib)
{
  CombinedBin c;
  double sw_wo = 0, swv_wo = 0, sw_d = 0, swv_d = 0;
  for (int s = 0; s < N_DPT_SAMPLES; ++s)
  {
    int thr = (int)Dpt_threshold[s];
    TProfile *pwo = (TProfile *)f->Get(Form("pv%d_woBDT_Dpt%d", n, thr));
    TProfile *pw = (TProfile *)f->Get(Form("pv%d_%s_Dpt%d", n, sel.Data(), thr));
    if (!pwo || !pw)
      continue;
    double nwo = pwo->GetBinEntries(ib), nw = pw->GetBinEntries(ib);
    if (nwo < MIN_ENTRIES_WO || nw < MIN_ENTRIES_W)
      continue;
    double vwo = pwo->GetBinContent(ib), ewo = pwo->GetBinError(ib);
    double ew = pw->GetBinError(ib);
    double d = pw->GetBinContent(ib) - vwo;
    double de = sqrt(fabs(ew * ew - ewo * ewo));
    if (ewo <= 0 || de <= 0)
      continue;
    sw_wo += 1. / (ewo * ewo);
    swv_wo += vwo / (ewo * ewo);
    sw_d += 1. / (de * de);
    swv_d += d / (de * de);
    c.nwo += nwo;
    c.nsamples++;
  }
  if (c.nsamples > 0)
  {
    c.vwo = swv_wo / sw_wo;
    c.ewo = 1. / sqrt(sw_wo);
    c.d = swv_d / sw_d;
    c.de = 1. / sqrt(sw_d);
  }
  return c;
}

// BDT-cut class g
inline CombinedBin combineSamples(TFile *f, int n, int g, int ib)
{
  return combineSamplesSel(f, n, Form("wBDT_%s", cen_name[g]), ib);
}
