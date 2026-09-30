// =============================================================================
//  check_BDT_sys.C
//
//  Quick health check of a hadd-ed flow_MC_BDT_sys.C output, optionally
//  compared with a reference run (e.g. before / after a code change).
//
//  Per pT bin it prints
//    * vn(noBDT) and its pull vn/sigma          -> MC closure (expect ~0)
//    * N_eff of the noBDT profile               -> statistical power
//    * dvn = vn(BDT) - vn(noBDT), Barlow sigma, pull dvn/sigma
//  and per centrality class
//    * chi2/ndf of vn(noBDT)/sigma and of the dvn pulls, over pT
//    * weighted mean of dvn over pT (constant fit) and its error
//  With a reference file, the same numbers are shown side by side with
//  N_eff(new)/N_eff(ref) and sigma_dvn(new)/sigma_dvn(ref).
//
//  Each file is read in one of two modes:
//     stitch : the cross-section-weighted profiles  pv{n}_woBDT / pv{n}_wBDT_<cen>
//     ivw    : per-Dpt-sample unweighted profiles, inverse-variance combined
//              (SampleCombination.h); N_eff is then the raw candidate count
//     auto   : ivw if the file has the per-sample profiles, else stitch
//  The same file can be passed twice with different modes to compare them.
//
//  Run:
//     root -l -b -q 'check_BDT_sys.C("new_combined.root")'
//     root -l -b -q 'check_BDT_sys.C("new.root","ref.root","auto","stitch")'
//  or compiled:
//     g++ check_BDT_sys.C $(root-config --cflags --libs) -O2 -o check_BDT_sys.exe
//     ./check_BDT_sys.exe new_combined.root [ref_combined.root [mode_new [mode_ref]]]
// =============================================================================

#include <cstdio>
#include <cmath>
#include <iostream>

#include "TFile.h"
#include "TProfile.h"
#include "TString.h"

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/BDT_sys/MC_based/Analysis_bin_BDT.h"
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/BDT_sys/MC_based/SampleCombination.h"

using namespace std;

struct BinStat
{
  double vwo, ewo, neff; // noBDT
  double d, de;          // BDT - noBDT, Barlow error
};

static bool getStats(TFile *f, int n, int g, int ib, BinStat &s, bool ivw)
{
  if (ivw)
  {
    CombinedBin c = combineSamples(f, n, g, ib);
    s.vwo = c.vwo;
    s.ewo = c.ewo;
    s.neff = c.nwo;
    s.d = c.d;
    s.de = c.de;
    return true;
  }
  TProfile *pwo = (TProfile *)f->Get(Form("pv%d_woBDT", n));
  TProfile *pw = (TProfile *)f->Get(Form("pv%d_wBDT_%s", n, cen_name[g]));
  if (!pwo || !pw)
    return false;
  s.vwo = pwo->GetBinContent(ib);
  s.ewo = pwo->GetBinError(ib);
  s.neff = pwo->GetBinEffectiveEntries(ib);
  double ew = pw->GetBinError(ib);
  s.d = pw->GetBinContent(ib) - s.vwo;
  s.de = sqrt(fabs(ew * ew - s.ewo * s.ewo));
  return true;
}

struct Summary
{
  double chi2_wo = 0, chi2_d = 0, sw = 0, swd = 0;
  int ndf = 0, nbig = 0; // nbig: bins with |dvn pull| > 2
};

static Summary summarize(TFile *f, int n, int g, int npt, bool ivw)
{
  Summary s;
  for (int ib = 1; ib <= npt; ++ib)
  {
    BinStat b;
    if (!getStats(f, n, g, ib, b, ivw) || b.ewo <= 0 || b.de <= 0)
      continue;
    s.chi2_wo += pow(b.vwo / b.ewo, 2);
    double p = b.d / b.de;
    s.chi2_d += p * p;
    if (fabs(p) > 2)
      s.nbig++;
    s.sw += 1. / (b.de * b.de);
    s.swd += b.d / (b.de * b.de);
    s.ndf++;
  }
  return s;
}

// "auto" -> ivw if the per-sample profiles exist
static bool resolveMode(TFile *f, TString mode)
{
  if (mode == "ivw" && !hasSampleProfiles(f))
  {
    cerr << "FATAL: ivw mode requested but " << f->GetName() << " has no per-sample profiles" << endl;
    exit(1);
  }
  return mode == "ivw" || (mode == "auto" && hasSampleProfiles(f));
}

// ---- random-cut control: pulls of dvn(random subset) with the same
// per-sample IVW machinery; if the Barlow errors are right they are ~N(0,1).
// Printed for all pT bins and for pT < 10 GeV, next to the BDT (cent0to10).
static void randomControl(TFile *f)
{
  if (!f->Get(Form("pv2_wRAND0_Dpt%d", (int)Dpt_threshold[0])))
    return;
  printf("\n==== random-cut control (%d replicas, per-sample IVW) ====\n", N_RAND);
  for (int n = 2; n <= 3; ++n)
  {
    const int npt = (n == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
    const double *edges = (n == 2) ? pt_edges_v2 : pt_edges_v3;
    for (int lowOnly = 0; lowOnly <= 1; ++lowOnly)
    {
      double sp = 0, sp2 = 0;
      int np = 0, nbig = 0, nworse = 0;
      // BDT chi2 in the same bins
      double chi2_bdt = 0;
      int ndf = 0;
      for (int ib = 1; ib <= npt; ++ib)
      {
        if (lowOnly && edges[ib] > 10)
          continue;
        CombinedBin b = combineSamples(f, n, 0, ib);
        if (b.de <= 0)
          continue;
        chi2_bdt += pow(b.d / b.de, 2);
        ndf++;
      }
      TString per;
      for (int r = 0; r < N_RAND; ++r)
      {
        double chi2 = 0;
        for (int ib = 1; ib <= npt; ++ib)
        {
          if (lowOnly && edges[ib] > 10)
            continue;
          CombinedBin b = combineSamplesSel(f, n, Form("wRAND%d", r), ib);
          if (b.de <= 0)
            continue;
          double p = b.d / b.de;
          sp += p;
          sp2 += p * p;
          np++;
          chi2 += p * p;
          if (fabs(p) > 2)
            nbig++;
        }
        if (chi2 >= chi2_bdt)
          nworse++;
        per += Form(" %.1f", chi2);
      }
      double mean = np ? sp / np : 0, rms = np ? sqrt(sp2 / np - mean * mean) : 0;
      printf("v%d %-10s pulls: N=%d mean=%+.3f RMS=%.3f (+-%.3f)  |p|>2: %.1f%% (expect 4.6%%)\n",
             n, lowOnly ? "pT<10" : "all pT", np, mean, rms, np ? rms / sqrt(2. * np) : 0,
             np ? 100. * nbig / np : 0);
      printf("   chi2/ndf per replica (ndf=%d):%s\n", ndf, per.Data());
      printf("   BDT %s chi2 = %.1f/%d ; replicas with chi2 >= BDT: %d/%d\n",
             cen_name[0], chi2_bdt, ndf, nworse, N_RAND);
    }
  }
}

int check_BDT_sys(TString fnew, TString fref = "", TString mode_new = "auto", TString mode_ref = "auto")
{
  TFile *fn = TFile::Open(fnew);
  if (!fn || fn->IsZombie())
  {
    cerr << "cannot open " << fnew << endl;
    return 1;
  }
  TFile *fr = nullptr;
  if (fref.Length())
  {
    fr = TFile::Open(fref);
    if (!fr || fr->IsZombie())
    {
      cerr << "cannot open " << fref << endl;
      return 1;
    }
  }

  const bool ivw_new = resolveMode(fn, mode_new);
  const bool ivw_ref = fr ? resolveMode(fr, mode_ref) : false;
  printf("NEW: %s  [%s]\n", fnew.Data(), ivw_new ? "ivw" : "stitch");
  if (fr)
    printf("REF: %s  [%s]\n", fref.Data(), ivw_ref ? "ivw" : "stitch");

  for (int n = 2; n <= 3; ++n)
  {
    const int npt = (n == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
    const double *edges = (n == 2) ? pt_edges_v2 : pt_edges_v3;

    // ---- per-pT detail for the first class (classes share the noBDT reference)
    const int g0 = 0;
    printf("\n==== v%d  per-pT  (dvn with %s cut) ====\n", n, cen_name[g0]);
    if (!fr)
      printf("%-8s %-20s %6s %8s | %-20s %6s\n", "pT", "vn(noBDT)", "pull", "Neff",
             "dvn +- sigma", "pull");
    else
      printf("%-8s | %-17s %6s %6s | %-17s %6s %6s | %8s %8s\n", "pT",
             "vn(noBDT) NEW", "pull", "dpull", "vn(noBDT) REF", "pull", "dpull",
             "Neff n/r", "sig_d n/r");

    for (int ib = 1; ib <= npt; ++ib)
    {
      BinStat a, r;
      if (!getStats(fn, n, g0, ib, a, ivw_new))
      {
        cerr << "missing profiles in " << fnew << endl;
        return 1;
      }
      TString pt = Form("%g-%g", edges[ib - 1], edges[ib]);
      double pa = a.ewo > 0 ? a.vwo / a.ewo : 0, da = a.de > 0 ? a.d / a.de : 0;
      if (!fr)
      {
        printf("%-8s %+.5f+-%.5f %+6.1f %8.0f | %+.5f+-%.5f %+6.1f\n", pt.Data(),
               a.vwo, a.ewo, pa, a.neff, a.d, a.de, da);
        continue;
      }
      if (!getStats(fr, n, g0, ib, r, ivw_ref))
      {
        cerr << "missing profiles in " << fref << endl;
        return 1;
      }
      double pr = r.ewo > 0 ? r.vwo / r.ewo : 0, dr = r.de > 0 ? r.d / r.de : 0;
      printf("%-8s | %+.4f+-%.4f %+6.1f %+6.1f | %+.4f+-%.4f %+6.1f %+6.1f | %8.2f %8.2f\n",
             pt.Data(), a.vwo, a.ewo, pa, da, r.vwo, r.ewo, pr, dr,
             r.neff > 0 ? a.neff / r.neff : 0, r.de > 0 ? a.de / r.de : 0);
    }

    // ---- per-class summary
    printf("\n---- v%d  summary per BDT-cut class ----\n", n);
    printf("%-12s %-10s %-12s %-6s %-24s\n", "class", "chi2_wo", "chi2_dvn", "|p|>2",
           "<dvn> (const fit)");
    for (int g = 0; g < N_CENTBINS; ++g)
    {
      for (int k = 0; k < (fr ? 2 : 1); ++k)
      {
        Summary s = summarize(k == 0 ? fn : fr, n, g, npt, k == 0 ? ivw_new : ivw_ref);
        double m = s.sw > 0 ? s.swd / s.sw : 0, me = s.sw > 0 ? 1. / sqrt(s.sw) : 0;
        printf("%-12s %4.1f/%-5d %5.1f/%-6d %-6d %+.5f +- %.5f  %s\n",
               cen_name[g], s.chi2_wo, s.ndf, s.chi2_d, s.ndf, s.nbig, m, me,
               fr ? (k == 0 ? "NEW" : "REF") : "");
      }
    }
  }
  if (ivw_new)
    randomControl(fn);
  return 0;
}

#ifndef __CLING__
int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    cout << "Usage: " << argv[0] << " new_combined.root [ref_combined.root [mode_new [mode_ref]]]"
         << "   (mode: auto | stitch | ivw)" << endl;
    return 1;
  }
  return check_BDT_sys(argv[1], argc > 2 ? argv[2] : "", argc > 3 ? argv[3] : "auto",
                       argc > 4 ? argv[4] : "auto");
}
#endif
