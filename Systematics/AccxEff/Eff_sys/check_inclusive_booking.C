// Usage: root -l -b -q 'check_inclusive_booking.C("/path/to/merged_output.root")'
//
// Verifies that entries in the cent0to50 (inclusive) histograms equal the
// sum of the five native 10%-wide bins, for the histogram groups that are
// filled via three different code paths:
//   - mass_fitted/h_mass_<cen>_<pt>_def                    (Step-2 Add() merge, v2 only)
//   - q_distributions/hq2_dist_<cen>_q2bin<n>               (direct per-event fill)
//   - q_distributions/hq3_dist_<cen>_q3bin<n>               (direct per-event fill)
//   - histograms_v2v3dist/h_v2dist_<cen>_q2bin<n>_<pt>      (Step-2 Add() merge, v2)
//   - histograms_v2v3dist/h_v3dist_<cen>_q3bin<n>_<pt>      (Step-2 Add() merge, v3)
//
// NOTE: h_mass_default_v3 / h_mass_q3 ("hist_mass_v3" in the analysis code) are
// never actually Write()'d anywhere in fit_mass_and_flow.C -- that's a pre-existing
// gap in the v3 branch (v2's counterparts do get written), present for every
// centrality bin, not just the new inclusive one. Don't use h_mass_*_pT*_def for a
// v3 check: most v3 pT-bin names are simply absent from the file, and the two that
// happen to share a label with a v2 pT bin (pT6to8, pT8to10) will silently resolve
// to the *v2* histogram instead, giving a false pass. Use h_v3dist_grouped instead,
// which is genuinely written for v3.

#include <vector>
#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"

void checkGroup(TFile* f, const char* dir, const char* prefix, const std::vector<TString>& suffixes)
{
    const char* natives[5] = {"cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"};
    int nBad = 0, nChecked = 0;

    for (auto& suf : suffixes)
    {
        double sumNative = 0;
        bool allFound = true;
        for (auto n : natives)
        {
            TString path = Form("%s/%s%s%s", dir, prefix, n, suf.Data());
            TH1D* h = (TH1D*)f->Get(path);
            if (!h) { std::cout << "MISSING: " << path << std::endl; allFound = false; continue; }
            sumNative += h->GetEntries();
        }

        TString inclPath = Form("%s/%s%s%s", dir, prefix, "cent0to50", suf.Data());
        TH1D* hIncl = (TH1D*)f->Get(inclPath);
        if (!hIncl) { std::cout << "MISSING: " << inclPath << std::endl; allFound = false; }
        if (!allFound) continue;

        nChecked++;
        double inclEntries = hIncl->GetEntries();
        bool ok = (std::abs(inclEntries - sumNative) < 1e-6);
        if (!ok) nBad++;
        std::cout << (ok ? "OK   " : "FAIL ")
                  << inclPath << "  incl=" << inclEntries
                  << "  sum(native)=" << sumNative << std::endl;
    }
    std::cout << prefix << ": " << (nChecked - nBad) << "/" << nChecked << " passed" << std::endl << std::endl;
}

void check_inclusive_booking(const char* fname)
{
    TFile* f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << fname << std::endl; return; }

    std::vector<TString> pt_v2 = {"pT2to3","pT3to4","pT4to5","pT5to6","pT6to8","pT8to10","pT10to15","pT15to30","pT30to100"};
    std::vector<TString> pt_v3 = {"pT2to4","pT4to6","pT6to8","pT8to10","pT10to20","pT20to50","pT50to100"};

    std::vector<TString> v2_def_suf, v3_def_suf;
    for (auto& pt : pt_v2) v2_def_suf.push_back(TString("_") + pt + "_def");
    for (auto& pt : pt_v3) v3_def_suf.push_back(TString("_") + pt + "_def");

    std::vector<TString> q2_suf, q3_suf;
    for (int iq = 0; iq < 12; iq++) q2_suf.push_back(Form("_q2bin%d", iq));
    for (int iq = 0; iq < 12; iq++) q3_suf.push_back(Form("_q3bin%d", iq));

    std::vector<TString> v2dist_suf, v3dist_suf;
    for (int iq = 0; iq < 12; iq++)
    {
        for (auto& pt : pt_v2) { v2dist_suf.push_back(TString(Form("_q2bin%d_", iq)) + pt); }
    }
    for (int iq = 0; iq < 12; iq++)
    {
        for (auto& pt : pt_v3) { v3dist_suf.push_back(TString(Form("_q3bin%d_", iq)) + pt); }
    }

    std::cout << "=== h_mass_default_v2 (mass_fitted/, Step-2 merge) ===" << std::endl;
    checkGroup(f, "mass_fitted", "h_mass_", v2_def_suf);

    std::cout << "=== hq2_dist (q_distributions/, direct per-event fill) ===" << std::endl;
    checkGroup(f, "q_distributions", "hq2_dist_", q2_suf);

    std::cout << "=== hq3_dist (q_distributions/, direct per-event fill) ===" << std::endl;
    checkGroup(f, "q_distributions", "hq3_dist_", q3_suf);

    std::cout << "=== h_v2dist_grouped (histograms_v2v3dist/, Step-2 merge) ===" << std::endl;
    checkGroup(f, "histograms_v2v3dist", "h_v2dist_", v2dist_suf);

    std::cout << "=== h_v3dist_grouped (histograms_v2v3dist/, Step-2 merge, v3-specific) ===" << std::endl;
    checkGroup(f, "histograms_v2v3dist", "h_v3dist_", v3dist_suf);

    f->Close();
}
