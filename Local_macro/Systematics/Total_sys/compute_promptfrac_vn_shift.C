// compute_promptfrac_vn_shift.C
//
// Standalone helper macro (independent of plot_ESE_systematics.C): measures the real
// (v_n,1 - v_n,2) prompt-fraction shift directly from data -- comparing the nominal
// (DCA < 0.0085 cm) D0 file against the loose-DCA (DCA < 0.135 cm, "SysPromFrac") D0
// file -- and combines it with the (centrality, pT)-binned "factor" already tabulated in
// PromptFraction_vn_factor.txt, producing a fresh table with a MEASURED vn_shift column
// (replacing that file's placeholder dvn = 0.01 assumption).
//
//   v_n,1 = nominal (tight-DCA-cut) v_n              -- prompt-enriched
//   v_n,2 = loose-DCA-cut ("SysPromFrac" file) v_n    -- closer to the fully inclusive sample
//   dvn   = v_n,1 - v_n,2
//
// The factor table's f_prompt values are (cent, pT) binned on v2's 9-bin scheme only (they
// come from the DCA/prompt-fraction fit, not from v_n itself). dvn, on the other hand, is
// measured directly from each harmonic's own v_n-vs-qbin graphs, so v2 uses its native 9 pT
// bins and v3 uses its native 6 pT bins. For v3 rows, the "factor" is borrowed from the
// nearest v2 pT bin in log(pT) (same convention as NearestV2PtBin() in
// plot_ESE_systematics.C) since the factor table itself has no v3-specific binning.
//
// dvn averaging: for each (harmonic, centrality, pT) bin, v_n,1 and v_n,2 are each averaged
// (unweighted mean) over their q2/q3 sub-bin values before differencing -- as requested,
// since the factor table has no q-bin resolution to match against.
//
// This macro only reads ROOT files and a text table and writes a text table; it does NOT
// modify plot_ESE_systematics.C or its (currently disabled) PromptFraction source. Point
// plot_ESE_systematics.C's prompt_frac_shift_file at the output of this macro, and flip
// source_is_table_shift[5] to true, only when you're ready to enable that systematic.
//
// Usage:  root -l -b -q 'compute_promptfrac_vn_shift.C()'

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>
#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TString.h"

// =========================================================================
// USER CONFIGURATION
// =========================================================================
static const char* nominal_d0_file    = "syst_inputs/D0_Flow_12NUQbin_diffq2q3_out_combined_Aug24.root";
static const char* promptfrac_d0_file = "syst_inputs/D0_Flow_SysPromFrac_combined.root";
static const char* factor_table_in    = "PromptFraction_vn_factor.txt";
static const char* factor_table_out   = "PromptFraction_vn_factor_measured.txt";

// =========================================================================
// Bin definitions -- MUST match plot_ESE_systematics.C
// =========================================================================
static const int N_CENTBINS = 5;
static const char* g_cen_name[N_CENTBINS] = {
    "cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"
};
static const int N_PTBINS_V2 = 9;
static const double g_pt_edges_v2[N_PTBINS_V2 + 1] = { 2, 3, 4, 5, 6, 8, 10, 15, 30, 100 };
static const int N_PTBINS_V3 = 6;
static const double g_pt_edges_v3[N_PTBINS_V3 + 1] = { 2, 4, 6, 8, 10, 20, 50 };
static const int N_QBINS = 12;

static double PtCenter(const double* edges, int ip) { return std::sqrt(edges[ip] * edges[ip + 1]); }

// Index into g_pt_edges_v2 whose bin center is nearest (in log pT) to pt_center -- used to
// borrow a v2-table factor for a v3 pT bin. Mirrors plot_ESE_systematics.C's own helper.
static int NearestV2PtBin(double pt_center)
{
    int best = 0;
    double bestd = std::numeric_limits<double>::max();
    for (int ip = 0; ip < N_PTBINS_V2; ++ip) {
        double d = std::fabs(std::log(PtCenter(g_pt_edges_v2, ip)) - std::log(pt_center));
        if (d < bestd) { bestd = d; best = ip; }
    }
    return best;
}

struct QbinAvg { double val = 0, err = 0; int npts = 0; bool ok = false; };

// Reads v{ivn}_vs_q{ivn}bin_<cent>_pT<lo>to<hi> from "vn_vs_qbin" and returns the unweighted
// mean of its q-bin y-values (statistical error propagated as sqrt(sum(e_i^2))/n).
static QbinAvg AverageOverQbins(TDirectory* dir, int ivn, const char* cen, double pt_lo, double pt_hi)
{
    QbinAvg out;
    TString name = Form("v%d_vs_q%dbin_%s_pT%dto%d", ivn, ivn, cen, (int)pt_lo, (int)pt_hi);
    auto* g = (TGraphErrors*)dir->Get(name.Data());
    if (!g || g->GetN() == 0) return out;

    double sum = 0, sumerr2 = 0;
    int n = 0;
    for (int iq = 0; iq < N_QBINS && iq < g->GetN(); ++iq) {
        double y = g->GetPointY(iq), e = g->GetErrorY(iq);
        if (e <= 0) continue;
        sum += y; sumerr2 += e * e; ++n;
    }
    if (n == 0) return out;
    out.val  = sum / n;
    out.err  = std::sqrt(sumerr2) / n;
    out.npts = n;
    out.ok   = true;
    return out;
}

struct FactorRow { double f1 = 0, f1e = 0, f2 = 0, f2e = 0, factor = 0; bool ok = false; };

static bool LoadFactorTable(const char* path, FactorRow table[N_CENTBINS][N_PTBINS_V2])
{
    std::ifstream fin(path);
    if (!fin.is_open()) {
        std::cerr << "ERROR: cannot open factor table: " << path << "\n";
        return false;
    }
    TString ptTag_v2[N_PTBINS_V2];
    for (int ip = 0; ip < N_PTBINS_V2; ++ip)
        ptTag_v2[ip] = Form("pT%dto%d", (int)g_pt_edges_v2[ip], (int)g_pt_edges_v2[ip + 1]);

    std::string line;
    int nfilled = 0;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::string centTok, ptTok;
        double f1, f1e, f2, f2e, factor, shift;
        if (!(iss >> centTok >> ptTok >> f1 >> f1e >> f2 >> f2e >> factor >> shift)) continue; // header/malformed

        int ic = -1;
        for (int k = 0; k < N_CENTBINS; ++k) if (centTok == g_cen_name[k]) { ic = k; break; }
        int ip = -1;
        for (int k = 0; k < N_PTBINS_V2; ++k) if (ptTok == ptTag_v2[k].Data()) { ip = k; break; }
        if (ic < 0 || ip < 0) continue;

        table[ic][ip] = FactorRow{f1, f1e, f2, f2e, factor, true};
        ++nfilled;
    }
    std::cout << "Loaded " << nfilled << " factor-table rows from " << path << "\n";
    return nfilled > 0;
}

void compute_promptfrac_vn_shift()
{
    TFile* fNom = TFile::Open(nominal_d0_file);
    TFile* fVar = TFile::Open(promptfrac_d0_file);
    if (!fNom || fNom->IsZombie() || !fVar || fVar->IsZombie()) {
        std::cerr << "ERROR: cannot open D0 files (" << nominal_d0_file << " / "
                   << promptfrac_d0_file << ").\n";
        return;
    }
    auto* dirNom = (TDirectory*)fNom->Get("vn_vs_qbin");
    auto* dirVar = (TDirectory*)fVar->Get("vn_vs_qbin");
    if (!dirNom || !dirVar) {
        std::cerr << "ERROR: vn_vs_qbin directory missing in one of the input files.\n";
        return;
    }

    FactorRow factor_v2[N_CENTBINS][N_PTBINS_V2];
    if (!LoadFactorTable(factor_table_in, factor_v2)) return;

    std::ofstream fout(factor_table_out);
    fout << "# f_prompt,1 = prompt fraction for DCA < 0.0085 cm (nominal)\n";
    fout << "# f_prompt,2 = prompt fraction over the whole DCA range (SysPromFrac, DCA < 0.135 cm)\n";
    fout << "# factor     = (1 - f_prompt,1) / (f_prompt,1 - f_prompt,2)   [carried over from "
         << factor_table_in << "]\n";
    fout << "# vn1_avg    = mean of v_n,1 (nominal) over the q2/q3 sub-bins of this (cent,pT) bin\n";
    fout << "# vn2_avg    = mean of v_n,2 (SysPromFrac) over the same q2/q3 sub-bins\n";
    fout << "# dvn_meas   = vn1_avg - vn2_avg   -- MEASURED (v_n,1 - v_n,2), replacing the\n";
    fout << "#              previous placeholder dvn = 0.01\n";
    fout << "# vn_shift   = factor * dvn_meas\n";
    fout << "# v3 rows borrow \"factor\" from the nearest v2 pT bin in log(pT) (the factor table\n";
    fout << "# has no v3-specific binning); dvn_meas for v3 is measured on v3's own pT bins.\n";
    fout << "#\n";
    fout << std::left
         << std::setw(5)  << "vn"
         << std::setw(13) << "centrality"
         << std::setw(12) << "pT"
         << std::right
         << std::setw(11) << "f_prompt,1" << std::setw(9) << "err"
         << std::setw(11) << "f_prompt,2" << std::setw(9) << "err"
         << std::setw(10) << "factor"
         << std::setw(12) << "vn1_avg" << std::setw(12) << "vn2_avg"
         << std::setw(12) << "dvn_meas"
         << std::setw(12) << "vn_shift" << "\n";

    int nwritten = 0;
    for (int ivn = 2; ivn <= 3; ++ivn) {
        int n_pt = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
        const double* edges = (ivn == 2) ? g_pt_edges_v2 : g_pt_edges_v3;
        for (int ic = 0; ic < N_CENTBINS; ++ic) {
            for (int ip = 0; ip < n_pt; ++ip) {
                QbinAvg a1 = AverageOverQbins(dirNom, ivn, g_cen_name[ic], edges[ip], edges[ip + 1]);
                QbinAvg a2 = AverageOverQbins(dirVar, ivn, g_cen_name[ic], edges[ip], edges[ip + 1]);
                if (!a1.ok || !a2.ok) {
                    std::cerr << "WARNING: missing v" << ivn << " graph for " << g_cen_name[ic]
                              << " pT" << (int)edges[ip] << "to" << (int)edges[ip + 1] << " -- skipped.\n";
                    continue;
                }
                double dvn = a1.val - a2.val;

                int ip_f = (ivn == 2) ? ip : NearestV2PtBin(PtCenter(edges, ip));
                const FactorRow& fr = factor_v2[ic][ip_f];
                if (!fr.ok) {
                    std::cerr << "WARNING: no factor-table entry for " << g_cen_name[ic]
                              << " (borrowed v2 pT bin " << ip_f << ") -- skipped.\n";
                    continue;
                }
                double vn_shift = fr.factor * dvn;

                fout << std::left
                     << std::setw(5)  << ivn
                     << std::setw(13) << g_cen_name[ic]
                     << std::setw(12) << Form("pT%dto%d", (int)edges[ip], (int)edges[ip + 1])
                     << std::right << std::fixed << std::setprecision(4)
                     << std::setw(11) << fr.f1 << std::setw(9) << fr.f1e
                     << std::setw(11) << fr.f2 << std::setw(9) << fr.f2e
                     << std::setw(10) << fr.factor
                     << std::setw(12) << a1.val << std::setw(12) << a2.val
                     << std::setw(12) << dvn
                     << std::setw(12) << vn_shift << "\n";
                ++nwritten;
            }
        }
    }
    fout.close();
    fNom->Close();
    fVar->Close();

    std::cout << "Wrote " << nwritten << " rows to " << factor_table_out << "\n";
}
