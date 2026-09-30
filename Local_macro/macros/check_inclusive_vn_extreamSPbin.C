// check_inclusive_vn_extreamSPbin.C
//
// Computes the inclusive D0 v2 and v3 from the yield-weighted SP bin formula:
//
//   v_n = SUM( v_n^i * Y_i ) / SUM( Y_i )                          (Eq. 4)
//
//   Delta(v_n) = sqrt( SUM( (v_n^i - v_n)^2 * (DY_i)^2 ) ) / SUM( Y_i )  (Eq. 5)
//
// where v_n^i is the centre of the i-th SP bin (x-axis of TGraph),
//       Y_i   is the D0 yield (y-value),  DY_i is its statistical error.
//
// The macro also checks robustness by removing the two extreme SP bins
// (one from each side: the bins with minimum yield = outermost active bins)
// and reports the shift in v_n and its error.
//
// Input:  ROOT/Flow_MB0to31_May18_out_combined.root
//         directory: yields_vN
//
// Usage:
//   root -l -b -q 'check_inclusive_vn_extreamSPbin.C("ROOT/Flow_MB0to31_May18_out_combined.root",2,"0to10","10to15")'   # v2, all q-bins
//   root -l -b -q 'check_inclusive_vn_extreamSPbin.C("ROOT/Flow_MB0to31_May18_out_combined.root",3,"0to10","10to20")'   # v3, all q-bins
//   root -l -b -q  check_inclusive_vn_extreamSPbin.C                                                                         # default below
//
// pT bin tags per harmonics:
//   v2: 2to3  3to4  4to5  5to6  6to8  8to10  10to15  15to20  20to50
//   v3: 2to4  4to6  6to8  8to10  10to20  20to50  50to100

#include "TFile.h"
#include "TDirectoryFile.h"
#include "TGraphErrors.h"
#include <cmath>
#include <vector>
#include <limits>

// ============================================================================
// Compute vn and its error from Eq. 4 & 5.
// If skip_lo or skip_hi index >= 0, those points are excluded.
// ============================================================================
struct VnResult   { double vn, evn, sumY; int nused; };
struct CheckResult { bool valid; VnResult rAll, rExcl; double dVn, relV, sig; };

VnResult calcVn(const std::vector<double>& vx,
                const std::vector<double>& vy,
                const std::vector<double>& vey,
                int skip_lo = -1, int skip_hi = -1)
{
    int N = (int)vx.size();
    double sumYv = 0, sumY = 0;
    int nused = 0;
    for (int i = 0; i < N; i++) {
        if (i == skip_lo || i == skip_hi) continue;
        if (vy[i] <= 0.0 || vey[i] <= 0.0) continue;
        sumYv += vx[i] * vy[i];
        sumY  += vy[i];
        nused++;
    }
    if (sumY <= 0) return {0, 0, 0, 0};
    double vn = sumYv / sumY;
    double var = 0;
    for (int i = 0; i < N; i++) {
        if (i == skip_lo || i == skip_hi) continue;
        if (vy[i] <= 0.0 || vey[i] <= 0.0) continue;
        double dv = vx[i] - vn;
        var += dv * dv * vey[i] * vey[i];
    }
    return {vn, std::sqrt(var) / sumY, sumY, nused};
}

// ============================================================================
// Run the check for one vn order — prints detail and returns summary values
// ============================================================================
CheckResult runCheck(TDirectoryFile* dir,
                     int vn_order, const char* pt_tag, const char* cent_tag, int qbin)
{
    CheckResult bad = {false, {0,0,0,0}, {0,0,0,0}, 0, 0, 0};
    TString gname = Form("v%d_graph_pT%s_cent%s_q%dbin%d",
                          vn_order, pt_tag, cent_tag, vn_order, qbin);

    TGraphErrors* gr = (TGraphErrors*)dir->Get(gname.Data());
    if (!gr) {
        Printf("  ERROR: TGraph '%s' not found — skipping.", gname.Data());
        return bad;
    }

    int N = gr->GetN();
    std::vector<double> vx(N), vy(N), vey(N);
    for (int i = 0; i < N; i++) {
        gr->GetPoint(i, vx[i], vy[i]);
        vey[i] = gr->GetErrorY(i);
    }

    // find active bins (Y > 0 and DY > 0)
    int first_active = -1, last_active = -1;
    for (int i = 0;   i < N; i++) { if (vy[i]>0&&vey[i]>0){first_active=i; break;} }
    for (int i = N-1; i >= 0; i--){ if (vy[i]>0&&vey[i]>0){last_active =i; break;} }

    if (first_active < 0) { Printf("  No active bins found for %s", gname.Data()); return bad; }

    // all-bins result
    VnResult rAll = calcVn(vx, vy, vey);

    // exclude the two extreme bins (first_active and last_active)
    VnResult rExcl = calcVn(vx, vy, vey, first_active, last_active);

    double dVn     = rExcl.vn  - rAll.vn;
    double dEvn    = rExcl.evn - rAll.evn;
    double relV    = (rAll.vn  != 0) ? dVn  / rAll.vn  * 100.0 : 0;
    double relE    = (rAll.evn != 0) ? dEvn / rAll.evn * 100.0 : 0;
    double sig     = (rAll.evn > 0 || rExcl.evn > 0)
                     ? std::abs(dVn) / std::sqrt(rAll.evn*rAll.evn + rExcl.evn*rExcl.evn)
                     : 0;

    Printf("\n  ──────────────────────────────────────────────────────────────");
    Printf("  v%d  |  TGraph: %s", vn_order, gname.Data());
    Printf("  ──────────────────────────────────────────────────────────────");
    Printf("  Active SP bins : %d used of %d total  (first=%d  last=%d)",
           rAll.nused, N, first_active, last_active);
    Printf("  Excluded bins  : index %d  (x=%.3f, Y=%.2f)   and   index %d  (x=%.3f, Y=%.2f)",
           first_active, vx[first_active], vy[first_active],
           last_active,  vx[last_active],  vy[last_active]);
    Printf("  Sum of yields  : all = %.2f   excl = %.2f   removed = %.2f",
           rAll.sumY, rExcl.sumY, rAll.sumY - rExcl.sumY);
    Printf("  ──────────────────────────────────────────────────────────────");
    Printf("  Result (all bins)     :  v%d = %+.6f  ±  %.6f",
           vn_order, rAll.vn, rAll.evn);
    Printf("  Result (excl. extr.)  :  v%d = %+.6f  ±  %.6f",
           vn_order, rExcl.vn, rExcl.evn);
    Printf("  ──────────────────────────────────────────────────────────────");
    Printf("  Shift in v%d          :  %+.6f  (%+.3f%%)",
           vn_order, dVn, relV);
    Printf("  Shift in Delta(v%d)   :  %+.6f  (%+.3f%%)",
           vn_order, dEvn, relE);
    Printf("  v%d / Delta           :  all = %.2f sigma   excl = %.2f sigma",
           vn_order,
           (rAll.evn  > 0) ? rAll.vn  / rAll.evn  : 0,
           (rExcl.evn > 0) ? rExcl.vn / rExcl.evn : 0);
    Printf("  Shift signif.  :  %.2f sigma  (stat.)", sig);

    return {true, rAll, rExcl, dVn, relV, sig};
}

// ============================================================================
// Main
// ============================================================================
void check_inclusive_vn_extreamSPbin(
    const char* in_file  = "ROOT/Flow_MB0to31_May19_out_combined_v3.root",
    int         vn_order = 2,          // 2 or 3
    const char* cent_tag = "0to10",    // e.g. "0to10", "10to30", "30to50"
    const char* pt_tag   = "10to15")   // pT tag for the chosen harmonic
{
    TFile* f = TFile::Open(in_file);
    if (!f || f->IsZombie()) { Printf("ERROR: cannot open %s", in_file); return; }

    TDirectoryFile* dir = (TDirectoryFile*)f->Get("yields_vN");
    if (!dir) { Printf("ERROR: 'yields_vN' not found"); f->Close(); return; }

    Printf("\n==============================================================");
    Printf("  Inclusive D^0 flow  —  extreme SP bin robustness check");
    Printf("==============================================================");
    Printf("  v%d pT bin : %s GeV/c", vn_order, pt_tag);
    Printf("  Cent bin  : %s%%", cent_tag);
    Printf("  Formula   : v_n = SUM(v_n^i * Y_i) / SUM(Y_i)");
    Printf("  Excluded  : outermost active SP bin on each side (min yield)");

    CheckResult res[10];
    for (int qbin = 0; qbin <= 9; qbin++) {
        Printf("\n  ==== q%d bin : %d ====", vn_order, qbin);
        res[qbin] = runCheck(dir, vn_order, pt_tag, cent_tag, qbin);
    }

    // ── Summary table ─────────────────────────────────────────────────────────
    Printf("\n");
    Printf("  ==========================================================================================================");
    Printf("  SUMMARY TABLE  --  v%d,  pT %s GeV/c,  Cent %s%%", vn_order, pt_tag, cent_tag);
    Printf("  ==========================================================================================================");
    Printf("  %4s |  %-20s |  %-20s |  %-19s |  %s",
           "qbin",
           Form("v%d (all bins)",   vn_order),
           Form("v%d (excl.extr.)", vn_order),
           "Shift",
           "Signif.");
    Printf("  -----|----------------------|----------------------|---------------------|--------");
    for (int qbin = 0; qbin <= 9; qbin++) {
        if (!res[qbin].valid) {
            Printf("  %4d |  %-20s |  %-20s |  %-19s |  %s",
                   qbin, "N/A", "N/A", "N/A", "N/A");
            continue;
        }
        const VnResult& a = res[qbin].rAll;
        const VnResult& e = res[qbin].rExcl;
        Printf("  %4d |  %+.5f +/- %.5f |  %+.5f +/- %.5f |  %+.5f (%+.1f%%) |  %.2f sig",
               qbin,
               a.vn, a.evn,
               e.vn, e.evn,
               res[qbin].dVn, res[qbin].relV,
               res[qbin].sig);
    }
    Printf("  ==========================================================================================================");

    Printf("\n==============================================================\n");
    f->Close();
}

// ---------------------------------------------------------------------------
#ifndef __CINT__
#ifndef __CLING__
int main() { check_inclusive_vn_extreamSPbin(); return 0; }
#endif
#endif

