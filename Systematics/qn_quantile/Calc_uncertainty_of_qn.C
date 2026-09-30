// Run3 (2023 PbPb) version of Calc_uncertainty_of_q2.C, generalized to q2 AND q3.
//
// Same method as the Run2 macro: for a set of reference runs, read the
// per-cent-bin resolution/EP quantity from RescorSave-style .dat files,
// take the largest pairwise spread across runs in each 1%-wide centrality
// bin, normalize by the smallest value in that bin -> relative uncertainty.
// That relative uncertainty is then applied as a symmetric +/- scale factor
// to every nominal q2/q3 quantile edge in that centrality bin, and the
// varied edges are written out as two standalone header files (UP/DOWN),
// ready to be swapped in for the centrality-up/down-style systematic study.
//
// Run2 used trackmid2_1.dat (tracker-based). Here we use the HF-based
// quantities instead, per the new 2023 EP calibration: HF2_1.dat for q2,
// HF3_1.dat for q3, read from Rescor_files/Rescor_<run>/ for 6 reference runs.
//
// The .dat files cover the full 0-100% centrality range in 1% steps, but the
// analysis (and the nominal quantile header) only needs 0-50%, so only the
// first N_CENTBINS_1 rows are read/used here -- unlike the Run2 macro, which
// looped over all 100 bins unnecessarily.
//
// Also writes a plain-text report (fine 1%-bin table + envelope summary over
// the 5 physics centrality bins used downstream) for the analysis note.
//
// Run via ROOT: root.exe -l -q Calc_uncertainty_of_qn.C

#include <fstream>
#include <iostream>
#include <cmath>

#include "TString.h"

using namespace std;

// Nominal Run3 quantile edges: q2_cuts[N_CENTBINS_1][N_QBINS+1], q3_cuts[N_CENTBINS_1][N_QBINS+1]
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/qn_quantile/quantile_12NUQbin_diffq2q3_cuts_2023_MB0to31_Aug23.h"

const int N_CENTBINS_1 = 50; // 1%-wide cent bins, 0-50%; matches q2_cuts/q3_cuts row count
const int N_QBINS       = 12; // number of quantile bins -> N_QBINS+1 edges per row

// Physics centrality bins used downstream (Systematics/centrality_sys/Analysis_bin.h,
// N_CENTBINS=5): the report envelopes the fine 1%-bin uncertainty over these.
const int N_CENTBINS_PHYS = 5;
const int cen_edges_phys[N_CENTBINS_PHYS + 1] = {0, 10, 20, 30, 40, 50};

const int N_RUNS = 6;
const int runNumber[N_RUNS] = {374681, 374719, 374810, 374970, 375064, 375391};

const TString RescorDir = "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/qn_quantile/Rescor_files";
const TString OutDir    = "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/qn_quantile";

// Computes the max relative run-to-run spread of the quantity in <datName>
// (one file per run, "cen_min cen_max value value_err" per line) and fills
// maximum_uncer[N_CENTBINS_1] with the relative uncertainty per 1%-wide cent
// bin, for cent bins 0 to N_CENTBINS_1-1 (0-50%) only.
void computeMaxRelUncertainty(const char* datName, Float_t maximum_uncer[N_CENTBINS_1])
{
    Float_t cen_min[N_CENTBINS_1], cen_max[N_CENTBINS_1], evt_plane[N_RUNS][N_CENTBINS_1], evt_plane_err;

    ifstream inf[N_RUNS];
    for (int r = 0; r < N_RUNS; r++) {
        TString path = TString::Format("%s/Rescor_%d/%s", RescorDir.Data(), runNumber[r], datName);
        inf[r].open(path.Data());
        if (!inf[r].is_open())
            cerr << "ERROR: could not open " << path << endl;
    }

    for (int i = 0; i < N_CENTBINS_1; i++)
        for (int r = 0; r < N_RUNS; r++)
            inf[r] >> cen_min[i] >> cen_max[i] >> evt_plane[r][i] >> evt_plane_err;

    for (int i = 0; i < N_CENTBINS_1; i++) {
        // Largest pairwise difference between the 6 runs in this cent bin
        Float_t max_diff = 0.0;
        for (int a = 0; a < N_RUNS; a++)
            for (int b = a + 1; b < N_RUNS; b++) {
                Float_t diff = fabs(evt_plane[a][i] - evt_plane[b][i]);
                if (diff > max_diff) max_diff = diff;
            }

        // Smallest mean value among the 6 runs in this cent bin
        Float_t min_evt_plane = 1.0e5;
        for (int r = 0; r < N_RUNS; r++)
            if (evt_plane[r][i] < min_evt_plane) min_evt_plane = evt_plane[r][i];

        maximum_uncer[i] = max_diff / min_evt_plane;
        cout << "  cent bin " << i << " (" << cen_min[i] << "-" << cen_max[i]
             << "%): relative uncertainty = " << 100.0 * maximum_uncer[i] << "%" << endl;
    }
}

// Writes a plain-text report for the analysis note: the fine 1%-bin relative
// uncertainty for q2/q3, plus the envelope (worst-case) uncertainty over each
// of the 5 physics centrality bins actually used in the flow measurement.
// Up/down variations are +unc%/-unc% by construction (symmetric scale factor).
void writeUncertaintyReport(const TString& outName,
                             const Float_t maxUncer_q2[N_CENTBINS_1],
                             const Float_t maxUncer_q3[N_CENTBINS_1])
{
    ofstream rf(outName.Data());
    if (!rf.is_open()) {
        cerr << "ERROR: could not create " << outName << endl;
        return;
    }

    rf << "qn quantile-edge systematic: relative uncertainty report\n";
    rf << "(up variation = +unc%, down variation = -unc%, symmetric by construction)\n\n";

    rf << "-- Fine 1%-wide centrality bins --\n";
    rf << "cent_low  cent_high  q2_unc[%]  q3_unc[%]\n";
    for (int i = 0; i < N_CENTBINS_1; i++) {
        char line[128];
        snprintf(line, sizeof(line), "%8d  %9d  %9.3f  %9.3f\n",
                 i, i + 1, 100.0 * maxUncer_q2[i], 100.0 * maxUncer_q3[i]);
        rf << line;
    }

    rf << "\n-- Envelope over physics centrality bins (worst-case 1%-subbin) --\n";
    rf << "cent_low  cent_high  q2_unc[%]  q3_unc[%]\n";
    cout << "\n=== Envelope (max) relative uncertainty per physics centrality bin ===\n";
    cout << "cent_low  cent_high  q2_unc[%]  q3_unc[%]\n";
    for (int c = 0; c < N_CENTBINS_PHYS; c++) {
        int lo = cen_edges_phys[c], hi = cen_edges_phys[c + 1];
        Float_t maxq2 = 0.0, maxq3 = 0.0;
        for (int i = lo; i < hi; i++) {
            if (maxUncer_q2[i] > maxq2) maxq2 = maxUncer_q2[i];
            if (maxUncer_q3[i] > maxq3) maxq3 = maxUncer_q3[i];
        }
        char line[128];
        snprintf(line, sizeof(line), "%8d  %9d  %9.3f  %9.3f\n", lo, hi, 100.0 * maxq2, 100.0 * maxq3);
        rf << line;
        cout << line;
    }

    rf.close();
    cout << "\nWrote " << outName << endl;
}

void writeVariedHeader(const TString& outName,
                        double q2_binning_var[N_CENTBINS_1][N_QBINS + 1],
                        double q3_binning_var[N_CENTBINS_1][N_QBINS + 1])
{
    ofstream hf(outName.Data());
    if (!hf.is_open()) {
        cerr << "ERROR: could not create " << outName << endl;
        return;
    }

    hf << "#ifndef Q_CUTS_2023_H\n#define Q_CUTS_2023_H\n\n";

    hf << "double q2_cuts[" << N_CENTBINS_1 << "][" << N_QBINS + 1 << "] = {\n";
    for (int i = 0; i < N_CENTBINS_1; i++) {
        hf << "  {";
        for (int j = 0; j <= N_QBINS; j++) {
            hf << q2_binning_var[i][j];
            if (j < N_QBINS) hf << ",";
        }
        hf << "}";
        if (i < N_CENTBINS_1 - 1) hf << ",";
        hf << "\n";
    }
    hf << "};\n\n";

    hf << "double q3_cuts[" << N_CENTBINS_1 << "][" << N_QBINS + 1 << "] = {\n";
    for (int i = 0; i < N_CENTBINS_1; i++) {
        hf << "  {";
        for (int j = 0; j <= N_QBINS; j++) {
            hf << q3_binning_var[i][j];
            if (j < N_QBINS) hf << ",";
        }
        hf << "}";
        if (i < N_CENTBINS_1 - 1) hf << ",";
        hf << "\n";
    }
    hf << "};\n\n#endif\n";

    hf.close();
    cout << "Wrote " << outName << endl;
}

void Calc_uncertainty_of_qn()
{
    Float_t maxUncer_q2[N_CENTBINS_1], maxUncer_q3[N_CENTBINS_1];

    cout << "=== q2 relative uncertainty (from HF2_1.dat, " << N_RUNS << " runs, cent 0-50%) ===" << endl;
    computeMaxRelUncertainty("HF2_1.dat", maxUncer_q2);

    cout << endl << "=== q3 relative uncertainty (from HF3_1.dat, " << N_RUNS << " runs, cent 0-50%) ===" << endl;
    computeMaxRelUncertainty("HF3_1.dat", maxUncer_q3);

    static double q2_binning_up[N_CENTBINS_1][N_QBINS + 1], q2_binning_down[N_CENTBINS_1][N_QBINS + 1];
    static double q3_binning_up[N_CENTBINS_1][N_QBINS + 1], q3_binning_down[N_CENTBINS_1][N_QBINS + 1];

    for (int i = 0; i < N_CENTBINS_1; i++) {
        for (int j = 0; j <= N_QBINS; j++) {
            q2_binning_up[i][j]   = q2_cuts[i][j] * (1.0 + maxUncer_q2[i]);
            q2_binning_down[i][j] = q2_cuts[i][j] * (1.0 - maxUncer_q2[i]);
            q3_binning_up[i][j]   = q3_cuts[i][j] * (1.0 + maxUncer_q3[i]);
            q3_binning_down[i][j] = q3_cuts[i][j] * (1.0 - maxUncer_q3[i]);
        }
    }

    cout << endl;
    writeVariedHeader(OutDir + "/quantile_12NUQbin_diffq2q3_cuts_2023_MB0to31_qnUP.h",   q2_binning_up,   q3_binning_up);
    writeVariedHeader(OutDir + "/quantile_12NUQbin_diffq2q3_cuts_2023_MB0to31_qnDOWN.h", q2_binning_down, q3_binning_down);

    writeUncertaintyReport(OutDir + "/qn_quantile_uncertainty_report.txt", maxUncer_q2, maxUncer_q3);
}
