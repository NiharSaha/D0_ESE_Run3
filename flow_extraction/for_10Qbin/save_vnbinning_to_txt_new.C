// save_vnbinning_to_txt_new.C
//
// New binning strategy:
//   - N_TAIL = 6 equal-count bins per side, built as 3-thirds x 2-halves:
//       1) split tail into 3 equal-count thirds (at p33, p66)
//       2) split each third into 2 equal-count halves (at median)
//       → 6 bins all with equal counts, far tail well covered
//       tail region: outside 1.5-sigma for v2, outside 1.0-sigma for v3
//   - N_CORE = 10 equal-population (quantile) bins per side
//       core region: inside 1.5-sigma for v2, inside 1.0-sigma for v3
//   - N_TOTAL = 2*(N_TAIL+N_CORE) = 32 bins total (fixed for all cen/pT)
//   - Mirror symmetry preserved

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFileCollection.h"
#include "TLegend.h"

#include "TMath.h"
#include "TComplex.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include <TNtuple.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <algorithm>
#include <iomanip>
#include <cmath>

using namespace std;

// -------------------------------------------------------
// ensure bin edges are strictly increasing
// -------------------------------------------------------
static void fix_bin_edges(Double_t *edges, int nBins, Double_t eps = 1e-6) {
    for (int i = 1; i <= nBins; ++i) {
        if (!(edges[i] > edges[i-1]))
            edges[i] = edges[i-1] + eps;
    }
}

static void fill_missing_edges(Double_t *edges, int nBins,
                                Double_t sentinel = -1e9,
                                Double_t eps = 1e-12) {
    std::vector<int> idx;
    for (int i = 0; i <= nBins; ++i)
        if (edges[i] != sentinel) idx.push_back(i);
    if (idx.empty()) {
        for (int i = 0; i <= nBins; ++i) edges[i] = (Double_t)i * eps;
        return;
    }
    int first = idx.front();
    for (int i = 0; i < first; ++i)
        edges[i] = edges[first] - (first - i) * eps;
    for (size_t k = 0; k + 1 < idx.size(); ++k) {
        int a = idx[k], b = idx[k+1];
        Double_t va = edges[a], vb = edges[b];
        int gap = b - a;
        for (int i = 1; i < gap; ++i)
            edges[a + i] = va + (vb - va) * (Double_t(i) / gap);
    }
    int last = idx.back();
    for (int i = last + 1; i <= nBins; ++i)
        edges[i] = edges[last] + (i - last) * eps;
    fix_bin_edges(edges, nBins, eps);
}

// -------------------------------------------------------
// Build tail edges: split tail [tail_min, sigma_lo] into
// n_thirds equal-count thirds, then split each third into
// 2 equal-count halves → n_tail = 2*n_thirds bins, all
// with approximately equal counts.
// N_TAIL must be even (6 = 3 thirds x 2 halves).
// -------------------------------------------------------
static void make_tail_edges(Double_t *edges,
                             int       n_tail,    // must be even (6)
                             Double_t  tail_min,
                             Double_t  sigma_lo,
                             TH1D     *h_full)
{
    const int n_thirds = n_tail / 2;   // 3

    // -------------------------------------------------------
    // Step 1: zero out histogram outside tail region
    // -------------------------------------------------------
    TH1D *h_tail = (TH1D*)h_full->Clone("h_tail_tmp");
    for (int ib = 1; ib <= h_tail->GetNbinsX(); ib++) {
        double bc = h_tail->GetBinCenter(ib);
        if (bc < tail_min || bc > sigma_lo)
            h_tail->SetBinContent(ib, 0);
    }
    h_tail->SetEntries(h_tail->Integral());

    double total_tail = h_tail->Integral();
    std::cout << Form("    Tail [%.4f, %.4f]: total_counts=%.0f\n",
                      tail_min, sigma_lo, total_tail);

    // -------------------------------------------------------
    // Step 2: find third boundaries (p33, p66)
    // -------------------------------------------------------
    Double_t third_probs[n_thirds+1], third_q[n_thirds+1];
    for (int i = 0; i <= n_thirds; i++)
        third_probs[i] = (double)i / n_thirds;   // 0, 0.333, 0.667, 1.0
    h_tail->GetQuantiles(n_thirds+1, third_q, third_probs);
    third_q[0]        = tail_min;   // force exact outer boundaries
    third_q[n_thirds] = sigma_lo;

    std::cout << "    Third boundaries:";
    for (int i = 0; i <= n_thirds; i++)
        std::cout << Form(" %.4f", third_q[i]);
    std::cout << "\n";

    // -------------------------------------------------------
    // Step 3: split each third into 2 equal-count halves
    // -------------------------------------------------------
    edges[0]      = tail_min;   // forced
    edges[n_tail] = sigma_lo;   // forced

    for (int t = 0; t < n_thirds; t++) {
        double lo = third_q[t];
        double hi = third_q[t+1];

        TH1D *h_third = (TH1D*)h_full->Clone(
                            Form("h_third_%d_tmp", t));
        for (int ib = 1; ib <= h_third->GetNbinsX(); ib++) {
            double bc = h_third->GetBinCenter(ib);
            if (bc < lo || bc > hi)
                h_third->SetBinContent(ib, 0);
        }
        h_third->SetEntries(h_third->Integral());

        Double_t med_prob[1] = {0.5};
        Double_t med_q[1]    = {0.0};
        h_third->GetQuantiles(1, med_q, med_prob);

        // clamp median within [lo, hi]
        if (med_q[0] <= lo) med_q[0] = lo + 1e-6;
        if (med_q[0] >= hi) med_q[0] = hi - 1e-6;

        edges[2*t]   = third_q[t];   // left boundary of third t
        edges[2*t+1] = med_q[0];     // median splits third t into 2

        std::cout << Form("    third %d: [%8.4f, %8.4f]  median=%.4f"
                          "  → bin%d:[%.4f,%.4f]  bin%d:[%.4f,%.4f]\n",
                          t, lo, hi, med_q[0],
                          2*t,   edges[2*t],   edges[2*t+1],
                          2*t+1, edges[2*t+1], third_q[t+1]);

        delete h_third;
    }

    // -------------------------------------------------------
    // Step 4: diagnostic — print all n_tail bins
    // -------------------------------------------------------
    std::cout << "    Final tail edges (3-thirds x 2-halves):\n";
    for (int i = 0; i < n_tail; i++)
        std::cout << Form("      bin %d: [%8.4f, %8.4f]  width=%7.4f  [third %d]\n",
                          i, edges[i], edges[i+1],
                          edges[i+1] - edges[i], i/2);

    delete h_tail;
}


int save_vnbinning_to_txt(int target_cen_group = -1, int target_pt = -1)
{
    // -------------------------------------------------------
    // Output filename
    // -------------------------------------------------------
    std::string outname;
    const char *env_job  = std::getenv("SLURM_ARRAY_JOB_ID");
    if (!env_job) env_job = std::getenv("SLURM_JOB_ID");
    const char *env_task = std::getenv("SLURM_ARRAY_TASK_ID");
    if (target_cen_group >= 0 && target_pt >= 0)
        outname = TString::Format("vnbinning_new_c%d_pt%d", target_cen_group, target_pt).Data();
    else if (target_cen_group >= 0)
        outname = TString::Format("vnbinning_new_c%d_all", target_cen_group).Data();
    else if (target_pt >= 0)
        outname = TString::Format("vnbinning_new_all_pt%d", target_pt).Data();
    else
        outname = "vnbinning_new_all";
    if (env_job)  outname += std::string("_job")  + env_job;
    if (env_task) outname += std::string("_task") + env_task;
    outname += ".txt";

    // -------------------------------------------------------
    // Fixed binning parameters — array sizes NEVER change
    // -------------------------------------------------------
    const int    N_TAIL      = 6;       // tail bins per side (3 thirds x 2 halves, equal-count)
    const int    N_CORE      = 10;      // core bins per side (equal-quantile)
    const int    N_NEG       = N_TAIL + N_CORE;      // 16 bins on negative side
    const int    N_VBINS     = 2 * N_NEG;            // 32 total bins (fixed)
    // array edge count = N_VBINS + 1 = 33

    const int    N_HBINS     = 5000;    // histogram resolution for quantile finding
    const double TAIL_MIN    = -15.0;   // absolute outer edge of SP distribution
    const double CORE_MAX    =  0.0;    // negative core ends at zero (centre)

    const int N_CENTBINS = 6;
    Int_t  cen_edges[N_CENTBINS+1] = {0,10,20,30,40,50,80};

    const int N_PTBINS = 12;
    const char *pt_label[N_PTBINS] = {"1to2","2to3","3to4","4to5","5to6","6to8",
                                       "8to10","10to15","15to20","20to40","40to60","60to100"};
    double pt_edges[N_PTBINS+1] = {1,2,3,4,5,6,8,10,15,20,40,60,100};

    // storage: always [N_CENTBINS][N_PTBINS][N_VBINS+1]
    Double_t vnbinning_v2[N_CENTBINS][N_PTBINS][N_VBINS+1];
    Double_t vnbinning_v3[N_CENTBINS][N_PTBINS][N_VBINS+1];

    // -------------------------------------------------------
    // Open input ntuple
    // -------------------------------------------------------
    TFile *f = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/SP_MB11to21_Feb20_v4/ROOT/Flow_combined.root");
    if (!f || f->IsZombie()) { std::cerr << "Cannot open ntuple file\n"; return 1; }
    TNtuple *nt = (TNtuple*)f->Get("nt");
    if (!nt) { std::cerr << "Cannot find TNtuple 'nt'\n"; return 1; }

    // branch variables
    Float_t cen_val, q2_val, q3_val, pT_val, mass_val, v2_val, v3_val, dca_val, y_val;
    nt->SetBranchAddress("cent",         &cen_val);
    nt->SetBranchAddress("pT",           &pT_val);
    nt->SetBranchAddress("mass",         &mass_val);
    nt->SetBranchAddress("dca",          &dca_val);
    nt->SetBranchAddress("y",            &y_val);
    nt->SetBranchAddress("q2_hf_total",  &q2_val);
    nt->SetBranchAddress("q3_hf_total",  &q3_val);
    nt->SetBranchAddress("v2",           &v2_val);
    nt->SetBranchAddress("v3",           &v3_val);

    // histograms reused per (cen, pT) bin
    // full range for 1-sigma detection
    TH1D *h_v2_full = new TH1D("h_v2_full","h_v2_full", N_HBINS, TAIL_MIN, -TAIL_MIN);
    TH1D *h_v3_full = new TH1D("h_v3_full","h_v3_full", N_HBINS, TAIL_MIN, -TAIL_MIN);
    // negative core (filled after sigma_lo is known)
    TH1D *h_v2_core = new TH1D("h_v2_core","h_v2_core", N_HBINS, TAIL_MIN, CORE_MAX);
    TH1D *h_v3_core = new TH1D("h_v3_core","h_v3_core", N_HBINS, TAIL_MIN, CORE_MAX);

    TString cuts, cuts_core_v2, cuts_core_v3;

    for (int ic = 0; ic < N_CENTBINS; ++ic) {
        if (target_cen_group >= 0 && ic != target_cen_group) continue;

        for (int ip = 0; ip < N_PTBINS; ++ip) {
            if (target_pt >= 0 && ip != target_pt) continue;

            std::cout << "\n====================================================\n";
            std::cout << Form("  cen=%d-%d  pT=%s\n",
                              cen_edges[ic], cen_edges[ic+1], pt_label[ip]);
            std::cout << "====================================================\n";

            // base selection cuts
            cuts = TString::Format(
                "dca<0.0085 && mass<2.0 && mass>1.74 && y> -1 && y < 1.0 "
                "&& cent>=%d && cent<%d && pT>=%g && pT<%g",
                cen_edges[ic], cen_edges[ic+1],
                pt_edges[ip],  pt_edges[ip+1]);

            // --------------------------------------------------
            // Step 1: Fill full-range histograms
            // --------------------------------------------------
            h_v2_full->Reset();
            h_v3_full->Reset();
            nt->Draw("v2>>h_v2_full", cuts, "goff");
            nt->Draw("v3>>h_v3_full", cuts, "goff");

            Int_t n_v2 = (Int_t)h_v2_full->GetEntries();
            Int_t n_v3 = (Int_t)h_v3_full->GetEntries();

            if (n_v2 == 0)
                std::cerr << "WARNING: zero v2 entries for ic=" << ic
                          << " ip=" << ip << std::endl;
            if (n_v3 == 0)
                std::cerr << "WARNING: zero v3 entries for ic=" << ic
                          << " ip=" << ip << std::endl;

            // initialise all edges to sentinel
            const Double_t SENT = -1e9;
            for (int ib = 0; ib <= N_VBINS; ++ib) {
                vnbinning_v2[ic][ip][ib] = SENT;
                vnbinning_v3[ic][ip][ib] = SENT;
            }

            // --------------------------------------------------
            // Step 2: Find core/tail boundary edges
            //   v2: 1.5-sigma (86.64%) — v2 is near-Gaussian, use wider core
            //   v3: 1.0-sigma (68.27%) — v3 is broader/non-Gaussian, use narrower core
            // --------------------------------------------------
            const double P_LOW_V2  = (1.0 - 0.8664) / 2.0;  // 1.5-sigma = 0.06681
            const double P_HIGH_V2 = 1.0 - P_LOW_V2;
            const double P_LOW_V3  = (1.0 - 0.6827) / 2.0;  // 1.0-sigma = 0.15865
            const double P_HIGH_V3 = 1.0 - P_LOW_V3;

            // default fallback if histogram empty
            double v2_sig_lo = -2.3, v3_sig_lo = -2.3;

            if (n_v2 > 0) {
                Double_t probs_v2[2] = {P_LOW_V2, P_HIGH_V2};
                Double_t qv2[2]      = {0, 0};
                h_v2_full->GetQuantiles(2, qv2, probs_v2);
                v2_sig_lo = qv2[0];
                if (v2_sig_lo < TAIL_MIN) v2_sig_lo = TAIL_MIN;
                if (v2_sig_lo > -0.01)   v2_sig_lo = -0.01;
                std::cout << Form("  v2 1.5-sigma edge: [%.4f, %.4f]\n",
                                  v2_sig_lo, -v2_sig_lo);
            }
            if (n_v3 > 0) {
                Double_t probs_v3[2] = {P_LOW_V3, P_HIGH_V3};
                Double_t qv3[2]      = {0, 0};
                h_v3_full->GetQuantiles(2, qv3, probs_v3);
                v3_sig_lo = qv3[0];
                if (v3_sig_lo < TAIL_MIN) v3_sig_lo = TAIL_MIN;
                if (v3_sig_lo > -0.01)   v3_sig_lo = -0.01;
                std::cout << Form("  v3 1.0-sigma edge: [%.4f, %.4f]\n",
                                  v3_sig_lo, -v3_sig_lo);
            }

            // --------------------------------------------------
            // Step 3: Build equal-count tail edges (3-thirds x 2-halves)
            // --------------------------------------------------
            Double_t tail_edges_v2[N_TAIL+1], tail_edges_v3[N_TAIL+1];

            std::cout << "  v2 tail:\n";
            if (n_v2 > 0)
                make_tail_edges(tail_edges_v2, N_TAIL, TAIL_MIN, v2_sig_lo, h_v2_full);
            std::cout << "  v3 tail:\n";
            if (n_v3 > 0)
                make_tail_edges(tail_edges_v3, N_TAIL, TAIL_MIN, v3_sig_lo, h_v3_full);

            // --------------------------------------------------
            // Step 5: Equal-quantile core edges on negative half
            //         [v2_sig_lo, 0]
            // --------------------------------------------------
            // Set histogram range dynamically
            h_v2_core->Reset();
            h_v2_core->SetBins(N_HBINS, v2_sig_lo - 1e-6, CORE_MAX + 1e-6);
            h_v3_core->Reset();
            h_v3_core->SetBins(N_HBINS, v3_sig_lo - 1e-6, CORE_MAX + 1e-6);

            cuts_core_v2 = cuts + TString::Format(" && v2>=%g && v2<=0", v2_sig_lo);
            cuts_core_v3 = cuts + TString::Format(" && v3>=%g && v3<=0", v3_sig_lo);

            nt->Draw("v2>>h_v2_core", cuts_core_v2, "goff");
            nt->Draw("v3>>h_v3_core", cuts_core_v3, "goff");

            // Adaptive N_CORE_EFF: reduce if stats are low (but N_TOTAL stays fixed)
            // We still fill N_CORE edges; if stats low, equal-quantile compresses them
            Double_t core_probs[N_CORE+1];
            Double_t core_q_v2[N_CORE+1], core_q_v3[N_CORE+1];
            for (int i = 0; i <= N_CORE; i++)
                core_probs[i] = (double)i / N_CORE;

            // initialise to uniformly-spaced fallback
            for (int i = 0; i <= N_CORE; i++) {
                core_q_v2[i] = v2_sig_lo + i * (-v2_sig_lo) / N_CORE;
                core_q_v3[i] = v3_sig_lo + i * (-v3_sig_lo) / N_CORE;
            }

            if (h_v2_core->GetEntries() > N_CORE) {
                h_v2_core->GetQuantiles(N_CORE+1, core_q_v2, core_probs);
                std::cout << Form("  v2 core quantiles (N=%d):",
                                  (int)h_v2_core->GetEntries());
                for (int i = 0; i <= N_CORE; i++)
                    std::cout << Form(" %.3f", core_q_v2[i]);
                std::cout << "\n";
            } else {
                std::cout << "  v2 core: low stats, using uniform fallback edges\n";
            }

            if (h_v3_core->GetEntries() > N_CORE) {
                h_v3_core->GetQuantiles(N_CORE+1, core_q_v3, core_probs);
                std::cout << Form("  v3 core quantiles (N=%d):",
                                  (int)h_v3_core->GetEntries());
                for (int i = 0; i <= N_CORE; i++)
                    std::cout << Form(" %.3f", core_q_v3[i]);
                std::cout << "\n";
            } else {
                std::cout << "  v3 core: low stats, using uniform fallback edges\n";
            }

            // --------------------------------------------------
            // Step 6: Assemble negative-half edges
            // Layout: [TAIL_MIN, ...N_TAIL..., sig_lo, ...N_CORE..., 0.0]
            //  index:     0      1 .. N_TAIL   N_TAIL  N_TAIL+1 .. N_NEG
            // --------------------------------------------------
            // negative half: indices 0 .. N_NEG  (N_NEG+1 = 16 values)
            // tail portion: tail_edges_v2[0..N_TAIL], skip last (= sig_lo, covered by core)
            // core portion: core_q_v2[0..N_CORE]

            Double_t neg_v2[N_NEG+1], neg_v3[N_NEG+1];

            for (int i = 0; i < N_TAIL; i++) {
                neg_v2[i] = tail_edges_v2[i];   // [TAIL_MIN .. sigma_lo)
                neg_v3[i] = tail_edges_v3[i];
            }
            for (int i = 0; i <= N_CORE; i++) {
                neg_v2[N_TAIL + i] = core_q_v2[i];  // [sigma_lo .. 0]
                neg_v3[N_TAIL + i] = core_q_v3[i];
            }
            // force exact boundaries
            neg_v2[0]     = TAIL_MIN;
            neg_v3[0]     = TAIL_MIN;
            neg_v2[N_NEG] = 0.0;
            neg_v3[N_NEG] = 0.0;

            // --------------------------------------------------
            // Step 7: Mirror to positive side
            // full array: [TAIL_MIN, ...neg..., 0, ...pos..., +15]
            //   index 0       .. N_NEG   = negative half (including 0 at centre)
            //   index N_NEG+1 .. N_VBINS = positive half (mirror)
            // --------------------------------------------------
            for (int i = 0; i <= N_NEG; i++) {
                vnbinning_v2[ic][ip][i] = neg_v2[i];
                vnbinning_v3[ic][ip][i] = neg_v3[i];
            }
            for (int i = 1; i <= N_NEG; i++) {
                vnbinning_v2[ic][ip][N_NEG + i] = -neg_v2[N_NEG - i];
                vnbinning_v3[ic][ip][N_NEG + i] = -neg_v3[N_NEG - i];
            }

            // --------------------------------------------------
            // Step 8: Safety — fix monotonicity
            // --------------------------------------------------
            fill_missing_edges(vnbinning_v2[ic][ip], N_VBINS);
            fill_missing_edges(vnbinning_v3[ic][ip], N_VBINS);

            // --------------------------------------------------
            // Step 9: Sentinel check
            // --------------------------------------------------
            for (int ib = 0; ib <= N_VBINS; ++ib) {
                if (vnbinning_v2[ic][ip][ib] < -999.0 ||
                    vnbinning_v3[ic][ip][ib] < -999.0) {
                    std::cerr << "ERROR: sentinel survived at ic=" << ic
                              << " ip=" << ip << " ib=" << ib << std::endl;
                    return 1;
                }
            }

            // diagnostic: print full edge array
            std::cout << Form("  Full v2 edges (%d bins):\n", N_VBINS);
            for (int i = 0; i <= N_VBINS; i++) {
                std::string label;
                if      (i < N_TAIL)                            label = "tail-";
                else if (i == N_NEG)                            label = "center";
                else if (i > N_NEG && i <= N_NEG + N_CORE)     label = "core+";
                else if (i > N_NEG + N_CORE)                   label = "tail+";
                else                                            label = "core-";
                std::cout << Form("    edge[%2d] = %9.4f  [%s]\n",
                                  i, vnbinning_v2[ic][ip][i], label.c_str());
            }

        }  // end ip loop
    }  // end ic loop

    // -------------------------------------------------------
    // Write output text file
    // -------------------------------------------------------
    std::ofstream out(outname.c_str());
    if (!out.is_open()) { std::cerr << "Cannot open output file\n"; return 1; }
    out << std::fixed << std::setprecision(6);

    auto print_array_line = [&](std::ofstream &os, const char *type,
                                int ic, int ip, Double_t arr[]) {
        os << type << "_cent" << cen_edges[ic] << "to" << cen_edges[ic+1]
           << "_pT" << pt_label[ip] << "={";
        for (int ib = 0; ib <= N_VBINS; ++ib) {
            os << arr[ib];
            if (ib < N_VBINS) os << ", ";
        }
        os << "}\n";
    };

    if (target_cen_group >= 0 && target_pt >= 0) {
        print_array_line(out, "v2", target_cen_group, target_pt,
                         vnbinning_v2[target_cen_group][target_pt]);
        print_array_line(out, "v3", target_cen_group, target_pt,
                         vnbinning_v3[target_cen_group][target_pt]);
    } else {
        for (int ic = 0; ic < N_CENTBINS; ++ic) {
            if (target_cen_group >= 0 && ic != target_cen_group) continue;
            for (int ip = 0; ip < N_PTBINS; ++ip) {
                if (target_pt >= 0 && ip != target_pt) continue;
                print_array_line(out, "v2", ic, ip, vnbinning_v2[ic][ip]);
                print_array_line(out, "v3", ic, ip, vnbinning_v3[ic][ip]);
            }
        }
    }
    out.close();
    std::cout << "\nSaved vnbinning to " << outname << "\n";
    std::cout << Form("Array structure: N_TAIL=%d  N_CORE=%d"
                      "  N_NEG=%d  N_VBINS=%d  edges_per_bin=%d\n",
                      N_TAIL, N_CORE, N_NEG, N_VBINS, N_VBINS+1);
    f->Close();
    return 0;
}


int main(int argc, char **argv)
{
    if (argc == 1) {
        save_vnbinning_to_txt();
    } else if (argc == 2) {
        int idx = atoi(argv[1]);
        if (idx < 0 || idx > 6) {
            std::cerr << "ERROR: target_cen_group must be 0-6, got " << idx << std::endl;
            return 1;
        }
        save_vnbinning_to_txt(idx);
    } else if (argc == 3) {
        int cen_idx = atoi(argv[1]);
        int pt_idx  = atoi(argv[2]);
        if (cen_idx < 0 || cen_idx > 6) {
            std::cerr << "ERROR: target_cen_group must be 0-6, got " << cen_idx << std::endl;
            return 1;
        }
        if (pt_idx < 0 || pt_idx > 11) {
            std::cerr << "ERROR: target_pt must be 0-11, got " << pt_idx << std::endl;
            return 1;
        }
        save_vnbinning_to_txt(cen_idx, pt_idx);
    } else {
        std::cout << "Usage: ./save_vnbinning_to_txt_new [cen_group] [pt_index]\n";
        return 1;
    }
    return 0;
}
