#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFileCollection.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TComplex.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "THnSparse.h"
#include "TNtuple.h"
#include <climits>
#include <iomanip>

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/charge_vn/check_setup/20quantile_cuts_2023_MB0to1_charge_check_Jun11.h"

using namespace std;

void flow_Analysis_chg(TString input_txt, TString output_path, int istart, int iend)
{

    const int N_QBINS = 20;

    const int N_CENTBINS = 10;
    const int min_centbin[N_CENTBINS] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45}; // in 1% units
    const int max_centbin[N_CENTBINS] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50};
    const TString cent_label[N_CENTBINS] = {"cent0to5", "cent5to10", "cent10to15", "cent15to20", "cent20to25", "cent25to30", "cent30to35", "cent35to40", "cent40to45", "cent45to50"};

    const double MAX_PT_ANA = 20.0;
    const double MIN_PT_ANA = 0.5;
    const double MAX_ETA_ANA = 1.0;

    // pT binning for ESE profiles (0.5–3.0 GeV, 25 x 0.1 GeV bins)
    const int N_PTBINS_ESE = 10;
    const double PT_MIN_ESE = 0.5;
    const double PT_MAX_ESE = 3.0;

    // v2/v3 vs pT per (1% cent, q-bin) — for ESE, comparison with D0 vn
    TProfile *hp_v2[N_CENTBINS][N_QBINS];
    TProfile *hp_v3[N_CENTBINS][N_QBINS];
    TH1D *hnevt_q2[N_CENTBINS][N_QBINS];
    TH1D *hnevt_q3[N_CENTBINS][N_QBINS];

    TH2D *hQ2_trk_vs_HFp[N_CENTBINS];
    TH2D *hQ2_trk_vs_HFm[N_CENTBINS];
    TH2D *hQ2_trk_vs_HFpHm[N_CENTBINS];
    TH2D *hQ3_trk_vs_HFp[N_CENTBINS];
    TH2D *hQ3_trk_vs_HFm[N_CENTBINS];
    TH2D *hQ3_trk_vs_HFpHm[N_CENTBINS];

    // Exact resolution denominators per cent vs q-bin (2D overview)
    TH2D *hres_v2_plus_cent_vs_q;
    TH2D *hres_v2_minus_cent_vs_q;
    TH2D *hres_v3_plus_cent_vs_q;
    TH2D *hres_v3_minus_cent_vs_q;
    TH2D *hres_v2_eff_cent_vs_q; // geometric mean: sqrt(plus * minus)
    TH2D *hres_v3_eff_cent_vs_q;
    // Combined (geometric mean) resolution vs q-bin, one TH1D per 1% cent bin
    TH1D *hres_v2_eff_vsq[N_CENTBINS];
    TH1D *hres_v3_eff_vsq[N_CENTBINS];

    // v2/v3 vs q2/q3 bin per broad cent bin — integrated over 0.5-3.0 GeV/c
    TProfile *hp_v2_vsq2_pt0p5to3[N_CENTBINS];
    TProfile *hp_v3_vsq3_pt0p5to3[N_CENTBINS];

    // v2/v3 vs q2/q3 bin per broad cent bin — integrated over 1.0-3.0 GeV/c
    TProfile *hp_v2_vsq2_pt1to3[N_CENTBINS];
    TProfile *hp_v3_vsq3_pt1to3[N_CENTBINS];

    // q2/q3 distributions per (cent bin, q-bin) — to extract <q> per quantile
    TH1D *hq2_dist[N_CENTBINS][N_QBINS];
    TH1D *hq3_dist[N_CENTBINS][N_QBINS];

    // Accumulators for mean-q TProfile filling (filled after file loop)
    /*Double_t sum_v2_pt0p5to3[N_CENTBINS][N_QBINS] = {};
    Double_t sum_v3_pt0p5to3[N_CENTBINS][N_QBINS] = {};
    Double_t sum_v2_pt1to3[N_CENTBINS][N_QBINS] = {};
    Double_t sum_v3_pt1to3[N_CENTBINS][N_QBINS] = {};
    Long64_t cnt_v2_pt0p5to3[N_CENTBINS][N_QBINS] = {};
    Long64_t cnt_v3_pt0p5to3[N_CENTBINS][N_QBINS] = {};
    Long64_t cnt_v2_pt1to3[N_CENTBINS][N_QBINS] = {};
    Long64_t cnt_v3_pt1to3[N_CENTBINS][N_QBINS] = {};*/

    auto file_res = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Charge_Resolution_MB0to1_Jun11_ForCheck/ROOT/Resolution_chg_qbin_out_combined.root");
    if (!file_res || file_res->IsZombie())
    {
        std::cerr << "ERROR: Cannot open resolution file" << std::endl;
        return;
    }
    Double_t v2_den_Dy_plus[N_CENTBINS][N_QBINS];
    Double_t v2_den_Dy_minus[N_CENTBINS][N_QBINS];
    Double_t v2_den_Dy_eff[N_CENTBINS][N_QBINS]; // geometric mean: sqrt(plus * minus)

    Double_t v3_den_Dy_plus[N_CENTBINS][N_QBINS];
    Double_t v3_den_Dy_minus[N_CENTBINS][N_QBINS];
    Double_t v3_den_Dy_eff[N_CENTBINS][N_QBINS]; // geometric mean: sqrt(plus * minus)

    Float_t q2_hf_total = 0.0 , q3_hf_total = 0.0;

    TH1D *hQ2Q2_HFmHFp[N_CENTBINS][N_QBINS];
    TH1D *hQ2Q2_HFmTrk[N_CENTBINS][N_QBINS];
    TH1D *hQ2Q2_HFpTrk[N_CENTBINS][N_QBINS];
    TH1D *hQ3Q3_HFmHFp[N_CENTBINS][N_QBINS];
    TH1D *hQ3Q3_HFmTrk[N_CENTBINS][N_QBINS];
    TH1D *hQ3Q3_HFpTrk[N_CENTBINS][N_QBINS];
    TString hname;

    for (Int_t i = 0; i < N_CENTBINS; i++)
    {
        for (Int_t j = 0; j < N_QBINS; j++)
        {
            hname = "Q2Q2_HFmHFp_Re_cen" + to_string(i) + "_" + to_string(i + 1) + "_q2bin_" + to_string(j);
            hQ2Q2_HFmHFp[i][j] = (TH1D *)file_res->Get(hname);
            hname = "Q2Q2_HFmTrk_Re_cen" + to_string(i) + "_" + to_string(i + 1) + "_q2bin_" + to_string(j);
            hQ2Q2_HFmTrk[i][j] = (TH1D *)file_res->Get(hname);
            hname = "Q2Q2_HFpTrk_Re_cen" + to_string(i) + "_" + to_string(i + 1) + "_q2bin_" + to_string(j);
            hQ2Q2_HFpTrk[i][j] = (TH1D *)file_res->Get(hname);

            hname = "Q3Q3_HFmHFp_Re_cen" + to_string(i) + "_" + to_string(i + 1) + "_q3bin_" + to_string(j);
            hQ3Q3_HFmHFp[i][j] = (TH1D *)file_res->Get(hname);
            hname = "Q3Q3_HFmTrk_Re_cen" + to_string(i) + "_" + to_string(i + 1) + "_q3bin_" + to_string(j);
            hQ3Q3_HFmTrk[i][j] = (TH1D *)file_res->Get(hname);
            hname = "Q3Q3_HFpTrk_Re_cen" + to_string(i) + "_" + to_string(i + 1) + "_q3bin_" + to_string(j);
            hQ3Q3_HFpTrk[i][j] = (TH1D *)file_res->Get(hname);

            v2_den_Dy_plus[i][j] = TMath::Sqrt((hQ2Q2_HFmHFp[i][j]->GetMean() * hQ2Q2_HFmTrk[i][j]->GetMean()) / (hQ2Q2_HFpTrk[i][j]->GetMean()));
            v2_den_Dy_minus[i][j] = TMath::Sqrt((hQ2Q2_HFmHFp[i][j]->GetMean() * hQ2Q2_HFpTrk[i][j]->GetMean()) / (hQ2Q2_HFmTrk[i][j]->GetMean()));
            v2_den_Dy_eff[i][j] = TMath::Sqrt(v2_den_Dy_plus[i][j] * v2_den_Dy_minus[i][j]);

            v3_den_Dy_plus[i][j] = TMath::Sqrt((hQ3Q3_HFmHFp[i][j]->GetMean() * hQ3Q3_HFmTrk[i][j]->GetMean()) / (hQ3Q3_HFpTrk[i][j]->GetMean()));
            v3_den_Dy_minus[i][j] = TMath::Sqrt((hQ3Q3_HFmHFp[i][j]->GetMean() * hQ3Q3_HFpTrk[i][j]->GetMean()) / (hQ3Q3_HFmTrk[i][j]->GetMean()));
            v3_den_Dy_eff[i][j] = TMath::Sqrt(v3_den_Dy_plus[i][j] * v3_den_Dy_minus[i][j]);
        }
    }
    file_res->Close();

    ifstream file_stream(input_txt.Data());
    TString outfile = TString::Format("%s/ROOT/flow_Analysis_chg_%d_%d.root", output_path.Data(), istart, iend);
    TFile *fout = new TFile(outfile, "RECREATE");

    // Create output subdirectories first, then book histograms inside them
    // so ROOT registers them in the correct directory (not at the file root).
    TDirectory *dir_ese_prof = fout->mkdir("ESE_TProfile");
    TDirectory *dir_ntrk = fout->mkdir("NTrack_QBin");
    TDirectory *dir_qdist = fout->mkdir("Qdist_perBin");
    TDirectory *dir_vsq = fout->mkdir("vsQbin_TProfile");
    TDirectory *dir_vsq_meanq = fout->mkdir("vsQmean_TProfile");
    TDirectory *dir_vsq_hdrmean = fout->mkdir("vsQ_HeaderMean_TProfile");

    dir_ese_prof->cd();
    for (int ic = 0; ic < N_CENTBINS; ++ic)
        for (int iq = 0; iq < N_QBINS; ++iq)
        {
            TString n2 = TString::Format("hp_v2_cen%d_q2bin%d", ic, iq);
            TString n3 = TString::Format("hp_v3_cen%d_q3bin%d", ic, iq);
            hp_v2[ic][iq] = new TProfile(n2, n2, N_PTBINS_ESE, PT_MIN_ESE, PT_MAX_ESE);
            hp_v3[ic][iq] = new TProfile(n3, n3, N_PTBINS_ESE, PT_MIN_ESE, PT_MAX_ESE);
            hp_v2[ic][iq]->Sumw2();
            hp_v3[ic][iq]->Sumw2();
        }

    dir_vsq->cd();
    for (int ib = 0; ib < N_CENTBINS; ++ib)
    {
        TString nq2 = TString::Format("hp_v2_vsq2_pt0p5to3_%s", cent_label[ib].Data());
        TString nq3 = TString::Format("hp_v3_vsq3_pt0p5to3_%s", cent_label[ib].Data());
        hp_v2_vsq2_pt0p5to3[ib] = new TProfile(nq2, nq2, N_QBINS, 0, N_QBINS);
        hp_v3_vsq3_pt0p5to3[ib] = new TProfile(nq3, nq3, N_QBINS, 0, N_QBINS);
        hp_v2_vsq2_pt0p5to3[ib]->Sumw2();
        hp_v3_vsq3_pt0p5to3[ib]->Sumw2();

        TString nq2_pt1 = TString::Format("hp_v2_vsq2_pt1to3_%s", cent_label[ib].Data());
        TString nq3_pt1 = TString::Format("hp_v3_vsq3_pt1to3_%s", cent_label[ib].Data());
        hp_v2_vsq2_pt1to3[ib] = new TProfile(nq2_pt1, nq2_pt1, N_QBINS, 0, N_QBINS);
        hp_v3_vsq3_pt1to3[ib] = new TProfile(nq3_pt1, nq3_pt1, N_QBINS, 0, N_QBINS);
        hp_v2_vsq2_pt1to3[ib]->Sumw2();
        hp_v3_vsq3_pt1to3[ib]->Sumw2();
    }

    dir_ntrk->cd();
    for (int ic = 0; ic < N_CENTBINS; ++ic)
        for (int iq = 0; iq < N_QBINS; ++iq)
        {
            TString nq2 = TString::Format("hNtrk_q2_cen%d_q2bin%d", ic, iq);
            TString nq3 = TString::Format("hNtrk_q3_cen%d_q3bin%d", ic, iq);
            hnevt_q2[ic][iq] = new TH1D(nq2, nq2, N_PTBINS_ESE, PT_MIN_ESE, PT_MAX_ESE);
            hnevt_q3[ic][iq] = new TH1D(nq3, nq3, N_PTBINS_ESE, PT_MIN_ESE, PT_MAX_ESE);
        }

    dir_qdist->cd();
    for (int ic = 0; ic < N_CENTBINS; ++ic)
        for (int iq = 0; iq < N_QBINS; ++iq)
        {
            // x-range: use the quantile edges for this (cent, q-bin)
            // so the histogram naturally captures the right range
            double xlo2 = q2_cuts[ic][iq];
            double xhi2 = q2_cuts[ic][iq + 1];
            double xlo3 = q3_cuts[ic][iq];
            double xhi3 = q3_cuts[ic][iq + 1];

            TString nd2 = TString::Format("hq2_dist_cen%d_q2bin%d", ic, iq);
            TString nd3 = TString::Format("hq3_dist_cen%d_q3bin%d", ic, iq);
            hq2_dist[ic][iq] = new TH1D(nd2, nd2, 100, xlo2, xhi2);
            hq3_dist[ic][iq] = new TH1D(nd3, nd3, 100, xlo3, xhi3);
            hq2_dist[ic][iq]->Sumw2();
            hq3_dist[ic][iq]->Sumw2();
        }

    dir_vsq_meanq->cd();

    TProfile *hp_v2_vsmeanq2_pt0p5to3[N_CENTBINS];
    TProfile *hp_v3_vsmeanq3_pt0p5to3[N_CENTBINS];
    TProfile *hp_v2_vsmeanq2_pt1to3[N_CENTBINS];
    TProfile *hp_v3_vsmeanq3_pt1to3[N_CENTBINS];

    for (int ic = 0; ic < N_CENTBINS; ++ic)
    {
        // x-axis spans full q range for this cent bin
        double xlo2 = q2_cuts[ic][0];
        double xhi2 = q2_cuts[ic][N_QBINS];
        double xlo3 = q3_cuts[ic][0];
        double xhi3 = q3_cuts[ic][N_QBINS];

        TString nm2a = TString::Format("hp_v2_vsmeanq2_pt0p5to3_%s", cent_label[ic].Data());
        TString nm3a = TString::Format("hp_v3_vsmeanq3_pt0p5to3_%s", cent_label[ic].Data());
        TString nm2b = TString::Format("hp_v2_vsmeanq2_pt1to3_%s", cent_label[ic].Data());
        TString nm3b = TString::Format("hp_v3_vsmeanq3_pt1to3_%s", cent_label[ic].Data());

        hp_v2_vsmeanq2_pt0p5to3[ic] = new TProfile(nm2a, nm2a, N_QBINS, q2_cuts[ic]);
        hp_v3_vsmeanq3_pt0p5to3[ic] = new TProfile(nm3a, nm3a, N_QBINS, q3_cuts[ic]);
        hp_v2_vsmeanq2_pt1to3[ic] = new TProfile(nm2b, nm2b, N_QBINS, q2_cuts[ic]);
        hp_v3_vsmeanq3_pt1to3[ic] = new TProfile(nm3b, nm3b, N_QBINS, q3_cuts[ic]);
        hp_v2_vsmeanq2_pt0p5to3[ic]->Sumw2();
        hp_v3_vsmeanq3_pt0p5to3[ic]->Sumw2();
        hp_v2_vsmeanq2_pt1to3[ic]->Sumw2();
        hp_v3_vsmeanq3_pt1to3[ic]->Sumw2();
    }

    dir_vsq_hdrmean->cd();

    TProfile *hp_v2_vshdrmeanq2_pt0p5to3[N_CENTBINS];
    TProfile *hp_v3_vshdrmeanq3_pt0p5to3[N_CENTBINS];
    TProfile *hp_v2_vshdrmeanq2_pt1to3[N_CENTBINS];
    TProfile *hp_v3_vshdrmeanq3_pt1to3[N_CENTBINS];

    for (int ic = 0; ic < N_CENTBINS; ++ic)
    {
        TString nm2a_hdr = TString::Format("hp_v2_vshdrmeanq2_pt0p5to3_%s", cent_label[ic].Data());
        TString nm3a_hdr = TString::Format("hp_v3_vshdrmeanq3_pt0p5to3_%s", cent_label[ic].Data());
        TString nm2b_hdr = TString::Format("hp_v2_vshdrmeanq2_pt1to3_%s", cent_label[ic].Data());
        TString nm3b_hdr = TString::Format("hp_v3_vshdrmeanq3_pt1to3_%s", cent_label[ic].Data());

        // Initialize using the variable bin edges from the header file
        hp_v2_vshdrmeanq2_pt0p5to3[ic] = new TProfile(nm2a_hdr, nm2a_hdr, N_QBINS, q2_cuts[ic]);
        hp_v3_vshdrmeanq3_pt0p5to3[ic] = new TProfile(nm3a_hdr, nm3a_hdr, N_QBINS, q3_cuts[ic]);
        hp_v2_vshdrmeanq2_pt1to3[ic] = new TProfile(nm2b_hdr, nm2b_hdr, N_QBINS, q2_cuts[ic]);
        hp_v3_vshdrmeanq3_pt1to3[ic] = new TProfile(nm3b_hdr, nm3b_hdr, N_QBINS, q3_cuts[ic]);
        hp_v2_vshdrmeanq2_pt0p5to3[ic]->Sumw2();
        hp_v3_vshdrmeanq3_pt0p5to3[ic]->Sumw2();
        hp_v2_vshdrmeanq2_pt1to3[ic]->Sumw2();
        hp_v3_vshdrmeanq3_pt1to3[ic]->Sumw2();
    }

    // Maximum number of tracks per event
    const int MAXTRK = 50000;

    string filename;
    int ifile = 0;

    while (file_stream >> filename)
    {
        if (ifile < istart)
        {
            ifile++;
            continue;
        }
        if (ifile >= iend)
            break;

        TFile *fin = TFile::Open(filename.c_str());

        if (!fin || fin->IsZombie())
        {
            std::cout << "Warning: Skipping bad file: " << filename << std::endl;
            if (fin)
            {
                fin->Close();
                delete fin;
            }
            ifile++;
            continue;
        }

        std::cout << ">>> Processing ifile=" << ifile << " : " << filename << std::endl;

        TTree *tree = (TTree *)fin->Get("Ana/ntEvtInfo");

        // Event-level variables
        Int_t centrality;
        Float_t ephfmQ[3], ephfpQ[3], ephfmSumW[3], ephfpSumW[3];
        Float_t ephfpAngle[3], ephfmAngle[3], eptkAngle[2], eptkQ[2];

        // Track-level variables (arrays, one entry per track per event)
        Int_t trk_mult;
        Float_t pT[MAXTRK], eta[MAXTRK], phi[MAXTRK];

        tree->SetBranchAddress("centrality", &centrality);
        tree->SetBranchAddress("ephfpAngle", ephfpAngle);
        tree->SetBranchAddress("ephfmAngle", ephfmAngle);
        tree->SetBranchAddress("ephfmQ", ephfmQ);
        tree->SetBranchAddress("ephfpQ", ephfpQ);
        tree->SetBranchAddress("ephfmSumW", ephfmSumW);
        tree->SetBranchAddress("ephfpSumW", ephfpSumW);
        tree->SetBranchAddress("eptkAngle", eptkAngle);
        tree->SetBranchAddress("eptkQ", eptkQ);
        tree->SetBranchAddress("trk_mult", &trk_mult);
        tree->SetBranchAddress("pT", pT);
        tree->SetBranchAddress("eta", eta);
        tree->SetBranchAddress("phi", phi);

        tree->SetBranchStatus("*", 0);
        for (const auto &p : {"centrality", "ephfpAngle", "ephfmAngle", "ephfpQ", "ephfmQ",
                              "eptkAngle", "eptkQ", "ephfmSumW", "ephfpSumW",
                              "trk_mult", "pT", "eta", "phi"})
            tree->SetBranchStatus(p, 1);

        Int_t n_entries = tree->GetEntries();

        for (Long64_t ii = 0; ii < n_entries; ii++)
        {
            tree->GetEntry(ii);

            if (ii % 100000 == 0)
                printf("Processing entry %lld of %lld : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);

            Int_t cent = centrality / 2;
            if (cent < 0 || cent >= 90)
                continue;

            TComplex aux_Q2_HFm(ephfmQ[1] * TMath::Cos(2.0 * ephfmAngle[1]), ephfmQ[1] * TMath::Sin(2.0 * ephfmAngle[1]), 0);
            TComplex aux_Q2_HFp(ephfpQ[1] * TMath::Cos(2.0 * ephfpAngle[1]), ephfpQ[1] * TMath::Sin(2.0 * ephfpAngle[1]), 0);
            TComplex aux_Q2_Trk(eptkQ[0] * TMath::Cos(2.0 * eptkAngle[0]), eptkQ[0] * TMath::Sin(2.0 * eptkAngle[0]), 0);

            TComplex aux_Q3_HFm(ephfmQ[2] * TMath::Cos(3.0 * ephfmAngle[2]), ephfmQ[2] * TMath::Sin(3.0 * ephfmAngle[2]), 0);
            TComplex aux_Q3_HFp(ephfpQ[2] * TMath::Cos(3.0 * ephfpAngle[2]), ephfpQ[2] * TMath::Sin(3.0 * ephfpAngle[2]), 0);
            TComplex aux_Q3_Trk(eptkQ[1] * TMath::Cos(3.0 * eptkAngle[1]), eptkQ[1] * TMath::Sin(3.0 * eptkAngle[1]), 0);

            int ese_cent_idx = -1;
            int q2_idx = -1, q3_idx = -1;
            Double_t res2_plus = 0, res2_minus = 0, res3_plus = 0, res3_minus = 0;

            for (int i_cen = 0; i_cen < N_CENTBINS; i_cen++)
            {
                if (cent >= min_centbin[i_cen] && cent < max_centbin[i_cen])
                {
                    ese_cent_idx = i_cen;

                    float sumW2 = ephfpSumW[1] + ephfmSumW[1];
                    float sumW3 = ephfpSumW[2] + ephfmSumW[2];
                    if (sumW2 <= 0 || sumW3 <= 0)
                        break;
                    q2_hf_total = (ephfpQ[1] + ephfmQ[1]) / sumW2;
                    q3_hf_total = (ephfpQ[2] + ephfmQ[2]) / sumW3;

                    q2_idx = TMath::BinarySearch(21, q2_cuts[i_cen], (double)q2_hf_total);
                    q3_idx = TMath::BinarySearch(21, q3_cuts[i_cen], (double)q3_hf_total);

                    // BinarySearch returns -1 if value < first edge: clamp to valid range
                    q2_idx = std::max(0, std::min(19, q2_idx));
                    q3_idx = std::max(0, std::min(19, q3_idx));

                    // Fill q distributions for <q> extraction
                    hq2_dist[ese_cent_idx][q2_idx]->Fill(q2_hf_total);
                    hq3_dist[ese_cent_idx][q3_idx]->Fill(q3_hf_total);

                    // Pre-fetch resolution factors for this event (ESE, q-binned)
                    res2_plus = v2_den_Dy_plus[i_cen][q2_idx];
                    res2_minus = v2_den_Dy_minus[i_cen][q2_idx];
                    res3_plus = v3_den_Dy_plus[i_cen][q3_idx];
                    res3_minus = v3_den_Dy_minus[i_cen][q3_idx];
                    break;
                }
            }

            double hdr_mq2 = 0.0;
            double hdr_mq3 = 0.0;
            if (ese_cent_idx >= 0 && q2_idx >= 0 && q3_idx >= 0) {
                hdr_mq2 = (q2_cuts[ese_cent_idx][q2_idx] + q2_cuts[ese_cent_idx][q2_idx + 1]) / 2.0;
                hdr_mq3 = (q3_cuts[ese_cent_idx][q3_idx] + q3_cuts[ese_cent_idx][q3_idx + 1]) / 2.0;
            }

            for (int itrk = 0; itrk < trk_mult; ++itrk)
            {
                float Eta = eta[itrk];
                float Pt = pT[itrk];
                float Phi = phi[itrk];

                if (Pt >= MAX_PT_ANA || Pt < MIN_PT_ANA || fabs(Eta) >= MAX_ETA_ANA)
                    continue;

                TComplex aux_Q2_trk(TMath::Cos(2.0 * Phi), TMath::Sin(2.0 * Phi), 0);
                TComplex aux_Q3_trk(TMath::Cos(3.0 * Phi), TMath::Sin(3.0 * Phi), 0);
                Double_t v2_obs, v3_obs;
                Double_t res2_final, res3_final;

                // Eta-gap: positive eta uses HFm reference, negative eta uses HFp reference
                if (Eta > 0)
                {
                    v2_obs = (aux_Q2_trk * TComplex::Conjugate(aux_Q2_HFm)).Re();
                    v3_obs = (aux_Q3_trk * TComplex::Conjugate(aux_Q3_HFm)).Re();
                    res2_final = res2_plus;
                    res3_final = res3_plus;
                }
                else
                {
                    v2_obs = (aux_Q2_trk * TComplex::Conjugate(aux_Q2_HFp)).Re();
                    v3_obs = (aux_Q3_trk * TComplex::Conjugate(aux_Q3_HFp)).Re();
                    res2_final = res2_minus;
                    res3_final = res3_minus;
                }

                if (ese_cent_idx >= 0 && q2_idx >= 0 && q3_idx >= 0)
                {
                    if (res2_final > 0)
                        hp_v2[ese_cent_idx][q2_idx]->Fill(Pt, v2_obs / res2_final);
                    if (res3_final > 0)
                        hp_v3[ese_cent_idx][q3_idx]->Fill(Pt, v3_obs / res3_final);

                    // Count tracks per (cent, q-bin) vs pT for diagnostics
                    hnevt_q2[ese_cent_idx][q2_idx]->Fill(Pt);
                    hnevt_q3[ese_cent_idx][q3_idx]->Fill(Pt);

                    // v2/v3 vs q-bin integrated over 0.5-3.0 GeV/c
                    if (Pt >= PT_MIN_ESE && Pt < PT_MAX_ESE)
                    {
                        

                        if (res2_final > 0)
                        {
                            hp_v2_vsq2_pt0p5to3[ese_cent_idx]->Fill(q2_idx, v2_obs / res2_final);
                            hp_v2_vsmeanq2_pt0p5to3[ese_cent_idx]->Fill(q2_hf_total, v2_obs / res2_final);
                            hp_v2_vshdrmeanq2_pt0p5to3[ese_cent_idx]->Fill(hdr_mq2, v2_obs / res2_final);
                            //sum_v2_pt0p5to3[ese_cent_idx][q2_idx] += v2_obs / res2_final;
                            //cnt_v2_pt0p5to3[ese_cent_idx][q2_idx]++;
                        }
                        if (res3_final > 0)
                        {
                            hp_v3_vsq3_pt0p5to3[ese_cent_idx]->Fill(q3_idx, v3_obs / res3_final);
                            hp_v3_vsmeanq3_pt0p5to3[ese_cent_idx]->Fill(q3_hf_total, v3_obs / res3_final);
                            hp_v3_vshdrmeanq3_pt0p5to3[ese_cent_idx]->Fill(hdr_mq3, v3_obs / res3_final);
                            //sum_v3_pt0p5to3[ese_cent_idx][q3_idx] += v3_obs / res3_final;
                            //cnt_v3_pt0p5to3[ese_cent_idx][q3_idx]++;
                        }
                    }

                    // v2/v3 vs q-bin integrated over 1.0-3.0 GeV/c
                    if (Pt >= 1.0 && Pt < PT_MAX_ESE)
                    {

                        if (res2_final > 0)
                        {
                            hp_v2_vsq2_pt1to3[ese_cent_idx]->Fill(q2_idx, v2_obs / res2_final);
                            hp_v2_vsmeanq2_pt1to3[ese_cent_idx]->Fill(q2_hf_total, v2_obs / res2_final);
                            hp_v2_vshdrmeanq2_pt1to3[ese_cent_idx]->Fill(hdr_mq2, v2_obs / res2_final);
                            //sum_v2_pt1to3[ese_cent_idx][q2_idx] += v2_obs / res2_final;
                            //cnt_v2_pt1to3[ese_cent_idx][q2_idx]++;
                        }
                        if (res3_final > 0)
                        {
                            hp_v3_vsq3_pt1to3[ese_cent_idx]->Fill(q3_idx, v3_obs / res3_final);
                            hp_v3_vsmeanq3_pt1to3[ese_cent_idx]->Fill(q3_hf_total, v3_obs / res3_final);
                            hp_v3_vshdrmeanq3_pt1to3[ese_cent_idx]->Fill(hdr_mq3, v3_obs / res3_final);
                            //sum_v3_pt1to3[ese_cent_idx][q3_idx] += v3_obs / res3_final;
                            //cnt_v3_pt1to3[ese_cent_idx][q3_idx]++;
                        }
                    }
                }
            } // -- track loop --
        } // -- Event loop ---

        fin->Close();
        delete fin;
        ifile++;
    } // -- file loop ---

    // Extract <q> per (cent, q-bin) from filled distributions
    /*Double_t q2_mean_data[N_CENTBINS][N_QBINS];
    Double_t q3_mean_data[N_CENTBINS][N_QBINS];

    for (int ic = 0; ic < N_CENTBINS; ++ic)
        for (int iq = 0; iq < N_QBINS; ++iq)
        {
            q2_mean_data[ic][iq] = hq2_dist[ic][iq]->GetMean();
            q3_mean_data[ic][iq] = hq3_dist[ic][iq]->GetMean();
        }*/

    /*dir_vsq_meanq->cd();

    TProfile *hp_v2_vsmeanq2_pt0p5to3[N_CENTBINS];
    TProfile *hp_v3_vsmeanq3_pt0p5to3[N_CENTBINS];
    TProfile *hp_v2_vsmeanq2_pt1to3[N_CENTBINS];
    TProfile *hp_v3_vsmeanq3_pt1to3[N_CENTBINS];

    for (int ic = 0; ic < N_CENTBINS; ++ic)
    {
        // x-axis spans full q range for this cent bin
        double xlo2 = q2_cuts[ic][0];
        double xhi2 = q2_cuts[ic][N_QBINS];
        double xlo3 = q3_cuts[ic][0];
        double xhi3 = q3_cuts[ic][N_QBINS];

        TString nm2a = TString::Format("hp_v2_vsmeanq2_pt0p5to3_%s", cent_label[ic].Data());
        TString nm3a = TString::Format("hp_v3_vsmeanq3_pt0p5to3_%s", cent_label[ic].Data());
        TString nm2b = TString::Format("hp_v2_vsmeanq2_pt1to3_%s", cent_label[ic].Data());
        TString nm3b = TString::Format("hp_v3_vsmeanq3_pt1to3_%s", cent_label[ic].Data());

        hp_v2_vsmeanq2_pt0p5to3[ic] = new TProfile(nm2a, nm2a, N_QBINS, q2_cuts[ic]);
        hp_v3_vsmeanq3_pt0p5to3[ic] = new TProfile(nm3a, nm3a, N_QBINS, q3_cuts[ic]);
        hp_v2_vsmeanq2_pt1to3[ic] = new TProfile(nm2b, nm2b, N_QBINS, q2_cuts[ic]);
        hp_v3_vsmeanq3_pt1to3[ic] = new TProfile(nm3b, nm3b, N_QBINS, q3_cuts[ic]);

        for (int iq = 0; iq < N_QBINS; ++iq)
        {
            double mq2 = q2_mean_data[ic][iq];
            double mq3 = q3_mean_data[ic][iq];

            if (cnt_v2_pt0p5to3[ic][iq] > 0)
                hp_v2_vsmeanq2_pt0p5to3[ic]->Fill(mq2, sum_v2_pt0p5to3[ic][iq] / cnt_v2_pt0p5to3[ic][iq]);
            if (cnt_v3_pt0p5to3[ic][iq] > 0)
                hp_v3_vsmeanq3_pt0p5to3[ic]->Fill(mq3, sum_v3_pt0p5to3[ic][iq] / cnt_v3_pt0p5to3[ic][iq]);
            if (cnt_v2_pt1to3[ic][iq] > 0)
                hp_v2_vsmeanq2_pt1to3[ic]->Fill(mq2, sum_v2_pt1to3[ic][iq] / cnt_v2_pt1to3[ic][iq]);
            if (cnt_v3_pt1to3[ic][iq] > 0)
                hp_v3_vsmeanq3_pt1to3[ic]->Fill(mq3, sum_v3_pt1to3[ic][iq] / cnt_v3_pt1to3[ic][iq]);
        }
    }*/

    //============================
    // Write all histograms + TGraphErrors
    //============================

    fout->Write(0, TObject::kOverwrite);
    fout->Close();
}

int main(int argc, char *argv[])
{
    if (argc == 5) // Expecting 4 arguments: input_txt, output_path, istart, iend
    {
        TString input_txt = argv[1];
        TString output_path = argv[2];
        int istart = std::stoi(argv[3]);
        int iend = std::stoi(argv[4]);

        flow_Analysis_chg(input_txt, output_path, istart, iend);
    }
    else
    {
        std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
        return 1;
    }
    return 0;
}
