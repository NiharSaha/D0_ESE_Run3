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

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/charge_vn/quantile_cuts_2023_MB0to1_charge.h"

using namespace std;

void flow_Analysis_chg(TString input_txt, TString output_path, int istart, int iend)
{

    const int N_QBINS = 10;
    const int N_CENTBINS_1 = 90;

    const int N_CENTBINS = 5;
    const int min_centbin[N_CENTBINS] = {0, 10, 20, 30, 40}; // in 1% units
    const int max_centbin[N_CENTBINS] = {10, 20, 30, 40, 50};
    const TString cent_label[N_CENTBINS] = {"cent0to10","cent10to20","cent20to30", "cent30to40", "cent40to50"};

    const double MAX_PT_ANA = 20.0;
    const double MIN_PT_ANA = 0.5;
    const double MAX_ETA_ANA = 1.0;

    // pT binning for ESE profiles (0.5–3.0 GeV, 25 x 0.1 GeV bins)
    const int    N_PTBINS_ESE = 25;
    const double PT_MIN_ESE   = 0.5;
    const double PT_MAX_ESE   = 3.0;

    // pT binning for inclusive profiles (variable bins, wider range)
    const int N_PTBINS_INCL = 18;
    const double pt_edges_incl[N_PTBINS_INCL+1] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 12.5, 15.0, 20.0};

    // v2/v3 vs pT per (1% cent, q-bin) — for ESE, comparison with D0 vn
    TProfile *hp_v2[N_CENTBINS_1][N_QBINS];
    TProfile *hp_v3[N_CENTBINS_1][N_QBINS];
    TH1D *hnevt_q2[N_CENTBINS_1][N_QBINS];
    TH1D *hnevt_q3[N_CENTBINS_1][N_QBINS];

    // v2/v3 vs pT per broad cent bin (q-integrated) — for comparison with published results
    TProfile *hp_v2_incl[N_CENTBINS];
    TProfile *hp_v3_incl[N_CENTBINS];
    TH1D *hnevt_incl[N_CENTBINS];

    // v2/v3 vs q2/q3 bin per broad cent bin — integrated over 0.5-3.0 GeV/c
    TProfile *hp_v2_vsq2_pt0p5to3[N_CENTBINS];
    TProfile *hp_v3_vsq3_pt0p5to3[N_CENTBINS];

    // v2/v3 vs q2/q3 bin per broad cent bin — integrated over 1.0-3.0 GeV/c
    TProfile *hp_v2_vsq2_pt1to3[N_CENTBINS];
    TProfile *hp_v3_vsq3_pt1to3[N_CENTBINS];

    auto file_res = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Resolution_MB0to1_Apr27_charge_pub_v0/ROOT/Resolution_chg_qbin_out_combined.root");

    Double_t v2_den_Dy_plus[N_CENTBINS_1][N_QBINS];
    Double_t v2_den_Dy_minus[N_CENTBINS_1][N_QBINS];

    Double_t v3_den_Dy_plus[N_CENTBINS_1][N_QBINS];
    Double_t v3_den_Dy_minus[N_CENTBINS_1][N_QBINS];

    Float_t q2_hf_total, q3_hf_total;

    TH1D *hQ2Q2_HFmHFp[N_CENTBINS_1][N_QBINS];
    TH1D *hQ2Q2_HFmTrk[N_CENTBINS_1][N_QBINS];
    TH1D *hQ2Q2_HFpTrk[N_CENTBINS_1][N_QBINS];
    TH1D *hQ3Q3_HFmHFp[N_CENTBINS_1][N_QBINS];
    TH1D *hQ3Q3_HFmTrk[N_CENTBINS_1][N_QBINS];
    TH1D *hQ3Q3_HFpTrk[N_CENTBINS_1][N_QBINS];
    TString hname;

    for (Int_t i = 0; i < N_CENTBINS_1; i++)
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

            v3_den_Dy_plus[i][j] = TMath::Sqrt((hQ3Q3_HFmHFp[i][j]->GetMean() * hQ3Q3_HFmTrk[i][j]->GetMean()) / (hQ3Q3_HFpTrk[i][j]->GetMean()));
            v3_den_Dy_minus[i][j] = TMath::Sqrt((hQ3Q3_HFmHFp[i][j]->GetMean() * hQ3Q3_HFpTrk[i][j]->GetMean()) / (hQ3Q3_HFmTrk[i][j]->GetMean()));
        }
    }
    file_res->Close();

    // ---------------------------------------------------------------
    // Load inclusive resolution file (no q-bin dimension)
    // Histogram naming: Q2Q2_HFmHFp_Re_cen{i}_{i+1}, etc.
    // ---------------------------------------------------------------
    auto file_res_incl = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Resolution_MB0to1_Apr15_charge_pub_v0/ROOT/Resolution_incl_out_combined.root");

    // Broad cent bin labels matching Calculate_Resolution_incl.C naming convention
    // e.g. "Q2Q2_HFmHFp_Re_cent0to10", "Q2Q2_HFmHFp_Re_cent10to30", ...
    Double_t v2_den_incl_plus[N_CENTBINS];
    Double_t v2_den_incl_minus[N_CENTBINS];
    Double_t v3_den_incl_plus[N_CENTBINS];
    Double_t v3_den_incl_minus[N_CENTBINS];

    for (Int_t ib = 0; ib < N_CENTBINS; ib++)
    {
        TString lbl = TString::Format("cent%dto%d", min_centbin[ib], max_centbin[ib]);
        TH1D *hQ2_mm = (TH1D *)file_res_incl->Get("Q2Q2_HFmHFp_Re_" + lbl);
        TH1D *hQ2_mt = (TH1D *)file_res_incl->Get("Q2Q2_HFmTrk_Re_" + lbl);
        TH1D *hQ2_pt = (TH1D *)file_res_incl->Get("Q2Q2_HFpTrk_Re_" + lbl);
        TH1D *hQ3_mm = (TH1D *)file_res_incl->Get("Q3Q3_HFmHFp_Re_" + lbl);
        TH1D *hQ3_mt = (TH1D *)file_res_incl->Get("Q3Q3_HFmTrk_Re_" + lbl);
        TH1D *hQ3_pt = (TH1D *)file_res_incl->Get("Q3Q3_HFpTrk_Re_" + lbl);

        v2_den_incl_plus[ib] = TMath::Sqrt((hQ2_mm->GetMean() * hQ2_mt->GetMean()) / hQ2_pt->GetMean());
        v2_den_incl_minus[ib] = TMath::Sqrt((hQ2_mm->GetMean() * hQ2_pt->GetMean()) / hQ2_mt->GetMean());
        v3_den_incl_plus[ib] = TMath::Sqrt((hQ3_mm->GetMean() * hQ3_mt->GetMean()) / hQ3_pt->GetMean());
        v3_den_incl_minus[ib] = TMath::Sqrt((hQ3_mm->GetMean() * hQ3_pt->GetMean()) / hQ3_mt->GetMean());
    }
    file_res_incl->Close();

    ifstream file_stream(input_txt.Data());
    TString outfile = TString::Format("%s/ROOT/flow_Analysis_chg_%d_%d.root", output_path.Data(), istart, iend);
    TFile *fout = new TFile(outfile, "RECREATE");

    

    // Create output subdirectories first, then book histograms inside them
    // so ROOT registers them in the correct directory (not at the file root).
    TDirectory *dir_ese_prof  = fout->mkdir("ESE_TProfile");
    TDirectory *dir_ese_gr    = fout->mkdir("ESE_TGraphErrors");
    TDirectory *dir_incl_prof = fout->mkdir("Inclusive_TProfile");
    TDirectory *dir_incl_gr   = fout->mkdir("Inclusive_TGraphErrors");
    TDirectory *dir_vsq       = fout->mkdir("vsQbin_TProfile");
    TDirectory *dir_ntrk      = fout->mkdir("NTrack_QBin");

    dir_ese_prof->cd();
    for (int ic = 0; ic < N_CENTBINS_1; ++ic)
        for (int iq = 0; iq < N_QBINS; ++iq)
        {
            TString n2 = TString::Format("hp_v2_cen%d_q2bin%d", ic, iq);
            TString n3 = TString::Format("hp_v3_cen%d_q3bin%d", ic, iq);
            hp_v2[ic][iq] = new TProfile(n2, n2, N_PTBINS_ESE, PT_MIN_ESE, PT_MAX_ESE);
            hp_v3[ic][iq] = new TProfile(n3, n3, N_PTBINS_ESE, PT_MIN_ESE, PT_MAX_ESE);
            hp_v2[ic][iq]->Sumw2();
            hp_v3[ic][iq]->Sumw2();
        }

    dir_incl_prof->cd();
    for (int ib = 0; ib < N_CENTBINS; ++ib)
    {
        TString ni2 = TString::Format("hp_v2_incl_cent%dto%d", min_centbin[ib], max_centbin[ib]);
        TString ni3 = TString::Format("hp_v3_incl_cent%dto%d", min_centbin[ib], max_centbin[ib]);
        hp_v2_incl[ib] = new TProfile(ni2, ni2, N_PTBINS_INCL, pt_edges_incl);
        hp_v3_incl[ib] = new TProfile(ni3, ni3, N_PTBINS_INCL, pt_edges_incl);
        hp_v2_incl[ib]->Sumw2();
        hp_v3_incl[ib]->Sumw2();
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
    for (int ic = 0; ic < N_CENTBINS_1; ++ic)
        for (int iq = 0; iq < N_QBINS; ++iq)
        {
            TString nq2 = TString::Format("hNtrk_q2_cen%d_q2bin%d", ic, iq);
            TString nq3 = TString::Format("hNtrk_q3_cen%d_q3bin%d", ic, iq);
            hnevt_q2[ic][iq] = new TH1D(nq2, nq2, N_PTBINS_ESE, PT_MIN_ESE, PT_MAX_ESE);
            hnevt_q3[ic][iq] = new TH1D(nq3, nq3, N_PTBINS_ESE, PT_MIN_ESE, PT_MAX_ESE);
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

            float sumW2 = ephfpSumW[1] + ephfmSumW[1];
            float sumW3 = ephfpSumW[2] + ephfmSumW[2];

            if (sumW2 <= 0 || sumW3 <= 0)
                continue;
            q2_hf_total = (ephfpQ[1] + ephfmQ[1]) / sumW2;
            q3_hf_total = (ephfpQ[2] + ephfmQ[2]) / sumW3;

            int q2_idx = TMath::BinarySearch(11, q2_cuts[cent], (double)q2_hf_total);
            int q3_idx = TMath::BinarySearch(11, q3_cuts[cent], (double)q3_hf_total);

            // BinarySearch returns -1 if value < first edge: clamp to valid range
            q2_idx = std::max(0, std::min(9, q2_idx));
            q3_idx = std::max(0, std::min(9, q3_idx));

            // Count events per (cent, q-bin)
            //hnevt_q2[cent][q2_idx]->Fill(0.5);
            //hnevt_q3[cent][q3_idx]->Fill(0.5);

            // Broad centrality bin index: 0=0-10%, 1=10-30%, 2=30-50%; -1 outside range
            int broad_cent_idx = -1;
            for (int ib = 0; ib < N_CENTBINS; ++ib)
                if (cent >= min_centbin[ib] && cent < max_centbin[ib])
                {
                    broad_cent_idx = ib;
                    break;
                }
            //if (broad_cent_idx >= 0)
            //    hnevt_incl[broad_cent_idx]->Fill(0.5);

            // Pre-fetch resolution factors for this event (ESE, q-binned)
            Double_t res2_plus = v2_den_Dy_plus[cent][q2_idx];
            Double_t res2_minus = v2_den_Dy_minus[cent][q2_idx];
            Double_t res3_plus = v3_den_Dy_plus[cent][q3_idx];
            Double_t res3_minus = v3_den_Dy_minus[cent][q3_idx];

            // Resolution factors for inclusive profiles (no q-bin), indexed by broad cent bin
            Double_t res2_incl_plus = (broad_cent_idx >= 0) ? v2_den_incl_plus[broad_cent_idx] : 0.0;
            Double_t res2_incl_minus = (broad_cent_idx >= 0) ? v2_den_incl_minus[broad_cent_idx] : 0.0;
            Double_t res3_incl_plus = (broad_cent_idx >= 0) ? v3_den_incl_plus[broad_cent_idx] : 0.0;
            Double_t res3_incl_minus = (broad_cent_idx >= 0) ? v3_den_incl_minus[broad_cent_idx] : 0.0;

            for (int itrk = 0; itrk < trk_mult; ++itrk)
            {
                float Eta = eta[itrk];
                float Pt = pT[itrk];
                float Phi = phi[itrk];

                if (Pt >= MAX_PT_ANA || Pt < MIN_PT_ANA || fabs(Eta) >= MAX_ETA_ANA) continue;

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

                if (res2_final > 0)
                    hp_v2[cent][q2_idx]->Fill(Pt, v2_obs / res2_final);
                if (res3_final > 0)
                    hp_v3[cent][q3_idx]->Fill(Pt, v3_obs / res3_final);

                // Count tracks per (cent, q-bin) vs pT for diagnostics
                hnevt_q2[cent][q2_idx]->Fill(Pt);
                hnevt_q3[cent][q3_idx]->Fill(Pt);

                // q-integrated profiles for comparison with published results (use inclusive resolution)
                if (broad_cent_idx >= 0)
                {
                    Double_t res2_incl_final = (Eta > 0) ? res2_incl_plus : res2_incl_minus;
                    Double_t res3_incl_final = (Eta > 0) ? res3_incl_plus : res3_incl_minus;
                    if (res2_incl_final > 0)
                        hp_v2_incl[broad_cent_idx]->Fill(Pt, v2_obs / res2_incl_final);
                    if (res3_incl_final > 0)
                        hp_v3_incl[broad_cent_idx]->Fill(Pt, v3_obs / res3_incl_final);

                    // v2/v3 vs q-bin integrated over 0.5-3.0 GeV/c
                    if (Pt >= PT_MIN_ESE && Pt < PT_MAX_ESE)
                    {
                        if (res2_final > 0)
                            hp_v2_vsq2_pt0p5to3[broad_cent_idx]->Fill(q2_idx + 0.5, v2_obs / res2_final);
                        if (res3_final > 0)
                            hp_v3_vsq3_pt0p5to3[broad_cent_idx]->Fill(q3_idx + 0.5, v3_obs / res3_final);
                    }

                    // v2/v3 vs q-bin integrated over 1.0-3.0 GeV/c
                    if (Pt >= 1.0 && Pt < PT_MAX_ESE)
                    {
                        if (res2_final > 0)
                            hp_v2_vsq2_pt1to3[broad_cent_idx]->Fill(q2_idx + 0.5, v2_obs / res2_final);
                        if (res3_final > 0)
                            hp_v3_vsq3_pt1to3[broad_cent_idx]->Fill(q3_idx + 0.5, v3_obs / res3_final);
                    }
                }

            } // -- track loop --
        } // -- Event loop ---

        fin->Close();
        delete fin;
        ifile++;
    } // -- file loop ---

    //============================
    // Write all histograms + TGraphErrors
    //============================

    // Helper: convert a TProfile to a TGraphErrors (skips empty bins)
    auto ProfileToGraph = [](TProfile *hp, const char *gname) -> TGraphErrors *
    {
        int n = hp->GetNbinsX();
        std::vector<double> vx, vy, vex, vey;
        for (int ib = 1; ib <= n; ++ib)
        {
            if (hp->GetBinEntries(ib) > 0)
            {
                vx.push_back(hp->GetBinCenter(ib));
                vy.push_back(hp->GetBinContent(ib));
                vex.push_back(0.0);
                vey.push_back(hp->GetBinError(ib));
            }
        }
        TGraphErrors *gr = new TGraphErrors((Int_t)vx.size(), vx.data(), vy.data(), vex.data(), vey.data());
        gr->SetName(gname);
        return gr;
    };

    // --- ESE TGraphErrors ---
    /*dir_ese_gr->cd();
    for (int ic = 0; ic < N_CENTBINS_1; ++ic)
        for (int iq = 0; iq < N_QBINS; ++iq)
        {
            TGraphErrors *gr2 = ProfileToGraph(hp_v2[ic][iq], TString::Format("gr_v2_cen%d_q2bin%d", ic, iq));
            TGraphErrors *gr3 = ProfileToGraph(hp_v3[ic][iq], TString::Format("gr_v3_cen%d_q3bin%d", ic, iq));
            gr2->Write(); delete gr2;
            gr3->Write(); delete gr3;
        }

    // --- Inclusive TGraphErrors ---
    dir_incl_gr->cd();
    for (int ib = 0; ib < N_CENTBINS; ++ib)
    {
        TGraphErrors *gr2i = ProfileToGraph(hp_v2_incl[ib], TString::Format("gr_v2_incl_cent%dto%d", min_centbin[ib], max_centbin[ib]));
        TGraphErrors *gr3i = ProfileToGraph(hp_v3_incl[ib], TString::Format("gr_v3_incl_cent%dto%d", min_centbin[ib], max_centbin[ib]));
        gr2i->Write(); delete gr2i;
        gr3i->Write(); delete gr3i;
    }*/

    // --- vsQbin TProfiles (already owned by dir_vsq, flushed by fout->Write) ---
    // No explicit Write() needed; fout->Write() below recurses into subdirectories.

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
