#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <chrono>
#include <ctime>
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

using namespace std;

void flow_Analysis_chg(TString input_txt, TString output_path, int istart, int iend)
{
    const int N_CENTBINS_1 = 50;
    const int N_CENTBINS = 6;
    const int min_centbin[N_CENTBINS] = {0, 5, 10, 20, 30, 40}; // in 1% units
    const int max_centbin[N_CENTBINS] = {5, 10, 20, 30, 40, 50};
    const TString cent_label[N_CENTBINS] = {"cent0to5", "cent5to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"};

    const double MAX_PT_ANA = 20.0;
    const double MIN_PT_ANA = 0.5;
    const double MAX_ETA_ANA = 1.0;

    // pT binning for inclusive profiles (variable bins, wider range)
    const int N_PTBINS_INCL = 18;
    const double pt_edges_incl[N_PTBINS_INCL + 1] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 12.5, 15.0, 20.0};

    // Combined inclusive resolution per broad cent bin
    TH1D *hres_v2_incl_plus;
    TH1D *hres_v2_incl_minus;
    TH1D *hres_v3_incl_plus;
    TH1D *hres_v3_incl_minus;
    TH1D *hres_v2_incl_eff;
    TH1D *hres_v3_incl_eff;

    // Q-correlation diagnostic histograms
    TH2D *hQ2_trk_vs_HFp[N_CENTBINS];
    TH2D *hQ2_trk_vs_HFm[N_CENTBINS];
    TH2D *hQ2_trk_vs_HFpHm[N_CENTBINS];
    TH2D *hQ3_trk_vs_HFp[N_CENTBINS];
    TH2D *hQ3_trk_vs_HFm[N_CENTBINS];
    TH2D *hQ3_trk_vs_HFpHm[N_CENTBINS];

    TH2D *hQ2_trk_vs_HFp_all;
    TH2D *hQ2_trk_vs_HFm_all;
    TH2D *hQ2_trk_vs_HFpHm_all;
    TH2D *hQ3_trk_vs_HFp_all;
    TH2D *hQ3_trk_vs_HFm_all;
    TH2D *hQ3_trk_vs_HFpHm_all;

    // v2/v3 vs pT per broad cent bin — for inclusive results
    TProfile *hp_v2_incl[N_CENTBINS];
    TProfile *hp_v3_incl[N_CENTBINS];

    // 1D distributions for v2/res and v3/res
    TH1D *hv2_dist_incl[N_CENTBINS];
    TH1D *hv3_dist_incl[N_CENTBINS];

    // ---------------------------------------------------------------
    // Load inclusive resolution file
    // ---------------------------------------------------------------
    //w/ Era 2023 corrections
    //auto file_res_incl = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Charge_Resolution_INCLUSIVE_MB0to1_July21_v2/ROOT/Resolution_incl_out_combined.root");
    
    //w/out Era 2023 corrections
    //auto file_res_incl = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/INCLUSIVE_Charge_Resolution_MB0to1_July22_WO_Era/ROOT/Resolution_incl_out_combined.root");

    //w/ Era data+ w/o Era calibration
    //auto file_res_incl = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/INCLUSIVE_Charge_Resolution_MB0to1_July22_W_AND_WO_Era/ROOT/Resolution_incl_out_combined.root");

    //w/ Runwise EP calibration
    auto file_res_incl = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/INCLUSIVE_Charge_Resolution_MB0to1_wEra_Aug6/ROOT/Resolution_incl_out_combined.root");
    
    if (!file_res_incl || file_res_incl->IsZombie()) {
        std::cerr << "Error: Could not open inclusive resolution file!" << std::endl;
        return;
    }

    Double_t v2_den_incl_plus[N_CENTBINS];
    Double_t v2_den_incl_minus[N_CENTBINS];
    Double_t v2_den_incl_eff[N_CENTBINS];
    Double_t v3_den_incl_plus[N_CENTBINS];
    Double_t v3_den_incl_minus[N_CENTBINS];
    Double_t v3_den_incl_eff[N_CENTBINS];

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
        v2_den_incl_eff[ib] = TMath::Sqrt(v2_den_incl_plus[ib] * v2_den_incl_minus[ib]);
        
        v3_den_incl_plus[ib] = TMath::Sqrt((hQ3_mm->GetMean() * hQ3_mt->GetMean()) / hQ3_pt->GetMean());
        v3_den_incl_minus[ib] = TMath::Sqrt((hQ3_mm->GetMean() * hQ3_pt->GetMean()) / hQ3_mt->GetMean());
        v3_den_incl_eff[ib] = TMath::Sqrt(v3_den_incl_plus[ib] * v3_den_incl_minus[ib]);
    }
    file_res_incl->Close();

    // ---------------------------------------------------------------
    // Output File and Directory Setup
    // ---------------------------------------------------------------
    ifstream file_stream(input_txt.Data());
    TString outfile = TString::Format("%s/ROOT/flow_Analysis_chg_incl_%d_%d.root", output_path.Data(), istart, iend);
    TFile *fout = new TFile(outfile, "RECREATE");

    TDirectory *dir_incl_prof = fout->mkdir("Inclusive_TProfile");
    TDirectory *dir_dist = fout->mkdir("SP_Distributions");
    TDirectory *dir_checks = fout->mkdir("Diagnostics");

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

    dir_dist->cd();
    for (int ib = 0; ib < N_CENTBINS; ++ib)
    {
        TString nd2 = TString::Format("hv2_dist_incl_cent%dto%d", min_centbin[ib], max_centbin[ib]);
        TString nd3 = TString::Format("hv3_dist_incl_cent%dto%d", min_centbin[ib], max_centbin[ib]);
        
        // 2000 bins from -10 to 10 for detailed shape tracking
        hv2_dist_incl[ib] = new TH1D(nd2, nd2 + ";v_{2}^{obs} / R_{2};Tracks", 2000, -10.0, 10.0);
        hv3_dist_incl[ib] = new TH1D(nd3, nd3 + ";v_{3}^{obs} / R_{3};Tracks", 2000, -10.0, 10.0);
        
        hv2_dist_incl[ib]->Sumw2();
        hv3_dist_incl[ib]->Sumw2();
    }

    dir_checks->cd();
    hres_v2_incl_plus  = new TH1D("hres_v2_incl_plus",  "v2 res+ vs broad cent;centrality;R_{2}^{+}", N_CENTBINS, 0, N_CENTBINS);
    hres_v2_incl_minus = new TH1D("hres_v2_incl_minus", "v2 res- vs broad cent;centrality;R_{2}^{-}", N_CENTBINS, 0, N_CENTBINS);
    hres_v3_incl_plus  = new TH1D("hres_v3_incl_plus",  "v3 res+ vs broad cent;centrality;R_{3}^{+}", N_CENTBINS, 0, N_CENTBINS);
    hres_v3_incl_minus = new TH1D("hres_v3_incl_minus", "v3 res- vs broad cent;centrality;R_{3}^{-}", N_CENTBINS, 0, N_CENTBINS);
    hres_v2_incl_eff   = new TH1D("hres_v2_incl_eff",   "v2 eff res (geom. mean) vs broad cent;centrality;R_{2}^{eff}", N_CENTBINS, 0, N_CENTBINS);
    hres_v3_incl_eff   = new TH1D("hres_v3_incl_eff",   "v3 eff res (geom. mean) vs broad cent;centrality;R_{3}^{eff}", N_CENTBINS, 0, N_CENTBINS);
    
    hres_v2_incl_plus->Sumw2(); hres_v2_incl_minus->Sumw2(); 
    hres_v3_incl_plus->Sumw2(); hres_v3_incl_minus->Sumw2();
    hres_v2_incl_eff->Sumw2(); hres_v3_incl_eff->Sumw2();

    for (int ib = 0; ib < N_CENTBINS; ++ib) {
        hres_v2_incl_plus->GetXaxis()->SetBinLabel(ib+1,  cent_label[ib].Data());
        hres_v2_incl_minus->GetXaxis()->SetBinLabel(ib+1, cent_label[ib].Data());
        hres_v3_incl_plus->GetXaxis()->SetBinLabel(ib+1,  cent_label[ib].Data());
        hres_v3_incl_minus->GetXaxis()->SetBinLabel(ib+1, cent_label[ib].Data());
        hres_v2_incl_eff->GetXaxis()->SetBinLabel(ib+1,   cent_label[ib].Data());
        hres_v3_incl_eff->GetXaxis()->SetBinLabel(ib+1,   cent_label[ib].Data());

        hres_v2_incl_plus->SetBinContent(ib+1,  v2_den_incl_plus[ib]);
        hres_v2_incl_minus->SetBinContent(ib+1, v2_den_incl_minus[ib]);
        hres_v3_incl_plus->SetBinContent(ib+1,  v3_den_incl_plus[ib]);
        hres_v3_incl_minus->SetBinContent(ib+1, v3_den_incl_minus[ib]);
        hres_v2_incl_eff->SetBinContent(ib+1,   v2_den_incl_eff[ib]);
        hres_v3_incl_eff->SetBinContent(ib+1,   v3_den_incl_eff[ib]);
    }

    for (int ib = 0; ib < N_CENTBINS; ++ib) {
        TString tag = cent_label[ib].Data();
        hQ3_trk_vs_HFp[ib] = new TH2D(TString::Format("hQ3_trk_vs_HFp_%s", tag.Data()), TString::Format("|Q3| Trk vs |Q3| HF+ %s;|Q3 Trk|;|Q3 HF+|", tag.Data()), 200, 0, 500, 200, 0, 500);
        hQ3_trk_vs_HFm[ib] = new TH2D(TString::Format("hQ3_trk_vs_HFm_%s", tag.Data()), TString::Format("|Q3| Trk vs |Q3| HF- %s;|Q3 Trk|;|Q3 HF-|", tag.Data()), 200, 0, 500, 200, 0, 500);
        hQ3_trk_vs_HFpHm[ib] = new TH2D(TString::Format("hQ3_trk_vs_HFpHm_%s", tag.Data()), TString::Format("|Q3| Trk vs |Q3(HF+ + HF-)| %s;|Q3 Trk|;|Q3(HF+ + HF-)|", tag.Data()), 200, 0, 500, 200, 0, 500);

        hQ2_trk_vs_HFp[ib] = new TH2D(TString::Format("hQ2_trk_vs_HFp_%s", tag.Data()), TString::Format("Q2 Trk vs Q2 HF+ %s;Re(Q2 Trk);Re(Q2 HF+)", tag.Data()), 200, 0, 500, 200, 0, 500);
        hQ2_trk_vs_HFm[ib] = new TH2D(TString::Format("hQ2_trk_vs_HFm_%s", tag.Data()), TString::Format("Q2 Trk vs Q2 HF- %s;Re(Q2 Trk);Re(Q2 HF-)", tag.Data()), 200, 0, 500, 200, 0, 500);
        hQ2_trk_vs_HFpHm[ib] = new TH2D(TString::Format("hQ2_trk_vs_HFpHm_%s", tag.Data()), TString::Format("Q2 Trk vs Q2 HF+HF- %s;Re(Q2 Trk);Re(Q2 (HF+ + HF-))", tag.Data()), 200, 0, 500, 200, 0, 500);
        
        hQ3_trk_vs_HFp[ib]->Sumw2(); hQ3_trk_vs_HFm[ib]->Sumw2(); hQ3_trk_vs_HFpHm[ib]->Sumw2();
        hQ2_trk_vs_HFp[ib]->Sumw2(); hQ2_trk_vs_HFm[ib]->Sumw2(); hQ2_trk_vs_HFpHm[ib]->Sumw2();
    }

    hQ3_trk_vs_HFp_all = new TH2D("hQ3_trk_vs_HFp_all","|Q3| Trk vs |Q3| HF+ (all cent);|Q3 Trk|;|Q3 HF+|",200,0,500,200,0,500);
    hQ3_trk_vs_HFm_all = new TH2D("hQ3_trk_vs_HFm_all","|Q3| Trk vs |Q3| HF- (all cent);|Q3 Trk|;|Q3 HF-|",200,0,500,200,0,500);
    hQ3_trk_vs_HFpHm_all = new TH2D("hQ3_trk_vs_HFpHm_all","|Q3| Trk vs |Q3(HF+ + HF-)| (all cent);|Q3 Trk|;|Q3(HF+ + HF-)|",200,0,500,200,0,500);
    hQ2_trk_vs_HFp_all = new TH2D("hQ2_trk_vs_HFp_all","Q2 Trk vs Q2 HF+ (all cent)",200,0,500,200,0,500);
    hQ2_trk_vs_HFm_all = new TH2D("hQ2_trk_vs_HFm_all","Q2 Trk vs Q2 HF- (all cent)",200,0,500,200,0,500);
    hQ2_trk_vs_HFpHm_all = new TH2D("hQ2_trk_vs_HFpHm_all","Q2 Trk vs Q2 HF+HF- (all cent)",200,0,500,200,0,500);
    
    hQ3_trk_vs_HFp_all->Sumw2(); hQ3_trk_vs_HFm_all->Sumw2(); hQ3_trk_vs_HFpHm_all->Sumw2();
    hQ2_trk_vs_HFp_all->Sumw2(); hQ2_trk_vs_HFm_all->Sumw2(); hQ2_trk_vs_HFpHm_all->Sumw2();

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

        TTree *tree = (TTree *)fin->Get("Ana/ntInfo");

        Int_t centrality;
        Float_t ephfmQ[3], ephfpQ[3];
        Float_t ephfpAngle[3], ephfmAngle[3], eptkAngle[2], eptkQ[2];
        Int_t trk_mult;
        Float_t pT[MAXTRK], eta[MAXTRK], phi[MAXTRK];

        tree->SetBranchStatus("*", 0);
        for (const auto &p : {"centrality", "ephfpAngle", "ephfmAngle", "ephfpQ", "ephfmQ",
                              "eptkmidAngle", "eptkmidQ", "trk_mult", "pT", "eta", "phi"})
            tree->SetBranchStatus(p, 1);


	tree->SetBranchAddress("centrality", &centrality);
        tree->SetBranchAddress("ephfpAngle", ephfpAngle);
        tree->SetBranchAddress("ephfmAngle", ephfmAngle);
        tree->SetBranchAddress("ephfmQ", ephfmQ);
        tree->SetBranchAddress("ephfpQ", ephfpQ);
        tree->SetBranchAddress("eptkmidAngle", eptkAngle);
        tree->SetBranchAddress("eptkmidQ", eptkQ);
        tree->SetBranchAddress("trk_mult", &trk_mult);
        tree->SetBranchAddress("pT", pT);
        tree->SetBranchAddress("eta", eta);
        tree->SetBranchAddress("phi", phi);


        Int_t n_entries = tree->GetEntries();

        for (Long64_t ii = 0; ii < n_entries; ii++)
        {
            tree->GetEntry(ii);

            if (ii % 100000 == 0)
                printf("Processing entry %lld of %lld : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);

            Int_t cent = centrality / 2;
            if (cent < 0 || cent >= N_CENTBINS_1)
                continue;

            int broad_cent_idx = -1;
            for (int ib = 0; ib < N_CENTBINS; ++ib)
            {
                if (cent >= min_centbin[ib] && cent < max_centbin[ib])
                {
                    broad_cent_idx = ib;
                    break;
                }
            }
            
            // Skip event if it doesn't fall into the broad centrality bins of interest
            if (broad_cent_idx < 0) continue; 

            TComplex aux_Q2_HFm(ephfmQ[1] * TMath::Cos(2.0 * ephfmAngle[1]), ephfmQ[1] * TMath::Sin(2.0 * ephfmAngle[1]), 0);
            TComplex aux_Q2_HFp(ephfpQ[1] * TMath::Cos(2.0 * ephfpAngle[1]), ephfpQ[1] * TMath::Sin(2.0 * ephfpAngle[1]), 0);
            TComplex aux_Q2_Trk(eptkQ[0] * TMath::Cos(2.0 * eptkAngle[0]), eptkQ[0] * TMath::Sin(2.0 * eptkAngle[0]), 0);

            TComplex aux_Q3_HFm(ephfmQ[2] * TMath::Cos(3.0 * ephfmAngle[2]), ephfmQ[2] * TMath::Sin(3.0 * ephfmAngle[2]), 0);
            TComplex aux_Q3_HFp(ephfpQ[2] * TMath::Cos(3.0 * ephfpAngle[2]), ephfpQ[2] * TMath::Sin(3.0 * ephfpAngle[2]), 0);
            TComplex aux_Q3_Trk(eptkQ[1] * TMath::Cos(3.0 * eptkAngle[1]), eptkQ[1] * TMath::Sin(3.0 * eptkAngle[1]), 0);

            Double_t res2_incl_plus = v2_den_incl_plus[broad_cent_idx];
            Double_t res2_incl_minus = v2_den_incl_minus[broad_cent_idx];
            Double_t res3_incl_plus = v3_den_incl_plus[broad_cent_idx];
            Double_t res3_incl_minus = v3_den_incl_minus[broad_cent_idx];

            double q3_trk_mod   = TComplex::Abs(aux_Q3_Trk);
            double q3_hfp_mod   = TComplex::Abs(aux_Q3_HFp);
            double q3_hfm_mod   = TComplex::Abs(aux_Q3_HFm);
            double q3_hfphm_mod = TComplex::Abs(aux_Q3_HFp + aux_Q3_HFm);
            
            double q2_trk_mod  = TComplex::Abs(aux_Q2_Trk);
            double q2_hfp_mod  = TComplex::Abs(aux_Q2_HFp);
            double q2_hfm_mod  = TComplex::Abs(aux_Q2_HFm);
            double q2_hfphm_mod = TComplex::Abs(aux_Q2_HFp + aux_Q2_HFm);

            // Populate global correlations
            hQ3_trk_vs_HFp_all->Fill(q3_trk_mod, q3_hfp_mod);
            hQ3_trk_vs_HFm_all->Fill(q3_trk_mod, q3_hfm_mod);
            hQ3_trk_vs_HFpHm_all->Fill(q3_trk_mod, q3_hfphm_mod);

            hQ2_trk_vs_HFp_all->Fill(q2_trk_mod, q2_hfp_mod);
            hQ2_trk_vs_HFm_all->Fill(q2_trk_mod, q2_hfm_mod);
            hQ2_trk_vs_HFpHm_all->Fill(q2_trk_mod, q2_hfphm_mod);

            // Populate per-centrality correlations
            hQ3_trk_vs_HFp[broad_cent_idx]->Fill(q3_trk_mod, q3_hfp_mod);
            hQ3_trk_vs_HFm[broad_cent_idx]->Fill(q3_trk_mod, q3_hfm_mod);
            hQ3_trk_vs_HFpHm[broad_cent_idx]->Fill(q3_trk_mod, q3_hfphm_mod);

            hQ2_trk_vs_HFp[broad_cent_idx]->Fill(q2_trk_mod, q2_hfp_mod);
            hQ2_trk_vs_HFm[broad_cent_idx]->Fill(q2_trk_mod, q2_hfm_mod);
            hQ2_trk_vs_HFpHm[broad_cent_idx]->Fill(q2_trk_mod, q2_hfphm_mod);

            for (int itrk = 0; itrk < trk_mult; ++itrk)
            {
                float Eta = eta[itrk];
                float Pt = pT[itrk];
                float Phi = phi[itrk];

                if (Pt >= MAX_PT_ANA || Pt < MIN_PT_ANA || fabs(Eta) >= MAX_ETA_ANA) continue;

                TComplex aux_Q2_trk(TMath::Cos(2.0 * Phi), TMath::Sin(2.0 * Phi), 0);
                TComplex aux_Q3_trk(TMath::Cos(3.0 * Phi), TMath::Sin(3.0 * Phi), 0);
                
                Double_t v2_obs, v3_obs;
                Double_t res2_incl_final, res3_incl_final;

                if (Eta > 0)
                {
                    v2_obs = (aux_Q2_trk * TComplex::Conjugate(aux_Q2_HFm)).Re();
                    v3_obs = (aux_Q3_trk * TComplex::Conjugate(aux_Q3_HFm)).Re();
                    res2_incl_final = res2_incl_plus;
                    res3_incl_final = res3_incl_plus;
                }
                else
                {
                    v2_obs = (aux_Q2_trk * TComplex::Conjugate(aux_Q2_HFp)).Re();
                    v3_obs = (aux_Q3_trk * TComplex::Conjugate(aux_Q3_HFp)).Re();
                    res2_incl_final = res2_incl_minus;
                    res3_incl_final = res3_incl_minus;
                }

                if (res2_incl_final > 0)
                {
                    double SP_v2 = v2_obs / res2_incl_final;
                    hp_v2_incl[broad_cent_idx]->Fill(Pt, SP_v2);
                    hv2_dist_incl[broad_cent_idx]->Fill(SP_v2); // Fill the 1D distribution
                }
                if (res3_incl_final > 0)
                {
                    double SP_v3 = v3_obs / res3_incl_final;
                    hp_v3_incl[broad_cent_idx]->Fill(Pt, SP_v3);
                    hv3_dist_incl[broad_cent_idx]->Fill(SP_v3); // Fill the 1D distribution
                }

            } 
        } 

        fin->Close();
        delete fin;
        ifile++;
    } 

    fout->Write(0, TObject::kOverwrite);
    fout->Close();
}

int main(int argc, char *argv[])
{
    if (argc == 5)
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
