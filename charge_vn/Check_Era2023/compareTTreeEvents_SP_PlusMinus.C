/*
 * compareTTreeEvents_SP_AllHarmonics.C
 * Compares Event Plane angles (HF+, HF-, Tracker) and Q-vector components (Qx, Qy, |Q|) for n=2 and n=3.
 * Computes SP Resolutions and saves all distributions directly to a ROOT file.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TComplex.h"
#include "TFile.h"
#include <map>        
#include <utility>


using namespace std;

// Helper to properly wrap angular differences based on harmonic 'n'
float getDeltaPsi(float angNew, float angOld, float n)
{
    float dPsi = angNew - angOld;
    float period = TMath::TwoPi() / n;
    float halfPeriod = period / 2.0;
    while (dPsi > halfPeriod)
        dPsi -= period;
    while (dPsi < -halfPeriod)
        dPsi += period;
    return dPsi;
}

void addFilesToChain(TChain *chain, const string &fileList)
{
    ifstream infile(fileList.c_str());
    string line;
    if (!infile.is_open())
    {
        cout << "Error: Could not open " << fileList << endl;
        return;
    }

    cout << "--- Reading files from: " << fileList << " ---" << endl;
    int fileCount = 0;
    while (getline(infile, line))
    {
        if (!line.empty())
        {
            chain->Add(line.c_str());
            fileCount++;
        }
    }
    cout << "Successfully chained " << fileCount << " files." << endl
         << endl;
    infile.close();
}

void compareTTreeEvents_SP_PlusMinus(
    string listOld = "inputFiles_charge_MB0_woRun3Era_Aug7.txt",
    string listNew = "inputFiles_charge_wEra_MB01_Aug6_small.txt",
    string treeName_old = "Ana/ntInfo",
    string treeName_new = "Ana/ntInfo",
    string outFile = "Comparison_Output_Aug7.root")
{
    TChain *ch_old = new TChain(treeName_old.c_str());
    TChain *ch_new = new TChain(treeName_new.c_str());

    addFilesToChain(ch_old, listOld);
    addFilesToChain(ch_new, listNew);

    long nEntriesOld = ch_old->GetEntries();
    long nEntriesNew = ch_new->GetEntries();

    cout << "=== TTree Event-by-Event Comparison (All Harmonics) ===" << endl;
    cout << "Entries Old (w/o Era): " << nEntriesOld << endl;
    cout << "Entries New (w/ Era) : " << nEntriesNew << endl;

    if (nEntriesOld == 0 || nEntriesNew == 0)
    {
        cout << "Error: One or both TTrees are empty!" << endl;
        return;
    }
    if (nEntriesOld != nEntriesNew)
    {
        cout << "WARNING: Entry counts do not match! The fast loop requires identical row counts." << endl;
    }

    // --- Branch Setup ---
    int centrality_old, centrality_new;
    int run_old, run_new;
    int lumi_old, lumi_new;
    int event_old, event_new;


    

    float ephfpAngle_old[3], ephfmAngle_old[3], eptkAngle_old[2];
    float ephfpQ_old[3], ephfmQ_old[3], eptkQ_old[2];

    float ephfpAngle_new[3], ephfmAngle_new[3], eptkAngle_new[2];
    float ephfpQ_new[3], ephfmQ_new[3], eptkQ_new[2];

    ch_old->SetBranchAddress("RunNo", &run_old);
    ch_old->SetBranchAddress("LumiNo", &lumi_old);
    ch_old->SetBranchAddress("EvtNo", &event_old);
    ch_old->SetBranchAddress("centrality", &centrality_old);
    ch_old->SetBranchAddress("ephfpAngle", ephfpAngle_old);
    ch_old->SetBranchAddress("ephfmAngle", ephfmAngle_old);
    ch_old->SetBranchAddress("eptkmidAngle", eptkAngle_old);
    ch_old->SetBranchAddress("ephfpQ", ephfpQ_old);
    ch_old->SetBranchAddress("ephfmQ", ephfmQ_old);
    ch_old->SetBranchAddress("eptkmidQ", eptkQ_old);

    ch_new->SetBranchAddress("RunNo", &run_new);
    ch_new->SetBranchAddress("LumiNo", &lumi_new);
    ch_new->SetBranchAddress("EvtNo", &event_new);
    ch_new->SetBranchAddress("centrality", &centrality_new);
    ch_new->SetBranchAddress("ephfpAngle", ephfpAngle_new);
    ch_new->SetBranchAddress("ephfmAngle", ephfmAngle_new);
    ch_new->SetBranchAddress("eptkmidAngle", eptkAngle_new);
    ch_new->SetBranchAddress("ephfpQ", ephfpQ_new);
    ch_new->SetBranchAddress("ephfmQ", ephfmQ_new);
    ch_new->SetBranchAddress("eptkmidQ", eptkQ_new);

    TFile *fOut = new TFile(outFile.c_str(), "RECREATE");
    fOut->cd();

    int idx_tk[2] = {0, 1}; // Tracker indices for n=2, n=3
    int idx_hf[2] = {1, 2}; // HF indices for n=2, n=3

    // Angle Histograms
    TH2D *h2_corr_tk[2], *h2_corr_hfp[2], *h2_corr_hfm[2];
    TH1D *h1_res_tk[2], *h1_res_hfp[2], *h1_res_hfm[2];

    // Q-Vector 2D Correlation Histograms (Qx, Qy, |Q|)
    TH2D *h2_Qx_tk[2], *h2_Qy_tk[2], *h2_Qmag_tk[2];
    TH2D *h2_Qx_hfp[2], *h2_Qy_hfp[2], *h2_Qmag_hfp[2];
    TH2D *h2_Qx_hfm[2], *h2_Qy_hfm[2], *h2_Qmag_hfm[2];

    // Q-Vector 1D Distributions (Old)
    TH1D *h1_Qx_tk_old[2], *h1_Qy_tk_old[2], *h1_Qmag_tk_old[2];
    TH1D *h1_Qx_hfp_old[2], *h1_Qy_hfp_old[2], *h1_Qmag_hfp_old[2];
    TH1D *h1_Qx_hfm_old[2], *h1_Qy_hfm_old[2], *h1_Qmag_hfm_old[2];

    // Q-Vector 1D Distributions (New)
    TH1D *h1_Qx_tk_new[2], *h1_Qy_tk_new[2], *h1_Qmag_tk_new[2];
    TH1D *h1_Qx_hfp_new[2], *h1_Qy_hfp_new[2], *h1_Qmag_hfp_new[2];
    TH1D *h1_Qx_hfm_new[2], *h1_Qy_hfm_new[2], *h1_Qmag_hfm_new[2];

    // Profiles & Resolutions
    TProfile *p_HFmHFp_old[2], *p_HFmTrk_old[2], *p_HFpTrk_old[2];
    TProfile *p_HFmHFp_new[2], *p_HFmTrk_new[2], *p_HFpTrk_new[2];

    TH1D *h_sp_plus_old[2], *h_sp_plus_new[2];
    TH1D *h_sp_minus_old[2], *h_sp_minus_new[2];

    TH1D *h_vn_plus_old[2], *h_vn_plus_new[2];
    TH1D *h_vn_minus_old[2], *h_vn_minus_new[2];

    const int N_CENTBIN = 50;
    const float CENT_MAX = 50.0;

    // Output limits for Q-vector correlation axes
    const float lim_tk = 150.0;
    const float lim_hf_old = 50.0;
    const float lim_hf_new = 500.0;

    // Initialize all histograms in a loop for n=2 and n=3
    for (int i = 0; i < 2; i++)
    {
        int n = i + 2;
        float lim = TMath::Pi() / n;

        // === 1. Angles ===
        h2_corr_tk[i] = new TH2D(Form("h2_corr_tk_n%d", n), Form("Tracker n=%d;Old Angle;New Angle", n), 100, -lim, lim, 100, -lim, lim);
        h1_res_tk[i] = new TH1D(Form("h1_res_tk_n%d", n), Form("Tracker n=%d Residuals;#Delta#Psi_{%d};Events", n, n), 100, -0.5, 0.5);

        h2_corr_hfp[i] = new TH2D(Form("h2_corr_hfp_n%d", n), Form("HF+ n=%d;Old Angle;New Angle", n), 100, -lim, lim, 100, -lim, lim);
        h1_res_hfp[i] = new TH1D(Form("h1_res_hfp_n%d", n), Form("HF+ n=%d Residuals;#Delta#Psi_{%d};Events", n, n), 100, -0.5, 0.5);

        h2_corr_hfm[i] = new TH2D(Form("h2_corr_hfm_n%d", n), Form("HF- n=%d;Old Angle;New Angle", n), 100, -lim, lim, 100, -lim, lim);
        h1_res_hfm[i] = new TH1D(Form("h1_res_hfm_n%d", n), Form("HF- n=%d Residuals;#Delta#Psi_{%d};Events", n, n), 100, -0.5, 0.5);

        // === 2. Q-Vector 2D Correlations ===
        h2_Qx_tk[i] = new TH2D(Form("h2_Qx_tk_n%d", n), Form("Tracker Q_{x,n=%d};Q_{x} (Old);Q_{x} (New)", n), 100, -lim_tk, lim_tk, 100, -lim_tk, lim_tk);
        h2_Qy_tk[i] = new TH2D(Form("h2_Qy_tk_n%d", n), Form("Tracker Q_{y,n=%d};Q_{y} (Old);Q_{y} (New)", n), 100, -lim_tk, lim_tk, 100, -lim_tk, lim_tk);
        h2_Qmag_tk[i] = new TH2D(Form("h2_Qmag_tk_n%d", n), Form("Tracker |Q_{n=%d}|;|Q| (Old);|Q| (New)", n), 100, 0, lim_tk, 100, 0, lim_tk);

        h2_Qx_hfp[i] = new TH2D(Form("h2_Qx_hfp_n%d", n), Form("HF+ Q_{x,n=%d};Q_{x} (Old);Q_{x} (New)", n), 100, -lim_hf_old, lim_hf_old, 100, -lim_hf_new, lim_hf_new);
        h2_Qy_hfp[i] = new TH2D(Form("h2_Qy_hfp_n%d", n), Form("HF+ Q_{y,n=%d};Q_{y} (Old);Q_{y} (New)", n), 100, -lim_hf_old, lim_hf_old, 100, -lim_hf_new, lim_hf_new);
        h2_Qmag_hfp[i] = new TH2D(Form("h2_Qmag_hfp_n%d", n), Form("HF+ |Q_{n=%d}|;|Q| (Old);|Q| (New)", n), 100, 0, lim_hf_old, 100, 0, lim_hf_new);

        h2_Qx_hfm[i] = new TH2D(Form("h2_Qx_hfm_n%d", n), Form("HF- Q_{x,n=%d};Q_{x} (Old);Q_{x} (New)", n), 100, -lim_hf_old, lim_hf_old, 100, -lim_hf_new, lim_hf_new);
        h2_Qy_hfm[i] = new TH2D(Form("h2_Qy_hfm_n%d", n), Form("HF- Q_{y,n=%d};Q_{y} (Old);Q_{y} (New)", n), 100, -lim_hf_old, lim_hf_old, 100, -lim_hf_new, lim_hf_new);
        h2_Qmag_hfm[i] = new TH2D(Form("h2_Qmag_hfm_n%d", n), Form("HF- |Q_{n=%d}|;|Q| (Old);|Q| (New)", n), 100, 0, lim_hf_old, 100, 0, lim_hf_new);

        // === 3. Q-Vector 1D Distributions (Old) ===
        h1_Qx_tk_old[i] = new TH1D(Form("h1_Qx_tk_old_n%d", n), Form("Tracker Q_{x,n=%d} (Old);Q_{x};Events", n), 100, -lim_tk, lim_tk);
        h1_Qy_tk_old[i] = new TH1D(Form("h1_Qy_tk_old_n%d", n), Form("Tracker Q_{y,n=%d} (Old);Q_{y};Events", n), 100, -lim_tk, lim_tk);
        h1_Qmag_tk_old[i] = new TH1D(Form("h1_Qmag_tk_old_n%d", n), Form("Tracker |Q_{n=%d}| (Old);|Q|;Events", n), 100, 0, lim_tk);

        h1_Qx_hfp_old[i] = new TH1D(Form("h1_Qx_hfp_old_n%d", n), Form("HF+ Q_{x,n=%d} (Old);Q_{x};Events", n), 100, -lim_hf_old, lim_hf_old);
        h1_Qy_hfp_old[i] = new TH1D(Form("h1_Qy_hfp_old_n%d", n), Form("HF+ Q_{y,n=%d} (Old);Q_{y};Events", n), 100, -lim_hf_old, lim_hf_old);
        h1_Qmag_hfp_old[i] = new TH1D(Form("h1_Qmag_hfp_old_n%d", n), Form("HF+ |Q_{n=%d}| (Old);|Q|;Events", n), 100, 0, lim_hf_old);

        h1_Qx_hfm_old[i] = new TH1D(Form("h1_Qx_hfm_old_n%d", n), Form("HF- Q_{x,n=%d} (Old);Q_{x};Events", n), 100, -lim_hf_old, lim_hf_old);
        h1_Qy_hfm_old[i] = new TH1D(Form("h1_Qy_hfm_old_n%d", n), Form("HF- Q_{y,n=%d} (Old);Q_{y};Events", n), 100, -lim_hf_old, lim_hf_old);
        h1_Qmag_hfm_old[i] = new TH1D(Form("h1_Qmag_hfm_old_n%d", n), Form("HF- |Q_{n=%d}| (Old);|Q|;Events", n), 100, 0, lim_hf_old);

        // === 4. Q-Vector 1D Distributions (New) ===
        h1_Qx_tk_new[i] = new TH1D(Form("h1_Qx_tk_new_n%d", n), Form("Tracker Q_{x,n=%d} (New);Q_{x};Events", n), 100, -lim_tk, lim_tk);
        h1_Qy_tk_new[i] = new TH1D(Form("h1_Qy_tk_new_n%d", n), Form("Tracker Q_{y,n=%d} (New);Q_{y};Events", n), 100, -lim_tk, lim_tk);
        h1_Qmag_tk_new[i] = new TH1D(Form("h1_Qmag_tk_new_n%d", n), Form("Tracker |Q_{n=%d}| (New);|Q|;Events", n), 100, 0, lim_tk);

        h1_Qx_hfp_new[i] = new TH1D(Form("h1_Qx_hfp_new_n%d", n), Form("HF+ Q_{x,n=%d} (New);Q_{x};Events", n), 100, -lim_hf_new, lim_hf_new);
        h1_Qy_hfp_new[i] = new TH1D(Form("h1_Qy_hfp_new_n%d", n), Form("HF+ Q_{y,n=%d} (New);Q_{y};Events", n), 100, -lim_hf_new, lim_hf_new);
        h1_Qmag_hfp_new[i] = new TH1D(Form("h1_Qmag_hfp_new_n%d", n), Form("HF+ |Q_{n=%d}| (New);|Q|;Events", n), 100, 0, lim_hf_new);

        h1_Qx_hfm_new[i] = new TH1D(Form("h1_Qx_hfm_new_n%d", n), Form("HF- Q_{x,n=%d} (New);Q_{x};Events", n), 100, -lim_hf_new, lim_hf_new);
        h1_Qy_hfm_new[i] = new TH1D(Form("h1_Qy_hfm_new_n%d", n), Form("HF- Q_{y,n=%d} (New);Q_{y};Events", n), 100, -lim_hf_new, lim_hf_new);
        h1_Qmag_hfm_new[i] = new TH1D(Form("h1_Qmag_hfm_new_n%d", n), Form("HF- |Q_{n=%d}| (New);|Q|;Events", n), 100, 0, lim_hf_new);

        // === 5. Profiles & SP Resolutions ===
        p_HFmHFp_old[i] = new TProfile(Form("p_HFmHFp_old_n%d", n), Form("<Q_{HF-}^{%d} Q_{HF+}^{*%d}> (Old);Centrality %%;", n, n), N_CENTBIN, 0, CENT_MAX);
        p_HFmTrk_old[i] = new TProfile(Form("p_HFmTrk_old_n%d", n), Form("<Q_{HF-}^{%d} Q_{Trk}^{*%d}> (Old);Centrality %%;", n, n), N_CENTBIN, 0, CENT_MAX);
        p_HFpTrk_old[i] = new TProfile(Form("p_HFpTrk_old_n%d", n), Form("<Q_{HF+}^{%d} Q_{Trk}^{*%d}> (Old);Centrality %%;", n, n), N_CENTBIN, 0, CENT_MAX);

        p_HFmHFp_new[i] = new TProfile(Form("p_HFmHFp_new_n%d", n), Form("<Q_{HF-}^{%d} Q_{HF+}^{*%d}> (New);Centrality %%;", n, n), N_CENTBIN, 0, CENT_MAX);
        p_HFmTrk_new[i] = new TProfile(Form("p_HFmTrk_new_n%d", n), Form("<Q_{HF-}^{%d} Q_{Trk}^{*%d}> (New);Centrality %%;", n, n), N_CENTBIN, 0, CENT_MAX);
        p_HFpTrk_new[i] = new TProfile(Form("p_HFpTrk_new_n%d", n), Form("<Q_{HF+}^{%d} Q_{Trk}^{*%d}> (New);Centrality %%;", n, n), N_CENTBIN, 0, CENT_MAX);

        h_sp_plus_old[i] = new TH1D(Form("h_sp_plus_old_n%d", n), Form("n=%d SP Resolution (Plus, Old);Centrality %%;R_{%d}^{plus}", n, n), N_CENTBIN, 0, CENT_MAX);
        h_sp_plus_new[i] = new TH1D(Form("h_sp_plus_new_n%d", n), Form("n=%d SP Resolution (Plus, New);Centrality %%;R_{%d}^{plus}", n, n), N_CENTBIN, 0, CENT_MAX);
        h_sp_minus_old[i] = new TH1D(Form("h_sp_minus_old_n%d", n), Form("n=%d SP Resolution (Minus, Old);Centrality %%;R_{%d}^{minus}", n, n), N_CENTBIN, 0, CENT_MAX);
        h_sp_minus_new[i] = new TH1D(Form("h_sp_minus_new_n%d", n), Form("n=%d SP Resolution (Minus, New);Centrality %%;R_{%d}^{minus}", n, n), N_CENTBIN, 0, CENT_MAX);

        h_vn_plus_old[i] = new TH1D(Form("h_vn_plus_old_n%d", n), Form("v_{%d}{SP, HF+} (Old);Centrality %%;v_{%d}", n, n), N_CENTBIN, 0, CENT_MAX);
        h_vn_plus_new[i] = new TH1D(Form("h_vn_plus_new_n%d", n), Form("v_{%d}{SP, HF+} (New);Centrality %%;v_{%d}", n, n), N_CENTBIN, 0, CENT_MAX);

        h_vn_minus_old[i] = new TH1D(Form("h_vn_minus_old_n%d", n), Form("v_{%d}{SP, HF-} (Old);Centrality %%;v_{%d}", n, n), N_CENTBIN, 0, CENT_MAX);
        h_vn_minus_new[i] = new TH1D(Form("h_vn_minus_new_n%d", n), Form("v_{%d}{SP, HF-} (New);Centrality %%;v_{%d}", n, n), N_CENTBIN, 0, CENT_MAX);
    }

    // --- Build Index for Synchronization ---
    cout << "Building robust memory index on new tree to bypass ROOT TChain limits..." << endl;
    std::map<std::pair<int, int>, long> syncIndex;
    
    for (long j = 0; j < nEntriesNew; j++) {
        ch_new->GetEntry(j);
        syncIndex[std::make_pair(run_new, event_new)] = j;

        if (j < 5) cout << "NEW TREE - Run: " << run_new << ", Evt: " << event_new << endl;
    }
    cout << "Successfully indexed " << syncIndex.size() << " unique events from the new tree." << endl;

    // --- Synchronized Event Loop ---
    cout << "Processing events..." << endl;
    long matched_events = 0;
    
    // Loop over the old tree and search for the matching event in the new tree
    for (long i = 0; i < nEntriesOld; i++)
    {
        ch_old->GetEntry(i);

        if (i < 5) cout << "OLD TREE - Run: " << run_old << ", Evt: " << event_old << endl;
        // Look up the old event in our custom index
        auto it = syncIndex.find(std::make_pair(run_old, event_old));

        
        
        // If the event doesn't exist in the new dataset, skip it
        if (it == syncIndex.end()) {
            continue;
        }

        // Fetch the exact matching event in the new tree using the stored row number
        ch_new->GetEntry(it->second);

        matched_events++;

        if (matched_events % 100000 == 0)
            cout << "Processed " << matched_events << " synchronized events..." << endl;

        int cent_old = centrality_old / 2;
        int cent_new = centrality_new / 2;

        for (int harm = 0; harm < 2; harm++)
        {
            float n = (float)(harm + 2);
            int t = idx_tk[harm];
            int h = idx_hf[harm];

            // 1. Angle Correlations & Residuals
            h2_corr_tk[harm]->Fill(eptkAngle_old[t], eptkAngle_new[t]);
            h1_res_tk[harm]->Fill(getDeltaPsi(eptkAngle_new[t], eptkAngle_old[t], n));

            h2_corr_hfp[harm]->Fill(ephfpAngle_old[h], ephfpAngle_new[h]);
            h1_res_hfp[harm]->Fill(getDeltaPsi(ephfpAngle_new[h], ephfpAngle_old[h], n));

            h2_corr_hfm[harm]->Fill(ephfmAngle_old[h], ephfmAngle_new[h]);
            h1_res_hfm[harm]->Fill(getDeltaPsi(ephfmAngle_new[h], ephfmAngle_old[h], n));

            // 2. Explicit Q-Vector Component Calculations (Old)
            double old_Qx_tk = eptkQ_old[t] * TMath::Cos(n * eptkAngle_old[t]);
            double old_Qy_tk = eptkQ_old[t] * TMath::Sin(n * eptkAngle_old[t]);
            double old_Qmag_tk = eptkQ_old[t];

            double old_Qx_hfp = ephfpQ_old[h] * TMath::Cos(n * ephfpAngle_old[h]);
            double old_Qy_hfp = ephfpQ_old[h] * TMath::Sin(n * ephfpAngle_old[h]);
            double old_Qmag_hfp = ephfpQ_old[h];

            double old_Qx_hfm = ephfmQ_old[h] * TMath::Cos(n * ephfmAngle_old[h]);
            double old_Qy_hfm = ephfmQ_old[h] * TMath::Sin(n * ephfmAngle_old[h]);
            double old_Qmag_hfm = ephfmQ_old[h];

            // 3. Explicit Q-Vector Component Calculations (New)
            double new_Qx_tk = eptkQ_new[t] * TMath::Cos(n * eptkAngle_new[t]);
            double new_Qy_tk = eptkQ_new[t] * TMath::Sin(n * eptkAngle_new[t]);
            double new_Qmag_tk = eptkQ_new[t];

            double new_Qx_hfp = ephfpQ_new[h] * TMath::Cos(n * ephfpAngle_new[h]);
            double new_Qy_hfp = ephfpQ_new[h] * TMath::Sin(n * ephfpAngle_new[h]);
            double new_Qmag_hfp = ephfpQ_new[h];

            double new_Qx_hfm = ephfmQ_new[h] * TMath::Cos(n * ephfmAngle_new[h]);
            double new_Qy_hfm = ephfmQ_new[h] * TMath::Sin(n * ephfmAngle_new[h]);
            double new_Qmag_hfm = ephfmQ_new[h];

            // 4. Fill Qx, Qy, and |Q| 2D Correlations
            h2_Qx_tk[harm]->Fill(old_Qx_tk, new_Qx_tk);
            h2_Qy_tk[harm]->Fill(old_Qy_tk, new_Qy_tk);
            h2_Qmag_tk[harm]->Fill(old_Qmag_tk, new_Qmag_tk);

            h2_Qx_hfp[harm]->Fill(old_Qx_hfp, new_Qx_hfp);
            h2_Qy_hfp[harm]->Fill(old_Qy_hfp, new_Qy_hfp);
            h2_Qmag_hfp[harm]->Fill(old_Qmag_hfp, new_Qmag_hfp);

            h2_Qx_hfm[harm]->Fill(old_Qx_hfm, new_Qx_hfm);
            h2_Qy_hfm[harm]->Fill(old_Qy_hfm, new_Qy_hfm);
            h2_Qmag_hfm[harm]->Fill(old_Qmag_hfm, new_Qmag_hfm);

            // 5. Fill Qx, Qy, and |Q| 1D Distributions
            h1_Qx_tk_old[harm]->Fill(old_Qx_tk);
            h1_Qy_tk_old[harm]->Fill(old_Qy_tk);
            h1_Qmag_tk_old[harm]->Fill(old_Qmag_tk);

            h1_Qx_tk_new[harm]->Fill(new_Qx_tk);
            h1_Qy_tk_new[harm]->Fill(new_Qy_tk);
            h1_Qmag_tk_new[harm]->Fill(new_Qmag_tk);

            h1_Qx_hfp_old[harm]->Fill(old_Qx_hfp);
            h1_Qy_hfp_old[harm]->Fill(old_Qy_hfp);
            h1_Qmag_hfp_old[harm]->Fill(old_Qmag_hfp);

            h1_Qx_hfp_new[harm]->Fill(new_Qx_hfp);
            h1_Qy_hfp_new[harm]->Fill(new_Qy_hfp);
            h1_Qmag_hfp_new[harm]->Fill(new_Qmag_hfp);

            h1_Qx_hfm_old[harm]->Fill(old_Qx_hfm);
            h1_Qy_hfm_old[harm]->Fill(old_Qy_hfm);
            h1_Qmag_hfm_old[harm]->Fill(old_Qmag_hfm);

            h1_Qx_hfm_new[harm]->Fill(new_Qx_hfm);
            h1_Qy_hfm_new[harm]->Fill(new_Qy_hfm);
            h1_Qmag_hfm_new[harm]->Fill(new_Qmag_hfm);

            // 6. Fill Event-Averaged Profiles
            if (cent_old >= 0 && cent_old < 50)
            {
                TComplex old_Q_HFm(old_Qx_hfm, old_Qy_hfm, 0);
                TComplex old_Q_HFp(old_Qx_hfp, old_Qy_hfp, 0);
                TComplex old_Q_Trk(old_Qx_tk, old_Qy_tk, 0);

                p_HFmHFp_old[harm]->Fill(cent_old, (old_Q_HFm * TComplex::Conjugate(old_Q_HFp)).Re());
                p_HFmTrk_old[harm]->Fill(cent_old, (old_Q_HFm * TComplex::Conjugate(old_Q_Trk)).Re());
                p_HFpTrk_old[harm]->Fill(cent_old, (old_Q_HFp * TComplex::Conjugate(old_Q_Trk)).Re());
            }

            if (cent_new >= 0 && cent_new < 50)
            {
                TComplex new_Q_HFm(new_Qx_hfm, new_Qy_hfm, 0);
                TComplex new_Q_HFp(new_Qx_hfp, new_Qy_hfp, 0);
                TComplex new_Q_Trk(new_Qx_tk, new_Qy_tk, 0);

                p_HFmHFp_new[harm]->Fill(cent_new, (new_Q_HFm * TComplex::Conjugate(new_Q_HFp)).Re());
                p_HFmTrk_new[harm]->Fill(cent_new, (new_Q_HFm * TComplex::Conjugate(new_Q_Trk)).Re());
                p_HFpTrk_new[harm]->Fill(cent_new, (new_Q_HFp * TComplex::Conjugate(new_Q_Trk)).Re());
            }
        }
    }

    // --- Compute Final SP Resolutions ---
    for (int harm = 0; harm < 2; harm++)
    {
        for (int b = 1; b <= N_CENTBIN; b++)
        {

            // Old Math
            double hmm_old = p_HFmHFp_old[harm]->GetBinContent(b);
            double hmt_old = p_HFmTrk_old[harm]->GetBinContent(b);
            double hpt_old = p_HFpTrk_old[harm]->GetBinContent(b);

            if (hpt_old != 0 && hmt_old != 0)
            {
                double plus_term = (hmm_old * hmt_old) / hpt_old;
                double minus_term = (hmm_old * hpt_old) / hmt_old;

                double res_plus_old = (plus_term > 0) ? TMath::Sqrt(plus_term) : 0;
                double res_minus_old = (minus_term > 0) ? TMath::Sqrt(minus_term) : 0;

                if (res_plus_old > 0)
                {
                    h_sp_plus_old[harm]->SetBinContent(b, res_plus_old);
                    h_vn_plus_old[harm]->SetBinContent(b, hpt_old / res_plus_old); // v_n calculation
                }
                if (res_minus_old > 0)
                {
                    h_sp_minus_old[harm]->SetBinContent(b, res_minus_old);
                    h_vn_minus_old[harm]->SetBinContent(b, hmt_old / res_minus_old); // v_n calculation
                }
            }

            // New Math
            double hmm_new = p_HFmHFp_new[harm]->GetBinContent(b);
            double hmt_new = p_HFmTrk_new[harm]->GetBinContent(b);
            double hpt_new = p_HFpTrk_new[harm]->GetBinContent(b);

            if (hpt_new != 0 && hmt_new != 0)
            {
                double plus_term_new = (hmm_new * hmt_new) / hpt_new;
                double minus_term_new = (hmm_new * hpt_new) / hmt_new;

                double res_plus_new = (plus_term_new > 0) ? TMath::Sqrt(plus_term_new) : 0;
                double res_minus_new = (minus_term_new > 0) ? TMath::Sqrt(minus_term_new) : 0;

                if (res_plus_new > 0)
                {
                    h_sp_plus_new[harm]->SetBinContent(b, res_plus_new);
                    h_vn_plus_new[harm]->SetBinContent(b, hpt_new / res_plus_new); // v_n calculation
                }
                if (res_minus_new > 0)
                {
                    h_sp_minus_new[harm]->SetBinContent(b, res_minus_new);
                    h_vn_minus_new[harm]->SetBinContent(b, hmt_new / res_minus_new); // v_n calculation
                }
            }
        }
    }

    // --- Save and Close ---
    fOut->Write();
    fOut->Close();

    cout << "Finished! All n=2 and n=3 distributions saved to: " << outFile << endl;
}
