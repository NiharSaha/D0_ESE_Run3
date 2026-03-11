#include <iostream>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"

void ExtractQuantiles_quantileFunc() {
    // --- 1. CONFIGURATION ---
    // Update this path to your exact merged file location on scratch
  TString inputPath = "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Q2quantile_Jan20_MB11to21_v2/ROOT/Quantiles_ESE_out_combined.root";
  TString outputHeader = "ESE_Cuts_Jan22.h";

    const int N_CENTBIN = 4;        // Your broad centrality groups
    const int N_CENTBIN_1PC = 90;   // 1% bins
    const int N_Q_BINS = 10;        // Deciles (10% steps)

    // Centrality binning definitions (matching your previous macro)
    int min_centbin[] = {0, 10, 30, 50};
    int max_centbin[] = {10, 30, 50, 90};

    // Define probabilities for GetQuantiles (0%, 10%, ..., 100%)
    double probs[N_Q_BINS + 1];
    for (int i = 0; i <= N_Q_BINS; i++) probs[i] = i * 0.1;

    // --- 2. OPEN INPUT ---
    TFile *fIn = TFile::Open(inputPath, "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "Error: Could not open input file: " << inputPath << std::endl;
        return;
    }

    // Storage for results
    double q2_cuts_table[90][11];
    double q3_cuts_table[90][11];

    Long64_t total_q2_events = 0;
    Long64_t total_q3_events = 0;
    
    // --- 3. EXTRACTION LOOP ---
    std::cout << ">>> Extracting quantiles from 1% histograms..." << std::endl;

    for (int i = 0; i < N_CENTBIN; i++) {
        for (int j = min_centbin[i]; j < max_centbin[i]; j++) {
            // Adjust the paths below if your histograms are inside subdirectories (e.g., "q2_raw_1pc/")
            TH1D *hQ2 = (TH1D*)fIn->Get(Form("q2_raw_1pc/hist_q2_tot_cent%d_c%d", i, j));
            TH1D *hQ3 = (TH1D*)fIn->Get(Form("q3_raw_1pc/hist_q3_tot_cent%d_c%d", i, j));

            if (hQ2 && hQ2->GetEntries() > 0) {
                hQ2->GetQuantiles(N_Q_BINS + 1, q2_cuts_table[j], probs);
		total_q2_events += hQ2->GetEntries(); // Track stats
	    } else {
                for(int k=0; k<=N_Q_BINS; k++) q2_cuts_table[j][k] = 0.0;
            }

            if (hQ3 && hQ3->GetEntries() > 0) {
                hQ3->GetQuantiles(N_Q_BINS + 1, q3_cuts_table[j], probs);
		total_q3_events += hQ3->GetEntries(); // Track stats
            } else {
                for(int k=0; k<=N_Q_BINS; k++) q3_cuts_table[j][k] = 0.0;
            }
        }
    }

    // --- 4. GENERATE HEADER FILE ---
    std::ofstream hFile(outputHeader);
    hFile << "#ifndef ESE_CUTS_H" << std::endl;
    hFile << "#define ESE_CUTS_H" << std::endl << std::endl;
    hFile << "/* Auto-generated ESE Quantile Table */" << std::endl;
    hFile << "/* Produced from " << total_q2_events << " total events */" << std::endl; // Statistics printout in header
    hFile << "/* Rows: 90 Centrality Bins (1% each) */" << std::endl;
    hFile << "/* Columns: 11 Quantile Boundaries (0%, 10%, ..., 100%) */" << std::endl << std::endl;

    // Write Q2 Array
    hFile << "double q2_cuts[90][11] = {" << std::endl;
    for (int j = 0; j < N_CENTBIN_1PC; j++) {
        hFile << "  {";
        for (int k = 0; k <= N_Q_BINS; k++) {
            hFile << Form("%.8f", q2_cuts_table[j][k]) << (k == N_Q_BINS ? "" : ", ");
        }
        hFile << "}" << (j == N_CENTBIN_1PC - 1 ? "" : ",") << " // Cent " << j << "%" << std::endl;
    }
    hFile << "};" << std::endl << std::endl;

    // Write Q3 Array
    hFile << "double q3_cuts[90][11] = {" << std::endl;
    for (int j = 0; j < N_CENTBIN_1PC; j++) {
        hFile << "  {";
        for (int k = 0; k <= N_Q_BINS; k++) {
            hFile << Form("%.8f", q3_cuts_table[j][k]) << (k == N_Q_BINS ? "" : ", ");
        }
        hFile << "}" << (j == N_CENTBIN_1PC - 1 ? "" : ",") << " // Cent " << j << "%" << std::endl;
    }
    hFile << "};" << std::endl << std::endl;

    hFile << "#endif" << std::endl;
    hFile.close();

    fIn->Close();

    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << ">>> SUCCESS: " << outputHeader << " generated." << std::endl;
    std::cout << ">>> Total events used for q2 quantiles: " << total_q2_events << std::endl;
    std::cout << ">>> Total events used for q3 quantiles: " << total_q3_events << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
}
