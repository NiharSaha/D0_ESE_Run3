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
#include "THnSparse.h"

// Include your auto-generated quantile cuts header
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/ESE_Cuts_MB11to21_Jan22.h"

using namespace std;

void Calculate_Resolution(TString input_txt, TString output_path, int istart, int iend) {
    
    // --- 1. CONFIGURATION ---
    const int N_CENTBIN = 4;
    Int_t min_centbin[N_CENTBIN] = {0, 10, 30, 50};
    Int_t max_centbin[N_CENTBIN] = {10, 30, 50, 90};
    const int N_Q_BINS = 10;

    // Output file setup
    TString outfile = TString::Format("%s/ROOT/Resolution_%d_%d.root", output_path.Data(), istart, iend);
    
    // TH2D to store resolution ingredients
    // X-axis: 0.5 (HF+ * HF-), 1.5 (HF+ * Trk), 2.5 (HF- * Trk)
    // Y-axis: Index = (CentralityGroup * 10) + QuantileIndex (0 to 39)
    TH2D *hRes_v2 = new TH2D("hRes_v2", "v2 Resolution Ingredients;Type;Cent_Quant_Idx", 3, 0, 3, 40, 0, 40);
    TH2D *hRes_v3 = new TH2D("hRes_v3", "v3 Resolution Ingredients;Type;Cent_Quant_Idx", 3, 0, 3, 40, 0, 40);
    hRes_v2->Sumw2();
    hRes_v3->Sumw2();

    // --- 2. FILE LOOP ---
    ifstream file_stream(input_txt.Data());
    string filename;
    int ifile = 0;

    while (file_stream >> filename) {
        if (ifile < istart) { ifile++; continue; }
        if (ifile >= iend) break;

        TFile *fin = TFile::Open(filename.c_str());
        if(!fin || fin->IsZombie()) {
            std::cout << "Warning: Skipping bad file: " << filename << std::endl;
            if (fin) { fin->Close(); delete fin; }
            ifile++; continue;
        }

        std::cout << ">>> Resolution Pass: Processing file " << ifile << " : " << filename << std::endl;

        TTree *tree = (TTree*)fin->Get("d0Analyzer/VCNtuple_D02kpi");
        TTree *t_eventinfoana = (TTree*)fin->Get("eventinfoana/EventInfoNtuple");
        tree->AddFriend(t_eventinfoana);


        Int_t centrality;
        Float_t ephfmQ[3], ephfpQ[3], ephfmSumW[3], ephfpSumW[3], ephfpAngle[3], ephfmAngle[3], eptkAngle[2], eptkQ[2];

        tree->SetBranchAddress("centrality", &centrality);
        tree->SetBranchAddress("ephfpAngle", ephfpAngle);
        tree->SetBranchAddress("ephfmAngle", ephfmAngle);
        tree->SetBranchAddress("ephfmQ", ephfmQ);
        tree->SetBranchAddress("ephfpQ", ephfpQ);
        tree->SetBranchAddress("ephfmSumW", ephfmSumW);
        tree->SetBranchAddress("ephfpSumW", ephfpSumW);
        tree->SetBranchAddress("eptkAngle", eptkAngle);
        tree->SetBranchAddress("eptkQ", eptkQ);

        tree->SetBranchStatus("*", 0);
        for (const auto& p : {"centrality", "ephfpAngle", "ephfmAngle", "ephfpQ", "ephfmQ", "eptkAngle", "eptkQ", "ephfmSumW", "ephfpSumW"})
            tree->SetBranchStatus(p, 1);

        // --- 4. EVENT LOOP ---
        Int_t n_entries = tree->GetEntries();
        for (Long64_t ii = 0; ii < n_entries; ii++) {
            tree->GetEntry(ii);

            int cent = centrality / 2;
            if (cent < 0 || cent >= 90) continue;

            // Determine Broad Centrality Group (0-3)
            int i_idx = -1;
            for (int i = 0; i < N_CENTBIN; i++) {
                if (cent >= min_centbin[i] && cent < max_centbin[i]) { i_idx = i; break; }
            }
            if (i_idx == -1) continue;

            // ESE Quantile Determination
            Double_t Q2_w = ephfmSumW[1] + ephfpSumW[1];
            Double_t Q3_w = ephfmSumW[2] + ephfpSumW[2];
            if (Q2_w <= 0 || Q3_w <= 0) continue;

            Double_t q2_val = (ephfmQ[1] + ephfpQ[1]) / Q2_w;
            Double_t q3_val = (ephfmQ[2] + ephfpQ[2]) / Q3_w;

            Int_t k2_idx = TMath::BinarySearch(11, q2_cuts[cent], q2_val);
            Int_t k3_idx = TMath::BinarySearch(11, q3_cuts[cent], q3_val);
            if (k2_idx < 0) k2_idx = 0;
	    if (k2_idx >= 10) k2_idx = 9;

	    if (k3_idx < 0) k3_idx = 0;
	    if (k3_idx >= 10) k3_idx = 9;

            Int_t k_indices[2] = {k2_idx, k3_idx};

            // Calculate Correlations for v2 (n=0) and v3 (n=1)
            for (int n = 0; n < 2; ++n) {
                int harm = n + 2;
                Int_t current_y_bin = i_idx * N_Q_BINS + k_indices[n];
                TH2D *hRes_current = (n == 0) ? hRes_v2 : hRes_v3;

                // Sub-event Q-vectors: A (HF+), B (HF-), C (Tracker)
                TComplex QA(ephfpQ[harm-1]*cos(harm*ephfpAngle[harm-1]), ephfpQ[harm-1]*sin(harm*ephfpAngle[harm-1]));
                TComplex QB(ephfmQ[harm-1]*cos(harm*ephfmAngle[harm-1]), ephfmQ[harm-1]*sin(harm*ephfmAngle[harm-1]));
                TComplex QC(eptkQ[n]*cos(harm*eptkAngle[n]), eptkQ[n]*sin(harm*eptkAngle[n]));

                Double_t WA = ephfpSumW[harm-1];
                Double_t WB = ephfmSumW[harm-1];

                if (WA > 0 && WB > 0)
                    hRes_current->Fill(0.5, (Double_t)current_y_bin, (QA * TComplex::Conjugate(QB)).Re() / (WA * WB));
                if (WA > 0)
                    hRes_current->Fill(1.5, (Double_t)current_y_bin, (QA * TComplex::Conjugate(QC)).Re() / WA);
                if (WB > 0)
                    hRes_current->Fill(2.5, (Double_t)current_y_bin, (QB * TComplex::Conjugate(QC)).Re() / WB);
            }
        }
        fin->Close(); delete fin;
        ifile++;
    }

    // --- 5. SAVE ---
    TFile *fout = new TFile(outfile, "RECREATE");
    hRes_v2->Write();
    hRes_v3->Write();
    fout->Close();
    std::cout << ">>> Resolution Ingredients saved to: " << outfile << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc == 5) // Expecting 4 arguments: input_txt, output_path, istart, iend

    {
        TString input_txt = argv[1];
        TString output_path = argv[2];
        int istart = std::stoi(argv[3]);
        int iend = std::stoi(argv[4]);

        Calculate_Resolution(input_txt, output_path, istart, iend);
    }
    else
    {
        std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
        return 1;
    }
    return 0;
}
