#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TSystem.h"

void ExtractQuantiles_sortFunc_new(int centBin) {
    // 1. Open input file
    TFile *fIn = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Quantiles_MB11to21_Feb2_weight_v3/ROOT/Quantiles_ESE_out_combined.root", "READ");
    if (!fIn || fIn->IsZombie()) return;

    // 2. Get the Ntuple
    TNtuple *nt = (TNtuple*)fIn->Get(Form("Q_ntuple/nt_qn_cent%i_%i", centBin, centBin+1));
    if (!nt || nt->GetEntries() == 0) {
        std::cout << "Bin " << centBin << " empty or missing." << std::endl;
        fIn->Close();
        return;
    }

    // 3. Set Branch Addresses (Including Weights)
    float q2_val, q3_val, q2_w, q3_w;
    nt->SetBranchAddress("q2_hf_total", &q2_val);
    nt->SetBranchAddress("q3_hf_total", &q3_val);
    nt->SetBranchAddress("q2_hf_total_w", &q2_w);
    nt->SetBranchAddress("q3_hf_total_w", &q3_w);

    std::vector<double> vQ2, vQ3;
    
    // Ensure parent directories exist
    if (gSystem->AccessPathName("temp_cuts"))  gSystem->mkdir("temp_cuts");
    if (gSystem->AccessPathName("temp_roots")) gSystem->mkdir("temp_roots");

    // 4. Create Output File
    TFile *fOut = new TFile(Form("temp_roots/hist_bin%d.root", centBin), "RECREATE");

    // Create Internal Directory Structure for better organization
    TDirectory *dirWeights = fOut->mkdir("Weights");
    TDirectory *dirSlices  = fOut->mkdir("Slices");

    // Define Weight Histograms
    dirWeights->cd();
    TH1D *hWeightQ2 = new TH1D(Form("hWeightQ2_c%d", centBin), "Total Weight q2;E_{T} Sum;Entries", 2000, 0, 10000);
    TH1D *hWeightQ3 = new TH1D(Form("hWeightQ3_c%d", centBin), "Total Weight q3;E_{T} Sum;Entries", 2000, 0, 10000);

    // 5. Extract Data & Fill Weights
    for (int i=0; i < nt->GetEntries(); i++) {
        nt->GetEntry(i);
        if (q2_val > -1.0) {
            vQ2.push_back((double)q2_val);
            hWeightQ2->Fill(q2_w);
        }
        if (q3_val > -1.0) {
            vQ3.push_back((double)q3_val);
            hWeightQ3->Fill(q3_w);
        }
    }
    std::sort(vQ2.begin(), vQ2.end());
    std::sort(vQ3.begin(), vQ3.end());

    // 6. Calculate Quantile Boundaries
    double cuts2[11], cuts3[11];
    int n2 = vQ2.size();
    int n3 = vQ3.size();
    for (int k = 0; k <= 10; k++) {
        int idx2 = (k == 10) ? n2 - 1 : (k * n2 / 10);
        int idx3 = (k == 10) ? n3 - 1 : (k * n3 / 10);
        cuts2[k] = (n2 > 0) ? vQ2[idx2] : 0;
        cuts3[k] = (n3 > 0) ? vQ3[idx3] : 0;
    }

    // Save cut text files
    std::ofstream outCut(Form("temp_cuts/cuts_bin%d.txt", centBin));
    for (int k=0; k<=10; k++) outCut << cuts2[k] << (k==10?"":",");
    outCut << std::endl;
    for (int k=0; k<=10; k++) outCut << cuts3[k] << (k==10?"":",");
    outCut.close();

    // 7. Define and Fill qn Histograms
    dirSlices->cd();
    TH1D *hTotalQ2 = new TH1D(Form("hRawQ2_c%d", centBin), "Total q2", 7000, 0, 0.7);
    TH1D *hTotalQ3 = new TH1D(Form("hRawQ3_c%d", centBin), "Total q3", 7000, 0, 0.7);

    TH1D *hSliceQ2[10], *hSliceQ3[10];
    for(int k=0; k<10; k++) {
        hSliceQ2[k] = new TH1D(Form("hSliceQ2_c%d_q%d", centBin, k+1), "", 7000, 0, 0.7);
        hSliceQ3[k] = new TH1D(Form("hSliceQ3_c%d_q%d", centBin, k+1), "", 7000, 0, 0.7);
    }

    for (double val : vQ2) {
        hTotalQ2->Fill(val);
        for (int k=0; k<10; k++) {
            if (val >= cuts2[k] && (k==9 ? val <= cuts2[k+1] : val < cuts2[k+1])) {
                hSliceQ2[k]->Fill(val); 
                break; 
            }
        }
    }
    for (double val : vQ3) {
        hTotalQ3->Fill(val);
        for (int k=0; k<10; k++) {
            if (val >= cuts3[k] && (k==9 ? val <= cuts3[k+1] : val < cuts3[k+1])) {
                hSliceQ3[k]->Fill(val); 
                break; 
            }
        }
    }

    // 8. Write and Close
    fOut->Write();
    fOut->Close();
    fIn->Close();
}
