#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TList.h"

void ExtractQuantiles_sortFunc(int centBin, 
    const char* inputDir = "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Charge_Quantiles_20Qbin_MB0to1_Jun15/ROOT"
) {
     const int N_QBINS = 20;
    // 1. Scan directory and build TChain
    TChain *chain = new TChain(Form("Q_ntuple/nt_qn_cent%i_%i", centBin, centBin+1));

    TSystemDirectory dir(inputDir, inputDir);
    TList *files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "Cannot open directory: " << inputDir << std::endl;
        return;
    }

    int nFilesAdded = 0;
    TIter next(files);
    TSystemFile *sFile;
    while ((sFile = (TSystemFile*)next())) {
        TString fname = sFile->GetName();
        if (!sFile->IsDirectory() && fname.EndsWith(".root") && fname.Contains("Quantiles_ESE_ntuples_")) {
            TString fullPath = TString::Format("%s/%s", inputDir, fname.Data());

            // Check file is not corrupted before adding
            TFile *ftest = TFile::Open(fullPath.Data(), "READ");
            if (!ftest || ftest->IsZombie() || ftest->TestBit(TFile::kRecovered)) {
                std::cerr << "Skipping bad/recovered file: " << fullPath << std::endl;
                if (ftest) { ftest->Close(); delete ftest; }
                continue;
            }
            // Check ntuple exists in this file
            TObject *obj = ftest->Get(Form("Q_ntuple/nt_qn_cent%i_%i", centBin, centBin+1));
            if (!obj) {
                std::cerr << "Skipping file (ntuple not found): " << fullPath << std::endl;
                ftest->Close();
                delete ftest;
                continue;
            }
            ftest->Close();

            chain->Add(fullPath.Data());
            nFilesAdded++;
            std::cout << "Added: " << fullPath << std::endl;
        }
    }
    delete files;

    std::cout << "Total files added: " << nFilesAdded << " for centBin " << centBin << std::endl;

    if (chain->GetEntries() == 0) {
        std::cout << "Bin " << centBin << " empty or missing." << std::endl;
        delete chain;
        return;
    }

    // 2. Set Branch Addresses - disable unused branches for speed
    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("q2_hf_total",   1);
    chain->SetBranchStatus("q3_hf_total",   1);
    chain->SetBranchStatus("q2_hf_total_w", 1);
    chain->SetBranchStatus("q3_hf_total_w", 1);

    float q2_val, q3_val, q2_w, q3_w;
    chain->SetBranchAddress("q2_hf_total",   &q2_val);
    chain->SetBranchAddress("q3_hf_total",   &q3_val);
    chain->SetBranchAddress("q2_hf_total_w", &q2_w);
    chain->SetBranchAddress("q3_hf_total_w", &q3_w);

    std::vector<double> vQ2, vQ3;
    vQ2.reserve(chain->GetEntries());
    vQ3.reserve(chain->GetEntries());

    // Ensure output directories exist
    if (gSystem->AccessPathName("temp_cuts_MB0to1_charge"))  gSystem->mkdir("temp_cuts_MB0to1_charge");
    if (gSystem->AccessPathName("temp_roots_MB0to1_charge")) gSystem->mkdir("temp_roots_MB0to1_charge");

    // 3. Create Output File
    TFile *fOut = new TFile(Form("temp_roots_MB0to1_charge/hist_bin%d.root", centBin), "RECREATE");
    TDirectory *dirWeights = fOut->mkdir("Weights");
    TDirectory *dirSlices  = fOut->mkdir("Slices");

    dirWeights->cd();
    TH1D *hWeightQ2 = new TH1D(Form("hWeightQ2_c%d", centBin), "Total Weight q2;E_{T} Sum;Entries", 2000, 0, 10000);
    TH1D *hWeightQ3 = new TH1D(Form("hWeightQ3_c%d", centBin), "Total Weight q3;E_{T} Sum;Entries", 2000, 0, 10000);

    // 4. Extract Data & Fill Weights
    Long64_t nEntries = chain->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        chain->GetEntry(i);
        if (q2_val > -1.0) {
            vQ2.push_back((double)q2_val);
            hWeightQ2->Fill(q2_w);
        }
        if (q3_val > -1.0) {
            vQ3.push_back((double)q3_val);
            hWeightQ3->Fill(q3_w);
        }
    }
    delete chain; // free memory before sort

    std::sort(vQ2.begin(), vQ2.end());
    std::sort(vQ3.begin(), vQ3.end());

    // 5. Calculate Quantile Boundaries
    double cuts2[N_QBINS+1], cuts3[N_QBINS+1];
    int n2 = vQ2.size();
    int n3 = vQ3.size();
    for (int k = 0; k <= N_QBINS; k++) {
        int idx2 = (k == N_QBINS) ? n2 - 1 : (k * n2 / N_QBINS);
        int idx3 = (k == N_QBINS) ? n3 - 1 : (k * n3 / N_QBINS);
        cuts2[k] = (n2 > 0) ? vQ2[idx2] : 0;
        cuts3[k] = (n3 > 0) ? vQ3[idx3] : 0;
    }

    // Save cut text files
    std::ofstream outCut(Form("temp_cuts_MB0to1_charge/cuts_bin%d.txt", centBin));
    for (int k = 0; k <= N_QBINS; k++) outCut << cuts2[k] << (k == N_QBINS ? "" : ",");
    outCut << std::endl;
    for (int k = 0; k <= N_QBINS; k++) outCut << cuts3[k] << (k == N_QBINS ? "" : ",");
    outCut.close();

    // 6. Fill Slice Histograms
    dirSlices->cd();
    TH1D *hTotalQ2 = new TH1D(Form("hRawQ2_c%d", centBin), "Total q2", 7000, 0, 0.7);
    TH1D *hTotalQ3 = new TH1D(Form("hRawQ3_c%d", centBin), "Total q3", 7000, 0, 0.7);

    TH1D *hSliceQ2[N_QBINS], *hSliceQ3[N_QBINS];
    for (int k = 0; k < N_QBINS; k++) {
        hSliceQ2[k] = new TH1D(Form("hSliceQ2_c%d_q%d", centBin, k+1), "", 7000, 0, 0.7);
        hSliceQ3[k] = new TH1D(Form("hSliceQ3_c%d_q%d", centBin, k+1), "", 7000, 0, 0.7);
    }

    for (double val : vQ2) {
        hTotalQ2->Fill(val);
        for (int k = 0; k < N_QBINS; k++) {
            if (val >= cuts2[k] && (k == N_QBINS-1 ? val <= cuts2[k+1] : val < cuts2[k+1])) {
                hSliceQ2[k]->Fill(val);
                break;
            }
        }
    }
    for (double val : vQ3) {
        hTotalQ3->Fill(val);
        for (int k = 0; k < N_QBINS; k++) {
            if (val >= cuts3[k] && (k == N_QBINS-1 ? val <= cuts3[k+1] : val < cuts3[k+1])) {
                hSliceQ3[k]->Fill(val);
                break;
            }
        }
    }

    // 7. Write and Close
    fOut->Write();
    fOut->Close();
    std::cout << "Done centBin " << centBin << std::endl;
}
