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
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

const int N_QBINS = 20; 

void ExtractQuantiles(int centBin, const char *inputDir, const char *outTag, bool USE_SORT_METHOD, bool USE_NON_UNIFORM)
{

    // Determine the active number of bins depending on the mode
    int activeQBins = USE_NON_UNIFORM ? 12 : N_QBINS;

    // 1. Scan directory and build TChain
    TChain *chain = new TChain(Form("Q_ntuple/nt_qn_cent%i_%i", centBin, centBin + 1));

    TSystemDirectory dir(inputDir, inputDir);
    TList *files = dir.GetListOfFiles();
    if (!files)
    {
        std::cerr << "Cannot open directory: " << inputDir << std::endl;
        return;
    }

    int nFilesAdded = 0;
    TIter next(files);
    TSystemFile *sFile;
    while ((sFile = (TSystemFile *)next()))
    {
        TString fname = sFile->GetName();
        if (!sFile->IsDirectory() && fname.EndsWith(".root") && fname.Contains("Quantiles_ESE_ntuples_"))
        {
            TString fullPath = TString::Format("%s/%s", inputDir, fname.Data());

            // Check file is not corrupted before adding
            TFile *ftest = TFile::Open(fullPath.Data(), "READ");
            if (!ftest || ftest->IsZombie() || ftest->TestBit(TFile::kRecovered))
            {
                std::cerr << "Skipping bad/recovered file: " << fullPath << std::endl;
                if (ftest)
                {
                    ftest->Close();
                    delete ftest;
                }
                continue;
            }
            // Check ntuple exists in this file
            TObject *obj = ftest->Get(Form("Q_ntuple/nt_qn_cent%i_%i", centBin, centBin + 1));
            if (!obj)
            {
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

    if (chain->GetEntries() == 0)
    {
        std::cout << "Bin " << centBin << " empty or missing." << std::endl;
        delete chain;
        return;
    }

    // 2. Set Branch Addresses - disable unused branches for speed
    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("q2_hf_total", 1);
    chain->SetBranchStatus("q3_hf_total", 1);
    chain->SetBranchStatus("q2_hf_total_w", 1);
    chain->SetBranchStatus("q3_hf_total_w", 1);

    std::vector<double> vQ2, vQ3;

    TH1D *hForQuantiles2 = nullptr;
    TH1D *hForQuantiles3 = nullptr;

    if (USE_SORT_METHOD)
    {
        vQ2.reserve(chain->GetEntries());
        vQ3.reserve(chain->GetEntries());
    }
    else
    {
        hForQuantiles2 = new TH1D("hTemp2", "", 35000, 0, 0.7);
        hForQuantiles3 = new TH1D("hTemp3", "", 35000, 0, 0.7);
    }

    // Ensure output directories exist
    TString cutDir = Form("temp_cuts_%s", outTag);
    TString rootDir = Form("temp_roots_%s", outTag);

    if (gSystem->AccessPathName(cutDir))
        gSystem->mkdir(cutDir);
    if (gSystem->AccessPathName(rootDir))
        gSystem->mkdir(rootDir);

    // 3. Create Output File
    TFile *fOut = new TFile(Form("%s/hist_bin%d.root", rootDir.Data(), centBin), "RECREATE");
    TDirectory *dirWeights = fOut->mkdir("Weights");
    TDirectory *dirSlices = fOut->mkdir("Slices");

    dirWeights->cd();
    TH1D *hWeightQ2 = new TH1D(Form("hWeightQ2_c%d", centBin), "Total Weight q2;E_{T} Sum;Entries", 2000, 0, 10000);
    TH1D *hWeightQ3 = new TH1D(Form("hWeightQ3_c%d", centBin), "Total Weight q3;E_{T} Sum;Entries", 2000, 0, 10000);

    TTreeReader reader(chain);
    TTreeReaderValue<float> q2(reader, "q2_hf_total");
    TTreeReaderValue<float> q3(reader, "q3_hf_total");
    TTreeReaderValue<float> q2_w(reader, "q2_hf_total_w");
    TTreeReaderValue<float> q3_w(reader, "q3_hf_total_w");

    while (reader.Next())
    {
        if (USE_SORT_METHOD)
        {
            if (*q2 > -1.0)
                vQ2.push_back((double)*q2);
            if (*q3 > -1.0)
                vQ3.push_back((double)*q3);
        }
        else
        {
            if (*q2 > -1.0)
                hForQuantiles2->Fill(*q2);
            if (*q3 > -1.0)
                hForQuantiles3->Fill(*q3);
        }
        // Weights are always filled
        hWeightQ2->Fill(*q2_w);
        hWeightQ3->Fill(*q3_w);
    }

    // 4. Setup Probability Distributions
    std::vector<double> probs2(activeQBins + 1);
    std::vector<double> probs3(activeQBins + 1);

    if (USE_NON_UNIFORM)
    {
        // 12 bins q2 distribution
        double nuProbs2[] = {0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1.0};
        // 12 bins q3 distribution
        double nuProbs3[] = {0.0, 0.075, 0.175, 0.275, 0.40, 0.60, 0.70, 0.80, 0.875, 0.925, 0.95, 0.975, 1.0};
        
        for (int k = 0; k <= activeQBins; k++)
        {
            probs2[k] = nuProbs2[k];
            probs3[k] = nuProbs3[k];
        }
    }
    else
    {
        // Uniform bins for both
        for (int k = 0; k <= activeQBins; k++)
        {
            probs2[k] = (double)k / activeQBins;
            probs3[k] = (double)k / activeQBins;
        }
    }

    // 5. Calculate Quantiles
    std::vector<double> cuts2(activeQBins + 1, 0.0);
    std::vector<double> cuts3(activeQBins + 1, 0.0);

    if (USE_SORT_METHOD)
    {
        std::sort(vQ2.begin(), vQ2.end());
        std::sort(vQ3.begin(), vQ3.end());
        Long64_t n2 = vQ2.size(), n3 = vQ3.size();
        for (int k = 0; k <= activeQBins; k++)
        {
            cuts2[k] = (n2 > 0) ? vQ2[std::min((Long64_t)(probs2[k] * n2), n2 - 1)] : 0;
            cuts3[k] = (n3 > 0) ? vQ3[std::min((Long64_t)(probs3[k] * n3), n3 - 1)] : 0;
        }
    }
    else
    {
        hForQuantiles2->GetQuantiles(activeQBins + 1, cuts2.data(), probs2.data());
        hForQuantiles3->GetQuantiles(activeQBins + 1, cuts3.data(), probs3.data());
    }

    if (USE_SORT_METHOD)
    {
        vQ2.clear();
        vQ2.shrink_to_fit();
        vQ3.clear();
        vQ3.shrink_to_fit();
        std::cout << "DEBUG: Vectors cleared and memory released." << std::endl;
    }

    // Save cut text files
    std::ofstream outCut(Form("%s/cuts_bin%d.txt", cutDir.Data(), centBin));
    if (!outCut.is_open())
    {
        std::cerr << "ERROR: Cannot open cut file for centBin " << centBin << std::endl;
        return;
    }
    for (int k = 0; k <= activeQBins; k++)
        outCut << cuts2[k] << (k == activeQBins ? "" : ",");
    outCut << std::endl;
    for (int k = 0; k <= activeQBins; k++)
        outCut << cuts3[k] << (k == activeQBins ? "" : ",");
    outCut.close();

    // 6. Fill Slice Histograms
    dirSlices->cd();
    TH1D *hTotalQ2 = new TH1D(Form("hRawQ2_c%d", centBin), "Total q2", 7000, 0, 0.7);
    TH1D *hTotalQ3 = new TH1D(Form("hRawQ3_c%d", centBin), "Total q3", 7000, 0, 0.7);

    // Dynamic allocation via vectors to handle changing activeQBins
    std::vector<TH1D *> hSliceQ2(activeQBins);
    std::vector<TH1D *> hSliceQ3(activeQBins);

    for (int k = 0; k < activeQBins; k++)
    {
        hSliceQ2[k] = new TH1D(Form("hSliceQ2_c%d_q%d", centBin, k + 1), "", 7000, 0, 0.7);
        hSliceQ3[k] = new TH1D(Form("hSliceQ3_c%d_q%d", centBin, k + 1), "", 7000, 0, 0.7);
    }

    // Iterate once more (or reset reader) to fill slices
    reader.Restart();
    while (reader.Next())
    {
        hTotalQ2->Fill(*q2);
        hTotalQ3->Fill(*q3);
        for (int k = 0; k < activeQBins; k++)
        {
            if (*q2 >= cuts2[k] && (k == activeQBins - 1 ? *q2 <= cuts2[k + 1] : *q2 < cuts2[k + 1]))
                hSliceQ2[k]->Fill(*q2);
            if (*q3 >= cuts3[k] && (k == activeQBins - 1 ? *q3 <= cuts3[k + 1] : *q3 < cuts3[k + 1]))
                hSliceQ3[k]->Fill(*q3);
        }
    }
    // 7. Write and Close
    fOut->Write();
    fOut->Close();
    delete chain;
    std::cout << "Done centBin " << centBin << std::endl;
}
