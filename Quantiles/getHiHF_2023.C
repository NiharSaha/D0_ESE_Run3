// getHiHF_2023.C
// Usage: root -l -b -q 'getHiHF_2023.C("filelist_2023PbPb.txt", "hiHF_2023PbPb.root", "2023")'
//
// For skimEDM ntuples where HLT trigger + event selections were already
// applied upstream during skimming -- so we just read hiHF directly, no
// friend trees / trigger OR / filter AND needed.
//
// filelist_2023PbPb.txt format: one file path per line, e.g.
//   root://eos.cms.rcac.purdue.edu//store/user/.../TTree_fromEDM_1.root
//   # comment lines and blank lines are skipped

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <vector>

void getHiHF_2023(TString fileList    = "filelist_2023PbPb.txt",
                   TString outputFile = "hiHF_PbPb2023_new.root",
                   TString tag        = "2023",
                   Long64_t maxFiles  = -1)   // -1 = process all files in the list
{
    // --- Read file list ---
    std::ifstream infile(fileList.Data());
    if (!infile.is_open()) {
        std::cerr << "ERROR: could not open file list " << fileList << std::endl;
        return;
    }

    std::vector<TString> filePaths;
    std::string line;
    while (std::getline(infile, line)) {
        TString tline(line);
        tline.ReplaceAll(" ", "");
        if (tline.IsNull()) continue;
        if (tline.BeginsWith("#")) continue;
        filePaths.push_back(tline);
    }
    infile.close();

    std::cout << "Found " << filePaths.size() << " files listed in " << fileList << std::endl;

    if (maxFiles > 0 && maxFiles < (Long64_t)filePaths.size()) {
        filePaths.resize(maxFiles);
        std::cout << "Limiting to first " << maxFiles << " files." << std::endl;
    }

    if (filePaths.empty()) {
        std::cerr << "ERROR: no valid file paths found in list!" << std::endl;
        return;
    }

    // --- Build chain: note different TDirectory/tree name vs standard HiForest ---
    TChain* evtChain = new TChain("eventinfoana/EventInfoNtuple");

    int nAdded = 0;
    for (auto& fp : filePaths) {
        int ok = evtChain->Add(fp);
        if (ok) nAdded++;
        else std::cerr << "WARNING: failed to add (or tree missing) for: " << fp << std::endl;
    }
    std::cout << "Successfully added " << nAdded << " / " << filePaths.size() << " files to chain." << std::endl;
    if (nAdded == 0) { std::cerr << "ERROR: nothing added, aborting." << std::endl; return; }

    Long64_t nEntries = evtChain->GetEntries();
    std::cout << "Total entries: " << nEntries << std::endl;
    if (nEntries == 0) { std::cerr << "ERROR: chain has zero entries, aborting." << std::endl; return; }

    // --- Only read what we need ---
    evtChain->SetBranchStatus("*", 0);
    evtChain->SetBranchStatus("HFsumET", 1);

    Float_t hiHF = 0;
    Int_t status_hiHF = evtChain->SetBranchAddress("HFsumET", &hiHF);
    if (status_hiHF < 0) {
        std::cerr << "ERROR setting address for hiHF, status=" << status_hiHF
                   << " -- branch name may differ in this ntuple. Run EventInfoNtuple->Print() "
                   << "on one file to confirm the exact branch name." << std::endl;
        return;
    }

    // --- Histogram ---
    TString hName = "h_hiHF_" + tag;
    TH1D* h_hiHF = new TH1D(hName, "hiHF distribution (selected);hiHF;Events", 5000, 0, 10000);
    h_hiHF->Sumw2();

    // --- Progress bookkeeping: print every 10% ---
    Long64_t printEvery = nEntries / 10;
    if (printEvery == 0) printEvery = 1;

    std::cout << "\nStarting event loop...\n" << std::endl;

    for (Long64_t i = 0; i < nEntries; i++) {
        evtChain->GetEntry(i);

        if (i % printEvery == 10) {
            double frac = 100.0 * i / nEntries;
            std::cout << "Processed: " << Form("%.0f", frac) << "%" << std::endl;
        }

        h_hiHF->Fill(hiHF);
    }
    std::cout << "Processed: 100%" << std::endl;

    std::cout << "\nFilled entries: " << h_hiHF->GetEntries() << std::endl;

    TFile* fout = new TFile(outputFile, "RECREATE");
    h_hiHF->Write();
    fout->Close();

    std::cout << "Saved histogram '" << hName << "' to " << outputFile << std::endl;

    delete evtChain;
}
