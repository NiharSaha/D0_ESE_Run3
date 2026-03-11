#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include "TFile.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void ExtractQuantiles_sortFunc() {
    TString dirPath = "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Q2quantile_Jan20_MB11to21/ROOT/";
    
    std::cout << "--- Starting Memory-Efficient ESE Calculation ---" << std::endl;
    TChain *nt = new TChain("nt_ese_global");
    nt->Add(dirPath + "*.root"); 

    Long64_t totalEntries = nt->GetEntries();
    std::cout << "Chained " << nt->GetListOfFiles()->GetEntries() << " files." << std::endl;
    std::cout << "Total entries to process: " << totalEntries << std::endl;

    // We use vectors to store qn values for each 1% centrality bin (0-89)
    std::vector<double> q2_vecs[90];
    std::vector<double> q3_vecs[90];

    // TTreeReader is the modern, memory-efficient way to process large TChains
    TTreeReader reader(nt);
    TTreeReaderValue<float> cent(reader, "cent");
    TTreeReaderValue<float> q2(reader, "q2");
    TTreeReaderValue<float> q3(reader, "q3");

    Long64_t counter = 0;
    while (reader.Next()) {
        int c = (int)(*cent);
        if (c >= 0 && c < 90) {
            q2_vecs[c].push_back(*q2);
            q3_vecs[c].push_back(*q3);
        }

        counter++;
        if (counter % 10000000 == 0) {
            std::cout << "Processed " << counter / 1000000 << "M entries..." << std::endl;
        }
    }

    // --- Calculate Quantiles ---
    double q2_table[10][90] = {{0}};
    double q3_table[10][90] = {{0}};

    for (int c = 0; c < 90; c++) {
        if (q2_vecs[c].empty()) continue;

        std::sort(q2_vecs[c].begin(), q2_vecs[c].end());
        std::sort(q3_vecs[c].begin(), q3_vecs[c].end());

        for (int q = 0; q < 10; q++) {
            double percentile = (q + 1) * 0.1;
            int idx = (int)(percentile * q2_vecs[c].size()) - 1;
            if (idx < 0) idx = 0;

            q2_table[q][c] = q2_vecs[c][idx];
            q3_table[q][c] = q3_vecs[c][idx];
        }
    }

    // --- Generate Header ---
    std::ofstream header("ESE_Quantile_Table.h");
    header << "#ifndef ESE_QUANTILE_TABLE_H\n#define ESE_QUANTILE_TABLE_H\n\n";

    auto writeArray = [&](const char* name, double table[10][90]) {
        header << "const double " << name << "[10][90] = {\n";
        for (int q = 0; q < 10; q++) {
            header << "  {";
            for (int c = 0; c < 90; c++) {
                header << table[q][c] << (c == 89 ? "" : ", ");
            }
            header << "}" << (q == 9 ? "" : ",\n");
        }
        header << "\n};\n\n";
    };

    writeArray("q2_binning", q2_table);
    writeArray("q3_binning", q3_table);
    header << "#endif";
    header.close();

    std::cout << "SUCCESS: ESE_Quantile_Table.h generated using streaming." << std::endl;
}
