#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TSystem.h"

const int N_CENT_BINS = 50;


void create_final_quantile_header(bool USE_NON_UNIFORM = true, const std::string& methodTag = "HIST_NU") {
    
    // Set bins based on uniform/non-uniform flag
    int activeQBins = USE_NON_UNIFORM ? 12 : 20;

    // 1. Dynamic Output Filename
    std::string filename = "quantile_12NUQbin_diffq2q3_cuts_2023_MB0to31.h";
    std::ofstream hf(filename);
    if (!hf.is_open()) {
        std::cerr << "Error: Could not create file " << filename << std::endl;
        return;
    }

    hf << "#ifndef Q_CUTS_2023_H\n#define Q_CUTS_2023_H\n\n";

    // 2. Dynamic Fallback Zeros String
    // Generates "{0,0,0...}" exactly matching activeQBins + 1
    std::string fallback_zeros = "  {0";
    for (int i = 0; i < activeQBins; i++) {
        fallback_zeros += ",0";
    }
    fallback_zeros += "}";

    std::vector<std::string> q2_lines, q3_lines;

    // 3. Dynamic Input Directory
    std::string inDir = Form("temp_cuts_%s", methodTag.c_str());

    for (int i = 0; i < N_CENT_BINS; i++) {
        std::ifstream fIn(Form("%s/cuts_bin%d.txt", inDir.c_str(), i));
        std::string l2, l3;
        if (fIn && std::getline(fIn, l2) && std::getline(fIn, l3)) {
            q2_lines.push_back("  {" + l2 + "}");
            q3_lines.push_back("  {" + l3 + "}");
        } else {
            q2_lines.push_back(fallback_zeros);
            q3_lines.push_back(fallback_zeros);
        }
    }

    // 4. Write arrays with correct dimensions
    hf << "double q2_cuts[" << N_CENT_BINS << "][" << activeQBins + 1 << "] = {\n" << q2_lines[0];
    for (size_t i = 1; i < q2_lines.size(); ++i) hf << ",\n" << q2_lines[i];
    hf << "\n};\n\n";

    hf << "double q3_cuts[" << N_CENT_BINS << "][" << activeQBins + 1 << "] = {\n" << q3_lines[0];
    for (size_t i = 1; i < q3_lines.size(); ++i) hf << ",\n" << q3_lines[i];
    hf << "\n};\n\n#endif";
    hf.close();

    std::cout << "Header file: " << filename << " created from directory: " << inDir << std::endl;
}
