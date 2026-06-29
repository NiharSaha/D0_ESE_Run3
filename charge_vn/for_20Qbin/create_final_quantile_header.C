#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TSystem.h"

void create_final_quantile_header() {
    std::ofstream hf("quantile_cuts_20Qbin_2023_MB0to1_charge.h");
    hf << "#ifndef Q_CUTS_2023_H\n#define Q_CUTS_2023_H\n\n";

    std::vector<std::string> q2_lines, q3_lines;

    for (int i = 0; i < 50; i++) {
        std::ifstream fIn(Form("temp_cuts_MB0to1_charge/cuts_bin%d.txt", i));
        std::string l2, l3;
        if (fIn && std::getline(fIn, l2) && std::getline(fIn, l3)) {
            q2_lines.push_back("  {" + l2 + "}");
            q3_lines.push_back("  {" + l3 + "}");
        } else {
            q2_lines.push_back("  {0,0,0,0,0,0,0,0,0,0,0}");
            q3_lines.push_back("  {0,0,0,0,0,0,0,0,0,0,0}");
        }
    }

    hf << "double q2_cuts[50][21] = {\n" << q2_lines[0];
    for (size_t i = 1; i < q2_lines.size(); ++i) hf << ",\n" << q2_lines[i];
    hf << "\n};\n\n";

    hf << "double q3_cuts[50][21] = {\n" << q3_lines[0];
    for (size_t i = 1; i < q3_lines.size(); ++i) hf << ",\n" << q3_lines[i];
    hf << "\n};\n\n#endif";
    hf.close();

    std::cout << "Header file created." << std::endl;
}
