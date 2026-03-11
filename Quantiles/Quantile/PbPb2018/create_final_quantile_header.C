#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TSystem.h"

void create_final_quantile_header() {
    std::ofstream hf("q_cuts_2018.h");
    hf << "#ifndef Q_CUTS_2023_H\n#define Q_CUTS_2023_H\n\n";

    std::vector<std::string> q2_lines, q3_lines;

    for (int i = 0; i < 90; i++) {
        std::ifstream fIn(Form("temp_cuts/cuts_bin%d.txt", i));
        std::string l2, l3;
        if (fIn && std::getline(fIn, l2) && std::getline(fIn, l3)) {
            q2_lines.push_back("  {" + l2 + "}");
            q3_lines.push_back("  {" + l3 + "}");
        } else {
            q2_lines.push_back("  {0,0,0,0,0,0,0,0,0,0,0}");
            q3_lines.push_back("  {0,0,0,0,0,0,0,0,0,0,0}");
        }
    }

    hf << "double q2_cuts[90][11] = {\n" << q2_lines[0];
    for (size_t i = 1; i < q2_lines.size(); ++i) hf << ",\n" << q2_lines[i];
    hf << "\n};\n\n";

    hf << "double q3_cuts[90][11] = {\n" << q3_lines[0];
    for (size_t i = 1; i < q3_lines.size(); ++i) hf << ",\n" << q3_lines[i];
    hf << "\n};\n\n#endif";
    hf.close();

    std::cout << "Header q_cuts_2018.h created." << std::endl;
}
