#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TSystem.h"

void create_final_quantile_header() {
    std::ofstream hf("20quantile_cuts_2023_MB0to1_charge_check_Jun11.h");
    hf << "#ifndef Q_CUTS_2023_H\n#define Q_CUTS_2023_H\n\n";

    std::vector<std::string> q2_lines, q3_lines;

    int cent[11] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50};     
    for (int i = 0; i < 10; i++) {
      std::ifstream fIn(Form("temp_cuts_MB0to1_charge/cuts_bin%d_%d.txt", cent[i], cent[i+1]));
        std::string l2, l3;
        if (fIn && std::getline(fIn, l2) && std::getline(fIn, l3)) {
            q2_lines.push_back("  {" + l2 + "}");
            q3_lines.push_back("  {" + l3 + "}");
        } else {
            q2_lines.push_back("  {0,0,0,0,0,0,0,0,0,0,0}");
            q3_lines.push_back("  {0,0,0,0,0,0,0,0,0,0,0}");
        }
    }

    hf << "double q2_cuts[10][21] = {\n" << q2_lines[0];
    for (size_t i = 1; i < q2_lines.size(); ++i) hf << ",\n" << q2_lines[i];
    hf << "\n};\n\n";

    hf << "double q3_cuts[10][21] = {\n" << q3_lines[0];
    for (size_t i = 1; i < q3_lines.size(); ++i) hf << ",\n" << q3_lines[i];
    hf << "\n};\n\n#endif";
    hf.close();

    std::cout << "Header file created." << std::endl;
}
