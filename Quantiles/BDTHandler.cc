#include "BDTHandler.h"

BDTHandler::BDTHandler() {
  string bdt_cuts_path = "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Quantiles/bdt_cuts.csv";
  loadCuts(bdt_cuts_path);
}

tuple<int, int, int> BDTHandler::get_bins(float y, int centrality, float pT) const {
    int y_bin = -1, cent_bin = -1, pT_bin = -1;

    if (centrality >= 0 && centrality < 20) cent_bin = 0;
    else if (centrality >= 20 && centrality < 60) cent_bin = 1;
    else if (centrality >= 60 && centrality < 100) cent_bin = 2;
    else if (centrality >= 100 && centrality < 180) cent_bin = 3;

    if (pT >= 1 && pT < 2) pT_bin = 0;
    else if (pT >= 2 && pT < 3) pT_bin = 1;
    else if (pT >= 3 && pT < 4) pT_bin = 2;
    else if (pT >= 4 && pT < 5) pT_bin = 3;
    else if (pT >= 5 && pT < 6) pT_bin = 4;
    else if (pT >= 6 && pT < 8) pT_bin = 5;
    else if (pT >= 8 && pT < 10) pT_bin = 6;
    else if (pT >= 10 && pT < 15) pT_bin = 7;
    else if (pT >= 15 && pT < 20) pT_bin = 8;
    else if (pT >= 20 && pT < 40) pT_bin = 9;
    else if (pT >= 40 && pT < 60) pT_bin = 10;
    else if (pT >= 60 && pT < 100) pT_bin = 11;
    else if (pT >= 100 && pT < 900) pT_bin = 12;

    if (abs(y) < 1) y_bin = 0;
    else if (abs(y) >= 1 && abs(y) < 3) y_bin = 1;

    return make_tuple(y_bin, cent_bin, pT_bin);
}

void BDTHandler::loadCuts(const string &cuts_filename) {
    ifstream file(cuts_filename);
    cout << "cuts filename: " << cuts_filename << endl;
    if (!file.is_open()) {
        cout << "Error opening file: " << cuts_filename << " BREAK!!!! " << endl;
        return;
    }

    string line;
    while (getline(file, line)) {
        if (line[0] == '#') continue;
        stringstream ss(line);
        int y_bin = -9, cent_bin = -9, pT_bin = -9;
        double cut_value;
        ss >> y_bin >> cent_bin >> pT_bin >> cut_value;
        bdtCuts[y_bin][cent_bin][pT_bin] = cut_value;
    }
    file.close();
}

double BDTHandler::getBDTCut(float y, int centrality, float pT) const {
    int y_bin, cent_bin, pT_bin;
    tie(y_bin, cent_bin, pT_bin) = get_bins(y, centrality, pT);
    if (y_bin < 0 || cent_bin < 0 || pT_bin < 0) return 999;
    return bdtCuts[y_bin][cent_bin][pT_bin];
}

