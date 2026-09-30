#ifndef BDTHANDLER_H
#define BDTHANDLER_H

#include <tuple>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

using namespace std;

constexpr int y_bins = 2;
constexpr int cent_bins = 4;
constexpr int pT_bins = 13;

class BDTHandler {
private:
    double bdtCuts[y_bins][cent_bins][pT_bins];

    tuple<int, int, int> get_bins(float y, int centrality, float pT) const;

public:
    BDTHandler();
    void loadCuts(const string &cuts_filename);
    double getBDTCut(float y, int centrality, float pT) const;
};

#endif

