// root-friendly macro: call .x generate_vnbinning_header.C() or generate_vnbinning_header(in_dir,out_header)
// Scans directory for lines like `v2_cent0to10_pT1to2 = { ... }` and writes simplified arrays
// Arrays are written as: static const double vnbinning_v2_cent0to10_pT1to2[45] = { ... };
// Numbers printed with two decimal places.
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <regex>
#include <string>
#include <map>
#include <vector>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;

static string sanitize_ident(const string &s) {
    string out;
    for (char c: s) {
        if (isalnum((unsigned char)c) || c=='_') out.push_back(c);
        else out.push_back('_');
    }
    if (!out.empty() && isdigit((unsigned char)out.front())) out = string("_") + out;
    return out;
}

// --- Define fixed ordering ---
static const vector<string> CEN_ORDER = {
    "cent0to10", "cent10to20", "cent20to30",
    "cent30to40", "cent40to50", "cent50to80"
};

static const vector<string> PT_ORDER = {
    "pT1to2",   "pT2to3",   "pT3to4",   "pT4to5",
    "pT5to6",   "pT6to8",   "pT8to10",  "pT10to15",
    "pT15to20", "pT20to40", "pT40to60", "pT60to100"
};

static const vector<string> VN_ORDER = {"v2", "v3"};

// Returns sort key (vn_idx * 1000 + cen_idx * 100 + pt_idx), or -1 if not matched
static int get_sort_key(const string &ident) {
    // ident format: vnbinning_v2_cent0to10_pT1to2
    for (int iv = 0; iv < (int)VN_ORDER.size(); ++iv) {
        for (int ic = 0; ic < (int)CEN_ORDER.size(); ++ic) {
            for (int ip = 0; ip < (int)PT_ORDER.size(); ++ip) {
                string expected = sanitize_ident(
                    string("vnbinning_") + VN_ORDER[iv]
                    + "_" + CEN_ORDER[ic]
                    + "_" + PT_ORDER[ip]
                );
                if (ident == expected)
                    return iv * 10000 + ic * 100 + ip;
            }
        }
    }
    return -1;  // unknown key → goes to end
}

void generate_vnbinning_header(
    const char *input_dir  = "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/save_vnbinning_outputs_Mar5_v2/output",
    const char *out_header = "/home/saha115/D0_ESE/CMSSW_13_2_11/src/flow_extraction/vnbinning_generated_Mar5.h")
{
    const size_t TARGET_N = 63; // N_VBINS+1 (32+1)
    string INPUT_DIR  = input_dir;
    string OUT_HEADER = out_header;

    if (!fs::is_directory(INPUT_DIR)) {
        cerr << "Input path not a directory: " << INPUT_DIR << "\n"; return;
    }

    regex line_re(R"(^\s*(v[23]_[^=\s]+)\s*=\s*\{([^}]*)\})");

    // map: simple_name -> (latest_time, values)
    map<string, pair<fs::file_time_type, vector<double>>> arrays;

    for (auto &entry : fs::directory_iterator(INPUT_DIR)) {
        if (!entry.is_regular_file()) continue;
        if (entry.path().extension() != ".txt") continue;
        ifstream ifs(entry.path());
        if (!ifs.is_open()) continue;
        string line;
        auto mtime = fs::last_write_time(entry.path());
        while (getline(ifs, line)) {
            smatch m;
            if (!regex_search(line, m, line_re)) continue;
            string name = m[1].str();
            string vals = m[2].str();
            vector<double> v;
            string token;
            stringstream ss(vals);
            while (getline(ss, token, ',')) {
                size_t a = token.find_first_not_of(" \t\r\n");
                if (a == string::npos) continue;
                size_t b = token.find_last_not_of(" \t\r\n");
                string tok = token.substr(a, b - a + 1);
                try { v.push_back(stod(tok)); } catch(...) {}
            }
            if (v.empty()) continue;
            string simple = sanitize_ident(string("vnbinning_") + name);
            auto it = arrays.find(simple);
            if (it == arrays.end() || mtime > it->second.first)
                arrays[simple] = make_pair(mtime, move(v));
        }
    }

    if (arrays.empty()) { cerr << "No arrays found in " << INPUT_DIR << "\n"; return; }

    // --- Build ordered list of keys ---
    // 1) separate known keys (with valid sort key) from unknown
    vector<pair<int,string>> ordered_keys;
    vector<string> unknown_keys;

    for (auto &kv : arrays) {
        int key = get_sort_key(kv.first);
        if (key >= 0)
            ordered_keys.push_back({key, kv.first});
        else
            unknown_keys.push_back(kv.first);
    }

    // 2) sort known keys by sort key (vn → cent → pT)
    sort(ordered_keys.begin(), ordered_keys.end(),
         [](const pair<int,string> &a, const pair<int,string> &b){
             return a.first < b.first;
         });

    // 3) sort unknown keys alphabetically and append at end
    sort(unknown_keys.begin(), unknown_keys.end());

    // --- Write header ---
    ofstream os(OUT_HEADER);
    if (!os.is_open()) { cerr << "Cannot write header: " << OUT_HEADER << "\n"; return; }

    os << "// Auto-generated vnbinning header\n";
    os << "#pragma once\n\n";
    os << "#include <cstddef>\n\n";

    // write known keys in physics order
    int written = 0;
    for (auto &kp : ordered_keys) {
        const string &ident = kp.second;
        vector<double> vals = arrays[ident].second;
        if (vals.size() < TARGET_N) vals.resize(TARGET_N, 0.0);
        else if (vals.size() > TARGET_N) vals.resize(TARGET_N);

        os << "static const double " << ident << "[" << TARGET_N << "] = {";
        for (size_t i = 0; i < TARGET_N; ++i) {
            os.setf(ios::fixed);
            os << setprecision(2) << vals[i];
            if (i + 1 < TARGET_N) os << ", ";
        }
        os << "};\n\n";
        ++written;
    }

    // write unknown keys at end with a warning comment
    if (!unknown_keys.empty()) {
        os << "// WARNING: unrecognized array names (not matched to cent/pT order):\n";
        for (auto &ident : unknown_keys) {
            vector<double> vals = arrays[ident].second;
            if (vals.size() < TARGET_N) vals.resize(TARGET_N, 0.0);
            else if (vals.size() > TARGET_N) vals.resize(TARGET_N);

            os << "static const double " << ident << "[" << TARGET_N << "] = {";
            for (size_t i = 0; i < TARGET_N; ++i) {
                os.setf(ios::fixed);
                os << setprecision(2) << vals[i];
                if (i + 1 < TARGET_N) os << ", ";
            }
            os << "};\n\n";
            ++written;
            cerr << "WARNING: unrecognized array: " << ident << "\n";
        }
    }

    os.close();
    cout << "Wrote header: " << OUT_HEADER
         << " (" << written << " arrays, "
         << unknown_keys.size() << " unrecognized)\n";


}
