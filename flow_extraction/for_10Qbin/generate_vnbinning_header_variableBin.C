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

// Returns sort key or -1 if not matched
static int get_sort_key(const string &ident) {
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
    return -1;
}

void generate_vnbinning_header_variableBin(
    const char *input_dir  = "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/save_vnbinning_outputs_Mar8_v1_adaptBin/output",
    const char *out_header = "/home/saha115/D0_ESE/CMSSW_13_2_11/src/flow_extraction/vnbinning_generated_variBin_Mar8.h")
{
    string INPUT_DIR  = input_dir;
    string OUT_HEADER = out_header;

    if (!fs::is_directory(INPUT_DIR)) {
        cerr << "Input path not a directory: " << INPUT_DIR << "\n"; return;
    }

    // ── updated regex: matches  v2_cent30to40_pT4to5[47]={ ... } ─────────────
    // group 1: name (e.g. v2_cent30to40_pT4to5)
    // group 2: array size (e.g. 47)
    // group 3: values inside { }
    regex line_re(R"(^\s*(v[23]_[^=\[\s]+)\s*\[(\d+)\]\s*=\s*\{([^}]*)\})");
    // ─────────────────────────────────────────────────────────────────────────

    // map: ident -> (latest_time, size, values)
    struct ArrayInfo {
        fs::file_time_type mtime;
        size_t array_size;          // size read from [N] in txt file
        vector<double> values;
    };
    map<string, ArrayInfo> arrays;

    for (auto &entry : fs::directory_iterator(INPUT_DIR)) {
        if (!entry.is_regular_file()) continue;
        if (entry.path().extension() != ".txt") continue;
        ifstream ifs(entry.path());
        if (!ifs.is_open()) continue;

        auto mtime = fs::last_write_time(entry.path());

        // accumulate full logical lines (handle backslash line continuation)
        string logical_line, raw_line;
        while (getline(ifs, raw_line)) {
            // strip trailing backslash if present (wrapped lines in txt)
            if (!raw_line.empty() && raw_line.back() == '\\')
                raw_line.pop_back();
            logical_line += raw_line;

            // if the logical line contains a closing }, process it
            if (logical_line.find('}') == string::npos) continue;

            smatch m;
            if (regex_search(logical_line, m, line_re)) {
                string name      = m[1].str();                  // e.g. v2_cent30to40_pT4to5
                size_t arr_size  = (size_t)stoul(m[2].str());   // e.g. 47
                string vals_str  = m[3].str();

                vector<double> v;
                string token;
                stringstream ss(vals_str);
                while (getline(ss, token, ',')) {
                    size_t a = token.find_first_not_of(" \t\r\n");
                    if (a == string::npos) continue;
                    size_t b = token.find_last_not_of(" \t\r\n");
                    string tok = token.substr(a, b - a + 1);
                    try { v.push_back(stod(tok)); } catch(...) {}
                }
                if (v.empty()) { logical_line.clear(); continue; }

                // warn if parsed count doesn't match declared size
                if (v.size() != arr_size) {
                    cerr << "WARNING: " << name
                         << " declared size=" << arr_size
                         << " but parsed " << v.size() << " values. Using parsed count.\n";
                    arr_size = v.size();
                }

                string ident = sanitize_ident(string("vnbinning_") + name);
                auto it = arrays.find(ident);
                if (it == arrays.end() || mtime > it->second.mtime)
                    arrays[ident] = {mtime, arr_size, move(v)};
            }
            logical_line.clear();
        }
    }

    if (arrays.empty()) { cerr << "No arrays found in " << INPUT_DIR << "\n"; return; }

    // --- Build ordered list of keys ---
    vector<pair<int,string>> ordered_keys;
    vector<string> unknown_keys;

    for (auto &kv : arrays) {
        int key = get_sort_key(kv.first);
        if (key >= 0) ordered_keys.push_back({key, kv.first});
        else          unknown_keys.push_back(kv.first);
    }

    sort(ordered_keys.begin(), ordered_keys.end(),
         [](const pair<int,string> &a, const pair<int,string> &b){
             return a.first < b.first; });
    sort(unknown_keys.begin(), unknown_keys.end());

    // --- Write header ---
    ofstream os(OUT_HEADER);
    if (!os.is_open()) { cerr << "Cannot write header: " << OUT_HEADER << "\n"; return; }

    os << "// Auto-generated vnbinning header (variable bin sizes)\n";
    os << "#pragma once\n\n";
    os << "#include <cstddef>\n\n";

    auto write_array = [&](const string &ident) {
        const ArrayInfo &info = arrays[ident];
        size_t N = info.array_size;
        const vector<double> &vals = info.values;

        // ── write size constant: e.g. static const int vnbinning_v2_cent30to40_pT4to5_N = 47;
        os << "static const int " << ident << "_N = " << N << ";\n";

        // ── write array with actual size from txt file ────────────────────────
        os << "static const double " << ident << "[" << N << "] = {";
        for (size_t i = 0; i < N; ++i) {
            os.setf(ios::fixed);
            os << setprecision(3) << vals[i];   // 3 decimal places to preserve precision
            if (i + 1 < N) os << ", ";
        }
        os << "};\n\n";
    };

    int written = 0;
    for (auto &kp : ordered_keys) { write_array(kp.second); ++written; }

    if (!unknown_keys.empty()) {
        os << "// WARNING: unrecognized array names (not matched to cent/pT order):\n";
        for (auto &ident : unknown_keys) {
            write_array(ident); ++written;
            cerr << "WARNING: unrecognized array: " << ident << "\n";
        }
    }

    os.close();
    cout << "Wrote header: " << OUT_HEADER
         << " (" << written << " arrays, "
         << unknown_keys.size() << " unrecognized)\n";

    // ── print size summary for quick sanity check ─────────────────────────────
    cout << "\nArray sizes found:\n";
    for (auto &kp : ordered_keys) {
        cout << "  " << kp.second
             << " [" << arrays[kp.second].array_size << "]\n";
    }
}
