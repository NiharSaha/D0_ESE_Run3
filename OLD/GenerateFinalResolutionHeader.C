#include <iostream>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "TSystem.h"
#include "TString.h"

using namespace std;

void GenerateFinalResolutionHeader(
				   TString input_dir= "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Resolution_MB11to21_Jan22_v2/ROOT",
				   TString final_merged_name = "Full_Resolution_Ingredients.root",
				   TString header_file = "ESE_Resolution_Factors_Jan22_v2.h"

				   ) {
    
    // --- 1. MERGE FILES ---
    // This command merges all individual job outputs into one file for global averaging
    cout << ">>> Merging parallel resolution files in " << input_dir << "..." << endl;
    TString hadd_cmd = TString::Format("hadd -f %s/%s %s/Resolution_*.root", input_dir.Data(), final_merged_name.Data(), input_dir.Data());
    gSystem->Exec(hadd_cmd);

    // --- 2. OPEN MERGED FILE ---
    TFile *fIn = TFile::Open(input_dir + "/" + final_merged_name, "READ");
    if (!fIn || fIn->IsZombie()) {
        cout << "Error: Could not open merged file!" << endl;
        return;
    }

    TH2D *hRes_v2 = (TH2D*)fIn->Get("hRes_v2");
    TH2D *hRes_v3 = (TH2D*)fIn->Get("hRes_v3");

    if (!hRes_v2 || !hRes_v3) {
        cout << "Error: Resolution histograms not found in merged file!" << endl;
        return;
    }

    // --- 3. PREPARE HEADER FILE ---
    ofstream hFile(header_file);
    hFile << "#ifndef ESE_RESOLUTION_FACTORS_H" << endl;
    hFile << "#define ESE_RESOLUTION_FACTORS_H" << endl << endl;
    hFile << "/* Auto-generated Resolution Factors for SP Method */" << endl;
    hFile << "/* Index = (CentGroup * 10) + QuantileIndex [0-39] */" << endl << endl;

    auto ProcessHarmonic = [&](TH2D* h, TString name, int harmonic) {
        double R_yPlus[40];  // Rn(B) for D0 y > 0
        double R_yMinus[40]; // Rn(A) for D0 y < 0

        for (int i = 0; i < 40; i++) {
            // Bin 1: <QA*QB>, Bin 2: <QA*QC>, Bin 3: <QB*QC>
            double AB = h->GetBinContent(1, i + 1);
            double AC = h->GetBinContent(2, i + 1);
            double BC = h->GetBinContent(3, i + 1);

            // Rn for HF- (B): Used when D0 y > 0 (opposite subevent)
            // Formula: sqrt( (<AB> * <BC>) / <AC> )
            R_yPlus[i] = (AC > 1e-10 && (AB * BC / AC) > 0) ? TMath::Sqrt((AB * BC) / AC) : 0;

            // Rn for HF+ (A): Used when D0 y < 0 (opposite subevent)
            // Formula: sqrt( (<AB> * <AC>) / <BC> )
            R_yMinus[i] = (BC > 1e-10 && (AB * AC / BC) > 0) ? TMath::Sqrt((AB * AC) / BC) : 0;
        }

        // Write arrays to header
        hFile << "double R" << harmonic << "_yPlus[40] = {";
        for(int i=0; i<40; i++) hFile << Form("%.8f", R_yPlus[i]) << (i==39 ? "" : ", ");
        hFile << "};" << endl;

        hFile << "double R" << harmonic << "_yMinus[40] = {";
        for(int i=0; i<40; i++) hFile << Form("%.8f", R_yMinus[i]) << (i==39 ? "" : ", ");
        hFile << "};" << endl << endl;
    };

    ProcessHarmonic(hRes_v2, "v2", 2);
    ProcessHarmonic(hRes_v3, "v3", 3);

    hFile << "#endif" << endl;
    hFile.close();

    cout << ">>> SUCCESS: " << header_file <<"  generated." << endl;
}
