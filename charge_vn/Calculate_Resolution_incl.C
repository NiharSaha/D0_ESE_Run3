#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFileCollection.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TComplex.h"
#include "TProfile.h"
#include "THnSparse.h"

using namespace std;

void Calculate_Resolution_incl(TString input_txt, TString output_path, int istart, int iend)
{

    // --- 1. CONFIGURATION ---
    // Broad centrality bins for comparison with published results
    /*const int N_CENTBINS = 3;
    const int min_centbin[N_CENTBINS] = {0, 10, 30}; // in 1% units
    const int max_centbin[N_CENTBINS] = {10, 30, 50};
    const TString cent_label[N_CENTBINS] = {"cent0to10", "cent10to30", "cent30to50"};
*/
    const int N_CENTBINS = 5;
    const int min_centbin[N_CENTBINS] = {0, 10, 20, 30, 40}; // in 1% units
    const int max_centbin[N_CENTBINS] = {10, 20, 30, 40, 50};
    const TString cent_label[N_CENTBINS] = {"cent0to10","cent10to20", "cent20to30", "cent30to40", "cent40to50"};

    // Output file
    TString outfile = TString::Format("%s/ROOT/Resolution_incl_%d_%d.root", output_path.Data(), istart, iend);

    // --- 2. HISTOGRAM BOOKING (one set per broad cent bin) ---
    TH1D *hist_Q2Q2_HFmHFp_Re[N_CENTBINS];
    TH1D *hist_Q2Q2_HFmTrk_Re[N_CENTBINS];
    TH1D *hist_Q2Q2_HFpTrk_Re[N_CENTBINS];
    TH1D *hist_Q3Q3_HFmHFp_Re[N_CENTBINS];
    TH1D *hist_Q3Q3_HFmTrk_Re[N_CENTBINS];
    TH1D *hist_Q3Q3_HFpTrk_Re[N_CENTBINS];

    TH1D *aux_Q2 = new TH1D("aux_Q2", "", 6000, -30000, 30000);
    TH1D *aux_Q3 = new TH1D("aux_Q3", "", 3000, -15000, 15000);

    TString hname;
    for (Int_t ib = 0; ib < N_CENTBINS; ib++)
    {
        hname = "Q2Q2_HFmHFp_Re_" + cent_label[ib];
        hist_Q2Q2_HFmHFp_Re[ib] = (TH1D *)aux_Q2->Clone(hname);
        hname = "Q2Q2_HFmTrk_Re_" + cent_label[ib];
        hist_Q2Q2_HFmTrk_Re[ib] = (TH1D *)aux_Q2->Clone(hname);
        hname = "Q2Q2_HFpTrk_Re_" + cent_label[ib];
        hist_Q2Q2_HFpTrk_Re[ib] = (TH1D *)aux_Q2->Clone(hname);

        hname = "Q3Q3_HFmHFp_Re_" + cent_label[ib];
        hist_Q3Q3_HFmHFp_Re[ib] = (TH1D *)aux_Q3->Clone(hname);
        hname = "Q3Q3_HFmTrk_Re_" + cent_label[ib];
        hist_Q3Q3_HFmTrk_Re[ib] = (TH1D *)aux_Q3->Clone(hname);
        hname = "Q3Q3_HFpTrk_Re_" + cent_label[ib];
        hist_Q3Q3_HFpTrk_Re[ib] = (TH1D *)aux_Q3->Clone(hname);
    }
    delete aux_Q2;
    delete aux_Q3;

    // --- 3. FILE LOOP ---
    ifstream file_stream(input_txt.Data());
    string filename;
    int ifile = 0;

    while (file_stream >> filename)
    {
        if (ifile < istart)
        {
            ifile++;
            continue;
        }
        if (ifile >= iend)
            break;

        TFile *fin = TFile::Open(filename.c_str());
        if (!fin || fin->IsZombie())
        {
            std::cout << "Warning: Skipping bad file: " << filename << std::endl;
            if (fin)
            {
                fin->Close();
                delete fin;
            }
            ifile++;
            continue;
        }

        std::cout << ">>> Resolution (incl) Pass: Processing file " << ifile << " : " << filename << std::endl;

        TTree *tree = (TTree *)fin->Get("Ana/ntEvtInfo");

        Int_t centrality;
        Float_t ephfmQ[3], ephfpQ[3], ephfmSumW[3], ephfpSumW[3];
        Float_t ephfpAngle[3], ephfmAngle[3], eptkAngle[2], eptkQ[2];

        tree->SetBranchAddress("centrality", &centrality);
        tree->SetBranchAddress("ephfpAngle", ephfpAngle);
        tree->SetBranchAddress("ephfmAngle", ephfmAngle);
        tree->SetBranchAddress("ephfmQ", ephfmQ);
        tree->SetBranchAddress("ephfpQ", ephfpQ);
        tree->SetBranchAddress("ephfmSumW", ephfmSumW);
        tree->SetBranchAddress("ephfpSumW", ephfpSumW);
        tree->SetBranchAddress("eptkAngle", eptkAngle);
        tree->SetBranchAddress("eptkQ", eptkQ);

        tree->SetBranchStatus("*", 0);
        for (const auto &p : {"centrality", "ephfpAngle", "ephfmAngle", "ephfpQ", "ephfmQ",
                              "eptkAngle", "eptkQ", "ephfmSumW", "ephfpSumW"})
            tree->SetBranchStatus(p, 1);

        // --- 4. EVENT LOOP ---
        Int_t n_entries = tree->GetEntries();
        for (Long64_t ii = 0; ii < n_entries; ii++)
        {
            tree->GetEntry(ii);

            if (ii % 100000 == 0)
                printf("Processing entry %lld of %lld : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);

            int cent = centrality / 2;
            if (cent < 0 || cent >= max_centbin[N_CENTBINS - 1])
                continue;

            TComplex aux_Q2_HFm(ephfmQ[1] * TMath::Cos(2.0 * ephfmAngle[1]), ephfmQ[1] * TMath::Sin(2.0 * ephfmAngle[1]), 0);
            TComplex aux_Q2_HFp(ephfpQ[1] * TMath::Cos(2.0 * ephfpAngle[1]), ephfpQ[1] * TMath::Sin(2.0 * ephfpAngle[1]), 0);
            TComplex aux_Q2_Trk(eptkQ[0] * TMath::Cos(2.0 * eptkAngle[0]), eptkQ[0] * TMath::Sin(2.0 * eptkAngle[0]), 0);

            TComplex aux_Q3_HFm(ephfmQ[2] * TMath::Cos(3.0 * ephfmAngle[2]), ephfmQ[2] * TMath::Sin(3.0 * ephfmAngle[2]), 0);
            TComplex aux_Q3_HFp(ephfpQ[2] * TMath::Cos(3.0 * ephfpAngle[2]), ephfpQ[2] * TMath::Sin(3.0 * ephfpAngle[2]), 0);
            TComplex aux_Q3_Trk(eptkQ[1] * TMath::Cos(3.0 * eptkAngle[1]), eptkQ[1] * TMath::Sin(3.0 * eptkAngle[1]), 0);

            for (int i_cen = 0; i_cen < N_CENTBINS; i_cen++)
            {
                if (cent >= min_centbin[i_cen] && cent < max_centbin[i_cen])
                {
                    // Fill per-broad-cent histograms
                    hist_Q2Q2_HFmHFp_Re[i_cen]->Fill((aux_Q2_HFm * TComplex::Conjugate(aux_Q2_HFp)).Re());
                    hist_Q2Q2_HFmTrk_Re[i_cen]->Fill((aux_Q2_HFm * TComplex::Conjugate(aux_Q2_Trk)).Re());
                    hist_Q2Q2_HFpTrk_Re[i_cen]->Fill((aux_Q2_HFp * TComplex::Conjugate(aux_Q2_Trk)).Re());

                    hist_Q3Q3_HFmHFp_Re[i_cen]->Fill((aux_Q3_HFm * TComplex::Conjugate(aux_Q3_HFp)).Re());
                    hist_Q3Q3_HFmTrk_Re[i_cen]->Fill((aux_Q3_HFm * TComplex::Conjugate(aux_Q3_Trk)).Re());
                    hist_Q3Q3_HFpTrk_Re[i_cen]->Fill((aux_Q3_HFp * TComplex::Conjugate(aux_Q3_Trk)).Re());
                }
            }
        } // -- Event loop --

        fin->Close();
        delete fin;
        ifile++;
    } // -- File loop --

    // --- 5. SAVE ---
    TFile *fout = new TFile(outfile, "RECREATE");
    fout->cd();

    for (Int_t ib = 0; ib < N_CENTBINS; ib++)
    {
        hist_Q2Q2_HFmHFp_Re[ib]->Write();
        hist_Q2Q2_HFmTrk_Re[ib]->Write();
        hist_Q2Q2_HFpTrk_Re[ib]->Write();
        hist_Q3Q3_HFmHFp_Re[ib]->Write();
        hist_Q3Q3_HFmTrk_Re[ib]->Write();
        hist_Q3Q3_HFpTrk_Re[ib]->Write();
    }

    fout->Write(0, TObject::kOverwrite);
    fout->Close();
    std::cout << ">>> Inclusive Resolution saved to: " << outfile << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc == 5)
    {
        TString input_txt = argv[1];
        TString output_path = argv[2];
        int istart = std::stoi(argv[3]);
        int iend = std::stoi(argv[4]);

        Calculate_Resolution_incl(input_txt, output_path, istart, iend);
    }
    else
    {
        std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
        return 1;
    }
    return 0;
}
