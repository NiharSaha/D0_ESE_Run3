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

// Include your auto-generated quantile cuts header
//#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/ESE_Cuts_MB11to21_Jan22.h" //For 1/3 stat
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Quantiles/quantile_cuts_2023_MB0to31.h" //For 99% stat

using namespace std;

void Calculate_Resolution(TString input_txt, TString output_path, int istart, int iend) {
    
    // --- 1. CONFIGURATION ---
    const int N_CENTBIN = 4;
    //Int_t min_centbin[N_CENTBIN] = {0, 10, 30, 50};
    //Int_t max_centbin[N_CENTBIN] = {10, 30, 50, 90};
    const int N_CENTBIN_1 = 90;
    const int N_QBINS = 10;

    // Output file setup
    TString outfile = TString::Format("%s/ROOT/Resolution_%d_%d.root", output_path.Data(), istart, iend);

    TH1D *hist_Q2Q2_HFmHFp_Re_centbin[N_CENTBIN_1][N_QBINS]; // centrality
    TH1D *hist_Q2Q2_HFmHFp_Im_centbin[N_CENTBIN_1][N_QBINS]; // centrality
    TH1D *aux_hist_Q2Q2_HFmHFp = new TH1D("aux_hist_Q2Q2_HFmHFp", "", 6000, -30000, 30000);
    
    TH1D *hist_Q2Q2_HFmTrk_Re_centbin[N_CENTBIN_1][N_QBINS]; // centrality
    TH1D *hist_Q2Q2_HFmTrk_Im_centbin[N_CENTBIN_1][N_QBINS]; // centrality
    TH1D *aux_hist_Q2Q2_HFmTrk = new TH1D("aux_hist_Q2Q2_HFmTrk", "", 6000, -30000, 30000);
    
    TH1D *hist_Q2Q2_HFpTrk_Re_centbin[N_CENTBIN_1][N_QBINS]; // centrality
    TH1D *hist_Q2Q2_HFpTrk_Im_centbin[N_CENTBIN_1][N_QBINS]; // centrality
    TH1D *aux_hist_Q2Q2_HFpTrk = new TH1D("aux_hist_Q2Q2_HFpTrk", "", 6000, -30000, 30000);
    
    TH1D *hist_Q3Q3_HFmHFp_Re_centbin[N_CENTBIN_1][N_QBINS]; // centrality
    TH1D *hist_Q3Q3_HFmHFp_Im_centbin[N_CENTBIN_1][N_QBINS]; // centrality
    TH1D *aux_hist_Q3Q3_HFmHFp = new TH1D("aux_hist_Q3Q3_HFmHFp", "", 3000, -15000, 15000);
    
    TH1D *hist_Q3Q3_HFmTrk_Re_centbin[N_CENTBIN_1][N_QBINS]; // centrality
    TH1D *hist_Q3Q3_HFmTrk_Im_centbin[N_CENTBIN_1][N_QBINS]; // centrality
    TH1D *aux_hist_Q3Q3_HFmTrk = new TH1D("aux_hist_Q3Q3_HFmTrk", "", 3000, -15000, 15000);
    
    TH1D *hist_Q3Q3_HFpTrk_Re_centbin[N_CENTBIN_1][N_QBINS]; // centrality
    TH1D *hist_Q3Q3_HFpTrk_Im_centbin[N_CENTBIN_1][N_QBINS]; // centrality
    TH1D *aux_hist_Q3Q3_HFpTrk = new TH1D("aux_hist_Q3Q3_HFpTrk", "", 3000, -15000, 15000);
    

    
    TString hname;
    
    for (Int_t icent = 0; icent < N_CENTBIN_1; icent++) {
      for(Int_t iq=0; iq<N_QBINS; iq++){
	// Q2Q2 Histograms
	hname = "Q2Q2_HFmHFp_Re_cen"+to_string(icent)+"_"+to_string(icent+1)+"_q2bin_"+to_string(iq);
	hist_Q2Q2_HFmHFp_Re_centbin[icent][iq] = (TH1D *)aux_hist_Q2Q2_HFmHFp->Clone(hname);
	hname = "Q2Q2_HFmTrk_Re_cen"+to_string(icent)+"_"+to_string(icent+1)+"_q2bin_"+to_string(iq);
	hist_Q2Q2_HFmTrk_Re_centbin[icent][iq] = (TH1D *)aux_hist_Q2Q2_HFmTrk->Clone(hname);
	hname = "Q2Q2_HFpTrk_Re_cen"+to_string(icent)+"_"+to_string(icent+1)+"_q2bin_"+to_string(iq);
	hist_Q2Q2_HFpTrk_Re_centbin[icent][iq] = (TH1D *)aux_hist_Q2Q2_HFpTrk->Clone(hname);



	hname = "Q3Q3_HFmHFp_Re_cen"+to_string(icent)+"_"+to_string(icent+1)+"_q3bin_"+to_string(iq);
	hist_Q3Q3_HFmHFp_Re_centbin[icent][iq] = (TH1D *)aux_hist_Q3Q3_HFmHFp->Clone(hname);
	hname = "Q3Q3_HFmTrk_Re_cen"+to_string(icent)+"_"+to_string(icent+1)+"_q3bin_"+to_string(iq);
	hist_Q3Q3_HFmTrk_Re_centbin[icent][iq] = (TH1D *)aux_hist_Q3Q3_HFmTrk->Clone(hname);
	hname = "Q3Q3_HFpTrk_Re_cen"+to_string(icent)+"_"+to_string(icent+1)+"_q3bin_"+to_string(iq);
	hist_Q3Q3_HFpTrk_Re_centbin[icent][iq] = (TH1D *)aux_hist_Q3Q3_HFpTrk->Clone(hname);

	
      }
    }

    // --- 2. FILE LOOP ---
    ifstream file_stream(input_txt.Data());
    string filename;
    int ifile = 0;

    while (file_stream >> filename) {
        if (ifile < istart) { ifile++; continue; }
        if (ifile >= iend) break;

        TFile *fin = TFile::Open(filename.c_str());
        if(!fin || fin->IsZombie()) {
            std::cout << "Warning: Skipping bad file: " << filename << std::endl;
            if (fin) { fin->Close(); delete fin; }
            ifile++; continue;
        }

        std::cout << ">>> Resolution Pass: Processing file " << ifile << " : " << filename << std::endl;

        TTree *tree = (TTree*)fin->Get("d0Analyzer/VCNtuple_D02kpi");
        TTree *t_eventinfoana = (TTree*)fin->Get("eventinfoana/EventInfoNtuple");
        tree->AddFriend(t_eventinfoana);

	Float_t q2_hf_total, q3_hf_total;
        Int_t centrality;
        Float_t ephfmQ[3], ephfpQ[3], ephfmSumW[3], ephfpSumW[3], ephfpAngle[3], ephfmAngle[3], eptkAngle[2], eptkQ[2];

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
        for (const auto& p : {"centrality", "ephfpAngle", "ephfmAngle", "ephfpQ", "ephfmQ", "eptkAngle", "eptkQ", "ephfmSumW", "ephfpSumW"})
            tree->SetBranchStatus(p, 1);

        // --- 4. EVENT LOOP ---
        Int_t n_entries = tree->GetEntries();
        for (Long64_t ii = 0; ii < n_entries; ii++) {
            tree->GetEntry(ii);

            int cent = centrality / 2;
            if (cent < 0 || cent >= 90) continue;


	    TComplex aux_Q2_HFm(ephfmQ[1]*TMath::Cos(2.0*ephfmAngle[1]),ephfmQ[1]*TMath::Sin(2.0*ephfmAngle[1]),0);
	    TComplex aux_Q2_HFp(ephfpQ[1]*TMath::Cos(2.0*ephfpAngle[1]),ephfpQ[1]*TMath::Sin(2.0*ephfpAngle[1]),0);
	    TComplex aux_Q2_Trk(eptkQ[0]*TMath::Cos(2.0*eptkAngle[0]),eptkQ[0]*TMath::Sin(2.0*eptkAngle[0]),0);
	    
	    TComplex aux_Q3_HFm(ephfmQ[2]*TMath::Cos(3.0*ephfmAngle[2]), ephfmQ[2]*TMath::Sin(3.0*ephfmAngle[2]), 0);
	    TComplex aux_Q3_HFp(ephfpQ[2]*TMath::Cos(3.0*ephfpAngle[2]), ephfpQ[2]*TMath::Sin(3.0*ephfpAngle[2]), 0);
	    TComplex aux_Q3_Trk(eptkQ[1]*TMath::Cos(3.0*eptkAngle[1]), eptkQ[1]*TMath::Sin(3.0*eptkAngle[1]), 0);

	    if ((ephfpSumW[1] + ephfmSumW[1]) <= 0 || (ephfpSumW[2] + ephfmSumW[2]) <= 0) continue;

	    q2_hf_total = (ephfpQ[1] + ephfmQ[1]) / (ephfpSumW[1] + ephfmSumW[1]);
	    q3_hf_total = (ephfpQ[2] + ephfmQ[2]) / (ephfpSumW[2] + ephfmSumW[2]);

	    int q2_bin = TMath::BinarySearch(11, q2_cuts[cent], (double)q2_hf_total);
	    if (q2_bin >= 0 && q2_bin < 10) {
	        hist_Q2Q2_HFmHFp_Re_centbin[cent][q2_bin]->Fill((aux_Q2_HFm * TComplex::Conjugate(aux_Q2_HFp)).Re());
	        hist_Q2Q2_HFmTrk_Re_centbin[cent][q2_bin]->Fill((aux_Q2_HFm * TComplex::Conjugate(aux_Q2_Trk)).Re());
	        hist_Q2Q2_HFpTrk_Re_centbin[cent][q2_bin]->Fill((aux_Q2_HFp * TComplex::Conjugate(aux_Q2_Trk)).Re());
	    }

	    int q3_bin = TMath::BinarySearch(11, q3_cuts[cent], (double)q3_hf_total);
	    if (q3_bin >= 0 && q3_bin < 10) {
	        hist_Q3Q3_HFmHFp_Re_centbin[cent][q3_bin]->Fill((aux_Q3_HFm * TComplex::Conjugate(aux_Q3_HFp)).Re());
	        hist_Q3Q3_HFmTrk_Re_centbin[cent][q3_bin]->Fill((aux_Q3_HFm * TComplex::Conjugate(aux_Q3_Trk)).Re());
	        hist_Q3Q3_HFpTrk_Re_centbin[cent][q3_bin]->Fill((aux_Q3_HFp * TComplex::Conjugate(aux_Q3_Trk)).Re());
	    }
	}// -- Event loop --
	
	fin->Close(); delete fin;
	ifile++;
    }
    
    // --- 5. SAVE ---
    TFile *fout = new TFile(outfile, "RECREATE");

    fout->cd();
    for(Int_t i_cen = 0; i_cen < N_CENTBIN_1; i_cen++)
      {
	for(Int_t i_q = 0; i_q < N_QBINS; i_q++) // Assuming 10 quantiles for both
	  {
	    // --- Write Q2 Histograms ---
	    hist_Q2Q2_HFmHFp_Re_centbin[i_cen][i_q]->Write();
	    hist_Q2Q2_HFmTrk_Re_centbin[i_cen][i_q]->Write();
	    hist_Q2Q2_HFpTrk_Re_centbin[i_cen][i_q]->Write();
	    
	    // --- Write Q3 Histograms ---
	    hist_Q3Q3_HFmHFp_Re_centbin[i_cen][i_q]->Write();
	    hist_Q3Q3_HFmTrk_Re_centbin[i_cen][i_q]->Write();
	    hist_Q3Q3_HFpTrk_Re_centbin[i_cen][i_q]->Write();
	  }
      }

    fout->Write(0,TObject::kOverwrite);
    fout->Close();
    std::cout << ">>> Resolution Ingredients saved to: " << outfile << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc == 5) // Expecting 4 arguments: input_txt, output_path, istart, iend

    {
        TString input_txt = argv[1];
        TString output_path = argv[2];
        int istart = std::stoi(argv[3]);
        int iend = std::stoi(argv[4]);

        Calculate_Resolution(input_txt, output_path, istart, iend);
    }
    else
    {
        std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
        return 1;
    }
    return 0;
}
