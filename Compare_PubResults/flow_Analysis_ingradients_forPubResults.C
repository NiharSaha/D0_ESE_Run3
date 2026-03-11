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
#include "TNtuple.h"
#include <climits>
#include <iomanip>


#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/BDTHandler.h"
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/BDTHandler.cc"
//#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/ESE_Cuts_MB11to21_Jan22.h"


using namespace std;



void flow_Analysis_ingradients_forPubResults(TString input_txt, TString output_path, int istart, int iend) {


  BDTHandler bdtHandler;

  //const int N_QBINS = 10;  
  //const int N_CENTBINS_1 = 90;

  
  const int N_CENTBIN = 3;
  Int_t min_centbin[N_CENTBIN] = {0, 10, 30};
  Int_t max_centbin[N_CENTBIN] = {10, 30, 50};

  const double MAX_PT_ANA = 100.0;
  const double MIN_PT_ANA = 1.0;
  const double MAX_Y_ANA = 2.4;
  
  auto file_res = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Resolution_MB11to21_Mar2_PubResults_v4/ROOT/Resolution_out_combined.root");
    
  
  Double_t v2_den_Dy_plus[N_CENTBIN];
  Double_t v2_den_Dy_minus[N_CENTBIN];
  Double_t v3_den_Dy_plus[N_CENTBIN];
  Double_t v3_den_Dy_minus[N_CENTBIN];
  
  Float_t q2_hf_total, q3_hf_total;

  
  TH1D *hQ2Q2_HFmHFp[N_CENTBIN];
  TH1D *hQ2Q2_HFmTrk[N_CENTBIN];
  TH1D *hQ2Q2_HFpTrk[N_CENTBIN];
  TH1D *hQ3Q3_HFmHFp[N_CENTBIN];
  TH1D *hQ3Q3_HFmTrk[N_CENTBIN];
  TH1D *hQ3Q3_HFpTrk[N_CENTBIN];
  TString hname;

  
  for(Int_t i=0; i<N_CENTBIN; i++){
    hname = "Q2Q2_HFmHFp_Re_cen"+to_string(i)+"_"+to_string(i+1);
    hQ2Q2_HFmHFp[i] = (TH1D*)file_res->Get(hname);
    hname = "Q2Q2_HFmTrk_Re_cen"+to_string(i)+"_"+to_string(i+1);
    hQ2Q2_HFmTrk[i] = (TH1D*)file_res->Get(hname);
    hname = "Q2Q2_HFpTrk_Re_cen"+to_string(i)+"_"+to_string(i+1);
    hQ2Q2_HFpTrk[i] = (TH1D*)file_res->Get(hname);

    hname = "Q3Q3_HFmHFp_Re_cen"+to_string(i)+"_"+to_string(i+1);
    hQ3Q3_HFmHFp[i] = (TH1D*)file_res->Get(hname);
    hname = "Q3Q3_HFmTrk_Re_cen"+to_string(i)+"_"+to_string(i+1);
    hQ3Q3_HFmTrk[i] = (TH1D*)file_res->Get(hname);
    hname = "Q3Q3_HFpTrk_Re_cen"+to_string(i)+"_"+to_string(i+1);
    hQ3Q3_HFpTrk[i] = (TH1D*)file_res->Get(hname);

    v2_den_Dy_plus[i]  = TMath::Sqrt( (hQ2Q2_HFmHFp[i]->GetMean() * hQ2Q2_HFmTrk[i]->GetMean()) / hQ2Q2_HFpTrk[i]->GetMean() );
    v2_den_Dy_minus[i] = TMath::Sqrt( (hQ2Q2_HFmHFp[i]->GetMean() * hQ2Q2_HFpTrk[i]->GetMean()) / hQ2Q2_HFmTrk[i]->GetMean() );
    v3_den_Dy_plus[i]  = TMath::Sqrt( (hQ3Q3_HFmHFp[i]->GetMean() * hQ3Q3_HFmTrk[i]->GetMean()) / hQ3Q3_HFpTrk[i]->GetMean() );
    v3_den_Dy_minus[i] = TMath::Sqrt( (hQ3Q3_HFmHFp[i]->GetMean() * hQ3Q3_HFpTrk[i]->GetMean()) / hQ3Q3_HFmTrk[i]->GetMean() );
  }
  file_res->Close();

 
        
    ifstream file_stream(input_txt.Data());
    TString outfile = TString::Format("%s/ROOT/flow_Analysis_%d_%d.root", output_path.Data(), istart, iend);
    TFile *fout = new TFile(outfile, "RECREATE");
    auto nt = new TNtuple("nt","nt","cent:q2_hf_total:q3_hf_total:pT:y:dca:mass:v2:v3");

  
    Long64_t n_evts[N_CENTBIN] = {};

    string filename;
    int ifile = 0;

    while (file_stream >> filename) {
        if (ifile < istart) { ifile++; continue; }
        if (ifile >= iend) break;

        TFile *fin = TFile::Open(filename.c_str());

	if(!fin || fin->IsZombie()) {
	  std::cout << "Warning: Skipping bad file: " << filename << std::endl;
	  if (fin) { fin->Close(); delete fin; }
	  ifile ++;
	  continue;
	}


	std::cout << ">>> Processing ifile=" << ifile << " : " << filename << std::endl;
    

	TTree *tree = (TTree*)fin->Get("d0Analyzer/VCNtuple_D02kpi");
        TTree *t_eventinfoana = (TTree*)fin->Get("eventinfoana/EventInfoNtuple");
	tree->AddFriend(t_eventinfoana);
	
        // Standard Variable setup
	Int_t centrality;
        Int_t candSize;
        Float_t ephfmQ[3], ephfpQ[3], ephfmSumW[3], ephfpSumW[3], ephfpAngle[3], ephfmAngle[3], eptkAngle[2], eptkQ[2];
        vector<float> *pT=0, *phi=0, *mass=0, *y=0, *eta=0, *mva=0, *dca=0;

	
        tree->SetBranchAddress("candSize", &candSize);
        tree->SetBranchAddress("centrality", &centrality);
	tree->SetBranchAddress("ephfpAngle", ephfpAngle);
	tree->SetBranchAddress("ephfmAngle", ephfmAngle);
        tree->SetBranchAddress("ephfmQ", ephfmQ);
	tree->SetBranchAddress("ephfpQ", ephfpQ);
        tree->SetBranchAddress("ephfmSumW", ephfmSumW);
	tree->SetBranchAddress("ephfpSumW", ephfpSumW);
	tree->SetBranchAddress("eptkAngle", eptkAngle);
	tree->SetBranchAddress("eptkQ", eptkQ);
	tree->SetBranchAddress("pT", &pT);
	tree->SetBranchAddress("phi", &phi);
	tree->SetBranchAddress("mass", &mass);
	tree->SetBranchAddress("y", &y);
	tree->SetBranchAddress("eta", &eta);
	tree->SetBranchAddress("mva", &mva);
	tree->SetBranchAddress("twoTrackDCA", &dca);


	tree->SetBranchStatus("*", 0);
	for (const auto& p : {"candSize", "centrality", "ephfpAngle", "ephfmAngle", "ephfpQ", "ephfmQ", "eptkAngle", "eptkQ", "ephfmSumW", "ephfpSumW", "y", "eta", "pT", "phi", "mass", "mva", "twoTrackDCA"})
	  tree->SetBranchStatus(p, 1);


	
	Int_t n_entries = tree->GetEntries();
	
	for (Long64_t ii = 0; ii < n_entries; ii++) {
            tree->GetEntry(ii);

        if(ii % 10000 == 0)
          printf("Processing entry %lld of %lld : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);
        
            Int_t cent = centrality/2;
            // --- CHANGE: skip events outside coarse cent range ---
            if (cent < 0 || cent >= max_centbin[N_CENTBIN-1]) continue;

	      TComplex aux_Q2_HFm(ephfmQ[1]*TMath::Cos(2.0*ephfmAngle[1]),ephfmQ[1]*TMath::Sin(2.0*ephfmAngle[1]),0);
	      TComplex aux_Q2_HFp(ephfpQ[1]*TMath::Cos(2.0*ephfpAngle[1]),ephfpQ[1]*TMath::Sin(2.0*ephfpAngle[1]),0);
	      TComplex aux_Q2_Trk(eptkQ[0]*TMath::Cos(2.0*eptkAngle[0]),eptkQ[0]*TMath::Sin(2.0*eptkAngle[0]),0);
	      
	      TComplex aux_Q3_HFm(ephfmQ[2]*TMath::Cos(3.0*ephfmAngle[2]), ephfmQ[2]*TMath::Sin(3.0*ephfmAngle[2]), 0);
	      TComplex aux_Q3_HFp(ephfpQ[2]*TMath::Cos(3.0*ephfpAngle[2]), ephfpQ[2]*TMath::Sin(3.0*ephfpAngle[2]), 0);
	      TComplex aux_Q3_Trk(eptkQ[1]*TMath::Cos(3.0*eptkAngle[1]), eptkQ[1]*TMath::Sin(3.0*eptkAngle[1]), 0);
	      

	      float sumW2 = ephfpSumW[1] + ephfmSumW[1];
	      float sumW3 = ephfpSumW[2] + ephfmSumW[2];

	      if (sumW2 <= 0 || sumW3 <= 0) continue;
	      q2_hf_total = (ephfpQ[1] + ephfmQ[1]) / sumW2;
	      q3_hf_total = (ephfpQ[2] + ephfmQ[2]) / sumW3;
	

	      //int q2_idx = -1;
	      //int q3_idx = -1;
	      int i_cen_match = -1;
	      

	      // --- CHANGE: Find matching coarse cent bin (same logic as Resolution code) ---
	      for(int ic = 0; ic < N_CENTBIN; ic++){
	      if(cent >= min_centbin[ic] && cent < max_centbin[ic]){
	        i_cen_match = ic;
	        break;
	      }
	      }
	      if(i_cen_match < 0) continue;

	      n_evts[i_cen_match]++;

	      // --- CHANGE: Resolution lookup uses only i_cen_match (no q-bin index) ---
	      Double_t res2_plus  = v2_den_Dy_plus[i_cen_match];
	      Double_t res2_minus = v2_den_Dy_minus[i_cen_match];
	      Double_t res3_plus  = v3_den_Dy_plus[i_cen_match];
	      Double_t res3_minus = v3_den_Dy_minus[i_cen_match];
	      
	      
	      for (int icand = 0; icand < candSize; ++icand) {


	      float Y = y->at(icand);	      
	      float Pt   = pT->at(icand);
	      float Mass = mass->at(icand);
	      float Phi  = phi->at(icand);
	      float Mva  = mva->at(icand);
	      float Dca  = dca->at(icand);
	      

	      if(Pt >= MAX_PT_ANA || Pt < MIN_PT_ANA  || fabs(Y) >= MAX_Y_ANA) continue;


	      //I need to use 2*cent because of format
	      double bdt_cut = bdtHandler.getBDTCut(Y, 2*cent, Pt);         
	      if (Mva <= bdt_cut) continue;


	      TComplex aux_Q2_D0(TMath::Cos(2.0*Phi), TMath::Sin(2.0*Phi), 0);
	      TComplex aux_Q3_D0(TMath::Cos(3.0*Phi), TMath::Sin(3.0*Phi), 0);
	      Double_t v2_obs, v3_obs;
	      Double_t res2_final, res3_final;
    
	      if(Y > 0) {
		v2_obs = (aux_Q2_D0 * TComplex::Conjugate(aux_Q2_HFm)).Re();
		v3_obs = (aux_Q3_D0 * TComplex::Conjugate(aux_Q3_HFm)).Re();
		res2_final = res2_plus;
		res3_final = res3_plus;
	      } else {
		v2_obs = (aux_Q2_D0 * TComplex::Conjugate(aux_Q2_HFp)).Re();
		v3_obs = (aux_Q3_D0 * TComplex::Conjugate(aux_Q3_HFp)).Re();
		res2_final = res2_minus;
		res3_final = res3_minus;
	      }

	      if(res2_final > 0 && res3_final > 0) {
		nt->Fill(cent, q2_hf_total, q3_hf_total, Pt, Y, Dca, Mass, v2_obs/res2_final, v3_obs/res3_final);
	      }
	      
	      
	      }// -- cand loop --
	}// -- Event loop ---
	
	
        fin->Close();
	delete fin;
        ifile++;
    }// -- file loop ---
    
    
    //============================
    // Write all histograms!!!
    //============================


    fout->cd();
    nt->Write();

    fout->cd();
    fout->Write(0,TObject::kOverwrite);
    fout->Close();
}


int main(int argc, char *argv[])
{
    if (argc == 5) // Expecting 4 arguments: input_txt, output_path, istart, iend                                                          
                                                                                                                                           
    {
        TString input_txt = argv[1];
        TString output_path = argv[2];
        int istart = std::stoi(argv[3]);
        int iend = std::stoi(argv[4]);
	
        flow_Analysis_ingradients_forPubResults(input_txt, output_path, istart, iend);
    }
    else
    {
        std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
        return 1;
    }
    return 0;
}
