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

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/ESE_Cuts_MB11to21_Jan22.h"
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/ESE_Resolution_Factors_Jan22_v2.h"

using namespace std;

// ======================== To get BDT cuts ========================

const TString BDT_CSV_PATH = "/home/saha115/D0_ESE/CMSSW_13_2_11/src/bdt_cuts.csv";
double bdt_cuts[4][9];

// Function to load optimized BDT cuts from CSV
void LoadBDTCuts() {
  ifstream file(BDT_CSV_PATH.Data());
  if (!file.is_open()) {
    cout << "!!! ERROR: Could not open BDT CSV at: " << BDT_CSV_PATH << endl;
    return;
  }
  
  string line;
  int y, c, p;
  double cut;
  while (getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;
    stringstream ss(line);
    if (!(ss >> y >> c >> p >> cut)) continue;
    
    // Map y=0 values to our internal table
    if (y == 0 && c < 4 && p < 9) {
      bdt_cuts[c][p] = cut;
    }
  }
  file.close();
  cout << ">>> Successfully loaded BDT cuts from: " << BDT_CSV_PATH << endl;
}
    
// ======================================================================





void flow_Analysis(TString input_txt, TString output_path, int istart, int iend) {

  LoadBDTCuts();


  const int N_Q_BINS = 10;  
  
  const int N_CENTBIN = 4;
  Int_t min_centbin[N_CENTBIN]={0, 10, 30, 50};
  Int_t max_centbin[N_CENTBIN]={10, 30, 50, 90};
  TString label_centbin[N_CENTBIN]={"cent0to10", "cent10to30", "cent30to50", "cent50to90"};
  TString centbin[N_CENTBIN]={"0 <= cent < 10 %", "10 <= cent < 30 %", "30 <= cent < 50 %", "50 <= cent < 90 %"};
  
  const int N_PTBIN = 10;
  Float_t min_pTbin[N_PTBIN]={1.0, 2.0, 3.0, 4.0, 5.0, 6.0,8.0,10.0,15.0, 20.0};
  Float_t max_pTbin[N_PTBIN]={2.0, 3.0, 4.0, 5.0, 6.0,8.0,10.0,15.0,20.0, 40.0};
  TString label_pTbin[N_PTBIN]={"pT1to2", "pT2to3", "pT3to4","pT4to5", "pT5to6","pT6to8","pT8to10","pT10to15","pT15to20", "pT20to40"};
  TString pTbin[N_PTBIN]={"1.0 < pT < 2.0 GeV", "2.0 < pT < 3.0 GeV", "3.0 < pT < 4.0 GeV", "4.0 < pT < 5.0 GeV", "5.0 < pT < 6.0 GeV","6.0 < pT < 8.0","8.0 < pT < 10.0 GeV","10.0 < pT < 15.0 GeV","15.0 < pT < 20.0 GeV", "20.0 < pT < 40.0 GeV"};

  
  // 2. DEFINE SPARSE (Exactly your structure)
  const int nDim = 6;
  Int_t bins6D[nDim] = {60, 8, 90, 10, 2, 200}; 
  Double_t mins6D[nDim] = {1.7, 0, 0, 0, 0, -10.0};
  Double_t maxs6D[nDim] = {2.0, 8.0, 90.0, 10.0, 2.0, 10.0};
  THnSparseD *hn_v2 = new THnSparseD("hn_v2", "D0 v2 ESE;Mass;pT;Cent;q2;ObsType;SP_Value", 6, bins6D, mins6D, maxs6D);
  THnSparseD *hn_v3 = new THnSparseD("hn_v3", "D0 v3 ESE;Mass;pT;Cent;q3;ObsType;SP_Value", 6, bins6D, mins6D, maxs6D);
  hn_v2->Sumw2();
  hn_v3->Sumw2();
  


        
    ifstream file_stream(input_txt.Data());
    TString outfile = TString::Format("%s/ROOT/flow_Analysis_%d_%d.root", output_path.Data(), istart, iend);
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
        vector<float> *pT=0, *phi=0, *mass=0, *y=0, *eta=0, *mva=0;

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


	tree->SetBranchStatus("*", 0);
	for (const auto& p : {"candSize", "centrality", "ephfpAngle", "ephfmAngle", "ephfpQ", "ephfmQ", "eptkAngle", "eptkQ", "ephfmSumW", "ephfpSumW", "y", "eta", "pT", "phi", "mass", "mva"})
	  tree->SetBranchStatus(p, 1);

	Int_t n_entries = tree->GetEntries();
	
	for (Long64_t ii = 0; ii < n_entries; ii++) {
            tree->GetEntry(ii);

	    if(ii % 10000 == 0)
	      printf("Processing entry %lld of %lld : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);
	    
            Int_t cent = centrality/2;
            if (cent < 0 || cent >= 90) continue;

	      // Centrality binning for BDT matches min_centbin[4] = {0, 10, 30, 50}
	      int c_bdt = -1;
	      for (int i = 0; i < 4; i++) {
		if (cent >= min_centbin[i] && cent < max_centbin[i]) { c_bdt = i; break; }
	      }
	      if (c_bdt == -1) continue;	    
	    	    
	    Double_t Q2_total = ephfmQ[1] + ephfpQ[1];
	    Double_t Q3_total = ephfmQ[2] + ephfpQ[2];
	    Double_t Q2_weight = ephfmSumW[1] + ephfpSumW[1];
	    Double_t Q3_weight = ephfmSumW[2] + ephfpSumW[2];
	    
	    Double_t q2_val = (Q2_weight > 0) ? Q2_total/Q2_weight : -999;
	    Double_t q3_val = (Q3_weight > 0) ? Q3_total/Q3_weight : -999;
    
	    // 3. Determine Quantile Bins using Global Header Arrays
	    // TMath::BinarySearch returns the index k such that cuts[k] <= val < cuts[k+1]
	    Int_t k2_idx = TMath::BinarySearch(11, q2_cuts[cent], (Double_t)q2_val);
	    Int_t k3_idx = TMath::BinarySearch(11, q3_cuts[cent], (Double_t)q3_val);
	    
	    // Clamping to [0, N_Q_BINS-1] to handle boundary edges
	    if (k2_idx < 0) k2_idx = 0;
	    if (k2_idx >= 10) k2_idx = 9;

	    if (k3_idx < 0) k3_idx = 0;
	    if (k3_idx >= 10) k3_idx = 9;
	    
            int i_idx = -1;
	    for (int i = 0; i < N_CENTBIN; ++i) {
	      if (cent >= min_centbin[i] && cent < max_centbin[i]) {
		i_idx = i;
		break;
	      }
	    }
	    if (i_idx == -1) continue;


	    // Pre-calculate sub-event Q-vectors (A=HF+, B=HF-)
            TComplex QA_v2(ephfpQ[1]*cos(2*ephfpAngle[1]), ephfpQ[1]*sin(2*ephfpAngle[1]));
            TComplex QB_v2(ephfmQ[1]*cos(2*ephfmAngle[1]), ephfmQ[1]*sin(2*ephfmAngle[1]));
            TComplex QA_v3(ephfpQ[2]*cos(3*ephfpAngle[2]), ephfpQ[2]*sin(3*ephfpAngle[2]));
            TComplex QB_v3(ephfmQ[2]*cos(3*ephfmAngle[2]), ephfmQ[2]*sin(3*ephfmAngle[2]));



	    //events_binned_successfully++;
	    
	    for (int icand = 0; icand < candSize; ++icand) {

	      float Rap = y->at(icand);
	      if (TMath::Abs(Rap) > 1.0) continue;
	      
	      float Pt   = pT->at(icand);
	      float Mass = mass->at(icand);
	      float Phi  = phi->at(icand);
	      float Mva  = mva->at(icand);
	      
	      //Find the Analysis pT bin (0-9)
	      int p_idx = -1;
	      for (int p = 0; p < N_PTBIN; ++p) {
		if (Pt >= min_pTbin[p] && Pt < max_pTbin[p]) { p_idx = p; break; }
	      }
	      if (p_idx == -1) continue;


	      if (p_idx < 9 && Mva <= bdt_cuts[c_bdt][p_idx]) continue;
	      
	      int res_idx = c_bdt * 10 + 0; 	      
	      int q_map[2] = {k2_idx, k3_idx};
	      
	      for (int n = 0; n < 2; n++) {
		int harm = n + 2;
		int q_idx = q_map[n];
		int global_res_idx = c_bdt * 10 + q_idx;
		
		double Rn = 0;
		double y_bin = 0;
		TComplex Q_ref;
		
		// SP Condition: Use opposite sub-event
		if (Rap > 0) {
		  Rn = (n == 0) ? R2_yPlus[global_res_idx] : R3_yPlus[global_res_idx];
		  Q_ref = (n == 0) ? QB_v2 : QB_v3; // Use HF-
		  y_bin = 0.5;
		} else {
		  Rn = (n == 0) ? R2_yMinus[global_res_idx] : R3_yMinus[global_res_idx];
		  Q_ref = (n == 0) ? QA_v2 : QA_v3; // Use HF+
		  y_bin = 1.5;
		}
		
		if (Rn < 1e-5) continue;
		
		TComplex un(cos(harm * Phi), sin(harm * Phi));
		double sp_val = (un * TComplex::Conjugate(Q_ref)).Re() / Rn;
		
		Double_t x6D[6] = { (Double_t)Mass, (Double_t)p_idx + 0.5, (Double_t)cent + 0.5, 
				    (Double_t)q_idx + 0.5, y_bin, sp_val };
		
		if (n == 0) hn_v2->Fill(x6D);
		else hn_v3->Fill(x6D);
	      }
	      
	    }// -- cand loop --
	}// -- Event loop ---

	
        //fin->Close();
	//delete fin;
        ifile++;
    }// -- file loop ---
    
    
    //============================
    // Write all histograms!!!
    //============================

    TFile *fout = new TFile(outfile, "RECREATE");
    fout->mkdir("D0_vn_THn");
    fout->cd("D0_vn_THn");
    hn_v2->Write();
    hn_v3->Write();

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
	
        flow_Analysis(input_txt, output_path, istart, iend);
    }
    else
    {
        std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
        return 1;
    }
    return 0;
}
