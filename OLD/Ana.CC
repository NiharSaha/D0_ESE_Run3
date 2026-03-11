#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <chrono>
#include <ctime>
#include <cstdlib>

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


using namespace std::chrono;
using namespace std;



void Ana(){

  TH1::StatOverflows(kTRUE);
  TH2::StatOverflows(kTRUE);
  TH1::SetDefaultSumw2(); 
  TH2::SetDefaultSumw2();

  auto start = high_resolution_clock::now();
  // =================================================================
  // ### 1. Input and outfiles
  // =================================================================
  const TString input_file = "TTree_reProduction.txt";
  TFile * fout = new TFile("D0_ESE_Nov11.root","RECREATE");
  


  // =================================================================
  // ### 1. Define all histograms!
  // =================================================================

  const int N_CENTBIN = 4;
  Int_t min_centbin[N_CENTBIN]={0, 10, 30, 50};
  Int_t max_centbin[N_CENTBIN]={10, 30, 50, 90};
  Double_t centbinning[N_CENTBIN+1]={0., 10., 30., 50., 90.};
  TString label_centbin[N_CENTBIN]={"cent0to10", "cent10to30", "cent30to50", "cent50to90"};
  TString centbin[N_CENTBIN]={"0 <= cent < 10 %", "10 <= cent < 30 %", "30 <= cent < 50 %", "50 <= cent < 90 %"};


  const int N_CENTBIN_1PC = 90; 
  
  Int_t    min_centbin_1pc[N_CENTBIN_1PC];
  Int_t    max_centbin_1pc[N_CENTBIN_1PC];
  Double_t centbinning_1pc[N_CENTBIN_1PC + 1]; 
  TString  label_centbin_1pc[N_CENTBIN_1PC];
  TString  centbin_1pc[N_CENTBIN_1PC];
  
  for (int i = 0; i < N_CENTBIN_1PC; ++i) {
    min_centbin_1pc[i] = i;
    max_centbin_1pc[i] = i + 1;
    label_centbin_1pc[i] = TString::Format("cent%dto%d", i, i + 1);
    centbin_1pc[i]       = TString::Format("%d <= cent < %d %%", i, i + 1);
    centbinning_1pc[i] = (Double_t)i;
  }
  
  centbinning_1pc[N_CENTBIN_1PC] = (Double_t)N_CENTBIN_1PC;

  

  
  const int N_PTBIN = 6;
  Float_t min_pTbin[N_PTBIN]={4.0,6.0,8.0,10.0,15.0, 20.0};
  Float_t max_pTbin[N_PTBIN]={6.0,8.0,10.0,15.0,20.0, 40.0};
  TString label_pTbin[N_PTBIN]={"pT4to6","pT6to8","pT8to10","pT10to15","pT15to20", "pT20to40"};
  TString pTbin[N_PTBIN]={"4.0 < pT < 6.0 GeV","6.0 < pT < 8.0","8.0 < pT < 10.0 GeV","10.0 < pT < 15.0 GeV","15.0 < pT < 20.0 GeV", "20.0 < pT < 40.0 GeV"};
  Double_t ptbinning[N_PTBIN+1]={4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 40.0};
  

 
  TH1D * hist_Q2_tot[N_CENTBIN][N_CENTBIN_1PC];
  TH1D * hist_Q3_tot[N_CENTBIN][N_CENTBIN_1PC];

  Int_t    nBins = 500;
  Double_t xMin  = 0.0;
  Double_t xMax  = 0.50;
  
  for (int i = 0; i < N_CENTBIN; ++i) {     
    for (int j = 0; j < N_CENTBIN_1PC; ++j) { 
      
      TString h_q2_name  = TString::Format("hist_Q2_tot_%s_%s", label_centbin[i].Data(), label_centbin_1pc[j].Data());
      TString h_q2_title = TString::Format("Q2 (%s, %s)", centbin[i].Data(), centbin_1pc[j].Data());
      hist_Q2_tot[i][j] = new TH1D(h_q2_name.Data(), h_q2_title.Data(), nBins, xMin, xMax);
     
      TString h_q3_name  = TString::Format("hist_Q3_tot_%s_%s", label_centbin[i].Data(), label_centbin_1pc[j].Data());
      TString h_q3_title = TString::Format("Q3 (%s, %s)", centbin[i].Data(), centbin_1pc[j].Data()); 
      hist_Q3_tot[i][j] = new TH1D(h_q3_name.Data(), h_q3_title.Data(), nBins, xMin, xMax);
    }
  }





// =================================================================
// ### 1. File Loading and TTree Setup
// =================================================================
 
  TString Str;
  ifstream fpr(Form("%s",input_file.Data()), ios::in);
  if(!fpr.is_open()){
    cout << "List of input files not found!" << endl;
    return;
  }
  
  std::vector<TString> file_name_vector;
  string file_chain;
  while(getline(fpr, file_chain))
    {
      file_name_vector.push_back(file_chain);
    }
  
  TChain *tree = new TChain("d0Analyzer/VCNtuple_D02kpi");
  TChain *t_eventinfoana = new TChain("eventinfoana/EventInfoNtuple");

  
  for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++)
    {
      TFile *file = TFile::Open(*listIterator);
      cout << "Adding file:--- " << *listIterator << "--- to the chains" << endl;
      
      tree->Add(*listIterator);
      t_eventinfoana->Add(*listIterator);
    }
  
  tree->AddFriend(t_eventinfoana);



  Int_t MAXDsCandSize=30000;
  
  Int_t candSize;
  Int_t centrality;
  Float_t ephfpAngle[3];
  Float_t ephfmAngle[3];
  Float_t ephfpQ[3];
  Float_t ephfmQ[3];
  Float_t eptkAngle[2];
  Float_t eptkQ[2];

  //Float_t ephfmSumW[3];                                                                                                                                                                 
  //Float_t ephfpSumW[3];                                                                                                                                                                 
  Float_t ephfmSumW;
  Float_t ephfpSumW;
  
  //Float_t mva[MAXDsCandSize];
  //Double_t mvaNp[MAXDsCandSize];
  
  Float_t y[MAXDsCandSize];
  Float_t eta[MAXDsCandSize];
  Float_t pT[MAXDsCandSize];
  Float_t phi[MAXDsCandSize];
  Float_t mass[MAXDsCandSize];
  Float_t mva[MAXDsCandSize];
  Float_t mva_xg[MAXDsCandSize];
  Float_t VtxProb[MAXDsCandSize];
  Float_t VtxChi2[MAXDsCandSize];
  //Float_t 3DPointingAngle[MAXDsCandSize];
  //Float_t 3DDecayLengthSignificance[MAXDsCandSize];
  Float_t twoTrackDCA[MAXDsCandSize];

  tree->SetBranchAddress("candSize", &candSize);
  tree->SetBranchAddress("centrality", &centrality);
  tree->SetBranchAddress("ephfpAngle", &ephfpAngle);
  tree->SetBranchAddress("ephfmAngle", &ephfmAngle);
  tree->SetBranchAddress("ephfpQ", &ephfpQ);
  tree->SetBranchAddress("ephfmQ", &ephfmQ);
  tree->SetBranchAddress("eptkAngle", &eptkAngle);
  tree->SetBranchAddress("eptkQ", &eptkQ);
  tree->SetBranchAddress("ephfmSumW", &ephfmSumW);
  tree->SetBranchAddress("ephfpSumW", &ephfpSumW);
  
  tree->SetBranchStatus("*", 0);
  for (const auto& p : {"candSize", "centrality", "ephfpAngle", "ephfmAngle", "ephfpQ", "ephfmQ", "eptkAngle", "eptkQ", "ephfmSumW", "ephfpSumW"})
    tree->SetBranchStatus(p, 1);

   Int_t n_entries = tree->GetEntries();
   
   std::cout<<"nevents of first loop = "<<n_entries<<std::endl;
   
   for(Int_t ii=0; ii<n_entries; ii++){//loop in events                                                                                                     
     tree->GetEntry(ii);

     Int_t cent = centrality/2;
     
     if(ii % 100000 == 0)
       printf("current entry of first loop = %d out of %d : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);

     if (cent >= max_centbin[N_CENTBIN-1]) continue;

     TComplex aux_Q2_HFm(ephfmQ[1]*TMath::Cos(2.0*ephfmAngle[1]),ephfmQ[1]*TMath::Sin(2.0*ephfmAngle[1]),0); //q2 HF-
     TComplex aux_Q2_HFp(ephfpQ[1]*TMath::Cos(2.0*ephfpAngle[1]),ephfpQ[1]*TMath::Sin(2.0*ephfpAngle[1]),0); //q2 HF+
     TComplex aux_Q2_Trk(eptkQ[0]*TMath::Cos(2.0*eptkAngle[0]),eptkQ[0]*TMath::Sin(2.0*eptkAngle[0]),0);

     TComplex aux_Q3_HFm(ephfmQ[2]*TMath::Cos(3.0*ephfmAngle[2]),ephfmQ[2]*TMath::Sin(3.0*ephfmAngle[2]),0); // q3 HF-
     TComplex aux_Q3_HFp(ephfpQ[2]*TMath::Cos(3.0*ephfpAngle[2]),ephfpQ[2]*TMath::Sin(3.0*ephfpAngle[2]),0); //q3 HF+
     TComplex aux_Q3_Trk(eptkQ[1]*TMath::Cos(3.0*eptkAngle[1]),eptkQ[1]*TMath::Sin(3.0*eptkAngle[1]),0);

     //TComplex Q2_total = aux_Q2_HFm + aux_Q2_HFp;
     //TComplex Q3_total = aux_Q3_HFm + aux_Q3_HFp;

     Double_t Q2_total = ephfmQ[1] + ephfpQ[1];
     Double_t Q3_total = ephfpQ[2] + ephfpQ[2];

     Double_t Q2_weight = ephfmSumW + ephfpSumW;
     Double_t Q3_weight = ephfmSumW + ephfpSumW;

     Double_t Q2_total_norm = Q2_total/Q2_weight;
     Double_t Q3_total_norm = Q3_total/Q3_weight;
     
     Int_t j_idx = cent; //For 1% centrality bin
     
     for (int i_idx = 0; i_idx < N_CENTBIN; ++i_idx) {
       
    
       if (cent >= min_centbin[i_idx] && cent < max_centbin[i_idx]) {
      
	 hist_Q2_tot[i_idx][j_idx]->Fill(Q2_total_norm);
	 hist_Q3_tot[i_idx][j_idx]->Fill(Q3_total_norm);
	 
	 break;
       }       

     }

   }// ---- First event loop ------

   // =================================================================
   // ### Divide q2 in 10 slices for each cent bin
   // =================================================================

   std::cout << "Merging 1% Q2 histograms for plotting..." << std::endl;
   TH1D *hQ2_merged[N_CENTBIN];
   TH1D *hQ3_merged[N_CENTBIN];
   
   for (int i_idx = 0; i_idx < N_CENTBIN; ++i_idx) {
     TString hname_q2 = TString::Format("hQ2_merged_%s", label_centbin[i_idx].Data());
     TString htitle_q2 = TString::Format("Total Q2 for %s", centbin[i_idx].Data());

     TString hname_q3 = TString::Format("hQ3_merged_%s", label_centbin[i_idx].Data());
     TString htitle_q3 = TString::Format("Total Q3 for %s", centbin[i_idx].Data());
    
     // Find the first valid histogram in this bin to clone
     TH1D* h_first_q2 = nullptr;
     for (int j = min_centbin[i_idx]; j < max_centbin[i_idx]; ++j) {
       if (hist_Q2_tot[i_idx][j] && hist_Q2_tot[i_idx][j]->GetEntries() > 0) {
	 h_first_q2 = hist_Q2_tot[i_idx][j];
	 break; // Found it, stop this loop
       }
     }
     
     TH1D* h_first_q3 = nullptr;
     for (int j = min_centbin[i_idx]; j < max_centbin[i_idx]; ++j) {
       if (hist_Q3_tot[i_idx][j] && hist_Q3_tot[i_idx][j]->GetEntries() > 0) {
	 h_first_q3 = hist_Q3_tot[i_idx][j];
	 break; // Found it, stop this loop
       }
    }
     
     
     // Clone the first histogram to get the binning right
     if (h_first_q2) {
       hQ2_merged[i_idx] = (TH1D*)h_first_q2->Clone(hname_q2);
       hQ2_merged[i_idx]->SetTitle(htitle_q2);
       
       // Add all OTHER histograms
       for (int j = min_centbin[i_idx]; j < max_centbin[i_idx]; ++j) {
	 if (hist_Q2_tot[i_idx][j] && hist_Q2_tot[i_idx][j] != h_first_q2 && hist_Q2_tot[i_idx][j]->GetEntries() > 0) {
	   hQ2_merged[i_idx]->Add(hist_Q2_tot[i_idx][j]);
	 }
       }
     } else {
       std::cout << "Warning: No Q2 entries found for coarse bin " << i_idx << std::endl;
       hQ2_merged[i_idx] = nullptr;
     }
     
     // --- 4. Clone and add for Q3 (if found) ---
     if (h_first_q3) {
       hQ3_merged[i_idx] = (TH1D*)h_first_q3->Clone(hname_q3);
       hQ3_merged[i_idx]->SetTitle(htitle_q3);
       
       // Add all OTHER histograms
       for (int j = min_centbin[i_idx]; j < max_centbin[i_idx]; ++j) {
	 if (hist_Q3_tot[i_idx][j] && hist_Q3_tot[i_idx][j] != h_first_q3 && hist_Q3_tot[i_idx][j]->GetEntries() > 0) {
	   hQ3_merged[i_idx]->Add(hist_Q3_tot[i_idx][j]);
	 }
       }
     } else {
       std::cout << "Warning: No Q3 entries found for coarse bin " << i_idx << std::endl;
       hQ3_merged[i_idx] = nullptr;
     }
   }
   
   
   // =================================================================
   // ### Divide q2 and q3 in 10 slices for each cent bin
   // =================================================================
   const int N_Q_BINS = 10; // Assuming 10 bins for both
   Double_t probabilities[N_Q_BINS + 1];
   for (int k = 0; k <= N_Q_BINS; ++k) {
     probabilities[k] = (Double_t)k / N_Q_BINS;
   }
   
   // --- Define arrays for both Q2 and Q3 cuts ---
   Double_t q2_cuts_1pc[N_CENTBIN_1PC][N_Q_BINS + 1];
   TAxis   *q2_axes_1pc[N_CENTBIN_1PC];
   
   Double_t q3_cuts_1pc[N_CENTBIN_1PC][N_Q_BINS + 1];
   TAxis   *q3_axes_1pc[N_CENTBIN_1PC];
   
   std::cout << "Calculating q2 and q3 quantiles for each 1% bin..." << std::endl;
   
   // --- Single Combined Loop ---
   for (int i = 0; i < N_CENTBIN; ++i) { // Coarse bin loop
     for (int j = min_centbin[i]; j < max_centbin[i]; ++j) { // 1% bin loop
       
       // --- Process Q2 ---
       TH1D* h_1pc_q2 = hist_Q2_tot[i][j];
       if (h_1pc_q2 && h_1pc_q2->GetEntries() > 0) {
	 h_1pc_q2->GetQuantiles(N_Q_BINS + 1, &q2_cuts_1pc[j][0], probabilities);
	 q2_axes_1pc[j] = new TAxis(N_Q_BINS, &q2_cuts_1pc[j][0]);
       } else {
	 q2_axes_1pc[j] = nullptr;
       }
       
       // --- Process Q3 ---
       TH1D* h_1pc_q3 = hist_Q3_tot[i][j];
       if (h_1pc_q3 && h_1pc_q3->GetEntries() > 0) {
	 h_1pc_q3->GetQuantiles(N_Q_BINS + 1, &q3_cuts_1pc[j][0], probabilities);
	 q3_axes_1pc[j] = new TAxis(N_Q_BINS, &q3_cuts_1pc[j][0]);
       } else {
	 q3_axes_1pc[j] = nullptr;
       }
     }
   }
   

   // =================================================================
   // ### final histograms for each coarse centrality bin
   // =================================================================

   TH1D *hist_Q2_final_slices[N_CENTBIN][N_Q_BINS];
   TH1D *hist_Q3_final_slices[N_CENTBIN][N_Q_BINS];
   
   Int_t    nBins_q2 = hist_Q2_tot[0][0]->GetNbinsX();
   Double_t xMin_q2  = hist_Q2_tot[0][0]->GetXaxis()->GetXmin();
   Double_t xMax_q2  = hist_Q2_tot[0][0]->GetXaxis()->GetXmax();
   

   Int_t    nBins_q3 = hist_Q3_tot[0][0]->GetNbinsX();
   Double_t xMin_q3  = hist_Q3_tot[0][0]->GetXaxis()->GetXmin();
   Double_t xMax_q3  = hist_Q3_tot[0][0]->GetXaxis()->GetXmax();

   
   for (int i = 0; i < N_CENTBIN; ++i) {
     for (int k = 0; k < N_Q_BINS; ++k) { 
       TString hname = TString::Format("hQ2_slice_%s_s%d", label_centbin[i].Data(), k);
       TString htitle = TString::Format("Q2 Slice %d (%d-%d%%) for %s", k, k*10, (k+1)*10, centbin[i].Data());
       hist_Q2_final_slices[i][k] = new TH1D(hname, htitle, nBins_q2, xMin_q2, xMax_q2);

       
     }
   }

   for (int i = 0; i < N_CENTBIN; ++i) {
    for (int k = 0; k < N_Q_BINS; ++k) {
        // Name the histogram based on the Q3 distribution,
        // but binned by the Q2 slice
        TString hname = TString::Format("hQ3_slice_%s_s%d", label_centbin[i].Data(), k);
        TString htitle = TString::Format("Q3 dist for %s, Q2 Slice %d (%d-%d%%)", centbin[i].Data(), k, k*10, (k+1)*10);
        
        hist_Q3_final_slices[i][k] = new TH1D(hname, htitle, nBins_q3, xMin_q3, xMax_q3);
    }
   }
   


   std::cout << "Starting second event loop to fill final sliced histograms..." << std::endl;
   
   for(Int_t ii=0; ii<n_entries; ii++){//loop in events
     tree->GetEntry(ii);

     Int_t cent = centrality/2;
     
     if(ii % 100000 == 0)
       printf("current entry of second loop = %d out of %d : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);
     
     // Centrality cut
     if (cent >= max_centbin[N_CENTBIN-1]) continue;
     
     // --- Recalculate Q vectors ---
     TComplex aux_Q2_HFm(ephfmQ[1]*TMath::Cos(2.0*ephfmAngle[1]),ephfmQ[1]*TMath::Sin(2.0*ephfmAngle[1]),0);
     TComplex aux_Q2_HFp(ephfpQ[1]*TMath::Cos(2.0*ephfpAngle[1]),ephfpQ[1]*TMath::Sin(2.0*ephfpAngle[1]),0);
     TComplex aux_Q3_HFm(ephfmQ[2]*TMath::Cos(3.0*ephfmAngle[2]),ephfmQ[2]*TMath::Sin(3.0*ephfmAngle[2]),0);
     TComplex aux_Q3_HFp(ephfpQ[2]*TMath::Cos(3.0*ephfpAngle[2]),ephfpQ[2]*TMath::Sin(3.0*ephfpAngle[2]),0);

     //TComplex Q2_total = aux_Q2_HFm + aux_Q2_HFp;
     //TComplex Q3_total = aux_Q3_HFm + aux_Q3_HFp;

     Double_t Q2_total = ephfmQ[1] + ephfpQ[1];
     Double_t Q3_total = ephfpQ[2] + ephfpQ[2];
     
     Double_t Q2_weight = ephfmSumW + ephfpSumW;
     Double_t Q3_weight = ephfmSumW + ephfpSumW;

     
     Double_t q2_val = (Q2_weight > 0) ? Q2_total/Q2_weight : -999;
     Double_t q3_val = (Q3_weight > 0) ? Q3_total/Q3_weight : -999;

     Int_t j_idx = cent; // e.g., 7
     int i_idx = -1;
     
     for (int i = 0; i < N_CENTBIN; ++i) {
       if (cent >= min_centbin[i] && cent < max_centbin[i]) {
	 i_idx = i; // e.g., 0 (for 0-10%)
	 break;
       }
     }
     
     // Check if event is valid
     if (i_idx != -1 && q2_axes_1pc[j_idx] != nullptr) {
       
       TAxis *current_q2_axis = q2_axes_1pc[j_idx];
       
       
       Int_t k_bin = current_q2_axis->FindBin(q2_val);
       Int_t k_idx = k_bin - 1; // Convert to 0-9 index
       
       
       if (k_idx < 0) k_idx = 0;
       if (k_idx >= N_Q_BINS) k_idx = N_Q_BINS - 1;
       
       hist_Q2_final_slices[i_idx][k_idx]->Fill( q2_val );
     }


     if (i_idx != -1 && q3_axes_1pc[j_idx] != nullptr) {
       
       TAxis *current_q3_axis = q3_axes_1pc[j_idx];
       
       Int_t k_bin = current_q3_axis->FindBin(q3_val);
       Int_t k_idx = k_bin - 1; // Convert to 0-9 index
       
       
       if (k_idx < 0) k_idx = 0;
       if (k_idx >= N_Q_BINS) k_idx = N_Q_BINS - 1;
       
       hist_Q3_final_slices[i_idx][k_idx]->Fill( q3_val ); 
     }
   } // --- End of second event loop ---
   
   
   


   // =================================================================
   // ### Save all histograms for analysis
   // =================================================================

   fout->cd();
   fout->mkdir("Q2_raw_1pc");
   fout->cd("Q2_raw_1pc"); 

   for (int i = 0; i < N_CENTBIN; ++i) {
     for (int j = min_centbin[i]; j < max_centbin[i]; ++j) {   
       if (hist_Q2_tot[i][j] && hist_Q2_tot[i][j]->GetEntries() > 0) {
	 hist_Q2_tot[i][j]->Write();
       }
     }
   }
   fout->mkdir("Q3_raw_1pc");
   fout->cd("Q3_raw_1pc"); 
   for (int i = 0; i < N_CENTBIN; ++i) {
    for (int j = min_centbin[i]; j < max_centbin[i]; ++j) {
      if (hist_Q3_tot[i][j] && hist_Q3_tot[i][j]->GetEntries() > 0) {
	hist_Q3_tot[i][j]->Write();
      }
    }
   }
   
   fout->mkdir("Q2_final_slices");
   fout->cd("Q2_final_slices");
   
   for (int i = 0; i < N_CENTBIN; ++i) {
     hQ2_merged[i]->Write();
     for (int k = 0; k < N_Q_BINS; ++k) {
       TH1D* hist_to_save_q2 = hist_Q2_final_slices[i][k];
       if (hist_to_save_q2 && hist_to_save_q2->GetEntries() > 0) {
	 hist_to_save_q2->Write();
       }
     }
   }
   fout->mkdir("Q3_final_slices");
   fout->cd("Q3_final_slices");
   
   for (int i = 0; i < N_CENTBIN; ++i) {
     hQ3_merged[i]->Write();
     for (int k = 0; k < N_Q_BINS; ++k) {
       TH1D* hist_to_save_q3 = hist_Q3_final_slices[i][k];
       if (hist_to_save_q3 && hist_to_save_q3->GetEntries() > 0) {
	 hist_to_save_q3->Write();
       }
     }
   }
   
 


   fout->Close();

   auto stop = high_resolution_clock::now();
   auto duration = duration_cast<minutes>(stop - start);
   cout << "Total time taken: "<< duration.count() << "minutes" << endl;
   
   std::cout<<" \n Done!!! "<<std::endl;

}














/*int main(int argc, char **argv)
{
  using namespace std;

  TString inputfile;
  int itxtout;
  TString outfile;

  if(argc == 1)
    {
      std::cout<<"You did not pass any argument to the code other than the program name"<<std::endl;
      std::cout<<"You need to pass 4 arguments including the program name"<<std::endl;
    }

  if(argc >=1 && argc <=3)
    {
      std::cout<<"Only "<<argc<<" arguments you have given including the program name"<<std::endl;
      std::cout<<"You need to pass 4 arguments including the program name"<<std::endl;
    }

  if(argc == 4)
    {
      std::cout<<std::endl;
      std::cout<<"You have given 4 arguments including the program name;  Your program will run"<<std::endl;
      std::cout<<std::endl;

      inputfile = argv[1];
      itxtout = atoi(argv[2]);
      outfile = argv[3];

      hist_vn_mass(inputfile, itxtout, outfile);
    }
  return 0;
  }*/




