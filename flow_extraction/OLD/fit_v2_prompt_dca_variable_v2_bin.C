#include <cstdlib> 
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFileCollection.h"
#include "TLegend.h"

#include "TMath.h"
#include "TComplex.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include <TNtuple.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TCanvas.h>

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/ESE_Cuts_MB11to21_Jan22.h"

using namespace std;

int fit_v2_prompt_dca(int target1pct = -1, int target_pt = -1) //updated code with q2HF (Nov 20th 2024)
{

  //We can take this bining for first try!
  //const int N_VBINS=44;
  //Double_t vnbinning[N_VBINS+1]={-10,-8,-7,-6,-5.4,-1.8,-1.6,-1.4,-1.3,-1.0,-0.90,-0.80,-0.70,-0.60,-0.50,-0.40,-0.30,-0.25,-0.20,-0.15,-0.10,-0.05,0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.0,1.3,1.4,1.6,1.8,5.4,6,7,8,10}; // label binning 2??

  const int N_QBINS = 10;
  const int N_CENTBINS_1 = 90;
  const int N_CENTBINS = 4;
  Int_t min_centbin[N_CENTBINS]={0, 10, 30, 50};
  Int_t max_centbin[N_CENTBINS]={10, 30, 50, 90};
  TString cen_label[N_CENTBINS]={"cent0to10", "cent10to30", "cent30to50", "cent50to90"};
  std::string cen_name[N_CENTBINS]={"cent0to10", "cent10to30", "cent30to50", "cent50to90"};
  TString centbin[N_CENTBINS]={"0 <= cent < 10 %", "10 <= cent < 30 %", "30 <= cent < 50 %", "50 <= cent < 90 %"};
  Int_t cen_edges[N_CENTBINS+1]={0, 10, 30, 50, 90};

  auto get_cent_group = [&](int cent1pct)->int {
    for (int g=0; g<N_CENTBINS; ++g) {
      if (cent1pct >= min_centbin[g] && cent1pct < max_centbin[g]) return g;
    }
    return -1; 
  };
  
  const int N_PTBINS = 10;
  Float_t min_pTbin[N_PTBINS]={1.0, 2.0, 3.0, 4.0, 5.0, 6.0,8.0,10.0,15.0, 20.0};
  Float_t max_pTbin[N_PTBINS]={2.0, 3.0, 4.0, 5.0, 6.0,8.0,10.0,15.0,20.0, 40.0};
  TString label_pTbin[N_PTBINS]={"pT1to2", "pT2to3", "pT3to4","pT4to5", "pT5to6","pT6to8","pT8to10","pT10to15","pT15to20", "pT20to40"};
  std::string pt_name[N_PTBINS]={"pT1to2", "pT2to3", "pT3to4","pT4to5", "pT5to6","pT6to8","pT8to10","pT10to15","pT15to20", "pT20to40"};
  TString pTbin[N_PTBINS]={"1.0 < pT < 2.0 GeV", "2.0 < pT < 3.0 GeV", "3.0 < pT < 4.0 GeV", "4.0 < pT < 5.0 GeV", "5.0 < pT < 6.0 GeV","6.0 < pT < 8.0","8.0 < pT < 10.0 GeV","10.0 < pT < 15.0 GeV","15.0 < pT < 20.0 GeV", "20.0 < pT < 40.0 GeV"};
  double pt_edges[N_PTBINS+1] = {1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 40};

  const int N_VBINS = 44;
  const int N_BINS = 2300; //splitting -2.3 to 0 into 2500 bins
  const int N_EDGES = 10;
  Double_t vnbinning[N_CENTBINS][N_PTBINS][N_VBINS+1]; // For v2
  Double_t vnbinning_v3[N_CENTBINS][N_PTBINS][N_VBINS+1]; // For v3
  Double_t vnbinning_side[16]={-10,-8.2,-6.6,-5.2,-4.0,-3.5,-3.0,-2.3,2.3,3.0,3.5,4.0,5.2,6.6,8.2,10};
  

  
  const int N_MASSBINS = 52;
  Float_t fit_range_low = 1.74;
  Float_t fit_range_high = 2.00;
  Float_t width = (fit_range_high-fit_range_low)/N_MASSBINS;
  Double_t D0_mass = 1.8648;


  auto nt_Data = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/SP_MB11to21_Jan27/ROOT/flow_Analysis_out_combined.root");
  auto inf_MC = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/D0mass_MC_5p36TeV_OfficialMC_Jan26/ROOT/D0_Mass_MC_out_combined.root");

  auto input_nt = (TNtuple*)nt_Data->Get("nt");

  Float_t cen_val, q2_val, q3_val, pT_val, mass_val, v2_val, v3_val, dca_val, y_val;

  input_nt->SetBranchAddress("centrality",&cen_val);
  input_nt->SetBranchAddress("pT",&pT_val);
  input_nt->SetBranchAddress("mass",&mass_val);
  input_nt->SetBranchAddress("dca",&dca_val);
  input_nt->SetBranchAddress("y",&y_val);
  input_nt->SetBranchAddress("q2_hf_total",&q2_val);
  input_nt->SetBranchAddress("q3_hf_total",&q3_val);
  input_nt->SetBranchAddress("v2",&v2_val);
  input_nt->SetBranchAddress("v3",&v3_val);
  
  auto h_v2_bin = new TH1D("h_v2_bin","h_v2_bin",N_BINS,-2.3,0);
  auto h_v3_bin = new TH1D("h_v3_bin","h_v3_bin",N_BINS,-2.3,0);

  TString cuts;

  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {

      cuts = "dca<0.0085 && mass<2.0 && mass>1.74 && centrality>="+to_string(cen_edges[i_cen])+" && centrality<"+to_string(cen_edges[i_cen+1])+" && pT>="+to_string(pt_edges[i_pt])+" && pT<"+to_string(pt_edges[i_pt+1]);


      cout << endl << endl << cuts << endl;  
      h_v2_bin->Reset();
      //input_nt->Draw("v2>>h_v2_bin",cuts);
      //Int_t n_entries_v2 = h_v2_bin->Integral();
      input_nt->Project(h_v2_bin->GetName(),"v2",cuts);
      Int_t n_entries_v2 = static_cast<Int_t>(h_v2_bin->GetEntries());
      cout << endl << "Total entries v2: " << n_entries_v2 << endl;

      h_v3_bin->Reset();
      //input_nt->Draw("v3>>h_v3_bin",cuts);
      //Int_t n_entries_v3 = h_v3_bin->Integral();
      input_nt->Project(h_v3_bin->GetName(),"v3",cuts);
      Int_t n_entries_v3 = static_cast<Int_t>(h_v3_bin->GetEntries());
      cout << endl << "Total entries v3: " << n_entries_v3 << endl;

      for (int ib=0; ib<=N_VBINS; ++ib) vnbinning[i_cen][i_pt][ib] = 0.0;
      for (int i_bin=0; i_bin<8; i_bin++) {
        vnbinning[i_cen][i_pt][i_bin] = vnbinning_side[i_bin];
        vnbinning[i_cen][i_pt][i_bin+27] = vnbinning_side[i_bin+8];
      }
      vnbinning[i_cen][i_pt][17]=0.0;
      
      cout << "bins defined before loop"<< endl;

      if (n_entries_v2 > 0) {
        std::vector<double> cum(N_BINS+1,0.0);
        for (int b=1; b<=N_BINS; ++b) {
          cum[b] = cum[b-1] + h_v2_bin->GetBinContent(b);
        }
        for (int i_edge=1; i_edge<=N_EDGES; ++i_edge) {
          double target = i_edge * (n_entries_v2 / (double)N_EDGES);
          // find first bin index j such that cum[j] > target
          int j = 1;
          while (j <= N_BINS && cum[j] <= target) ++j;
          if (j <= N_BINS) {
            double edge = h_v2_bin->GetBinLowEdge(j);
            cout << "edge " << i_edge << " defined: +=" << edge << "  -=" << -edge << endl;
            vnbinning[i_cen][i_pt][i_edge+7] = edge;
            vnbinning[i_cen][i_pt][27-i_edge] = -edge;
          } else {
            // fallback: use max bin edge
            double edge = h_v2_bin->GetBinLowEdge(N_BINS+1);
            vnbinning[i_cen][i_pt][i_edge+7] = edge;
            vnbinning[i_cen][i_pt][27-i_edge] = -edge;
          }
        }
      } else {
        cout << "No entries for v2 quantiles, using default symmetric vnbinning for cen " << i_cen << " pt " << i_pt << endl;
      }

      if (n_entries_v3 > 0) {
        std::vector<double> cum3(N_BINS+1,0.0);
        for (int b=1; b<=N_BINS; ++b) {
          cum3[b] = cum3[b-1] + h_v3_bin->GetBinContent(b);
        }
        for (int i_edge=1; i_edge<=N_EDGES; ++i_edge) {
          double target3 = i_edge * (n_entries_v3 / (double)N_EDGES);
          int j3 = 1;
          while (j3 <= N_BINS && cum3[j3] <= target3) ++j3;
          if (j3 <= N_BINS) {
            double edge3 = h_v3_bin->GetBinLowEdge(j3);
            cout << "v3 edge " << i_edge << " defined: +=" << edge3 << "  -=" << -edge3 << endl;
            vnbinning_v3[i_cen][i_pt][i_edge+7] = edge3;
            vnbinning_v3[i_cen][i_pt][27-i_edge] = -edge3;
          } else {
            double edge3 = h_v3_bin->GetBinLowEdge(N_BINS+1);
            vnbinning_v3[i_cen][i_pt][i_edge+7] = edge3;
            vnbinning_v3[i_cen][i_pt][27-i_edge] = -edge3;
          }
        }
      } else {
        cout << "No entries for v3 quantiles, using default symmetric v3 binning for cen " << i_cen << " pt " << i_pt << endl;
      }
            
    } // --- PT BIN ---
  } // -- CENTBIN ---


  int ymax, a, b;
  Double_t yield, yield_error;
  Double_t chi2_ndf_for_q2_v2[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS], sigma_v2[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS];
  Double_t chi2_ndf_for_q2_v3[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS], sigma_v3[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS];
  TGraph* sigma_v2_fit[N_CENTBINS][N_QBINS][N_PTBINS];
  TGraph* sigma_v3_fit[N_CENTBINS][N_QBINS][N_PTBINS];

  TH1D *h_mass_default[N_CENTBINS][N_PTBINS];
  TH1D *h_mass_default_1percent_cen[N_CENTBINS_1][N_PTBINS];
  TH1D *h_mass[N_CENTBINS][N_QBINS][N_PTBINS];
  TH1D *h_mass_1pecent_cen[N_CENTBINS_1][N_QBINS][N_PTBINS];
  TH1D *h_v2[N_CENTBINS][N_QBINS][N_PTBINS];
  TH1D *h_v3[N_CENTBINS][N_QBINS][N_PTBINS];

  TH1D *h_mass_v2_fit[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS];
  TH1D *h_mass_v2_fit_1pecent_cen[N_CENTBINS_1][N_QBINS][N_PTBINS][N_VBINS];

  
  TH1D *h_mass_v3_fit[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS];
  TH1D *h_mass_v3_fit_1pecent_cen[N_CENTBINS_1][N_QBINS][N_PTBINS][N_VBINS];


  
  TString h_mass_name, h_v2_name, h_mass_v2_name, h_name_0, h_name_1;
  //histograms for 1% centrlity bins
  for (int i_cen=0; i_cen<N_CENTBINS_1; i_cen++) {
    
    if (target1pct >= 0 && i_cen != target1pct) continue;

    int cen_group_for_binning = get_cent_group(i_cen);
    if (cen_group_for_binning < 0) continue;
 
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {

      if (target_pt >= 0 && i_pt != target_pt) continue;
      
      h_name_0 = "h_mass_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_pt_"+pt_name[i_pt]+"_def";
      h_mass_default_1percent_cen[i_cen][i_pt] = new TH1D(h_name_0,h_name_0,N_MASSBINS,fit_range_low,fit_range_high);
      for (int i_q2=0; i_q2<N_QBINS; i_q2++) {
        h_name_0 = "h_mass_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_q2bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt];
        h_mass_1pecent_cen[i_cen][i_q2][i_pt] = new TH1D(h_name_0,h_name_0,N_MASSBINS,fit_range_low,fit_range_high);

        for(int i_v2=0; i_v2<N_VBINS; i_v2++)
        {
          h_name_1 = "hist_mass_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_q2bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt]+"_in_v2bin_"+to_string(vnbinning[i_cen/10][i_pt][i_v2])+"_"+to_string(vnbinning[i_cen/10][i_pt][i_v2+1]);
          h_mass_v2_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2] = new TH1D(h_name_1,h_name_1,N_MASSBINS,fit_range_low,fit_range_high);

          // v3 1% hist creation (mirror v2)
          h_name_1 = "hist_mass_v3_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_q2bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt]+"_in_v3bin_"+to_string(vnbinning_v3[i_cen/10][i_pt][i_v2])+"_"+to_string(vnbinning_v3[i_cen/10][i_pt][i_v2+1]);
          h_mass_v3_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2] = new TH1D(h_name_1,h_name_1,N_MASSBINS,fit_range_low,fit_range_high);
        }
      }
    }
  }

  //histograms for analysis bins
  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
      h_mass_name = "h_mass_cen_"+cen_name[i_cen]+"_pt_"+pt_name[i_pt]+"_def";
      h_mass_default[i_cen][i_pt] = new TH1D(h_mass_name, h_mass_name, N_MASSBINS, fit_range_low, fit_range_high); 
      for (int i_q2=0; i_q2<N_QBINS; i_q2++) {
        h_mass_name = "h_mass_cen_"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt];
        h_mass[i_cen][i_q2][i_pt] = new TH1D(h_mass_name, h_mass_name, N_MASSBINS, fit_range_low, fit_range_high);
        h_v2_name = "hist_v2_cen_"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt];
        h_v2[i_cen][i_q2][i_pt] = new TH1D(h_v2_name, h_v2_name, N_VBINS, vnbinning[i_cen][i_pt]);
        // create v3 histogram using vnbinning_v3
        TString h_v3_name = "hist_v3_cen_"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt];
        h_v3[i_cen][i_q2][i_pt] = new TH1D(h_v3_name, h_v3_name, N_VBINS, vnbinning_v3[i_cen][i_pt]);

        for(int i_v2=0; i_v2<N_VBINS; i_v2++)
        {
          h_mass_v2_name = "hist_mass_cen_"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt]+"_in_v2bin_"+to_string(vnbinning[i_cen][i_pt][i_v2])+"_"+to_string(vnbinning[i_cen][i_pt][i_v2+1]);
          h_mass_v2_fit[i_cen][i_q2][i_pt][i_v2]=new TH1D(h_mass_v2_name, h_mass_v2_name, N_MASSBINS, fit_range_low, fit_range_high);

          // create v3 analysis histograms (mirror v2)
          TString h_mass_v3_name = "hist_mass_v3_cen_"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt]+"_in_v3bin_"+to_string(vnbinning_v3[i_cen][i_pt][i_v2])+"_"+to_string(vnbinning_v3[i_cen][i_pt][i_v2+1]);
          h_mass_v3_fit[i_cen][i_q2][i_pt][i_v2] = new TH1D(h_mass_v3_name, h_mass_v3_name, N_MASSBINS, fit_range_low, fit_range_high);
        }
      }
    }
  }

  //Variables and graphs for results
  Double_t yield_v2[N_VBINS], yield_error_v2[N_VBINS], v2_x[N_VBINS], v2_x_err[N_VBINS], mean_val[N_PTBINS], mean_error[N_PTBINS];
  Double_t yield_v3[N_VBINS], yield_error_v3[N_VBINS], v3_x[N_VBINS], v3_x_err[N_VBINS];
  TGraphErrors* yield_v2_graph[N_CENTBINS][N_QBINS][N_PTBINS];
  TGraphErrors* yield_v3_graph[N_CENTBINS][N_QBINS][N_PTBINS];
  TGraphErrors* v2_cen[N_CENTBINS][N_QBINS];
  TGraphErrors* v3_cen[N_CENTBINS][N_QBINS];
  TGraph* chi2_ndf_for_q2_v2_graph[N_CENTBINS][N_QBINS][N_PTBINS];
  TGraph* chi2_ndf_for_q2_v3_graph[N_CENTBINS][N_QBINS][N_PTBINS];

  //Read and fill histograms from ntuples
  Long64_t N_ENTRIES = input_nt->GetEntriesFast();
  cout << "Reading ntuple to generate histograms" << endl;


  Long64_t maxEvents = 100; // <-- set how many events you want
  if (maxEvents > 0 && N_ENTRIES > maxEvents) N_ENTRIES = maxEvents;


  
  for(Long64_t i_entry=0; i_entry<N_ENTRIES; i_entry++) {
    input_nt->GetEntry(i_entry);

    if(i_entry%400000==0) cout <<i_entry<<" / "<<N_ENTRIES<< "  "<<100*i_entry/N_ENTRIES<<"%"<<endl;

    if(dca_val>=0.0085) continue;

    for (int i_cen=0; i_cen<N_CENTBINS_1; i_cen++) {

      if (target1pct >= 0 && i_cen != target1pct) continue;
      if(cen_val<(i_cen) || cen_val>=(i_cen+1)) continue;

      int cen_group = get_cent_group(i_cen);
      if (cen_group < 0) continue;
      
      for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {

	if (target_pt >= 0 && i_pt != target_pt) continue;

    if (pT_val<pt_edges[i_pt] || pT_val>=pt_edges[i_pt+1]) continue;
        h_mass_default_1percent_cen[i_cen][i_pt]->Fill(mass_val);

    for (int i_q2=0; i_q2<N_QBINS; i_q2++) {
          if (q2_val>=q2_cuts[i_cen][i_q2] && q2_val<q2_cuts[i_cen][i_q2+1]) {
            h_mass_1pecent_cen[i_cen][i_q2][i_pt]->Fill(mass_val);
            for(int i_v2=0; i_v2<N_VBINS; i_v2++)
            {
              if(v2_val>=vnbinning[cen_group][i_pt][i_v2] && v2_val<vnbinning[cen_group][i_pt][i_v2+1])
          {
                h_mass_v2_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2]->Fill(mass_val);
              }
              // fill v3 1% histograms
              if(v3_val>=vnbinning_v3[cen_group][i_pt][i_v2] && v3_val<vnbinning_v3[cen_group][i_pt][i_v2+1])
              {
                h_mass_v3_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2]->Fill(mass_val);
              }
            }
          }
        }
      }
    }
  }

  //Add the histograms of 1% centrality bins to get the 10% centrality bins
  if (target1pct < 0 && target_pt < 0 ) {
    for(int i_pt=0; i_pt<N_PTBINS; i_pt++) {
      for(int i_cen=0; i_cen<N_CENTBINS_1; i_cen++) {
        int cen_group = get_cent_group(i_cen);
        if (cen_group < 0) continue;
        
        h_mass_default[cen_group][i_pt]->Add(h_mass_default_1percent_cen[i_cen][i_pt]);
        for(int i_q2=0; i_q2<N_QBINS; i_q2++) {
          h_mass[cen_group][i_q2][i_pt]->Add(h_mass_1pecent_cen[i_cen][i_q2][i_pt]);
          for(int i_v2=0; i_v2<N_VBINS; i_v2++) {
            h_mass_v2_fit[cen_group][i_q2][i_pt][i_v2]->Add(h_mass_v2_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2]);
            h_mass_v3_fit[cen_group][i_q2][i_pt][i_v2]->Add(h_mass_v3_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2]);
          }
        }
      }
    }
  }

  TLatex* tex = new TLatex;
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);

  auto c1= new TCanvas("c1","c1",600,600);

  //New fitting macro (same as the default fit mechanics)

  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
      
      TH1D* h_mc_match_signal;
      TH1D* h_mc_match_all;

      h_mc_match_signal = (TH1D*)inf_MC->Get("hMass_Signal_"+label_pTbin[i_pt]);
      h_mc_match_all = (TH1D*)inf_MC->Get("hMass_Swap_"+label_pTbin[i_pt]);

      //for massfit defined by fit range from each q2 bin instead of full q2 bin
      for (int i_q2=0; i_q2<N_QBINS; i_q2++) {
        ymax = h_mass[i_cen][i_q2][i_pt]->GetBinContent(h_mass[i_cen][i_q2][i_pt]->GetMaximumBin());

        h_mass[i_cen][i_q2][i_pt]->SetMinimum(0);
        h_mass[i_cen][i_q2][i_pt]->SetMarkerSize(0.5);
        h_mass[i_cen][i_q2][i_pt]->SetTitle("");
        h_mass[i_cen][i_q2][i_pt]->SetMarkerStyle(20);
        h_mass[i_cen][i_q2][i_pt]->SetLineWidth(1);
        h_mass[i_cen][i_q2][i_pt]->SetOption("e");
        h_mass[i_cen][i_q2][i_pt]->GetXaxis()->SetRangeUser(fit_range_low,fit_range_high);
        h_mass[i_cen][i_q2][i_pt]->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
        h_mass[i_cen][i_q2][i_pt]->GetYaxis()->SetTitle("Entries / 2.5 MeV");
        h_mass[i_cen][i_q2][i_pt]->GetXaxis()->CenterTitle();
        h_mass[i_cen][i_q2][i_pt]->GetYaxis()->CenterTitle();
        h_mass[i_cen][i_q2][i_pt]->GetXaxis()->SetTitleOffset(1.3);
        h_mass[i_cen][i_q2][i_pt]->GetYaxis()->SetTitleOffset(2);
        h_mass[i_cen][i_q2][i_pt]->GetXaxis()->SetLabelOffset(0.007);
        h_mass[i_cen][i_q2][i_pt]->GetYaxis()->SetLabelOffset(0.007);
        h_mass[i_cen][i_q2][i_pt]->GetXaxis()->SetTitleSize(0.045);
        h_mass[i_cen][i_q2][i_pt]->GetYaxis()->SetTitleSize(0.045);
        h_mass[i_cen][i_q2][i_pt]->GetXaxis()->SetTitleFont(42);
        h_mass[i_cen][i_q2][i_pt]->GetYaxis()->SetTitleFont(42);
        h_mass[i_cen][i_q2][i_pt]->GetXaxis()->SetLabelFont(42);
        h_mass[i_cen][i_q2][i_pt]->GetYaxis()->SetLabelFont(42);
        h_mass[i_cen][i_q2][i_pt]->GetXaxis()->SetLabelSize(0.04);
        h_mass[i_cen][i_q2][i_pt]->GetYaxis()->SetLabelSize(0.04);
        h_mass[i_cen][i_q2][i_pt]->SetStats(kFALSE);

        h_mass[i_cen][i_q2][i_pt]->GetXaxis()->SetNoExponent(true);

        TF1* f = new TF1(Form("f_ptbin_%d",i_pt),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])) ) + [8] + [9]*x + [10]*x*x + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12])) + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])) )",fit_range_low,fit_range_high);

        f->SetLineColor(2);
        f->SetLineWidth(1);

        f->SetParameter(0,100.); //Prefactor to signal+swap
        f->SetParameter(1,D0_mass); //D0 mass
        f->SetParameter(2,0.03); //Sigma for 1st guassian in signal
        f->SetParameter(3,0.005); //Sigma for 2nd gaussian in signal
        f->SetParameter(4,0.1); // Prefactor for only 1st gaussian

        f->FixParameter(5,1); //No swap component for fitting double gaussian
        f->FixParameter(6,0); //always 0 in MC (scaling factor for width)
        f->FixParameter(7,0.1); //does not really mater here as yield is fix to 0 (swap component width)
        f->FixParameter(8,0); //Background
        f->FixParameter(9,0); //Background
        f->FixParameter(10,0); //Background
        f->FixParameter(11,0); //Prefactor for KK and PiPi
        f->FixParameter(12,0); //Mean of PiPi
        f->FixParameter(13,0); //Mean of KK


        f->SetParLimits(2,0.01,0.1);
        f->SetParLimits(3,0.001,0.05);
        f->SetParLimits(4,0,1);
        f->SetParLimits(5,0,1);
        

        f->FixParameter(1,1.8648); //for first few attempt fix mean of gaussian to get reasonable estimation of other pars; later open it up
        h_mc_match_signal->Fit(Form("f_ptbin_%d",i_pt),"q","",fit_range_low,fit_range_high);
        h_mc_match_signal->Fit(Form("f_ptbin_%d",i_pt),"q","",fit_range_low,fit_range_high);
        f->ReleaseParameter(1); //now let gaussian mean float
        h_mc_match_signal->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
        h_mc_match_signal->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
        h_mc_match_signal->Fit(Form("f_ptbin_%d",i_pt),"L q m","",fit_range_low,fit_range_high);

        //now fix signal double gaussian mean, sigma and gaus1,gaus2 yield ratio
        f->FixParameter(1,f->GetParameter(1));
        f->FixParameter(2,f->GetParameter(2));
        f->FixParameter(3,f->GetParameter(3));
        f->FixParameter(4,f->GetParameter(4));
        
        //now release swap bkg parameters to fit signal+swap MC
        f->ReleaseParameter(5);
        f->ReleaseParameter(7);
        //f->ReleaseParameter(8);

        f->SetParameter(7,0.1);

        //fit signal+swap MC
        h_mc_match_all->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
        h_mc_match_all->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
        h_mc_match_all->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
        //h_mc_match_all->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
        //h_mc_match_all->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);

        //now fix swap bkg parameters to fit data
        f->FixParameter(5,f->GetParameter(5));
        f->FixParameter(7,f->GetParameter(7));

        //Make sure all are fixed except background and scaling
        f->FixParameter(1,f->GetParameter(1));
        f->FixParameter(2,f->GetParameter(2));
        f->FixParameter(3,f->GetParameter(3));
        f->FixParameter(4,f->GetParameter(4));
        f->FixParameter(6,0);

        //now release poly bkg pars
        f->ReleaseParameter(8);
        f->ReleaseParameter(9);
        f->ReleaseParameter(10);

        f->SetParameter(8,1);
        f->SetParameter(9,1);
        f->SetParameter(10,1);

        h_mass[i_cen][i_q2][i_pt]->Fit(Form("f_ptbin_%d",i_pt),"q","",fit_range_low,fit_range_high);
        h_mass[i_cen][i_q2][i_pt]->Fit(Form("f_ptbin_%d",i_pt),"q","",fit_range_low,fit_range_high);
        f->ReleaseParameter(1); //allow data to have different mass peak mean than MC
        f->ReleaseParameter(6); //allow data to have different peak width than MC
        f->SetParameter(6,0);
        f->SetParLimits(6,-1,1);
        h_mass[i_cen][i_q2][i_pt]->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
        h_mass[i_cen][i_q2][i_pt]->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
        h_mass[i_cen][i_q2][i_pt]->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
        h_mass[i_cen][i_q2][i_pt]->Fit(Form("f_ptbin_%d",i_pt),"L q m","",fit_range_low,fit_range_high);

        if (i_pt!=0) {
          f->ReleaseParameter(11);
          f->SetParLimits(11,0,f->GetParameter(0)*(1-f->GetParameter(5)));
        }
        else {
          f->FixParameter(11,0);
        }

        h_mass[i_cen][i_q2][i_pt]->Fit(Form("f_ptbin_%d",i_pt),"L q m","",fit_range_low,fit_range_high);
        h_mass[i_cen][i_q2][i_pt]->Fit(Form("f_ptbin_%d",i_pt),"L q m","",fit_range_low,fit_range_high);
        f->SetParLimits(12,-0.03,0.03);
        f->SetParLimits(13,-0.03,0.03);
        
        f->ReleaseParameter(12);
        f->ReleaseParameter(13);
        h_mass[i_cen][i_q2][i_pt]->Fit(Form("f_ptbin_%d",i_pt),"L q m","",fit_range_low,fit_range_high);
        h_mass[i_cen][i_q2][i_pt]->Fit(Form("f_ptbin_%d",i_pt),"L q m","",fit_range_low,fit_range_high);

        h_mass[i_cen][i_q2][i_pt]->GetYaxis()->SetRangeUser(0,1.35*ymax);
        h_mass[i_cen][i_q2][i_pt]->GetYaxis()->SetMaxDigits(2);
        h_mass[i_cen][i_q2][i_pt]->Draw("ep");

        //draw D0 signal separately
        TF1* f1 = new TF1(Form("f_sig_%d",i_pt),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
        f1->SetLineColor(kOrange-3);
        f1->SetLineWidth(1);
        f1->SetLineStyle(2);
        f1->SetFillColorAlpha(kOrange-3,0.3);
        f1->SetFillStyle(1001);
        f1->FixParameter(0,f->GetParameter(0));
        f1->FixParameter(1,f->GetParameter(1));
        f1->FixParameter(2,f->GetParameter(2));
        f1->FixParameter(3,f->GetParameter(3));
        f1->FixParameter(4,f->GetParameter(4));
        f1->FixParameter(5,f->GetParameter(5));
        f1->FixParameter(6,f->GetParameter(6));

        f1->Draw("LSAME");

        //write Yield uncertainty         
        TF1* f1e = new TF1(Form("f_e_%d",i_pt),"[0]", fit_range_low, fit_range_high);
        f1e->SetLineColor(kOrange-3);
        f1e->SetLineWidth(1);
        f1e->SetLineStyle(2);
        f1e->SetFillColorAlpha(kOrange-3,0.3);
        f1e->SetFillStyle(1001);
        f1e->FixParameter(0,(f->GetParError(0)*f->GetParameter(5)+f->GetParError(5)*f->GetParameter(0))/(f->GetParameter(0)*f->GetParameter(5)));

        //draw swap separately
        TF1* f2 = new TF1(Form("f_swap_%d",i_pt),"[0]*((1-[1])*TMath::Gaus(x,[4],[3]*(1.0 +[2]))/(sqrt(2*3.14159)*[3]*(1.0 +[2])))", fit_range_low, fit_range_high);
        f2->SetLineColor(kGreen+4);
        f2->SetLineWidth(1);
        f2->SetLineStyle(1);
        f2->SetFillColorAlpha(kGreen+4,0.3);
        f2->SetFillStyle(1001);
        f2->FixParameter(0,f->GetParameter(0));
        f2->FixParameter(1,f->GetParameter(5));
        f2->FixParameter(2,f->GetParameter(6));
        f2->FixParameter(3,f->GetParameter(7));
        f2->FixParameter(4,f->GetParameter(1));

        f2->Draw("LSAME");

        //draw KK branch separately 
        TF1* f4 = new TF1(Form("f_swap_%d",i_pt),"[0]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[1]), 1.96*(1+[2])) + 4*[0]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[1]), 1.7734*(1+[3])) )", fit_range_low, fit_range_high);
        f4->SetLineColor(kViolet-4);
        f4->SetLineWidth(1);
        f4->SetLineStyle(1);
        f4->SetFillColorAlpha(kViolet-4,0.3);
        f4->SetFillStyle(1001);
        f4->FixParameter(0,f->GetParameter(11));
        f4->FixParameter(1,f->GetParameter(6));
        f4->FixParameter(2,f->GetParameter(12));
        f4->FixParameter(3,f->GetParameter(13));

        f4->Draw("LSAME");

        //draw poly bkg separately
        TF1* f3 = new TF1(Form("f_bkg_%d",i_pt),"[0] + [1]*x + [2]*x*x ", fit_range_low, fit_range_high);
        f3->SetLineColor(4);
        f3->SetLineWidth(1);
        f3->SetLineStyle(2);
        f3->FixParameter(0,f->GetParameter(8));
        f3->FixParameter(1,f->GetParameter(9));
        f3->FixParameter(2,f->GetParameter(10));
        //f3->FixParameter(3,0);

        f3->Draw("LSAME");

        TLegend* leg = new TLegend(0.65,0.6,0.81,0.9,NULL,"brNDC");
        leg->SetBorderSize(0);
        leg->SetTextSize(0.035);
        leg->SetTextFont(42);
        leg->SetFillStyle(0);
        leg->AddEntry(h_mass[i_cen][i_q2][i_pt]," Data","lep");
        leg->AddEntry(f," Fit","L");
        leg->AddEntry(f1," D^{0}+#bar{D^{#lower[0.2]{0}}} Signal","f");
        leg->AddEntry(f2," K-#pi swap","f");
        leg->AddEntry(f4," K-#bar{K}, #pi-#bar{#pi}","f");
        leg->AddEntry(f3," Combinatorial","l");
        leg->Draw("SAME");

        int Yield = f->GetParameter(0)*f->GetParameter(5)/width;
        int YieldEr = (f->GetParError(0)*f->GetParameter(5)+f->GetParError(5)*f->GetParameter(0))/width;

        tex->DrawLatex(0.14,0.75,Form("Yield = %d #pm %d",Yield,YieldEr));
        tex->DrawLatex(0.14,0.80,Form("Significance: %.2f",fabs(1.0*Yield/YieldEr)));

        // per-v2-bin fits & graphs (unchanged)
        for (int i_v2=0; i_v2<N_VBINS; i_v2++) {
          if (h_mass_v2_fit[i_cen][i_q2][i_pt][i_v2]->GetEntries()<10)
          {
            yield_v2[i_v2] = 0; 
            yield_error_v2[i_v2] = 0;
            chi2_ndf_for_q2_v2[i_cen][i_q2][i_pt][i_v2] = -10;
            sigma_v2[i_cen][i_q2][i_pt][i_v2] = -10;
            continue;
          }
          TF1* fitFcn_v2 = new TF1(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v2),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])) ) + [8] + [9]*x + [10]*x*x + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12])) + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])) )",fit_range_low,fit_range_high);

          fitFcn_v2->SetParameter(0,100);
          fitFcn_v2->FixParameter(1,f->GetParameter(1));
          fitFcn_v2->FixParameter(2,f->GetParameter(2));
          fitFcn_v2->FixParameter(3,f->GetParameter(3));
          fitFcn_v2->FixParameter(4,f->GetParameter(4));
          fitFcn_v2->FixParameter(5,f->GetParameter(5));
          fitFcn_v2->FixParameter(6,f->GetParameter(6));
          fitFcn_v2->FixParameter(7,f->GetParameter(7));
          fitFcn_v2->SetParameter(8,1);
          fitFcn_v2->SetParameter(9,1);
          fitFcn_v2->SetParameter(10,1);
          fitFcn_v2->FixParameter(11,0);

          h_mass_v2_fit[i_cen][i_q2][i_pt][i_v2]->Draw("AEP");
          h_mass_v2_fit[i_cen][i_q2][i_pt][i_v2]->Fit(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v2),"M","",fit_range_low,fit_range_high);
          h_mass_v2_fit[i_cen][i_q2][i_pt][i_v2]->Fit(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v2),"L q","",fit_range_low,fit_range_high);
          h_mass_v2_fit[i_cen][i_q2][i_pt][i_v2]->Fit(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v2),"L q","",fit_range_low,fit_range_high);
	  h_mass_v2_fit[i_cen][i_q2][i_pt][i_v2]->Fit(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v2),"L q","",fit_range_low,fit_range_high);
          h_mass_v2_fit[i_cen][i_q2][i_pt][i_v2]->Fit(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v2),"L m","",fit_range_low,fit_range_high);

          chi2_ndf_for_q2_v2[i_cen][i_q2][i_pt][i_v2] = 1.0*fitFcn_v2->GetChisquare()/fitFcn_v2->GetNDF();

          yield_v2[i_v2] = fitFcn_v2->GetParameter(0)*fitFcn_v2->GetParameter(5)/width;
          yield_error_v2[i_v2] = fitFcn_v2->GetParError(0)*fitFcn_v2->GetParameter(5)/width;
          sigma_v2[i_cen][i_q2][i_pt][i_v2] = yield_v2[i_v2]/yield_error_v2[i_v2];

          if (yield_v2[i_v2] <= 0) {
            yield_v2[i_v2] = 0; 
            yield_error_v2[i_v2] = 0;
            chi2_ndf_for_q2_v2[i_cen][i_q2][i_pt][i_v2] = -10;
            sigma_v2[i_cen][i_q2][i_pt][i_v2] = -10;
          }
          v2_x[i_v2] = (vnbinning[i_cen][i_pt][i_v2]+vnbinning[i_cen][i_pt][i_v2+1])/2;
          v2_x_err[i_v2] = fabs(vnbinning[i_cen][i_pt][i_v2+1]-vnbinning[i_cen][i_pt][i_v2])/2;

          h_v2[i_cen][i_q2][i_pt]->SetBinContent(i_v2+1,yield_v2[i_v2]);
          h_v2[i_cen][i_q2][i_pt]->SetBinError(i_v2+1,yield_error_v2[i_v2]);
          delete fitFcn_v2;
        }// loop over v2 bins

        yield_v2_graph[i_cen][i_q2][i_pt] = new TGraphErrors(N_VBINS,v2_x,yield_v2,v2_x_err,yield_error_v2);
        yield_v2_graph[i_cen][i_q2][i_pt]->SetNameTitle(("v2_graph_" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str(), ("v2_graph_" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str());
    
	//yield_v2_graph[i_cen][i_q2][i_pt]->SetNameTitle(("v2_graph_" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str(), ("v2_graph_" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str());


	// --- per-v3-bin fits & graphs (mirror v2) ---
        for (int i_v3=0; i_v3<N_VBINS; i_v3++) {
          if (h_mass_v3_fit[i_cen][i_q2][i_pt][i_v3]->GetEntries()<10)
          {
            yield_v3[i_v3] = 0; 
            yield_error_v3[i_v3] = 0;
            chi2_ndf_for_q2_v3[i_cen][i_q2][i_pt][i_v3] = -10;
            sigma_v3[i_cen][i_q2][i_pt][i_v3] = -10;
            continue;
          }
          TF1* fitFcn_v3 = new TF1(Form("f_ptbin_%d_v3bin_%d",i_pt,i_v3),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])) ) + [8] + [9]*x + [10]*x*x + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12])) + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])) )",fit_range_low,fit_range_high);
	  fitFcn_v3->SetParameter(0,100);
          fitFcn_v3->FixParameter(1,f->GetParameter(1));
          fitFcn_v3->FixParameter(2,f->GetParameter(2));
          fitFcn_v3->FixParameter(3,f->GetParameter(3));
          fitFcn_v3->FixParameter(4,f->GetParameter(4));
          fitFcn_v3->FixParameter(5,f->GetParameter(5));
          fitFcn_v3->FixParameter(6,f->GetParameter(6));
          fitFcn_v3->FixParameter(7,f->GetParameter(7));
          fitFcn_v3->SetParameter(8,1);
          fitFcn_v3->SetParameter(9,1);
          fitFcn_v3->SetParameter(10,1);
          fitFcn_v3->FixParameter(11,0);

          h_mass_v3_fit[i_cen][i_q2][i_pt][i_v3]->Draw("AEP");
          h_mass_v3_fit[i_cen][i_q2][i_pt][i_v3]->Fit(Form("f_ptbin_%d_v3bin_%d",i_pt,i_v3),"M","",fit_range_low,fit_range_high);
          h_mass_v3_fit[i_cen][i_q2][i_pt][i_v3]->Fit(Form("f_ptbin_%d_v3bin_%d",i_pt,i_v3),"L q","",fit_range_low,fit_range_high);
	  h_mass_v3_fit[i_cen][i_q2][i_pt][i_v3]->Fit(Form("f_ptbin_%d_v3bin_%d",i_pt,i_v3),"L q","",fit_range_low,fit_range_high);
          h_mass_v3_fit[i_cen][i_q2][i_pt][i_v3]->Fit(Form("f_ptbin_%d_v3bin_%d",i_pt,i_v3),"L m","",fit_range_low,fit_range_high);

          chi2_ndf_for_q2_v3[i_cen][i_q2][i_pt][i_v3] = 1.0*fitFcn_v3->GetChisquare()/fitFcn_v3->GetNDF();

          yield_v3[i_v3] = fitFcn_v3->GetParameter(0)*fitFcn_v3->GetParameter(5)/width;
          yield_error_v3[i_v3] = fitFcn_v3->GetParError(0)*fitFcn_v3->GetParameter(5)/width;
          sigma_v3[i_cen][i_q2][i_pt][i_v3] = yield_v3[i_v3]/yield_error_v3[i_v3];

          if (yield_v3[i_v3] <= 0) {
            yield_v3[i_v3] = 0; 
            yield_error_v3[i_v3] = 0;
            chi2_ndf_for_q2_v3[i_cen][i_q2][i_pt][i_v3] = -10;
            sigma_v3[i_cen][i_q2][i_pt][i_v3] = -10;
          }
          v3_x[i_v3] = (vnbinning_v3[i_cen][i_pt][i_v3]+vnbinning_v3[i_cen][i_pt][i_v3+1])/2;
          v3_x_err[i_v3] = fabs(vnbinning_v3[i_cen][i_pt][i_v3+1]-vnbinning_v3[i_cen][i_pt][i_v3])/2;

          h_v3[i_cen][i_q2][i_pt]->SetBinContent(i_v3+1,yield_v3[i_v3]);
          h_v3[i_cen][i_q2][i_pt]->SetBinError(i_v3+1,yield_error_v3[i_v3]);
          delete fitFcn_v3;
        } // loop over v3 bins

        yield_v3_graph[i_cen][i_q2][i_pt] = new TGraphErrors(N_VBINS,v3_x,yield_v3,v3_x_err,yield_error_v3);
        yield_v3_graph[i_cen][i_q2][i_pt]->SetNameTitle(("v3_graph_" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str(), ("v3_graph_" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str());

        // end v3
        //yield_v2_graph[i_cen][i_q2][i_pt]->SetNameTitle("v2_graph_"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2),"v2_graph_"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2));



    
      }// loop over q2 bins
    }//loop over pT bins
  }//loop over centrality bins

  Double_t pt_x[N_PTBINS], pt_x_error[N_PTBINS];
  TFile *outf = new TFile("prompt_vn_in_q2HF_total_graph_fullrange_withchi2_sigma_vs_v2v3.root","recreate");
  outf->cd();

  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    for (int i_q2=0; i_q2<N_QBINS; i_q2++){
      for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
        pt_x[i_pt]=(pt_edges[i_pt]+pt_edges[i_pt+1])/2;
        pt_x_error[i_pt]=(-1.0*pt_edges[i_pt]+pt_edges[i_pt+1])/2;
        mean_val[i_pt] = h_v2[i_cen][i_q2][i_pt]->GetMean();
        mean_error[i_pt] = h_v2[i_cen][i_q2][i_pt]->GetMeanError();
        h_v2[i_cen][i_q2][i_pt]->Write();
        h_v3[i_cen][i_q2][i_pt]->Write();

        chi2_ndf_for_q2_v2_graph[i_cen][i_q2][i_pt] = new TGraph(N_VBINS,v2_x,chi2_ndf_for_q2_v2[i_cen][i_q2][i_pt]);
	chi2_ndf_for_q2_v2_graph[i_cen][i_q2][i_pt]->SetNameTitle(
								  ("chi2ndf_graph_v2_pt" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str(),
								  ("chi2ndf_graph_v2_pt" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str()
								  );
	chi2_ndf_for_q2_v2_graph[i_cen][i_q2][i_pt]->Write();
        chi2_ndf_for_q2_v2_graph[i_cen][i_q2][i_pt]->GetXaxis()->SetTitle("v2");
        chi2_ndf_for_q2_v2_graph[i_cen][i_q2][i_pt]->GetYaxis()->SetTitle("chi2/ndf");
        chi2_ndf_for_q2_v2_graph[i_cen][i_q2][i_pt]->Draw("ACP");
        sigma_v2_fit[i_cen][i_q2][i_pt] = new TGraph(N_VBINS,v2_x,sigma_v2[i_cen][i_q2][i_pt]);
	sigma_v2_fit[i_cen][i_q2][i_pt]->SetNameTitle(
						      ("sigma_v2_graph_pt" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str(),
						      ("sigma_v2_graph_pt" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str()
						      );
	sigma_v2_fit[i_cen][i_q2][i_pt]->Write();
        sigma_v2_fit[i_cen][i_q2][i_pt]->GetXaxis()->SetTitle("v2");
        sigma_v2_fit[i_cen][i_q2][i_pt]->GetYaxis()->SetTitle("sigma");
        sigma_v2_fit[i_cen][i_q2][i_pt]->Draw("ACP");
	
        // v3 chi2 & sigma graphs
        chi2_ndf_for_q2_v3_graph[i_cen][i_q2][i_pt] = new TGraph(N_VBINS,v3_x,chi2_ndf_for_q2_v3[i_cen][i_q2][i_pt]);
    chi2_ndf_for_q2_v3_graph[i_cen][i_q2][i_pt]->SetNameTitle(
							      ("chi2ndf_graph_v3_pt" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str(),
							      ("chi2ndf_graph_v3_pt" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str()
							      );

    chi2_ndf_for_q2_v3_graph[i_cen][i_q2][i_pt]->Write();
        chi2_ndf_for_q2_v3_graph[i_cen][i_q2][i_pt]->GetXaxis()->SetTitle("v3");
        chi2_ndf_for_q2_v3_graph[i_cen][i_q2][i_pt]->GetYaxis()->SetTitle("chi2/ndf");
        chi2_ndf_for_q2_v3_graph[i_cen][i_q2][i_pt]->Draw("ACP");
        sigma_v3_fit[i_cen][i_q2][i_pt] = new TGraph(N_VBINS,v3_x,sigma_v3[i_cen][i_q2][i_pt]);
 sigma_v3_fit[i_cen][i_q2][i_pt]->SetNameTitle(
    ("sigma_v3_graph_pt" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str(),
    ("sigma_v3_graph_pt" + pt_name[i_pt] + "_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str()
);
    sigma_v3_fit[i_cen][i_q2][i_pt]->Write();
        sigma_v3_fit[i_cen][i_q2][i_pt]->GetXaxis()->SetTitle("v3");
        sigma_v3_fit[i_cen][i_q2][i_pt]->GetYaxis()->SetTitle("sigma");
        sigma_v3_fit[i_cen][i_q2][i_pt]->Draw("ACP");
      }
      v2_cen[i_cen][i_q2] = new TGraphErrors(N_PTBINS,pt_x,mean_val,pt_x_error,mean_error);
v2_cen[i_cen][i_q2]->SetNameTitle(
    ("v2_graph_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str(),
    ("v2_graph_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str()
);

      v2_cen[i_cen][i_q2]->Write();

      // write v3 pt graphs (mean of yield vs pT) - use h_v3 mean
      for (int ip=0; ip<N_PTBINS; ++ip) {
        mean_val[ip] = h_v3[i_cen][i_q2][ip]->GetMean();
        mean_error[ip] = h_v3[i_cen][i_q2][ip]->GetMeanError();
      }
      v3_cen[i_cen][i_q2] = new TGraphErrors(N_PTBINS,pt_x,mean_val,pt_x_error,mean_error);
      v3_cen[i_cen][i_q2]->SetNameTitle(
          ("v3_graph_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str(),
          ("v3_graph_cen" + cen_name[i_cen] + "_q2bin_" + to_string(i_q2)).c_str()
      );
      v3_cen[i_cen][i_q2]->Write();
    }
  }

  outf->Write(0,TObject::kOverwrite);
  outf->Close();

  return 0;
}





int main(int argc, char *argv[])
{
  if (argc==1) {
    fit_v2_prompt_dca();
  } else if (argc==2) {
    int idx = atoi(argv[1]);
    fit_v2_prompt_dca(idx);
  } else if (argc==3) {
    int cen_idx = atoi(argv[1]);
    int pt_idx  = atoi(argv[2]);
    fit_v2_prompt_dca(cen_idx, pt_idx);
  } else {
    std::cout << "Usage: ./fit_v2_prompt_dca [target1pct_index] [target_pt_index]" << std::endl;
    return 1;
  }
  return 0;
}
