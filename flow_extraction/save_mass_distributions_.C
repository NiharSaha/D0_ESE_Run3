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
#include <algorithm> //
#include <iomanip>

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/ESE_Cuts_MB11to21_Jan22.h"

using namespace std;




// add this helper to print non-increasing bin edges (move to global scope)
static bool check_bin_edges(const Double_t *edges, int nBins, const char* tag, int ic, int iq, int ip) {
  for (int i=1; i<=nBins; ++i) {
    if (!(edges[i] > edges[i-1])) {
      std::cout << "[DEBUG] Non-increasing edges for " << tag
                << " cen=" << ic << " q2=" << iq << " pt=" << ip
                << " at index " << (i-1) << "->" << i << " : "
                << edges[i-1] << " , " << edges[i] << std::endl;
      int s = std::max(0, i-3), e = std::min(nBins, i+3);
      std::cout << "[DEBUG] edges[" << s << ":" << e << "] =";
      for (int j=s; j<=e; ++j) std::cout << " " << edges[j];
      std::cout << std::endl;
      return false;
    }
  }
  return true;
}

// ensure bin edges are strictly increasing (fix tiny duplicates/zeros)
static void fix_bin_edges(Double_t *edges, int nBins, Double_t eps = 1e-12) {
  for (int i = 1; i <= nBins; ++i) {
    if (!(edges[i] > edges[i-1])) {
      edges[i] = edges[i-1] + eps;
      //std::cout << "[DEBUG] Fixed bin edge at idx " << i << " to " << edges[i] << std::endl;
    }
  }
}


int save_mass_distributions(int target1pct = -1, int target_pt = -1)
{
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
  const int N_BINS = 2300;
  const int N_EDGES = 10;
  Double_t vnbinning_v2[N_CENTBINS][N_PTBINS][N_VBINS+1];
  Double_t vnbinning_v3[N_CENTBINS][N_PTBINS][N_VBINS+1];
  Double_t vnbinning_side[16]={-10,-8.2,-6.6,-5.2,-4.0,-3.5,-3.0,-2.3,2.3,3.0,3.5,4.0,5.2,6.6,8.2,10};
  // readiness flags: true if vnbinning for (cent,pT) was computed in the first loop
  bool vnb_ready[N_CENTBINS][N_PTBINS] = {};

  const int N_MASSBINS = 52;
  Float_t fit_range_low = 1.74;
  Float_t fit_range_high = 2.00;
  Float_t width = (fit_range_high-fit_range_low)/N_MASSBINS;
  Double_t D0_mass = 1.8648;







  auto nt_Data = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/SP_MB11to21_Jan27/ROOT/flow_Analysis_out_combined.root");
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


  // determine which 10% centrality group we need (if any)
  int needed_cen_group = -1;
  if (target1pct >= 0) needed_cen_group = get_cent_group(target1pct);
  
  
  TString cuts;

  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    if (needed_cen_group >= 0 && i_cen != needed_cen_group) continue; // skip other 10% groups
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
      if (target_pt >= 0 && i_pt != target_pt) continue; // skip other pT bins


      /*cuts = "dca<0.0085 && mass<2.0 && mass>1.74 && centrality>="+to_string(cen_edges[i_cen])+" && centrality<"+to_string(cen_edges[i_cen+1])+" && pT>="+to_string(pt_edges[i_pt])+" && pT<"+to_string(pt_edges[i_pt+1]);

      cout << endl << endl << cuts << endl;  
      h_v2_bin->Reset();
      input_nt->Project(h_v2_bin->GetName(),"v2",cuts);
      Int_t n_entries_v2 = static_cast<Int_t>(h_v2_bin->GetEntries());
      cout << endl << "Total entries v2: " << n_entries_v2 << endl;

      h_v3_bin->Reset();
      input_nt->Project(h_v3_bin->GetName(),"v3",cuts);
      Int_t n_entries_v3 = static_cast<Int_t>(h_v3_bin->GetEntries());
      cout << endl << "Total entries v3: " << n_entries_v3 << endl;

      for (int ib=0; ib<=N_VBINS; ++ib){
	vnbinning_v2[i_cen][i_pt][ib] = 0.0;
	vnbinning_v3[i_cen][i_pt][ib] = 0.0;
      }
      for (int i_bin=0; i_bin<8; i_bin++) {
        vnbinning_v2[i_cen][i_pt][i_bin] = vnbinning_side[i_bin];
        vnbinning_v2[i_cen][i_pt][i_bin+27] = vnbinning_side[i_bin+8];
	vnbinning_v3[i_cen][i_pt][i_bin] = vnbinning_side[i_bin];
        vnbinning_v3[i_cen][i_pt][i_bin+27] = vnbinning_side[i_bin+8];
      }
      vnbinning_v2[i_cen][i_pt][17]=0.0;
      vnbinning_v3[i_cen][i_pt][17]=0.0;
      
      if (n_entries_v2 > 0) {
        std::vector<double> cum(N_BINS+1,0.0);
        for (int b=1; b<=N_BINS; ++b) cum[b] = cum[b-1] + h_v2_bin->GetBinContent(b);
        for (int i_edge=1; i_edge<=N_EDGES; ++i_edge) {
          double target = i_edge * (n_entries_v2 / (double)N_EDGES);
          int j = 1;
          while (j <= N_BINS && cum[j] <= target) ++j;
          if (j <= N_BINS) {
            double edge = h_v2_bin->GetBinLowEdge(j);
            vnbinning_v2[i_cen][i_pt][i_edge+7] = edge;
            vnbinning_v2[i_cen][i_pt][27-i_edge] = -edge;
          } else {
            double edge = h_v2_bin->GetBinLowEdge(N_BINS+1);
            vnbinning_v2[i_cen][i_pt][i_edge+7] = edge;
            vnbinning_v2[i_cen][i_pt][27-i_edge] = -edge;
          }
        }
      }
      if (n_entries_v3 > 0) {
        std::vector<double> cum3(N_BINS+1,0.0);
        for (int b=1; b<=N_BINS; ++b) cum3[b] = cum3[b-1] + h_v3_bin->GetBinContent(b);
        for (int i_edge=1; i_edge<=N_EDGES; ++i_edge) {
          double target3 = i_edge * (n_entries_v3 / (double)N_EDGES);
          int j3 = 1;
          while (j3 <= N_BINS && cum3[j3] <= target3) ++j3;
          if (j3 <= N_BINS) {
            double edge3 = h_v3_bin->GetBinLowEdge(j3);
            vnbinning_v3[i_cen][i_pt][i_edge+7] = edge3;
            vnbinning_v3[i_cen][i_pt][27-i_edge] = -edge3;
          } else {
            double edge3 = h_v3_bin->GetBinLowEdge(N_BINS+1);
            vnbinning_v3[i_cen][i_pt][i_edge+7] = edge3;
            vnbinning_v3[i_cen][i_pt][27-i_edge] = -edge3;
          }
        }
	}*/



      static const Double_t manual_vnb[45] = {
	-10.00, -8.50, -7.00, -6.00, -5.00, -4.00, -3.20, -2.50, -2.00, -1.80,
	-1.60, -1.40, -1.20, -1.00, -0.80, -0.60, -0.40, -0.30, -0.20, -0.15,
	-0.10, -0.05,  0.00,  0.05,  0.10,  0.15,  0.20,  0.30,  0.40,  0.60,
	0.80,  1.00,  1.20,  1.40,  1.60,  1.80,  2.00,  2.50,  3.20,  4.00,
	5.00,  6.00,  7.00,  8.50, 10.00
      };
      for (int ib=0; ib<=N_VBINS; ++ib) {
	vnbinning_v2[i_cen][i_pt][ib] = manual_vnb[ib];
	vnbinning_v3[i_cen][i_pt][ib] = manual_vnb[ib];
      }
      
      // mark this (cent,pT) as ready (we computed/initialized vnbinning for it)
      vnb_ready[i_cen][i_pt] = true;

      // fix possible non-increasing or duplicated edges around zero
      //fix_bin_edges(vnbinning_v2[i_cen][i_pt], N_VBINS);
      //fix_bin_edges(vnbinning_v3[i_cen][i_pt], N_VBINS);


      /*if ((needed_cen_group < 0 || i_cen == needed_cen_group) && (target_pt < 0 || i_pt == target_pt)) {
	cout << "vnbinning_v2 cen=" << i_cen << " pt=" << i_pt << " :";
	cout.setf(std::ios::fixed); cout << setprecision(8);
	for (int ib = 0; ib <= N_VBINS; ++ib) cout << " " << vnbinning_v2[i_cen][i_pt][ib];
	cout << endl;
	cout << "vnbinning_v3 cen=" << i_cen << " pt=" << i_pt << " :";
	for (int ib = 0; ib <= N_VBINS; ++ib) cout << " " << vnbinning_v3[i_cen][i_pt][ib];
	cout << endl;
	}*/
      


    }
  }

  // histograms declarations
  TH1D *h_mass_default[N_CENTBINS][N_PTBINS]={};
  TH1D *h_mass_default_1percent_cen[N_CENTBINS_1][N_PTBINS]={};

  //For q2/v2
  TH1D *h_mass_q2[N_CENTBINS][N_QBINS][N_PTBINS]={};
  TH1D *h_mass_1pecent_cen_q2[N_CENTBINS_1][N_QBINS][N_PTBINS]={};
  TH1D *h_mass_v2_fit[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS]={};
  TH1D *h_mass_v2_fit_1pecent_cen[N_CENTBINS_1][N_QBINS][N_PTBINS][N_VBINS]={};
  TH1D *h_v2[N_CENTBINS][N_QBINS][N_PTBINS]={};
  //For q3/v3
  TH1D *h_mass_q3[N_CENTBINS][N_QBINS][N_PTBINS]={};
  TH1D *h_mass_1pecent_cen_q3[N_CENTBINS_1][N_QBINS][N_PTBINS]={};
  TH1D *h_mass_v3_fit[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS]={};
  TH1D *h_mass_v3_fit_1pecent_cen[N_CENTBINS_1][N_QBINS][N_PTBINS][N_VBINS]={};
  TH1D *h_v3[N_CENTBINS][N_QBINS][N_PTBINS]={};

  TString h_mass_name, h_v2_name, h_v3_name,  h_mass_v2_name, h_mass_v3_name,  h_name_0, h_name_1;

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
        h_mass_1pecent_cen_q2[i_cen][i_q2][i_pt] = new TH1D(h_name_0,h_name_0,N_MASSBINS,fit_range_low,fit_range_high);

	h_name_0 = "h_mass_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_q3bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt];
	h_mass_1pecent_cen_q3[i_cen][i_q2][i_pt] = new TH1D(h_name_0,h_name_0,N_MASSBINS,fit_range_low,fit_range_high);

        for(int i_v2=0; i_v2<N_VBINS; i_v2++)
        {
          // use cen_group_for_binning (correct mapping) instead of i_cen/10
          h_name_1 = "hist_mass_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_q2bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt]+"_in_v2bin_idx_"+to_string(i_v2);
          h_mass_v2_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2] = new TH1D(h_name_1,h_name_1,N_MASSBINS,fit_range_low,fit_range_high);

	  h_name_1 = "hist_mass_v3_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_q2bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt]+"_in_v3bin_idx_"+to_string(i_v2);
          h_mass_v3_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2] = new TH1D(h_name_1,h_name_1,N_MASSBINS,fit_range_low,fit_range_high);
        }
      }
    }
  }

  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
      h_mass_name = "h_mass_"+cen_name[i_cen]+"_"+pt_name[i_pt]+"_def";
      h_mass_default[i_cen][i_pt] = new TH1D(h_mass_name, h_mass_name, N_MASSBINS, fit_range_low, fit_range_high);

      for (int i_q2=0; i_q2<N_QBINS; i_q2++) {
        h_mass_name = "h_mass_"+cen_name[i_cen]+"_q2bin"+to_string(i_q2)+"_"+pt_name[i_pt];
        h_mass_q2[i_cen][i_q2][i_pt] = new TH1D(h_mass_name, h_mass_name, N_MASSBINS, fit_range_low, fit_range_high);

	h_mass_name = "h_mass_"+cen_name[i_cen]+"_q3bin"+to_string(i_q2)+"_"+pt_name[i_pt];
	h_mass_q3[i_cen][i_q2][i_pt] = new TH1D(h_mass_name, h_mass_name, N_MASSBINS, fit_range_low, fit_range_high);

	h_v2_name = "hist_v2_"+cen_name[i_cen]+"_q2bin"+to_string(i_q2)+"_"+pt_name[i_pt];
        if (vnb_ready[i_cen][i_pt]) {
          if (!check_bin_edges(vnbinning_v2[i_cen][i_pt], N_VBINS, "v2", i_cen, i_q2, i_pt)) {
            h_v2[i_cen][i_q2][i_pt] = nullptr;
          } else {
            h_v2[i_cen][i_q2][i_pt] = new TH1D(h_v2_name, h_v2_name, N_VBINS, vnbinning_v2[i_cen][i_pt]);
          }
        } else {
          h_v2[i_cen][i_q2][i_pt] = nullptr;
        }
        h_v3_name = "hist_v3_"+cen_name[i_cen]+"_q2bin"+to_string(i_q2)+"_"+pt_name[i_pt];
        if (vnb_ready[i_cen][i_pt]) {
          if (!check_bin_edges(vnbinning_v3[i_cen][i_pt], N_VBINS, "v3", i_cen, i_q2, i_pt)) {
            h_v3[i_cen][i_q2][i_pt] = nullptr;
          } else {
            h_v3[i_cen][i_q2][i_pt] = new TH1D(h_v3_name, h_v3_name, N_VBINS, vnbinning_v3[i_cen][i_pt]);
          }
        } else {
          h_v3[i_cen][i_q2][i_pt] = nullptr;
        }
        for(int i_v2=0; i_v2<N_VBINS; i_v2++)
        {
          // use indices in the name to avoid duplicate/NaN edge strings
          h_mass_v2_name = "hist_mass_"+cen_name[i_cen]+"_q2bin"+to_string(i_q2)+"_"+pt_name[i_pt]+"_in_v2bin_idx_"+to_string(i_v2);
          h_mass_v2_fit[i_cen][i_q2][i_pt][i_v2]=new TH1D(h_mass_v2_name, h_mass_v2_name, N_MASSBINS, fit_range_low, fit_range_high);

          h_mass_v3_name = "hist_mass_v3_"+cen_name[i_cen]+"_q2bin"+to_string(i_q2)+"_"+pt_name[i_pt]+"_in_v3bin_idx_"+to_string(i_v2);
          h_mass_v3_fit[i_cen][i_q2][i_pt][i_v2] = new TH1D(h_mass_v3_name, h_mass_v3_name, N_MASSBINS, fit_range_low, fit_range_high);
        }
      }
    }
  }

  // Read and fill histograms from ntuples
  Long64_t N_ENTRIES = input_nt->GetEntriesFast();
  cout << "Reading ntuple to generate histograms" << endl;
  //Long64_t maxEvents = 1000; // small default; user can change after editing
  //if (maxEvents > 0 && N_ENTRIES > maxEvents) N_ENTRIES = maxEvents;


  for(Long64_t i_entry=0; i_entry<N_ENTRIES; i_entry++) {
    input_nt->GetEntry(i_entry);
    if(i_entry%100==0) cout <<i_entry<<" / "<<N_ENTRIES<< "  "<<100*i_entry/N_ENTRIES<<"%"<<endl;

    if(dca_val>=0.0085) continue;

    for (int i_cen=0; i_cen<N_CENTBINS_1; i_cen++) {

      if (target1pct >= 0 && i_cen != target1pct) continue;
      if(cen_val<(i_cen) || cen_val>=(i_cen+1)) continue;
      int cen_group = get_cent_group(i_cen);
      if (cen_group < 0) continue;

      //std::cout<<"CHECK-1"<<std::endl;
      for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
        if (target_pt >= 0 && i_pt != target_pt) continue;
        if (pT_val<pt_edges[i_pt] || pT_val>=pt_edges[i_pt+1]) continue;

	if (h_mass_default_1percent_cen[i_cen][i_pt]) {
          h_mass_default_1percent_cen[i_cen][i_pt]->Fill(mass_val);
        }
	
	//h_mass_default_1percent_cen[i_cen][i_pt]->Fill(mass_val);


	
	for (int iq=0; iq<N_QBINS; iq++) {

	  if (q2_val>=q2_cuts[i_cen][iq] && q2_val<q2_cuts[i_cen][iq+1]) {
            if (h_mass_1pecent_cen_q2[i_cen][iq][i_pt]){h_mass_1pecent_cen_q2[i_cen][iq][i_pt]->Fill(mass_val);}

	    for(int i_v2=0; i_v2<N_VBINS; i_v2++) {
              if (vnb_ready[cen_group][i_pt]) {
                if(v2_val>=vnbinning_v2[cen_group][i_pt][i_v2] && v2_val<vnbinning_v2[cen_group][i_pt][i_v2+1]) {
                  if (h_mass_v2_fit_1pecent_cen[i_cen][iq][i_pt][i_v2]) h_mass_v2_fit_1pecent_cen[i_cen][iq][i_pt][i_v2]->Fill(mass_val);
                }
              }
            }
          }
	  //std::cout<<"CHECK-2"<<std::endl;
	  if (q3_val>=q3_cuts[i_cen][iq] && q3_val<q3_cuts[i_cen][iq+1]) {
            if (h_mass_1pecent_cen_q3[i_cen][iq][i_pt]){h_mass_1pecent_cen_q3[i_cen][iq][i_pt]->Fill(mass_val);}

	    for(int i_v3=0; i_v3<N_VBINS; i_v3++) {
              if (vnb_ready[cen_group][i_pt]) {
                if(v3_val>=vnbinning_v3[cen_group][i_pt][i_v3] && v3_val<vnbinning_v3[cen_group][i_pt][i_v3+1]) {
                  if (h_mass_v3_fit_1pecent_cen[i_cen][iq][i_pt][i_v3]) h_mass_v3_fit_1pecent_cen[i_cen][iq][i_pt][i_v3]->Fill(mass_val);
                }
              }
            }
          } 
        }// -- iq loop --
	//std::cout<<"CHECK-3"<<std::endl;
	
      }//- -- pT loop ---
    }// --- CENT loop ---
  }// --- EVT loop ---

  //Add the histograms of 1% centrality bins to get the 10% centrality bins
  if (target1pct < 0 && target_pt < 0 ) {
    for(int i_pt=0; i_pt<N_PTBINS; i_pt++) {
      for(int i_cen=0; i_cen<N_CENTBINS_1; i_cen++) {
        int cen_group = get_cent_group(i_cen);
        if (cen_group < 0) continue;
        h_mass_default[cen_group][i_pt]->Add(h_mass_default_1percent_cen[i_cen][i_pt]);
        for(int iq=0; iq<N_QBINS; iq++) {
          h_mass_q2[cen_group][iq][i_pt]->Add(h_mass_1pecent_cen_q2[i_cen][iq][i_pt]);
	  h_mass_q3[cen_group][iq][i_pt]->Add(h_mass_1pecent_cen_q3[i_cen][iq][i_pt]);
          for(int i_v2=0; i_v2<N_VBINS; i_v2++) {
            h_mass_v2_fit[cen_group][iq][i_pt][i_v2]->Add(h_mass_v2_fit_1pecent_cen[i_cen][iq][i_pt][i_v2]);
            h_mass_v3_fit[cen_group][iq][i_pt][i_v2]->Add(h_mass_v3_fit_1pecent_cen[i_cen][iq][i_pt][i_v2]);
          }
        }
      }
    }
  }

  // Save histograms to file
  // choose output filename per-request to allow many jobs in parallel
  std::string outname = "mass_distributions";
  if (target1pct >= 0) outname += "_1pct_"+to_string(target1pct);
  if (target_pt >= 0)   outname += "_pt_"+to_string(target_pt);
  outname += ".root";

  // small diagnostic when saving: print entries for each histogram and warn if empty
  const bool VERBOSE_SAVE = true;
  TFile *outf = new TFile(outname.c_str(),"recreate");

  // write aggregated 10% histograms organized as /cent_<name>/pt_<name>/
  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    TString cen_dir = TString::Format("cent_%s", cen_name[i_cen].c_str());
    TDirectory *d_cen = outf->mkdir(cen_dir);
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
      TString pt_dir = TString::Format("%s/pt_%s", cen_dir.Data(), pt_name[i_pt].c_str());
      TDirectory *d_pt = outf->mkdir(pt_dir);
      d_pt->cd();
      if (h_mass_default[i_cen][i_pt]) {
        if (VERBOSE_SAVE) cout<<"WRITE: "<<h_mass_default[i_cen][i_pt]->GetName()<<" entries="<<h_mass_default[i_cen][i_pt]->GetEntries()<<endl;
        if (h_mass_default[i_cen][i_pt]->GetEntries()==0) cout<<"[WARN] empty: "<<h_mass_default[i_cen][i_pt]->GetName()<<endl;
        h_mass_default[i_cen][i_pt]->Write();
      }
      for (int iq=0; iq<N_QBINS; iq++) {
        if (h_mass_q2[i_cen][iq][i_pt]) {
          if (VERBOSE_SAVE) cout<<"WRITE: "<<h_mass_q2[i_cen][iq][i_pt]->GetName()<<" entries="<<h_mass_q2[i_cen][iq][i_pt]->GetEntries()<<endl;
          if (h_mass_q2[i_cen][iq][i_pt]->GetEntries()==0) cout<<"[WARN] empty: "<<h_mass_q2[i_cen][iq][i_pt]->GetName()<<endl;
          h_mass_q2[i_cen][iq][i_pt]->Write();
        }
        if (h_mass_q3[i_cen][iq][i_pt]) {
          if (VERBOSE_SAVE) cout<<"WRITE: "<<h_mass_q3[i_cen][iq][i_pt]->GetName()<<" entries="<<h_mass_q3[i_cen][iq][i_pt]->GetEntries()<<endl;
          if (h_mass_q3[i_cen][iq][i_pt]->GetEntries()==0) cout<<"[WARN] empty: "<<h_mass_q3[i_cen][iq][i_pt]->GetName()<<endl;
          h_mass_q3[i_cen][iq][i_pt]->Write();
        }
        if (h_v2[i_cen][iq][i_pt]) {
          if (VERBOSE_SAVE) cout<<"WRITE: "<<h_v2[i_cen][iq][i_pt]->GetName()<<" entries="<<h_v2[i_cen][iq][i_pt]->GetEntries()<<endl;
          if (h_v2[i_cen][iq][i_pt]->GetEntries()==0) cout<<"[WARN] empty: "<<h_v2[i_cen][iq][i_pt]->GetName()<<endl;
          h_v2[i_cen][iq][i_pt]->Write();
        }
        if (h_v3[i_cen][iq][i_pt]) {
          if (VERBOSE_SAVE) cout<<"WRITE: "<<h_v3[i_cen][iq][i_pt]->GetName()<<" entries="<<h_v3[i_cen][iq][i_pt]->GetEntries()<<endl;
          if (h_v3[i_cen][iq][i_pt]->GetEntries()==0) cout<<"[WARN] empty: "<<h_v3[i_cen][iq][i_pt]->GetName()<<endl;
          h_v3[i_cen][iq][i_pt]->Write();
        }
        // write per-v bin mass histograms inside the pt directory
        for (int iv=0; iv<N_VBINS; ++iv) {
          if (h_mass_v2_fit[i_cen][iq][i_pt][iv]) {
            if (VERBOSE_SAVE) cout<<"WRITE: "<<h_mass_v2_fit[i_cen][iq][i_pt][iv]->GetName()<<" entries="<<h_mass_v2_fit[i_cen][iq][i_pt][iv]->GetEntries()<<endl;
            if (h_mass_v2_fit[i_cen][iq][i_pt][iv]->GetEntries()==0) cout<<"[WARN] empty: "<<h_mass_v2_fit[i_cen][iq][i_pt][iv]->GetName()<<endl;
            h_mass_v2_fit[i_cen][iq][i_pt][iv]->Write();
          }
          if (h_mass_v3_fit[i_cen][iq][i_pt][iv]) {
            if (VERBOSE_SAVE) cout<<"WRITE: "<<h_mass_v3_fit[i_cen][iq][i_pt][iv]->GetName()<<" entries="<<h_mass_v3_fit[i_cen][iq][i_pt][iv]->GetEntries()<<endl;
            if (h_mass_v3_fit[i_cen][iq][i_pt][iv]->GetEntries()==0) cout<<"[WARN] empty: "<<h_mass_v3_fit[i_cen][iq][i_pt][iv]->GetName()<<endl;
            h_mass_v3_fit[i_cen][iq][i_pt][iv]->Write();
          }
        }
      }
      outf->cd();
    }
  }

  // write 1% histograms under /1pct/cent_<i>/pt_<j>/
  TDirectory *d_1pct = outf->mkdir("1pct");
  for (int i_cen=0; i_cen<N_CENTBINS_1; ++i_cen) {
    TString cen1_dir = TString::Format("1pct/cent_%02d", i_cen);
    for (int i_pt=0; i_pt<N_PTBINS; ++i_pt) {
      TString pt1_dir = TString::Format("%s/pt_%s", cen1_dir.Data(), pt_name[i_pt].c_str());
      TDirectory *d_pt1 = outf->mkdir(pt1_dir);
      d_pt1->cd();
      if (h_mass_default_1percent_cen[i_cen][i_pt]) {
        if (VERBOSE_SAVE) cout<<"WRITE: "<<h_mass_default_1percent_cen[i_cen][i_pt]->GetName()<<" entries="<<h_mass_default_1percent_cen[i_cen][i_pt]->GetEntries()<<endl;
        if (h_mass_default_1percent_cen[i_cen][i_pt]->GetEntries()==0) cout<<"[WARN] empty: "<<h_mass_default_1percent_cen[i_cen][i_pt]->GetName()<<endl;
        h_mass_default_1percent_cen[i_cen][i_pt]->Write();
      }
      for (int iq=0; iq<N_QBINS; ++iq) {
        if (h_mass_1pecent_cen_q2[i_cen][iq][i_pt]) {
          if (VERBOSE_SAVE) cout<<"WRITE: "<<h_mass_1pecent_cen_q2[i_cen][iq][i_pt]->GetName()<<" entries="<<h_mass_1pecent_cen_q2[i_cen][iq][i_pt]->GetEntries()<<endl;
          if (h_mass_1pecent_cen_q2[i_cen][iq][i_pt]->GetEntries()==0) cout<<"[WARN] empty: "<<h_mass_1pecent_cen_q2[i_cen][iq][i_pt]->GetName()<<endl;
          h_mass_1pecent_cen_q2[i_cen][iq][i_pt]->Write();
        }
        if (h_mass_1pecent_cen_q3[i_cen][iq][i_pt]) {
          if (VERBOSE_SAVE) cout<<"WRITE: "<<h_mass_1pecent_cen_q3[i_cen][iq][i_pt]->GetName()<<" entries="<<h_mass_1pecent_cen_q3[i_cen][iq][i_pt]->GetEntries()<<endl;
          if (h_mass_1pecent_cen_q3[i_cen][iq][i_pt]->GetEntries()==0) cout<<"[WARN] empty: "<<h_mass_1pecent_cen_q3[i_cen][iq][i_pt]->GetName()<<endl;
          h_mass_1pecent_cen_q3[i_cen][iq][i_pt]->Write();
        }
        for (int iv=0; iv<N_VBINS; ++iv) {
          if (h_mass_v2_fit_1pecent_cen[i_cen][iq][i_pt][iv]) {
            if (VERBOSE_SAVE) cout<<"WRITE: "<<h_mass_v2_fit_1pecent_cen[i_cen][iq][i_pt][iv]->GetName()<<" entries="<<h_mass_v2_fit_1pecent_cen[i_cen][iq][i_pt][iv]->GetEntries()<<endl;
            if (h_mass_v2_fit_1pecent_cen[i_cen][iq][i_pt][iv]->GetEntries()==0) cout<<"[WARN] empty: "<<h_mass_v2_fit_1pecent_cen[i_cen][iq][i_pt][iv]->GetName()<<endl;
            h_mass_v2_fit_1pecent_cen[i_cen][iq][i_pt][iv]->Write();
          }
          if (h_mass_v3_fit_1pecent_cen[i_cen][iq][i_pt][iv]) {
            if (VERBOSE_SAVE) cout<<"WRITE: "<<h_mass_v3_fit_1pecent_cen[i_cen][iq][i_pt][iv]->GetName()<<" entries="<<h_mass_v3_fit_1pecent_cen[i_cen][iq][i_pt][iv]->GetEntries()<<endl;
            if (h_mass_v3_fit_1pecent_cen[i_cen][iq][i_pt][iv]->GetEntries()==0) cout<<"[WARN] empty: "<<h_mass_v3_fit_1pecent_cen[i_cen][iq][i_pt][iv]->GetName()<<endl;
            h_mass_v3_fit_1pecent_cen[i_cen][iq][i_pt][iv]->Write();
          }
        }
      }
      outf->cd();
    }
  }
  outf->Write();
  outf->Close();

  
  cout<<"Saved  "<<outname<<endl;
  return 0;
}

int main(int argc, char *argv[])
{
  if (argc==1) {
    save_mass_distributions();
  } else if (argc==2) {
    int idx = atoi(argv[1]);
    save_mass_distributions(idx);
  } else if (argc==3) {
    int cen_idx = atoi(argv[1]);
    int pt_idx  = atoi(argv[2]);
    save_mass_distributions(cen_idx, pt_idx);
  } else {
    std::cout << "Usage: ./save_mass_distributions [target1pct_index] [target_pt_index]" << std::endl;
    return 1;
  }
  return 0;
}

