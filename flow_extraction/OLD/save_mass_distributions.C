// Macro to build and save mass distribution histograms
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

  // helper: formatted edge string (fixed precision)
  auto fmt_edge = [](double x)->std::string {
    char buf[32];
    snprintf(buf, sizeof(buf), "%.3f", x);
    return std::string(buf);
  };
  const bool debug = false; // set true for verbose checks

  // determine which 10% centrality group we need (if any)
  int needed_cen_group = -1;
  if (target1pct >= 0) needed_cen_group = get_cent_group(target1pct);

  TString cuts;

  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    if (needed_cen_group >= 0 && i_cen != needed_cen_group) continue; // skip other 10% groups
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
      if (target_pt >= 0 && i_pt != target_pt) continue; // skip other pT bins

      cuts = "dca<0.0085 && mass<2.0 && mass>1.74 && centrality>="+to_string(cen_edges[i_cen])+" && centrality<"+to_string(cen_edges[i_cen+1])+" && pT>="+to_string(pt_edges[i_pt])+" && pT<"+to_string(pt_edges[i_pt+1]);

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
        // use ROOT's GetQuantiles to compute the N_EDGES quantile boundaries
        std::vector<Double_t> probs(N_EDGES), qout(N_EDGES);
        for (int ie=1; ie<=N_EDGES; ++ie) probs[ie-1] = ie / (Double_t)N_EDGES; // 0.1,0.2,...,1.0
        h_v2_bin->GetQuantiles(N_EDGES, qout.data(), probs.data());
        for (int ie=1; ie<=N_EDGES; ++ie) {
          Double_t edge = qout[ie-1];
          vnbinning_v2[i_cen][i_pt][ie+7] = edge;
          vnbinning_v2[i_cen][i_pt][27-ie] = -edge;
          if (debug) cout << "v2 quantile " << ie << " (cen " << i_cen << " pt " << i_pt << "): +=" << edge << "  -=" << -edge << endl;
        }
      } else {
        if (debug) cout << "No entries for v2 quantiles, using default symmetric vnbinning_v2 for cen " << i_cen << " pt " << i_pt << endl;
      }

      if (n_entries_v3 > 0) {
        std::vector<Double_t> probs3(N_EDGES), qout3(N_EDGES);
        for (int ie=1; ie<=N_EDGES; ++ie) probs3[ie-1] = ie / (Double_t)N_EDGES;
        h_v3_bin->GetQuantiles(N_EDGES, qout3.data(), probs3.data());
        for (int ie=1; ie<=N_EDGES; ++ie) {
          Double_t edge3 = qout3[ie-1];
          vnbinning_v3[i_cen][i_pt][ie+7] = edge3;
          vnbinning_v3[i_cen][i_pt][27-ie] = -edge3;
          if (debug) cout << "v3 quantile " << ie << " (cen " << i_cen << " pt " << i_pt << "): +=" << edge3 << "  -=" << -edge3 << endl;
        }
      } else {
        if (debug) cout << "No entries for v3 quantiles, using default symmetric v3 binning for cen " << i_cen << " pt " << i_pt << endl;
      }

      // print full edge arrays for quick check
      if (debug) {
        cout << "vnbinning_v2[" << i_cen << "][" << i_pt << "] =";
        for (int ib=0; ib<=N_VBINS; ++ib) cout << " " << vnbinning_v2[i_cen][i_pt][ib];
        cout << endl;
        cout << "vnbinning_v3[" << i_cen << "][" << i_pt << "] =";
        for (int ib=0; ib<=N_VBINS; ++ib) cout << " " << vnbinning_v3[i_cen][i_pt][ib];
        cout << endl;
      }
      
    }// --PTBIN--
  }// -- CENTBIN --

  // histograms declarations (use static to avoid large stack usage)
  static TH1D *h_mass_default[N_CENTBINS][N_PTBINS];
  static TH1D *h_mass_default_1percent_cen[N_CENTBINS_1][N_PTBINS];
  static TH1D *h_mass[N_CENTBINS][N_QBINS][N_PTBINS];
  static TH1D *h_mass_1pecent_cen[N_CENTBINS_1][N_QBINS][N_PTBINS];
  static TH1D *h_mass_v2_fit[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS];
  static TH1D *h_mass_v2_fit_1pecent_cen[N_CENTBINS_1][N_QBINS][N_PTBINS][N_VBINS];
  static TH1D *h_mass_v3_fit[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS];
  static TH1D *h_mass_v3_fit_1pecent_cen[N_CENTBINS_1][N_QBINS][N_PTBINS][N_VBINS];
  static TH1D *h_v2[N_CENTBINS][N_QBINS][N_PTBINS];
  static TH1D *h_v3[N_CENTBINS][N_QBINS][N_PTBINS];
  // explicit init (defensive, static are zeroed but explicit makes intent clear)
  for (int a=0;a<N_CENTBINS;++a) for(int b=0;b<N_QBINS;++b) for(int c=0;c<N_PTBINS;++c)
    for(int d=0; d<N_VBINS; ++d) { h_mass_v2_fit[a][b][c][d]=nullptr; h_mass_v3_fit[a][b][c][d]=nullptr; }
  for (int a=0;a<N_CENTBINS;++a) for(int b=0;b<N_PTBINS;++b) { h_mass_default[a][b]=nullptr; h_v2[a][0][b]=nullptr; h_v3[a][0][b]=nullptr; }

  TString h_mass_name, h_v2_name, h_mass_v2_name, h_name_0, h_name_1;

  // create 1% histograms (loop order: i_cen, i_q2, i_pt, i_v2)
  for (int i_cen=0; i_cen<N_CENTBINS_1; ++i_cen) {
    if (target1pct >= 0 && i_cen != target1pct) continue;
    int cen_group_for_binning = get_cent_group(i_cen);
    if (cen_group_for_binning < 0) continue;
    for (int i_q2=0; i_q2<N_QBINS; ++i_q2) {
      for (int i_pt=0; i_pt<N_PTBINS; ++i_pt) {
        if (target_pt >= 0 && i_pt != target_pt) continue;
        h_name_0 = "h_mass_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_q2bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt];
        h_mass_1pecent_cen[i_cen][i_q2][i_pt] = new TH1D(h_name_0,h_name_0,N_MASSBINS,fit_range_low,fit_range_high);
        h_mass_1pecent_cen[i_cen][i_q2][i_pt]->SetDirectory(0);
        h_name_0 = "h_mass_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_pt_"+pt_name[i_pt]+"_def";
        if (!h_mass_default_1percent_cen[i_cen][i_pt]) { // create once per i_cen,i_pt
          h_mass_default_1percent_cen[i_cen][i_pt] = new TH1D(h_name_0,h_name_0,N_MASSBINS,fit_range_low,fit_range_high);
          h_mass_default_1percent_cen[i_cen][i_pt]->SetDirectory(0);
        }
        for (int i_v2=0; i_v2<N_VBINS; ++i_v2) {
          Double_t v2_e1 = vnbinning_v2[cen_group_for_binning][i_pt][i_v2];
          Double_t v2_e2 = vnbinning_v2[cen_group_for_binning][i_pt][i_v2+1];
          if (std::isnan(v2_e1) || std::isnan(v2_e2) || fabs(v2_e1 - v2_e2) < 1e-9) continue;
          h_name_1 = "hist_mass_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_q2bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt]+"_in_v2bin_"+fmt_edge(v2_e1)+"_"+fmt_edge(v2_e2);
          h_mass_v2_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2] = new TH1D(h_name_1,h_name_1,N_MASSBINS,fit_range_low,fit_range_high);
          h_mass_v2_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2]->SetDirectory(0);
          Double_t v3_e1 = vnbinning_v3[cen_group_for_binning][i_pt][i_v2];
          Double_t v3_e2 = vnbinning_v3[cen_group_for_binning][i_pt][i_v2+1];
          if (std::isnan(v3_e1) || std::isnan(v3_e2) || fabs(v3_e1 - v3_e2) < 1e-9) continue;
          h_name_1 = "hist_mass_v3_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_q2bin_"+to_string(i_q2)+"_pt_"+pt_name[i_pt]+"_in_v3bin_"+fmt_edge(v3_e1)+"_"+fmt_edge(v3_e2);
          h_mass_v3_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2] = new TH1D(h_name_1,h_name_1,N_MASSBINS,fit_range_low,fit_range_high);
          h_mass_v3_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2]->SetDirectory(0);
        }
      }
    }
  }

  // default (10% centrality) histograms: mass distributions without v2,v3 selection
  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
      h_name_0 = "h_mass_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_pt_"+pt_name[i_pt]+"_def";
      h_mass_default[i_cen][i_pt] = new TH1D(h_name_0,h_name_0,N_MASSBINS,fit_range_low,fit_range_high);
      h_mass_default[i_cen][i_pt]->SetDirectory(0);
    }
  }

  // fill default mass histograms (no selection) for all events
  {
    h_mass_name = "h_mass_cen_all_pt_all_def";
    auto h_mass_all_def = new TH1D(h_mass_name,h_mass_name,N_MASSBINS,fit_range_low,fit_range_high);
    h_mass_all_def->SetDirectory(0);
    input_nt->Project(h_mass_all_def->GetName(),"mass","");
  }

  // per-centrality mass histograms (no v2,v3 selection)
  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    int cen_min = cen_edges[i_cen];
    int cen_max = cen_edges[i_cen+1];
    TString cut_cen = "centrality>="+to_string(cen_min)+" && centrality<"+to_string(cen_max);
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
      h_mass_name = "h_mass_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_pt_"+pt_name[i_pt];
      auto h_mass_cen_pt = new TH1D(h_mass_name,h_mass_name,N_MASSBINS,fit_range_low,fit_range_high);
      h_mass_cen_pt->SetDirectory(0);
      input_nt->Project(h_mass_cen_pt->GetName(),"mass",cut_cen+" && pT>="+to_string(pt_edges[i_pt])+" && pT<"+to_string(pt_edges[i_pt+1]));
    }
  }

  // per-centrality, per-pt mass histograms (no v2,v3 selection)
  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    int cen_min = cen_edges[i_cen];
    int cen_max = cen_edges[i_cen+1];
    TString cut_cen = "centrality>="+to_string(cen_min)+" && centrality<"+to_string(cen_max);
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
      h_mass_name = "h_mass_cen_"+to_string(i_cen)+"_"+to_string(i_cen+1)+"_pt_"+pt_name[i_pt];
      auto h_mass_cen_pt = new TH1D(h_mass_name,h_mass_name,N_MASSBINS,fit_range_low,fit_range_high);
      h_mass_cen_pt->SetDirectory(0);
      input_nt->Project(h_mass_cen_pt->GetName(),"mass",cut_cen+" && pT>="+to_string(pt_edges[i_pt])+" && pT<"+to_string(pt_edges[i_pt+1]));
    }
  }

  // per-centrality, per-pt, per-v2 mass histograms (no v3 selection)
  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    int cen_min = cen_edges[i_cen];
    int cen_max = cen_edges[i_cen+1];
    TString cut_cen = "centrality>="+to_string(cen_min)+" && centrality<"+to_string(cen_max);
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
      for (int i_v2=0; i_v2<N_VBINS; i_v2++) {
        Double_t e1 = vnbinning_v2[i_cen][i_pt][i_v2];
        Double_t e2 = vnbinning_v2[i_cen][i_pt][i_v2+1];
        if (fabs(e1 - e2) < 1e-9) continue;
        h_mass_name = "hist_mass_cen_"+cen_name[i_cen]+"_q2bin_-1_pt_"+pt_name[i_pt]+"_in_v2bin_"+fmt_edge(e1)+"_"+fmt_edge(e2);
        h_mass_v2_fit[i_cen][i_q2][i_pt][i_v2]=new TH1D(h_mass_v2_name, h_mass_v2_name, N_MASSBINS, fit_range_low, fit_range_high);
        h_mass_v2_fit[i_cen][i_q2][i_pt][i_v2]->SetDirectory(0);

        Double_t f1 = vnbinning_v3[i_cen][i_pt][i_v2];
        Double_t f2 = vnbinning_v3[i_cen][i_pt][i_v2+1];
        if (fabs(f1 - f2) < 1e-9) continue;
        TString h_mass_v3_name = "hist_mass_v3_cen_"+cen_name[i_cen]+"_q2bin_-1_pt_"+pt_name[i_pt]+"_in_v3bin_"+fmt_edge(f1)+"_"+fmt_edge(f2);
        h_mass_v3_fit[i_cen][i_q2][i_pt][i_v2] = new TH1D(h_mass_v3_name, h_mass_v3_name, N_MASSBINS, fit_range_low, fit_range_high);
        h_mass_v3_fit[i_cen][i_q2][i_pt][i_v2]->SetDirectory(0);
      }// -- VBIN --
    }// -- PTBIN --
  }// -- CENTBIN --

  // per-centrality, per-pt, per-v3 mass histograms (no v2 selection)
  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    int cen_min = cen_edges[i_cen];
    int cen_max = cen_edges[i_cen+1];
    TString cut_cen = "centrality>="+to_string(cen_min)+" && centrality<"+to_string(cen_max);
    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
      for (int i_v3=0; i_v3<N_VBINS; i_v3++) {
        Double_t e1 = vnbinning_v3[i_cen][i_pt][i_v3];
        Double_t e2 = vnbinning_v3[i_cen][i_pt][i_v3+1];
        if (fabs(e1 - e2) < 1e-9) continue;
        h_mass_name = "hist_mass_cen_"+cen_name[i_cen]+"_q2bin_-1_pt_"+pt_name[i_pt]+"_in_v3bin_"+fmt_edge(e1)+"_"+fmt_edge(e2);
        h_mass_v3_fit[i_cen][i_q2][i_pt][i_v3] = new TH1D(h_mass_v3_name, h_mass_v3_name, N_MASSBINS, fit_range_low, fit_range_high);
        h_mass_v3_fit[i_cen][i_q2][i_pt][i_v3]->SetDirectory(0);
      }// -- VBIN --
    }// -- PTBIN --
  }// -- CENTBIN --

  // save histograms: ensure order i_cen, i_q2, i_pt, i_v2 / i_v3 when writing per-vn
  {
    TFile *outf = new TFile("mass_distributions.root","recreate");
    outf->cd();
    for (int i_cen=0; i_cen<N_CENTBINS; ++i_cen) {
      for (int i_q2=0; i_q2<N_QBINS; ++i_q2) {
        for (int i_pt=0; i_pt<N_PTBINS; ++i_pt) {
          if (h_mass[i_cen][i_q2][i_pt]) h_mass[i_cen][i_q2][i_pt]->Write();
          if (h_v2[i_cen][i_q2][i_pt]) h_v2[i_cen][i_q2][i_pt]->Write();
          if (h_v3[i_cen][i_q2][i_pt]) h_v3[i_cen][i_q2][i_pt]->Write();
          for (int i_v2=0; i_v2<N_VBINS; ++i_v2) {
            if (h_mass_v2_fit[i_cen][i_q2][i_pt][i_v2]) h_mass_v2_fit[i_cen][i_q2][i_pt][i_v2]->Write();
            if (h_mass_v3_fit[i_cen][i_q2][i_pt][i_v2]) h_mass_v3_fit[i_cen][i_q2][i_pt][i_v2]->Write();
          }
        }
      }
      // write defaults per i_cen (optional)
      for (int i_pt=0; i_pt<N_PTBINS; ++i_pt) if (h_mass_default[i_cen][i_pt]) h_mass_default[i_cen][i_pt]->Write();
    }
    // write 1% histograms (same order)
    for (int i_cen=0; i_cen<N_CENTBINS_1; ++i_cen) {
      for (int i_q2=0; i_q2<N_QBINS; ++i_q2) {
        for (int i_pt=0; i_pt<N_PTBINS; ++i_pt) {
          if (h_mass_1pecent_cen[i_cen][i_q2][i_pt]) h_mass_1pecent_cen[i_cen][i_q2][i_pt]->Write();
          if (h_mass_default_1percent_cen[i_cen][i_pt]) h_mass_default_1percent_cen[i_cen][i_pt]->Write();
          for (int i_v2=0; i_v2<N_VBINS; ++i_v2) {
            if (h_mass_v2_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2]) h_mass_v2_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2]->Write();
            if (h_mass_v3_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2]) h_mass_v3_fit_1pecent_cen[i_cen][i_q2][i_pt][i_v2]->Write();
          }
        }
      }
    }
    outf->Write();
    outf->Close();
   }
 
   return 0;
 }
