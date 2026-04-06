#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
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
#include <TParameter.h>

#include "Analysis_bin.h"

using namespace std;

int fit_mass_and_flow(int target_cent = -1)
{

  // --- Open ntuple data file (replaces pre-saved mass distributions) ---
  const char *inputNtuple = "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/SP_MB11to21_Feb20_v4/ROOT/Flow_combined.root";
  TFile *nt_Data_file = TFile::Open(inputNtuple);
  if (!nt_Data_file || nt_Data_file->IsZombie())
  {
    std::cerr << "Cannot open ntuple Data file." << std::endl;
    return 1;
  }
  TNtuple *input_nt = (TNtuple *)nt_Data_file->Get("nt");
  if (!input_nt)
  {
    std::cerr << "Cannot find ntuple 'nt' in data file." << std::endl;
    return 1;
  }
  const char *inputFile_MC = "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/MC_template_MC2023_Feb16_v2/ROOT/D0_MCtemplate_out_combined.root";
  TFile *inf_MC = TFile::Open(inputFile_MC);
  if (!inf_MC || inf_MC->IsZombie())
  {
    std::cerr << "Cannot open MC file." << std::endl;
    return 1;
  }

  // --- Get the v2 and v3 SP bining table from .h file ----
  static Double_t vnbinning_v2[N_CENTBINS][N_PTBINS][N_VBINS_V2 + 1];
  static Double_t vnbinning_v3[N_CENTBINS][N_PTBINS][N_VBINS_V3 + 1];
  for (int ic = 0; ic < N_CENTBINS; ++ic)
  {
    for (int ip = 0; ip < N_PTBINS; ++ip)
    {
      for (int iv = 0; iv <= N_VBINS_V2; ++iv)
        vnbinning_v2[ic][ip][iv] = vnb_v2_table[ic][ip][iv];
      for (int iv = 0; iv <= N_VBINS_V3; ++iv)
        vnbinning_v3[ic][ip][iv] = vnb_v3_table[ic][ip][iv];
    }
  }

  TH1::AddDirectory(kFALSE);
  gROOT->cd();

  // index of the requested 1%-bin (0-based) OR coarse centrality group (0..N_CENTBINS-1)
  int target_cent_idx = -1;
  int target_cen_group = -1;
  auto map_1pct_index_to_group = [](int idx) -> int
  {
    if (idx < 0)
      return -1;
    if (idx <= 9)
      return 0; // 0..9   -> 0-10%
    if (idx <= 19)
      return 1; // 10..19 -> 10-20%
    if (idx <= 29)
      return 2; // 20..29 -> 20-30%
    if (idx <= 39)
      return 3; // 30..40%
    if (idx <= 49)
      return 4; // 40..50%
    if (idx <= 79)
      return 5; // 50..80%
    return -1;
  };
  // Treat explicit coarse indices 0..N_CENTBINS-1 first; otherwise interpret 1..N_CENTBINS_1 as 1%-bin
  if (target_cent >= 0 && target_cent < N_CENTBINS)
  {
    // user supplied a coarse centrality index (0..N_CENTBINS-1)
    target_cen_group = target_cent;
    target_cent_idx = -1;
  }
  else if (target_cent >= 1 && target_cent <= N_CENTBINS_1)
  {
    // user supplied a 1%-bin (1..90) -> convert to 0-based and map to coarse group
    target_cent_idx = target_cent - 1;
    target_cen_group = map_1pct_index_to_group(target_cent_idx);
  }
  else
  {
    target_cent_idx = -1;
    target_cen_group = -1;
  }

  // ---- Coarse-bin (N_CENTBINS) histograms — filled by merging 1%-bin histograms ----
  TH1D *h_mass_default[N_CENTBINS][N_PTBINS] = {};
  TH1D *h_mass_q2[N_CENTBINS][N_QBINS][N_PTBINS] = {}; // q2-selected, used for v2 shape fits
  TH1D *h_mass_q3[N_CENTBINS][N_QBINS][N_PTBINS] = {}; // q3-selected, used for v3 shape fits
  TH1D *h_v2dist_grouped[N_CENTBINS][N_QBINS][N_PTBINS] = {};
  TH1D *h_v3dist_grouped[N_CENTBINS][N_QBINS][N_PTBINS] = {};
  TH1D *h_mass_v2_fit[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS_V2] = {};
  TH1D *h_mass_v3_fit[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS_V3] = {};

  // Create coarse-bin histograms
  for (int ic = 0; ic < N_CENTBINS; ++ic)
  {
    for (int ip = 0; ip < N_PTBINS; ++ip)
    {
      TString hdef = Form("h_mass_%s_%s_def", cen_name[ic], pt_name[ip]);
      h_mass_default[ic][ip] = new TH1D(hdef, hdef, N_MASSBINS, hist_range_low, hist_range_high);
      h_mass_default[ic][ip]->Sumw2();
      h_mass_default[ic][ip]->SetDirectory(nullptr);
    }
    for (int iq = 0; iq < N_QBINS; ++iq)
    {
      for (int ip = 0; ip < N_PTBINS; ++ip)
      {
        TString hm_q2 = Form("h_mass_%s_q2bin%d_%s", cen_name[ic], iq, pt_name[ip]);
        h_mass_q2[ic][iq][ip] = new TH1D(hm_q2, hm_q2, N_MASSBINS, hist_range_low, hist_range_high);
        h_mass_q2[ic][iq][ip]->Sumw2();
        h_mass_q2[ic][iq][ip]->SetDirectory(nullptr);

        TString hm_q3 = Form("h_mass_%s_q3bin%d_%s", cen_name[ic], iq, pt_name[ip]);
        h_mass_q3[ic][iq][ip] = new TH1D(hm_q3, hm_q3, N_MASSBINS, hist_range_low, hist_range_high);
        h_mass_q3[ic][iq][ip]->Sumw2();
        h_mass_q3[ic][iq][ip]->SetDirectory(nullptr);

        TString hv2d = Form("h_v2dist_%s_q2bin%d_%s", cen_name[ic], iq, pt_name[ip]);
        h_v2dist_grouped[ic][iq][ip] = new TH1D(hv2d, hv2d, 250, -25.0, 25.0); // This is just to see distribution, so const bin width is not a problem!
        h_v2dist_grouped[ic][iq][ip]->Sumw2();
        h_v2dist_grouped[ic][iq][ip]->SetDirectory(nullptr);

        TString hv3d = Form("h_v3dist_%s_q3bin%d_%s", cen_name[ic], iq, pt_name[ip]);
        h_v3dist_grouped[ic][iq][ip] = new TH1D(hv3d, hv3d, 500, -50.0, 50.0); // Same here!
        h_v3dist_grouped[ic][iq][ip]->Sumw2();
        h_v3dist_grouped[ic][iq][ip]->SetDirectory(nullptr);

        for (int iv = 0; iv < N_VBINS_V2; ++iv)
        {
          TString hv2 = Form("hist_mass_v2_%s_q2bin%d_pt_%s_in_v2bin_idx_%d", cen_name[ic], iq, pt_name[ip], iv);
          h_mass_v2_fit[ic][iq][ip][iv] = new TH1D(hv2, hv2, N_MASSBINS, hist_range_low, hist_range_high);
          h_mass_v2_fit[ic][iq][ip][iv]->Sumw2();
          h_mass_v2_fit[ic][iq][ip][iv]->SetDirectory(nullptr);
        }
        for (int iv = 0; iv < N_VBINS_V3; ++iv)
        {
          TString hv3 = Form("hist_mass_v3_%s_q3bin%d_pt_%s_in_v3bin_idx_%d", cen_name[ic], iq, pt_name[ip], iv);
          h_mass_v3_fit[ic][iq][ip][iv] = new TH1D(hv3, hv3, N_MASSBINS, hist_range_low, hist_range_high);
          h_mass_v3_fit[ic][iq][ip][iv]->Sumw2();
          h_mass_v3_fit[ic][iq][ip][iv]->SetDirectory(nullptr);
        }
      }
    }
  }

  // ---- 1%-bin (N_CENTBINS_1) histograms — filled directly from ntuple ----
  // Use static to avoid stack overflow for large pointer arrays
  static TH1D *h_mass_default_1pct[N_CENTBINS_1][N_PTBINS];
  static TH1D *h_mass_1pct_q2[N_CENTBINS_1][N_QBINS][N_PTBINS];
  static TH1D *h_mass_1pct_q3[N_CENTBINS_1][N_QBINS][N_PTBINS];
  static TH1D *h_v2dist_1pct[N_CENTBINS_1][N_QBINS][N_PTBINS];
  static TH1D *h_v3dist_1pct[N_CENTBINS_1][N_QBINS][N_PTBINS];
  static TH1D *h_mass_v2_fit_1pct[N_CENTBINS_1][N_QBINS][N_PTBINS][N_VBINS_V2];
  static TH1D *h_mass_v3_fit_1pct[N_CENTBINS_1][N_QBINS][N_PTBINS][N_VBINS_V3];
  // FIX: must zero-initialize so unused 1%-bins are guaranteed nullptr
  memset(h_mass_default_1pct, 0, sizeof(h_mass_default_1pct));
  memset(h_mass_1pct_q2, 0, sizeof(h_mass_1pct_q2));
  memset(h_mass_1pct_q3, 0, sizeof(h_mass_1pct_q3));
  memset(h_v2dist_1pct, 0, sizeof(h_v2dist_1pct));
  memset(h_v3dist_1pct, 0, sizeof(h_v3dist_1pct));
  memset(h_mass_v2_fit_1pct, 0, sizeof(h_mass_v2_fit_1pct));
  memset(h_mass_v3_fit_1pct, 0, sizeof(h_mass_v3_fit_1pct));

  // Create 1%-bin histograms (only for the centrality range requested)
  for (int i_cen1 = 0; i_cen1 < N_CENTBINS_1; ++i_cen1)
  {
    int cen_group = map_1pct_index_to_group(i_cen1);
    if (cen_group < 0)
      continue;
    if (target_cen_group >= 0 && cen_group != target_cen_group)
      continue;
    if (target_cent_idx >= 0 && i_cen1 != target_cent_idx)
      continue;

    for (int ip = 0; ip < N_PTBINS; ++ip)
    {
      TString hdef = Form("h_mass_cen%d_%s_def", i_cen1, pt_name[ip]);
      h_mass_default_1pct[i_cen1][ip] = new TH1D(hdef, hdef, N_MASSBINS, hist_range_low, hist_range_high);
      h_mass_default_1pct[i_cen1][ip]->Sumw2();
      h_mass_default_1pct[i_cen1][ip]->SetDirectory(nullptr);
    }
    for (int iq = 0; iq < N_QBINS; ++iq)
    {
      for (int ip = 0; ip < N_PTBINS; ++ip)
      {
        TString hq2 = Form("h_mass_cen%d_q2bin%d_%s", i_cen1, iq, pt_name[ip]);
        h_mass_1pct_q2[i_cen1][iq][ip] = new TH1D(hq2, hq2, N_MASSBINS, hist_range_low, hist_range_high);
        h_mass_1pct_q2[i_cen1][iq][ip]->Sumw2();
        h_mass_1pct_q2[i_cen1][iq][ip]->SetDirectory(nullptr);

        TString hq3 = Form("h_mass_cen%d_q3bin%d_%s", i_cen1, iq, pt_name[ip]);
        h_mass_1pct_q3[i_cen1][iq][ip] = new TH1D(hq3, hq3, N_MASSBINS, hist_range_low, hist_range_high);
        h_mass_1pct_q3[i_cen1][iq][ip]->Sumw2();
        h_mass_1pct_q3[i_cen1][iq][ip]->SetDirectory(nullptr);

        TString hv2d = Form("h_v2dist_cen%d_q2bin%d_%s", i_cen1, iq, pt_name[ip]);
        h_v2dist_1pct[i_cen1][iq][ip] = new TH1D(hv2d, hv2d, 250, -25.0, 25.0);
        h_v2dist_1pct[i_cen1][iq][ip]->Sumw2();
        h_v2dist_1pct[i_cen1][iq][ip]->SetDirectory(nullptr);

        TString hv3d = Form("h_v3dist_cen%d_q3bin%d_%s", i_cen1, iq, pt_name[ip]);
        h_v3dist_1pct[i_cen1][iq][ip] = new TH1D(hv3d, hv3d, 500, -50.0, 50.0);
        h_v3dist_1pct[i_cen1][iq][ip]->Sumw2();
        h_v3dist_1pct[i_cen1][iq][ip]->SetDirectory(nullptr);

        for (int iv = 0; iv < N_VBINS_V2; ++iv)
        {
          TString hv2 = Form("hist_mass_v2_cen%d_q2bin%d_%s_in_v2bin_idx_%d", i_cen1, iq, pt_name[ip], iv);
          h_mass_v2_fit_1pct[i_cen1][iq][ip][iv] = new TH1D(hv2, hv2, N_MASSBINS, hist_range_low, hist_range_high);
          h_mass_v2_fit_1pct[i_cen1][iq][ip][iv]->Sumw2();
          h_mass_v2_fit_1pct[i_cen1][iq][ip][iv]->SetDirectory(nullptr);
        }
        for (int iv = 0; iv < N_VBINS_V3; ++iv)
        {
          TString hv3 = Form("hist_mass_v3_cen%d_q3bin%d_%s_in_v3bin_idx_%d", i_cen1, iq, pt_name[ip], iv);
          h_mass_v3_fit_1pct[i_cen1][iq][ip][iv] = new TH1D(hv3, hv3, N_MASSBINS, hist_range_low, hist_range_high);
          h_mass_v3_fit_1pct[i_cen1][iq][ip][iv]->Sumw2();
          h_mass_v3_fit_1pct[i_cen1][iq][ip][iv]->SetDirectory(nullptr);
        }
      }
    }
  }

  // --- Step 1: Fill 1%-bin histograms from ntuple (exactly like save_mass_distributions.C) ---
  {
    Float_t cen_val, q2_val, q3_val, pT_val, mass_val, v2_val, v3_val, dca_val, y_val;
    input_nt->SetBranchAddress("cent", &cen_val); // cent stored as 0-based 0..79
    input_nt->SetBranchAddress("pT", &pT_val);
    input_nt->SetBranchAddress("mass", &mass_val);
    input_nt->SetBranchAddress("dca", &dca_val);
    input_nt->SetBranchAddress("y", &y_val);
    input_nt->SetBranchAddress("q2_hf_total", &q2_val);
    input_nt->SetBranchAddress("q3_hf_total", &q3_val);
    input_nt->SetBranchAddress("v2", &v2_val);
    input_nt->SetBranchAddress("v3", &v3_val);

    Long64_t N_ENTRIES = input_nt->GetEntriesFast();
    for (Long64_t i_entry = 0; i_entry < N_ENTRIES; i_entry++)
    {
      input_nt->GetEntry(i_entry);
      if (i_entry % 1000000 == 0)
        std::cout << i_entry << " / " << N_ENTRIES
                  << "  " << 100 * i_entry / N_ENTRIES << "%" << std::endl;

      if (dca_val >= 0.0085)
        continue;
      if (TMath::Abs(y_val) >= 1.0)
        continue;

      int i_cen1 = (int)cen_val; // 0-based 1%-bin
      if (i_cen1 < 0 || i_cen1 >= N_CENTBINS_1)
        continue;
      int cen_group = map_1pct_index_to_group(i_cen1);
      if (cen_group < 0 || cen_group >= N_CENTBINS)
        continue;
      // Apply centrality filter
      if (target_cen_group >= 0 && cen_group != target_cen_group)
        continue;
      if (target_cent_idx >= 0 && i_cen1 != target_cent_idx)
        continue;

      for (int i_pt = 0; i_pt < N_PTBINS; i_pt++)
      {
        if (pT_val < pt_edges[i_pt] || pT_val >= pt_edges[i_pt + 1])
          continue;

        if (h_mass_default_1pct[i_cen1][i_pt])
          h_mass_default_1pct[i_cen1][i_pt]->Fill(mass_val);

        for (int iq = 0; iq < N_QBINS; iq++)
        {
          // --- q2-binned: mass (q2), v2 distribution, v2-slice mass ---
          if (q2_val >= q2_cuts[i_cen1][iq] && q2_val < q2_cuts[i_cen1][iq + 1])
          {
            if (h_mass_1pct_q2[i_cen1][iq][i_pt])
              h_mass_1pct_q2[i_cen1][iq][i_pt]->Fill(mass_val);
            if (h_v2dist_1pct[i_cen1][iq][i_pt])
              h_v2dist_1pct[i_cen1][iq][i_pt]->Fill(v2_val);
            for (int i_v2 = 0; i_v2 < N_VBINS_V2; i_v2++)
            {
              if (v2_val >= vnbinning_v2[cen_group][i_pt][i_v2] && v2_val < vnbinning_v2[cen_group][i_pt][i_v2 + 1])
                if (h_mass_v2_fit_1pct[i_cen1][iq][i_pt][i_v2])
                  h_mass_v2_fit_1pct[i_cen1][iq][i_pt][i_v2]->Fill(mass_val);
            } // -- v2 slice --
          } // -- q2 cut --

          // --- q3-binned: mass (q3), v3 distribution, v3-slice mass ---
          if (q3_val >= q3_cuts[i_cen1][iq] && q3_val < q3_cuts[i_cen1][iq + 1])
          {
            if (h_mass_1pct_q3[i_cen1][iq][i_pt])
              h_mass_1pct_q3[i_cen1][iq][i_pt]->Fill(mass_val);
            if (h_v3dist_1pct[i_cen1][iq][i_pt])
              h_v3dist_1pct[i_cen1][iq][i_pt]->Fill(v3_val);
            for (int i_v3 = 0; i_v3 < N_VBINS_V3; i_v3++)
            {
              if (v3_val >= vnbinning_v3[cen_group][i_pt][i_v3] && v3_val < vnbinning_v3[cen_group][i_pt][i_v3 + 1])
                if (h_mass_v3_fit_1pct[i_cen1][iq][i_pt][i_v3])
                  h_mass_v3_fit_1pct[i_cen1][iq][i_pt][i_v3]->Fill(mass_val);
            } // -- v3 slice --
          } // -- q3 cut --
        } // ---- iq ----
      } // ---- i_pt ----
    } // ---- event loop ----
  } // end ntuple filling block

  // --- Step 2: Merge 1%-bin histograms into coarse-bin histograms ---
  for (int i_cen1 = 0; i_cen1 < N_CENTBINS_1; ++i_cen1)
  {
    int cen_group = map_1pct_index_to_group(i_cen1);
    if (cen_group < 0 || cen_group >= N_CENTBINS)
      continue;
    if (target_cen_group >= 0 && cen_group != target_cen_group)
      continue;
    if (target_cent_idx >= 0 && i_cen1 != target_cent_idx)
      continue;

    for (int ip = 0; ip < N_PTBINS; ++ip)
    {
      if (h_mass_default_1pct[i_cen1][ip] && h_mass_default[cen_group][ip])
        h_mass_default[cen_group][ip]->Add(h_mass_default_1pct[i_cen1][ip]);
    }
    for (int iq = 0; iq < N_QBINS; ++iq)
    {
      for (int ip = 0; ip < N_PTBINS; ++ip)
      {
        if (h_mass_1pct_q2[i_cen1][iq][ip] && h_mass_q2[cen_group][iq][ip])
          h_mass_q2[cen_group][iq][ip]->Add(h_mass_1pct_q2[i_cen1][iq][ip]);
        if (h_mass_1pct_q3[i_cen1][iq][ip] && h_mass_q3[cen_group][iq][ip])
          h_mass_q3[cen_group][iq][ip]->Add(h_mass_1pct_q3[i_cen1][iq][ip]);
        if (h_v2dist_1pct[i_cen1][iq][ip] && h_v2dist_grouped[cen_group][iq][ip])
          h_v2dist_grouped[cen_group][iq][ip]->Add(h_v2dist_1pct[i_cen1][iq][ip]);
        if (h_v3dist_1pct[i_cen1][iq][ip] && h_v3dist_grouped[cen_group][iq][ip])
          h_v3dist_grouped[cen_group][iq][ip]->Add(h_v3dist_1pct[i_cen1][iq][ip]);
        for (int iv = 0; iv < N_VBINS_V2; ++iv)
          if (h_mass_v2_fit_1pct[i_cen1][iq][ip][iv] && h_mass_v2_fit[cen_group][iq][ip][iv])
            h_mass_v2_fit[cen_group][iq][ip][iv]->Add(h_mass_v2_fit_1pct[i_cen1][iq][ip][iv]);
        for (int iv = 0; iv < N_VBINS_V3; ++iv)
          if (h_mass_v3_fit_1pct[i_cen1][iq][ip][iv] && h_mass_v3_fit[cen_group][iq][ip][iv])
            h_mass_v3_fit[cen_group][iq][ip][iv]->Add(h_mass_v3_fit_1pct[i_cen1][iq][ip][iv]);
      }
    }
  } // end merging 1%-bin -> coarse-bin

  // Prepare output structures
  Double_t yield_v2[N_VBINS_V2] = {}, yield_error_v2[N_VBINS_V2] = {}, v2_x[N_VBINS_V2] = {}, v2_x_err[N_VBINS_V2] = {}, mean_val_v2[N_PTBINS] = {}, mean_error_v2[N_PTBINS] = {};
  Double_t yield_v3[N_VBINS_V3] = {}, yield_error_v3[N_VBINS_V3] = {}, v3_x[N_VBINS_V3] = {}, v3_x_err[N_VBINS_V3] = {}, mean_val_v3[N_PTBINS] = {}, mean_error_v3[N_PTBINS] = {};
  Double_t pt_x[N_PTBINS] = {}, pt_x_error[N_PTBINS] = {};

  Double_t mean_v2_incl[N_PTBINS] = {}, mean_v2_incl_err[N_PTBINS] = {};
  Double_t mean_v3_incl[N_PTBINS] = {}, mean_v3_incl_err[N_PTBINS] = {};
  Double_t mean_v2_incl_simple[N_PTBINS] = {}, mean_v2_incl_err_simple[N_PTBINS] = {};
  Double_t mean_v3_incl_simple[N_PTBINS] = {}, mean_v3_incl_err_simple[N_PTBINS] = {};

  // Store v2_x and v3_x per (cen, pt) so summary loop uses correct x-values
  Double_t v2_x_store[N_CENTBINS][N_PTBINS][N_VBINS_V2] = {};
  Double_t v2_x_err_store[N_CENTBINS][N_PTBINS][N_VBINS_V2] = {};
  Double_t v3_x_store[N_CENTBINS][N_PTBINS][N_VBINS_V3] = {};
  Double_t v3_x_err_store[N_CENTBINS][N_PTBINS][N_VBINS_V3] = {};

  double chi2_ndf_for_q2_v2[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS_V2] = {};
  double chi2_ndf_for_q3_v3[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS_V3] = {};
  double sigma_v2[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS_V2] = {};
  double sigma_v3[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS_V3] = {};

  TGraphErrors *yield_v2_graph[N_CENTBINS][N_QBINS][N_PTBINS] = {};
  TGraphErrors *yield_v3_graph[N_CENTBINS][N_QBINS][N_PTBINS] = {};
  TGraph *sigma_v2_fit[N_CENTBINS][N_QBINS][N_PTBINS] = {};
  TGraph *sigma_v3_fit[N_CENTBINS][N_QBINS][N_PTBINS] = {};
  TGraph *chi2_ndf_for_q2_v2_graph[N_CENTBINS][N_QBINS][N_PTBINS] = {};
  TGraph *chi2_ndf_for_q3_v3_graph[N_CENTBINS][N_QBINS][N_PTBINS] = {};
  TGraphErrors *v2_cen[N_CENTBINS][N_QBINS] = {};
  TGraphErrors *v3_cen[N_CENTBINS][N_QBINS] = {};
  TGraphErrors *v2_cen_inclusive[N_CENTBINS] = {};
  TGraphErrors *v3_cen_inclusive[N_CENTBINS] = {};
  TGraphErrors *v2_cen_inclusive_simple[N_CENTBINS] = {};
  TGraphErrors *v3_cen_inclusive_simple[N_CENTBINS] = {};
  TH1D *h_v2_hist[N_CENTBINS][N_QBINS][N_PTBINS] = {};
  TH1D *h_v3_hist[N_CENTBINS][N_QBINS][N_PTBINS] = {};

  // create per-(cen,q,pt) histograms to hold yields per v-bin (prevent null deref)

  for (int ic = 0; ic < N_CENTBINS; ++ic)
  {
    for (int iq = 0; iq < N_QBINS; ++iq)
    {
      for (int ip = 0; ip < N_PTBINS; ++ip)
      {
        TString name_v2h = Form("h_v2_hist_%s_q2bin%d_%s", cen_name[ic], iq, pt_name[ip]);
        // h_v2_hist[ic][iq][ip] = new TH1D(name_v2h, name_v2h,N_VBINS,vnbinning_v2[ic][ip][0],vnbinning_v2[ic][ip][N_VBINS]);
        h_v2_hist[ic][iq][ip] = new TH1D(name_v2h, name_v2h, N_VBINS_V2, vnbinning_v2[ic][ip]);
        h_v2_hist[ic][iq][ip]->Sumw2();
        h_v2_hist[ic][iq][ip]->SetDirectory(nullptr);

        TString name_v3h = Form("h_v3_hist_%s_q3bin%d_%s", cen_name[ic], iq, pt_name[ip]);
        // h_v3_hist[ic][iq][ip] = new TH1D(name_v3h, name_v3h, N_VBINS,vnbinning_v3[ic][ip][0],vnbinning_v3[ic][ip][N_VBINS]);
        h_v3_hist[ic][iq][ip] = new TH1D(name_v3h, name_v3h, N_VBINS_V3, vnbinning_v3[ic][ip]);
        h_v3_hist[ic][iq][ip]->Sumw2();
        h_v3_hist[ic][iq][ip]->SetDirectory(nullptr);
      }
    }
  }

  TH1D *hist_mass = nullptr;

  TLatex *tex = new TLatex;
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  auto c1 = new TCanvas("c1", "c1", 600, 600);

  // Ensure the output directory for PDF mass plots exists
  const std::string plotDir = "prompt_mass_plot_withchi2_sigma";
  gSystem->mkdir(plotDir.c_str(), true);

  std::string outname = "Flow_output";
  if (target_cent >= 0)
    outname += std::string("_cen") + std::to_string(target_cent);
  // if (target_pt    >= 0) outname += std::string("_pt")  + std::to_string(target_pt);
  const char *env_job = std::getenv("SLURM_ARRAY_JOB_ID");
  if (!env_job)
    env_job = std::getenv("SLURM_JOB_ID");
  const char *env_task = std::getenv("SLURM_ARRAY_TASK_ID");
  if (env_job)
    outname += std::string("_job") + env_job;
  if (env_task)
    outname += std::string("_task") + env_task;
  outname += ".root";
  TFile *outf = new TFile(outname.c_str(), "recreate");
  // --- per-centrality sigma=-10 bad-bin log file (one per sbatch job/centrality) ---
  std::string logName = "sigma_bad_bins_";
  logName += (target_cen_group >= 0) ? cen_name[target_cen_group] : "allcent";
  logName += ".log";
  std::ofstream bad_bin_log(logName);
  // tag helper: is this v-bin in the central (mid) region?
  auto v2_is_mid = [](int iv) { return iv >= (N_SIDE_EDGES_V2-1) && iv < N_VBINS_V2-(N_SIDE_EDGES_V2-1); };
  auto v3_is_mid = [](int iv) { return iv >= (N_SIDE_EDGES_V3-1) && iv < N_VBINS_V3-(N_SIDE_EDGES_V3-1); };
  bad_bin_log << "# Bins where sigma_v2 or sigma_v3 was forced to -10 (problematic fits)\n";
  bad_bin_log << "# type | i_cen | cen_name | i_q | i_pt | pt_name | i_v | reason\n";

  TDirectory *dir_yield_vn = outf->mkdir("yields_vN");
  TDirectory *dir_hist_sp = outf->mkdir("histograms_SP");
  TDirectory *dir_hist_v2v3dist = outf->mkdir("histograms_v2v3dist"); // per-v TH1s
  TDirectory *dir_summary = outf->mkdir("summary_pt");                // v2/v3 vs pt graphs
  TDirectory *dir_diag_chi2_sigma = outf->mkdir("diagnostics_chi2_sigma");
  TDirectory *dir_mass_fitted = outf->mkdir("mass_fitted");

  outf->cd();

  for (int i_cen = 0; i_cen < N_CENTBINS; i_cen++)
  {
    if (target_cent >= 0 && i_cen != target_cen_group)
      continue;

    for (int i_pt = 0; i_pt < N_PTBINS; i_pt++)
    {

      TH1D *h_mc_match_signal;
      TH1D *h_mc_match_all;

      h_mc_match_signal = (TH1D *)inf_MC->Get(Form("hMass_Signal_%s", pt_name[i_pt]));
      h_mc_match_all = (TH1D *)inf_MC->Get(Form("hMass_Signal_Plus_Swap_%s", pt_name[i_pt]));

      if (h_mc_match_signal)
        h_mc_match_signal->SetTitle(Form("hMass_Signal_%s", pt_name[i_pt]));
      if (h_mc_match_all)
        h_mc_match_all->SetTitle(Form("hMass_Signal_Plus_Swap_%s", pt_name[i_pt]));

      if (i_pt == 10 || i_pt == 11)
      {
        h_mc_match_signal = (TH1D *)inf_MC->Get(Form("hMass_Signal_%s", pt_name[9]));
        h_mc_match_all = (TH1D *)inf_MC->Get(Form("hMass_Signal_Plus_Swap_%s", pt_name[9]));
        if (h_mc_match_signal)
          h_mc_match_signal->SetTitle(Form("hMass_Signal_%s", pt_name[9]));
        if (h_mc_match_all)
          h_mc_match_all->SetTitle(Form("hMass_Signal_Plus_Swap_%s", pt_name[9]));
      }

      TH1D *h_def = h_mass_default[i_cen][i_pt];

      /*TF1 *f_def = new TF1(Form("f_def_ptbin_%d", i_pt), "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])) ) + [8] + [9]*x + [10]*x*x + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12])) + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])) )", fit_range_low, fit_range_high);

      // Fix the smearing factor from inclusive q bins.
      if (h_def)
      {
        // f->ReleaseParameter(1);//mean
        f_def->FixParameter(1, 1.8648);
        f_def->FixParameter(6, 0);
        h_def->Fit(Form("f_def_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
        h_def->Fit(Form("f_def_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
        f_def->ReleaseParameter(1);
        f_def->ReleaseParameter(6);
        h_def->Fit(Form("f_def_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
        h_def->Fit(Form("f_def_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
        h_def->Fit(Form("f_def_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);

        double smearing_fixed = f_def->GetParameter(6);
        double mean_fixed = f_def->GetParameter(1);
      }*/

      for (int i_q = 0; i_q < N_QBINS; i_q++)
      {

        hist_mass = h_mass_q2[i_cen][i_q][i_pt];

        if (!hist_mass)
        {
          std::cerr << "[WARN] Missing q2 grouped histogram for cen=" << cen_name[i_cen] << ", q2=" << i_q << ", pt=" << pt_name[i_pt] << ")\n";
          continue;
        }

        int ymax = hist_mass->GetBinContent(hist_mass->GetMaximumBin());

        hist_mass->SetMinimum(0);
        hist_mass->SetMarkerSize(0.5);
        // hist_mass->SetTitle(Form("Mass_%s_%s_qbin%i", cen_name[i_cen].c_str(), pt_name[i_pt].c_str(),i_q ));
        hist_mass->SetTitle(Form("Mass_%s_%s_qbin%i", cen_name[i_cen], pt_name[i_pt], i_q));

        hist_mass->SetMarkerStyle(20);
        hist_mass->SetLineWidth(1);
        hist_mass->SetOption("e");
        hist_mass->GetXaxis()->SetRangeUser(fit_range_low, fit_range_high);
        hist_mass->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
        hist_mass->GetYaxis()->SetTitle("Entries / 5 MeV");
        hist_mass->SetStats(kFALSE);

        TF1 *f = new TF1(Form("f_ptbin_%d", i_pt), "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])) ) + [8] + [9]*x + [10]*x*x + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12])) + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])) )", fit_range_low, fit_range_high);
        // gROOT->GetListOfFunctions()->Remove(f);

        f->SetLineColor(2);
        f->SetLineWidth(1);

        f->SetParameter(0, 100.);   // [0]  overall normalization (total yield)
        f->SetParameter(1, 1.8648); // [1]  D0 mass mean (GeV/c^2)
        f->SetParameter(2, 0.03);   // [2]  sigma of wide Gaussian (core signal width)
        f->SetParameter(3, 0.005);  // [3]  sigma of narrow Gaussian (tail signal width)
        f->SetParameter(4, 0.1);    // [4]  fraction of narrow Gaussian in signal double-Gaussian
        f->FixParameter(5, 1);      // [5]  signal fraction (1 - swap fraction): signal/(signal+swap)
        f->FixParameter(6, 0);      // [6]  global width smearing factor (sigma scale shift, shared by sig/swap/KK)
        f->FixParameter(7, 0.1);    // [7]  sigma of K-pi swap Gaussian
        f->FixParameter(8, 0);      // [8]  polynomial background p0 (constant term)
        f->FixParameter(9, 0);      // [9]  polynomial background p1 (linear term)
        f->FixParameter(10, 0);     // [10] polynomial background p2 (quadratic term)
        f->FixParameter(11, 0);     // [11] KK/pipi Crystal Ball normalization
        f->FixParameter(12, 0);     // [12] CB mean shift for KK peak (1.96 GeV region)
        f->FixParameter(13, 0);     // [13] CB mean shift for pipi peak (1.7734 GeV region)

        /*double smearing_fixed = 0.0;
        double mean_fixed = 1.8648;
        if (h_def)
        {
          smearing_fixed = f_def->GetParameter(6);
          mean_fixed = f_def->GetParameter(1);
        }*/

        if (i_pt < 5)
        {
          f->SetParLimits(2, 0.01, 0.5);
          f->SetParLimits(3, 0.001, 0.25);
        }
        else
        {
          f->SetParLimits(2, 0.005, 0.15);
          f->SetParLimits(3, 0.001, 0.08);
        }

        f->SetParLimits(4, 0, 1);
        f->SetParLimits(5, 0, 1);

        // fit MC templates to get fixed params

        if (h_mc_match_signal)
        {
          f->FixParameter(1, 1.8648); // for first few attempt fix mean of gaussian to get reasonable estimation of other pars; later open it up

          h_mc_match_signal->Fit(Form("f_ptbin_%d", i_pt), "q", "", fit_range_low, fit_range_high);
          h_mc_match_signal->Fit(Form("f_ptbin_%d", i_pt), "q", "", fit_range_low, fit_range_high);
          // h_mc_match_signal->Fit(Form("f_ptbin_%d",i_pt),"q","",fit_range_low,fit_range_high);
          f->ReleaseParameter(1); // now let gaussian mean float
          h_mc_match_signal->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
          h_mc_match_signal->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
          h_mc_match_signal->Fit(Form("f_ptbin_%d", i_pt), "L q m", "", fit_range_low, fit_range_high);

          // now fix signal double gaussian mean, sigma and gaus1,gaus2 yield ratio
          f->FixParameter(1, f->GetParameter(1));
          f->FixParameter(2, f->GetParameter(2));
          f->FixParameter(3, f->GetParameter(3));
          f->FixParameter(4, f->GetParameter(4));
        }
        if (h_mc_match_all)
        {
          // now release swap bkg parameters to fit signal+swap MC
          f->ReleaseParameter(5);
          f->ReleaseParameter(7);
          // f->ReleaseParameter(8);
          f->SetParameter(7, 0.1);
          h_mc_match_all->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
          h_mc_match_all->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
          h_mc_match_all->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
          h_mc_match_all->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
          h_mc_match_all->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
        }

        // now fix swap bkg parameters to fit data
        f->FixParameter(5, f->GetParameter(5));
        f->FixParameter(7, f->GetParameter(7));
        // Make sure all are fixed except background and scaling
        f->FixParameter(1, f->GetParameter(1));
        f->FixParameter(2, f->GetParameter(2));
        f->FixParameter(3, f->GetParameter(3));
        f->FixParameter(4, f->GetParameter(4));
        f->FixParameter(6, 0);

        // now release poly bkg pars
        f->ReleaseParameter(8);
        f->ReleaseParameter(9);
        f->ReleaseParameter(10);

        if (h_def)
        {
          // f->ReleaseParameter(1);//mean
          f->FixParameter(1, 1.8648);
          f->FixParameter(6, 0);
          h_def->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
          h_def->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
          f->ReleaseParameter(1);
          f->ReleaseParameter(6);
          h_def->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
          h_def->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
          h_def->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);

          f->FixParameter(6, f->GetParameter(6)); // newly added: fix from mass fit inclusive q bin
          f->FixParameter(1, f->GetParameter(1));
        }

        // Fit data
        hist_mass->Fit(Form("f_ptbin_%d", i_pt), "q", "", fit_range_low, fit_range_high);
        hist_mass->Fit(Form("f_ptbin_%d", i_pt), "q", "", fit_range_low, fit_range_high);
        f->ReleaseParameter(1);
        hist_mass->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
        hist_mass->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
        hist_mass->Fit(Form("f_ptbin_%d", i_pt), "L q", "", fit_range_low, fit_range_high);
        hist_mass->Fit(Form("f_ptbin_%d", i_pt), "L q m", "", fit_range_low, fit_range_high);

        if (i_pt != 0)
        {
          f->ReleaseParameter(11);
          f->SetParLimits(11, 0, f->GetParameter(0) * (1 - f->GetParameter(5)));
        }
        else
        {
          f->FixParameter(11, 0);
        }

        hist_mass->Fit(Form("f_ptbin_%d", i_pt), "L q m", "", fit_range_low, fit_range_high);
        hist_mass->Fit(Form("f_ptbin_%d", i_pt), "L q m", "", fit_range_low, fit_range_high);

        f->SetParLimits(12, -0.03, 0.03);
        f->SetParLimits(13, -0.03, 0.03);

        f->ReleaseParameter(12);
        f->ReleaseParameter(13);

        hist_mass->Fit(Form("f_ptbin_%d", i_pt), "L q m", "", fit_range_low, fit_range_high);
        hist_mass->Fit(Form("f_ptbin_%d", i_pt), "L q m", "", fit_range_low, fit_range_high);
        hist_mass->Fit(Form("f_ptbin_%d", i_pt), "L q m", "", fit_range_low, fit_range_high);

        hist_mass->GetYaxis()->SetRangeUser(0, 1.35 * ymax);
        // hist_mass->GetXaxis()->SetRangeUser(fit_range_low, fit_range_high);
        hist_mass->Draw("ep");

        // draw D0 signal separately
        TF1 *f1 = new TF1(Form("f_sig_%d", i_pt), "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
        // gROOT->GetListOfFunctions()->Remove(f1);
        f1->SetLineColor(kOrange - 3);
        f1->SetLineWidth(1);
        f1->SetLineStyle(2);
        f1->SetFillColorAlpha(kOrange - 3, 0.3);
        f1->SetFillStyle(1001);
        f1->FixParameter(0, f->GetParameter(0));
        f1->FixParameter(1, f->GetParameter(1));
        f1->FixParameter(2, f->GetParameter(2));
        f1->FixParameter(3, f->GetParameter(3));
        f1->FixParameter(4, f->GetParameter(4));
        f1->FixParameter(5, f->GetParameter(5));
        f1->FixParameter(6, f->GetParameter(6));

        f1->Draw("LSAME");

        // write Yield uncertainty
        /*TF1* f1e = new TF1(Form("f_e_%d",i_pt),"[0]", fit_range_low, fit_range_high);
        f1e->SetLineColor(kOrange-3);
        f1e->SetLineWidth(1);
        f1e->SetLineStyle(2);
        f1e->SetFillColorAlpha(kOrange-3,0.3);
        f1e->SetFillStyle(1001);
        f1e->FixParameter(0,(f->GetParError(0)*f->GetParameter(5)+f->GetParError(5)*f->GetParameter(0))/(f->GetParameter(0)*f->GetParameter(5)));
        */
        // draw swap separately
        TF1 *f2 = new TF1(Form("f_swap_%d", i_pt), "[0]*((1-[1])*TMath::Gaus(x,[4],[3]*(1.0 +[2]))/(sqrt(2*3.14159)*[3]*(1.0 +[2])))", fit_range_low, fit_range_high);

        f2->SetLineColor(kGreen + 4);
        f2->SetLineWidth(1);
        f2->SetLineStyle(1);
        f2->SetFillColorAlpha(kGreen + 4, 0.3);
        f2->SetFillStyle(1001);
        f2->FixParameter(0, f->GetParameter(0));
        f2->FixParameter(1, f->GetParameter(5));
        f2->FixParameter(2, f->GetParameter(6));
        f2->FixParameter(3, f->GetParameter(7));
        f2->FixParameter(4, f->GetParameter(1));
        f2->Draw("LSAME");

        // draw poly bkg separately
        TF1 *f3 = new TF1(Form("f_bkg_%d", i_pt), "[0] + [1]*x + [2]*x*x ", fit_range_low, fit_range_high);
        f3->SetLineColor(4);
        f3->SetLineWidth(1);
        f3->SetLineStyle(2);
        f3->FixParameter(0, f->GetParameter(8));
        f3->FixParameter(1, f->GetParameter(9));
        f3->FixParameter(2, f->GetParameter(10));
        // f3->FixParameter(3,0);
        f3->Draw("LSAME");

        // draw KK branch separately
        TF1 *f4 = new TF1(Form("f_kk_%d", i_pt), "[0]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[1]), 1.96*(1+[2])) + 4*[0]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[1]), 1.7734*(1+[3])) )", fit_range_low, fit_range_high);
        f4->SetLineColor(kViolet - 4);
        f4->SetLineWidth(1);
        f4->SetLineStyle(1);
        f4->SetFillColorAlpha(kViolet - 4, 0.3);
        f4->SetFillStyle(1001);
        f4->FixParameter(0, f->GetParameter(11));
        f4->FixParameter(1, f->GetParameter(6));
        f4->FixParameter(2, f->GetParameter(12));
        f4->FixParameter(3, f->GetParameter(13));
        f4->Draw("LSAME");

        tex->DrawLatex(0.14, 0.86, Form("cent %s pT %s GeV/c q2 %d", cen_name[i_cen], pt_name[i_pt], i_q));
        // tex->DrawLatex(0.37,0.86,Form("%s GeV/c",pt_label[i_pt].Data()));
        // texCMS->DrawLatex(.18,.97,"#font[61]{CMS} #it{Preliminary}");
        // texCMS->DrawLatex(0.62,0.97, "#scale[0.8]{PbPb #sqrt{s_{NN}} = 5.02 TeV}");

        TLegend *leg = new TLegend(0.65, 0.6, 0.81, 0.9, NULL, "brNDC");
        leg->SetBorderSize(0);
        leg->SetTextSize(0.035);
        leg->SetTextFont(42);
        leg->SetFillStyle(0);
        // leg->AddEntry(h_mass_q2[i_cen][i_q][i_pt]," Data","lep"); //for full q2
        leg->AddEntry(h_mass_q2[i_cen][i_q][i_pt], " Data", "lep");
        leg->AddEntry(f, " Fit", "L");
        leg->AddEntry(f1, " D^{0}+#bar{D^{#lower[0.2]{0}}} Signal", "l");
        leg->AddEntry(f2, " K-#pi swap", "l");
        leg->AddEntry(f4, " K-#bar{K}, #pi-#bar{#pi}", "l");
        leg->AddEntry(f3, " Combinatorial", "l");
        leg->Draw("SAME");

        TString legname = Form("leg_mass_%s_q2bin%d_pt_%s", cen_name[i_cen], i_q, pt_name[i_pt]);
        leg->SetName(legname);

        double yield_val = f->GetParameter(0) * f->GetParameter(5) / width;
        double yield_err = (f->GetParError(0) * f->GetParameter(5) + f->GetParError(5) * f->GetParameter(0)) / width;
        double significance = (yield_err > 0) ? fabs(yield_val / yield_err) : 0.;
        double chi2ndf_mass = (f->GetNDF() > 0) ? f->GetChisquare() / f->GetNDF() : -1.;

        TLatex *texYield = new TLatex(0.14, 0.80, Form("Yield = %.0f #pm %.0f", yield_val, yield_err));
        texYield->SetNDC();
        texYield->SetTextFont(42);
        texYield->SetTextSize(0.035);
        TLatex *texSign = new TLatex(0.14, 0.75, Form("Y/#DeltaY = %.2f", significance));
        texSign->SetNDC();
        texSign->SetTextFont(42);
        texSign->SetTextSize(0.035);
        TLatex *texChi = new TLatex(0.14, 0.70, Form("#chi^{2}/ndf = %.2f", chi2ndf_mass));
        texChi->SetNDC();
        texChi->SetTextFont(42);
        texChi->SetTextSize(0.035);

        hist_mass->GetListOfFunctions()->Add(f);
        hist_mass->GetListOfFunctions()->Add(f1);
        hist_mass->GetListOfFunctions()->Add(f2);
        hist_mass->GetListOfFunctions()->Add(f3);
        hist_mass->GetListOfFunctions()->Add(f4);
        hist_mass->GetListOfFunctions()->Add(leg);
        hist_mass->GetListOfFunctions()->Add(texYield);
        hist_mass->GetListOfFunctions()->Add(texSign);
        hist_mass->GetListOfFunctions()->Add(texChi);

        // Write Mass histograms!
        dir_mass_fitted->cd();
        TH1D *h_mc_sig_clone = nullptr;
        TH1D *h_mc_all_clone = nullptr;

        if (h_mc_match_signal)
        {
          h_mc_sig_clone = (TH1D *)h_mc_match_signal->Clone(Form("%s_%s", h_mc_match_signal->GetName(), cen_name[i_cen]));
          h_mc_sig_clone->SetTitle(Form("%s_%s", h_mc_match_signal->GetTitle(), cen_name[i_cen]));
          h_mc_sig_clone->SetDirectory(nullptr);
        }
        if (h_mc_match_all)
        {
          h_mc_all_clone = (TH1D *)h_mc_match_all->Clone(Form("%s_%s", h_mc_match_all->GetName(), cen_name[i_cen]));
          h_mc_all_clone->SetTitle(Form("%s_%s", h_mc_match_all->GetTitle(), cen_name[i_cen]));
          h_mc_all_clone->SetDirectory(nullptr);
        }

        if (h_mc_sig_clone)
        {
          TString name_sig_c = h_mc_sig_clone->GetName();
          if (!dir_mass_fitted->Get(name_sig_c))
            h_mc_sig_clone->Write();
        }
        if (h_mc_all_clone)
        {
          TString name_all_c = h_mc_all_clone->GetName();
          if (!dir_mass_fitted->Get(name_all_c))
            h_mc_all_clone->Write();
        }

        /*if (h_def)
        {
          dir_mass_fitted->cd();
          TString name_def = h_def->GetName();
          if (!dir_mass_fitted->Get(name_def))
            h_def->Write();
        }*/

        if (h_def)
        {
          TString name_def = h_def->GetName();
          if (!dir_mass_fitted->Get(name_def))
            h_def->Write();
        }

        hist_mass->Write();

        c1->Clear();

        for (int i_v = 0; i_v < N_VBINS_V2; i_v++)
        {

          v2_x[i_v] = 0.5 * (vnbinning_v2[i_cen][i_pt][i_v] + vnbinning_v2[i_cen][i_pt][i_v + 1]);
          v2_x_err[i_v] = 0.5 * (fabs(vnbinning_v2[i_cen][i_pt][i_v + 1] - vnbinning_v2[i_cen][i_pt][i_v]));
          v2_x_store[i_cen][i_pt][i_v] = v2_x[i_v];
          v2_x_err_store[i_cen][i_pt][i_v] = v2_x_err[i_v];

          //=====================
          // For v2 calculations
          //=====================

          TH1D *h_v2 = h_mass_v2_fit[i_cen][i_q][i_pt][i_v];
          if (!h_v2 || h_v2->GetEntries() < 10)
          {
            yield_v2[i_v] = 0;
            yield_error_v2[i_v] = 0;
            chi2_ndf_for_q2_v2[i_cen][i_q][i_pt][i_v] = -10;
            sigma_v2[i_cen][i_q][i_pt][i_v] = -10;
            bad_bin_log << "v2 | " << i_cen << " | " << cen_name[i_cen]
                        << " | " << i_q << " | " << i_pt << " | " << pt_name[i_pt]
                        << " | " << i_v << " | low_stats\n";
            continue;
          }

          TF1 *fitFcn_v2 = new TF1(Form("fit_v2_%d_%d_%d", i_pt, i_q, i_v), "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])) ) + [8] + [9]*x + [10]*x*x + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12])) + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])) )", fit_range_low, fit_range_high);

          fitFcn_v2->SetParameter(0, 100);
          fitFcn_v2->FixParameter(1, f->GetParameter(1));
          fitFcn_v2->FixParameter(2, f->GetParameter(2));
          fitFcn_v2->FixParameter(3, f->GetParameter(3));
          fitFcn_v2->FixParameter(4, f->GetParameter(4));
          fitFcn_v2->FixParameter(5, f->GetParameter(5));
          fitFcn_v2->FixParameter(6, f->GetParameter(6));
          fitFcn_v2->FixParameter(7, f->GetParameter(7));
          fitFcn_v2->SetParameter(8, 1);
          fitFcn_v2->SetParameter(9, 1);
          fitFcn_v2->SetParameter(10, 1);
          fitFcn_v2->FixParameter(11, 0);

          h_v2->Fit(fitFcn_v2, "M", "", fit_range_low, fit_range_high);
          h_v2->Fit(fitFcn_v2, "L Q", "", fit_range_low, fit_range_high);
          h_v2->Fit(fitFcn_v2, "L Q", "", fit_range_low, fit_range_high);
          h_v2->Fit(fitFcn_v2, "L M", "", fit_range_low, fit_range_high);

          c1->cd();
          h_v2->SetStats(kFALSE);
          h_v2->GetXaxis()->SetRangeUser(fit_range_low, fit_range_high);
          h_v2->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
          h_v2->GetYaxis()->SetTitle("Entries / 5 MeV");
          h_v2->SetMarkerStyle(20);
          h_v2->SetMarkerSize(0.5);
          h_v2->SetLineWidth(1);
          h_v2->Draw("ep");

          // Build signal/swap/bkg/KK functions using fixed params from f
          TF1 *fv2_sig = new TF1(Form("fv2_sig_%d_%d_%d", i_pt, i_q, i_v), "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0+[6]))/(sqrt(2*3.14159)*[2]*(1.0+[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0+[6]))/(sqrt(2*3.14159)*[3]*(1.0+[6]))))", fit_range_low, fit_range_high);
          fv2_sig->FixParameter(0, fitFcn_v2->GetParameter(0));
          fv2_sig->FixParameter(1, f->GetParameter(1));
          fv2_sig->FixParameter(2, f->GetParameter(2));
          fv2_sig->FixParameter(3, f->GetParameter(3));
          fv2_sig->FixParameter(4, f->GetParameter(4));
          fv2_sig->FixParameter(5, f->GetParameter(5));
          fv2_sig->FixParameter(6, f->GetParameter(6));
          fv2_sig->SetLineColor(kOrange - 3);
          fv2_sig->SetLineWidth(1);
          fv2_sig->SetLineStyle(2);

          TF1 *fv2_swap = new TF1(Form("fv2_swap_%d_%d_%d", i_pt, i_q, i_v), "[0]*((1-[1])*TMath::Gaus(x,[4],[3]*(1.0+[2]))/(sqrt(2*3.14159)*[3]*(1.0+[2])))", fit_range_low, fit_range_high);
          fv2_swap->FixParameter(0, fitFcn_v2->GetParameter(0));
          fv2_swap->FixParameter(1, f->GetParameter(5));
          fv2_swap->FixParameter(2, f->GetParameter(6));
          fv2_swap->FixParameter(3, f->GetParameter(7));
          fv2_swap->FixParameter(4, f->GetParameter(1));
          fv2_swap->SetLineColor(kGreen + 4);
          fv2_swap->SetLineWidth(1);
          fv2_swap->SetLineStyle(1);

          TF1 *fv2_bkg = new TF1(Form("fv2_bkg_%d_%d_%d", i_pt, i_q, i_v), "[0]+[1]*x+[2]*x*x", fit_range_low, fit_range_high);
          fv2_bkg->FixParameter(0, fitFcn_v2->GetParameter(8));
          fv2_bkg->FixParameter(1, fitFcn_v2->GetParameter(9));
          fv2_bkg->FixParameter(2, fitFcn_v2->GetParameter(10));
          fv2_bkg->SetLineColor(4);
          fv2_bkg->SetLineWidth(1);
          fv2_bkg->SetLineStyle(2);

          TF1 *fv2_kk = new TF1(Form("fv2_kk_%d_%d_%d", i_pt, i_q, i_v), "[0]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[1]), 1.96*(1+[2])) + 4*[0]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[1]), 1.7734*(1+[3])))", fit_range_low, fit_range_high);
          fv2_kk->FixParameter(0, fitFcn_v2->GetParameter(11)); // normalization from v2 fit (fixed=0 if i_pt==0)
          fv2_kk->FixParameter(1, f->GetParameter(6));          // smearing factor from inclusive
          fv2_kk->FixParameter(2, f->GetParameter(12));         // CB mean shift
          fv2_kk->FixParameter(3, f->GetParameter(13));         // CB mean shift pipi
          fv2_kk->SetLineColor(kViolet - 4);
          fv2_kk->SetLineWidth(1);
          fv2_kk->SetLineStyle(1);
          fv2_kk->SetFillColorAlpha(kViolet - 4, 0.3);
          fv2_kk->SetFillStyle(1001);

          fitFcn_v2->SetLineColor(2);
          fitFcn_v2->SetLineWidth(1);
          fitFcn_v2->Draw("LSAME");
          fv2_sig->Draw("LSAME");
          fv2_swap->Draw("LSAME");
          fv2_bkg->Draw("LSAME");
          fv2_kk->Draw("LSAME");

          TLegend *leg_v2 = new TLegend(0.65, 0.6, 0.81, 0.9, NULL, "brNDC");
          leg_v2->SetBorderSize(0);
          leg_v2->SetTextSize(0.035);
          leg_v2->SetTextFont(42);
          leg_v2->SetFillStyle(0);
          leg_v2->AddEntry(h_v2, " Data", "lep");
          leg_v2->AddEntry(fitFcn_v2, " Fit", "L");
          leg_v2->AddEntry(fv2_sig, " D^{0} Signal", "l");
          leg_v2->AddEntry(fv2_swap, " K-#pi swap", "l");
          leg_v2->AddEntry(fv2_kk, " K-#bar{K}, #pi-#bar{#pi}", "l");
          leg_v2->AddEntry(fv2_bkg, " Combinatorial", "l");
          leg_v2->Draw("SAME");

          chi2_ndf_for_q2_v2[i_cen][i_q][i_pt][i_v] = (fitFcn_v2->GetNDF() > 0) ? (1.0 * fitFcn_v2->GetChisquare() / fitFcn_v2->GetNDF()) : -1.0;

          // tex->DrawLatex(0.14, 0.86, Form("cent %s pT %s q2bin %d v2bin %d", cen_name[i_cen], pt_name[i_pt], i_q, i_v));
          double yv2_val = fitFcn_v2->GetParameter(0) * fitFcn_v2->GetParameter(5) / width;
          double yv2_err = fitFcn_v2->GetParError(0) * fitFcn_v2->GetParameter(5) / width;
          TLatex *texCENT = new TLatex(0.14, 0.85, Form("%i < Cent < %i", cen_edges[i_cen], cen_edges[i_cen + 1]));
          texCENT->SetNDC();
          texCENT->SetTextFont(42);
          texCENT->SetTextSize(0.035);
          texCENT->Draw();
          TLatex *texPT = new TLatex(0.14, 0.80, Form("%.1f < p_{T} < %.1f GeV/c", pt_edges[i_pt], pt_edges[i_pt + 1]));
          texPT->SetNDC();
          texPT->SetTextFont(42);
          texPT->SetTextSize(0.035);
          texPT->Draw();
          TLatex *texVbin = new TLatex(0.14, 0.75, Form("%.3f < v_{2}^{i} < %.3f", vnbinning_v2[i_cen][i_pt][i_v], vnbinning_v2[i_cen][i_pt][i_v + 1]));
          texVbin->SetNDC();
          texVbin->SetTextFont(42);
          texVbin->SetTextSize(0.035);
          texVbin->Draw();
          TLatex *texY = new TLatex(0.14, 0.70, Form("Yield = %.0f #pm %.0f", yv2_val, yv2_err));
          texY->SetNDC();
          texY->SetTextFont(42);
          texY->SetTextSize(0.035);
          texY->Draw();
          TLatex *texC = new TLatex(0.14, 0.65, Form("#chi^{2}/ndf = %.2f", chi2_ndf_for_q2_v2[i_cen][i_q][i_pt][i_v]));
          texC->SetNDC();
          texC->SetTextFont(42);
          texC->SetTextSize(0.035);
          texC->Draw();

          c1->SaveAs(Form("prompt_mass_plot_withchi2_sigma/hmassfit_%s_%s_q2bin_%d_v2bin_%d.pdf", pt_name[i_pt], cen_name[i_cen], i_q, i_v));

          yield_v2[i_v] = fitFcn_v2->GetParameter(0) * fitFcn_v2->GetParameter(5) / width;
          yield_error_v2[i_v] = fitFcn_v2->GetParError(0) * fitFcn_v2->GetParameter(5) / width;
          // sigma_v2[i_cen][i_q][i_pt][i_v] = yield_v2[i_v]/yield_error_v2[i_v];
          sigma_v2[i_cen][i_q][i_pt][i_v] = (yield_error_v2[i_v] > 0) ? yield_v2[i_v] / yield_error_v2[i_v] : -10;
          bool v2_zero_error = !(yield_error_v2[i_v] > 0);

          if (yield_v2[i_v] <= 0)
          {
            yield_v2[i_v] = 0;
            yield_error_v2[i_v] = 0;
            chi2_ndf_for_q2_v2[i_cen][i_q][i_pt][i_v] = -10;
            sigma_v2[i_cen][i_q][i_pt][i_v] = -10;
          }
          if (sigma_v2[i_cen][i_q][i_pt][i_v] < 0)
          {
            const char *reason_v2 = v2_zero_error ? "zero_error" : "neg_yield";
            bad_bin_log << "v2 | " << i_cen << " | " << cen_name[i_cen]
            << " | " << i_q << " | " << i_pt << " | " << pt_name[i_pt]
            << " | " << i_v << " | " << reason_v2 << " | "
            << (v2_is_mid(i_v) ? "MID" : "TAIL") << "\n";
          }

          h_v2_hist[i_cen][i_q][i_pt]->SetBinContent(i_v + 1, yield_v2[i_v]);
          h_v2_hist[i_cen][i_q][i_pt]->SetBinError(i_v + 1, yield_error_v2[i_v]);

          c1->Clear();
          delete fitFcn_v2;
          delete fv2_sig;
          delete fv2_swap;
          delete fv2_bkg;
          delete fv2_kk;
          delete leg_v2;
          delete texCENT;
          delete texPT;
          delete texVbin;
          delete texY;
          delete texC;

        } // -- v2 loop --

        // --- Fit h_mass_q3 to derive q3-specific shape parameters for v3 fits ---
        TH1D *hist_mass_q3 = h_mass_q3[i_cen][i_q][i_pt];
        TF1 *f_q3 = (TF1 *)f->Clone(Form("f_q3_%d_%d_%d", i_cen, i_q, i_pt));
        if (hist_mass_q3 && hist_mass_q3->GetEntries() >= 10)
        {
          // Release mean [1], smearing [6], and background [8,9,10] to re-fit from q3-selected mass
          f_q3->ReleaseParameter(1);
          f_q3->ReleaseParameter(6);
          f_q3->ReleaseParameter(8);
          f_q3->ReleaseParameter(9);
          f_q3->ReleaseParameter(10);
          hist_mass_q3->Fit(f_q3, "L Q", "", fit_range_low, fit_range_high);
          hist_mass_q3->Fit(f_q3, "L Q", "", fit_range_low, fit_range_high);
          hist_mass_q3->Fit(f_q3, "L Q m", "", fit_range_low, fit_range_high);
          f_q3->FixParameter(1, f_q3->GetParameter(1));
          f_q3->FixParameter(6, f_q3->GetParameter(6));
        }
        // else: f_q3 retains all params from q2 fit as fallback

        //=====================
        // For v3 calculations
        //=====================
        for (int i_v = 0; i_v < N_VBINS_V3; i_v++)
        {

          v3_x[i_v] = 0.5 * (vnbinning_v3[i_cen][i_pt][i_v] + vnbinning_v3[i_cen][i_pt][i_v + 1]);
          v3_x_err[i_v] = 0.5 * (fabs(vnbinning_v3[i_cen][i_pt][i_v + 1] - vnbinning_v3[i_cen][i_pt][i_v]));
          v3_x_store[i_cen][i_pt][i_v] = v3_x[i_v];
          v3_x_err_store[i_cen][i_pt][i_v] = v3_x_err[i_v];

          TH1D *h_v3 = h_mass_v3_fit[i_cen][i_q][i_pt][i_v];
          if (!h_v3 || h_v3->GetEntries() < 10)
          {
            yield_v3[i_v] = 0;
            yield_error_v3[i_v] = 0;
            chi2_ndf_for_q3_v3[i_cen][i_q][i_pt][i_v] = -10;
            sigma_v3[i_cen][i_q][i_pt][i_v] = -10;
            bad_bin_log << "v3 | " << i_cen << " | " << cen_name[i_cen]
                        << " | " << i_q << " | " << i_pt << " | " << pt_name[i_pt]
                        << " | " << i_v << " | low_stats\n";
            continue;
          }

          TF1 *fitFcn_v3 = new TF1(Form("fit_v3_%d_%d_%d", i_pt, i_q, i_v), "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])) ) + [8] + [9]*x + [10]*x*x + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12])) + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])) )", fit_range_low, fit_range_high);

          fitFcn_v3->SetParameter(0, 100);
          fitFcn_v3->FixParameter(1, f_q3->GetParameter(1)); // from q3 fit
          fitFcn_v3->FixParameter(2, f_q3->GetParameter(2));
          fitFcn_v3->FixParameter(3, f_q3->GetParameter(3));
          fitFcn_v3->FixParameter(4, f_q3->GetParameter(4));
          fitFcn_v3->FixParameter(5, f_q3->GetParameter(5));
          fitFcn_v3->FixParameter(6, f_q3->GetParameter(6)); // from q3 fit
          fitFcn_v3->FixParameter(7, f_q3->GetParameter(7));
          fitFcn_v3->SetParameter(8, 1);
          fitFcn_v3->SetParameter(9, 1);
          fitFcn_v3->SetParameter(10, 1);
          fitFcn_v3->FixParameter(11, 0);

          // h_v3->Draw("AEP");
          h_v3->Fit(fitFcn_v3, "M", "", fit_range_low, fit_range_high);
          h_v3->Fit(fitFcn_v3, "L Q", "", fit_range_low, fit_range_high);
          h_v3->Fit(fitFcn_v3, "L Q", "", fit_range_low, fit_range_high);
          h_v3->Fit(fitFcn_v3, "L M", "", fit_range_low, fit_range_high);

          c1->cd();
          h_v3->SetStats(kFALSE);
          h_v3->GetXaxis()->SetRangeUser(fit_range_low, fit_range_high);
          h_v3->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
          h_v3->GetYaxis()->SetTitle("Entries / 5 MeV");
          h_v3->SetMarkerStyle(20);
          h_v3->SetMarkerSize(0.5);
          h_v3->SetLineWidth(1);
          h_v3->Draw("ep");

          // Signal
          TF1 *fv3_sig = new TF1(Form("fv3_sig_%d_%d_%d", i_pt, i_q, i_v), "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0+[6]))/(sqrt(2*3.14159)*[2]*(1.0+[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0+[6]))/(sqrt(2*3.14159)*[3]*(1.0+[6]))))", fit_range_low, fit_range_high);
          fv3_sig->FixParameter(0, fitFcn_v3->GetParameter(0));
          fv3_sig->FixParameter(1, f_q3->GetParameter(1));
          fv3_sig->FixParameter(2, f_q3->GetParameter(2));
          fv3_sig->FixParameter(3, f_q3->GetParameter(3));
          fv3_sig->FixParameter(4, f_q3->GetParameter(4));
          fv3_sig->FixParameter(5, f_q3->GetParameter(5));
          fv3_sig->FixParameter(6, f_q3->GetParameter(6));
          fv3_sig->SetLineColor(kOrange - 3);
          fv3_sig->SetLineWidth(1);
          fv3_sig->SetLineStyle(2);

          // Swap
          TF1 *fv3_swap = new TF1(Form("fv3_swap_%d_%d_%d", i_pt, i_q, i_v), "[0]*((1-[1])*TMath::Gaus(x,[4],[3]*(1.0+[2]))/(sqrt(2*3.14159)*[3]*(1.0+[2])))", fit_range_low, fit_range_high);
          fv3_swap->FixParameter(0, fitFcn_v3->GetParameter(0));
          fv3_swap->FixParameter(1, f_q3->GetParameter(5));
          fv3_swap->FixParameter(2, f_q3->GetParameter(6));
          fv3_swap->FixParameter(3, f_q3->GetParameter(7));
          fv3_swap->FixParameter(4, f_q3->GetParameter(1));
          fv3_swap->SetLineColor(kGreen + 4);
          fv3_swap->SetLineWidth(1);
          fv3_swap->SetLineStyle(1);

          // Combinatorial background
          TF1 *fv3_bkg = new TF1(Form("fv3_bkg_%d_%d_%d", i_pt, i_q, i_v), "[0]+[1]*x+[2]*x*x", fit_range_low, fit_range_high);
          fv3_bkg->FixParameter(0, fitFcn_v3->GetParameter(8));
          fv3_bkg->FixParameter(1, fitFcn_v3->GetParameter(9));
          fv3_bkg->FixParameter(2, fitFcn_v3->GetParameter(10));
          fv3_bkg->SetLineColor(4);
          fv3_bkg->SetLineWidth(1);
          fv3_bkg->SetLineStyle(2);

          // KK/pipi
          TF1 *fv3_kk = new TF1(Form("fv3_kk_%d_%d_%d", i_pt, i_q, i_v), "[0]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[1]), 1.96*(1+[2])) + 4*[0]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[1]), 1.7734*(1+[3])))", fit_range_low, fit_range_high);
          fv3_kk->FixParameter(0, fitFcn_v3->GetParameter(11));
          fv3_kk->FixParameter(1, f_q3->GetParameter(6));
          fv3_kk->FixParameter(2, f_q3->GetParameter(12));
          fv3_kk->FixParameter(3, f_q3->GetParameter(13));
          fv3_kk->SetLineColor(kViolet - 4);
          fv3_kk->SetLineWidth(1);
          fv3_kk->SetLineStyle(1);
          fv3_kk->SetFillColorAlpha(kViolet - 4, 0.3);
          fv3_kk->SetFillStyle(1001);

          fitFcn_v3->SetLineColor(2);
          fitFcn_v3->SetLineWidth(1);
          fitFcn_v3->Draw("LSAME");
          fv3_sig->Draw("LSAME");
          fv3_swap->Draw("LSAME");
          fv3_bkg->Draw("LSAME");
          fv3_kk->Draw("LSAME");

          TLegend *leg_v3 = new TLegend(0.65, 0.6, 0.81, 0.9, NULL, "brNDC");
          leg_v3->SetBorderSize(0);
          leg_v3->SetTextSize(0.035);
          leg_v3->SetTextFont(42);
          leg_v3->SetFillStyle(0);
          leg_v3->AddEntry(h_v3, " Data", "lep");
          leg_v3->AddEntry(fitFcn_v3, " Fit", "L");
          leg_v3->AddEntry(fv3_sig, " D^{0} Signal", "l");
          leg_v3->AddEntry(fv3_swap, " K-#pi swap", "l");
          leg_v3->AddEntry(fv3_kk, " K-#bar{K}, #pi-#bar{#pi}", "l");
          leg_v3->AddEntry(fv3_bkg, " Combinatorial", "l");
          leg_v3->Draw("SAME");

          chi2_ndf_for_q3_v3[i_cen][i_q][i_pt][i_v] = (fitFcn_v3->GetNDF() > 0) ? (1.0 * fitFcn_v3->GetChisquare() / fitFcn_v3->GetNDF()) : -1.0;

          double yv3_val = fitFcn_v3->GetParameter(0) * fitFcn_v3->GetParameter(5) / width;
          double yv3_err = fitFcn_v3->GetParError(0) * fitFcn_v3->GetParameter(5) / width;

          TLatex *texCENT_v3 = new TLatex(0.14, 0.85, Form("%i < Cent < %i", cen_edges[i_cen], cen_edges[i_cen + 1]));
          texCENT_v3->SetNDC();
          texCENT_v3->SetTextFont(42);
          texCENT_v3->SetTextSize(0.035);
          texCENT_v3->Draw();
          TLatex *texPT_v3 = new TLatex(0.14, 0.80, Form("%.1f < p_{T} < %.1f GeV/c", pt_edges[i_pt], pt_edges[i_pt + 1]));
          texPT_v3->SetNDC();
          texPT_v3->SetTextFont(42);
          texPT_v3->SetTextSize(0.035);
          texPT_v3->Draw();
          TLatex *texVbin_v3 = new TLatex(0.14, 0.75, Form("%.3f < v_{3}^{i} < %.3f", vnbinning_v3[i_cen][i_pt][i_v], vnbinning_v3[i_cen][i_pt][i_v + 1]));
          texVbin_v3->SetNDC();
          texVbin_v3->SetTextFont(42);
          texVbin_v3->SetTextSize(0.035);
          texVbin_v3->Draw();
          TLatex *texY_v3 = new TLatex(0.14, 0.70, Form("Yield = %.0f #pm %.0f", yv3_val, yv3_err));
          texY_v3->SetNDC();
          texY_v3->SetTextFont(42);
          texY_v3->SetTextSize(0.035);
          texY_v3->Draw();
          TLatex *texC_v3 = new TLatex(0.14, 0.65, Form("#chi^{2}/ndf = %.2f", chi2_ndf_for_q3_v3[i_cen][i_q][i_pt][i_v]));
          texC_v3->SetNDC();
          texC_v3->SetTextFont(42);
          texC_v3->SetTextSize(0.035);
          texC_v3->Draw();

          c1->SaveAs(Form("prompt_mass_plot_withchi2_sigma/hmassfit_%s_%s_q3bin_%d_v3bin_%d.pdf", pt_name[i_pt], cen_name[i_cen], i_q, i_v));

          // c1->SaveAs(Form("prompt_mass_plot_withchi2_sigma/hmassfit_pt_%s_cen%s_q2bin_%d_v3bin_%d.pdf", pt_name[i_pt].c_str(), cen_name[i_cen].c_str(), i_q, i_v));

          yield_v3[i_v] = fitFcn_v3->GetParameter(0) * fitFcn_v3->GetParameter(5) / width;
          yield_error_v3[i_v] = fitFcn_v3->GetParError(0) * fitFcn_v3->GetParameter(5) / width;
          // sigma_v3[i_cen][i_q][i_pt][i_v] = yield_v3[i_v]/yield_error_v3[i_v];
          sigma_v3[i_cen][i_q][i_pt][i_v] = (yield_error_v3[i_v] > 0) ? yield_v3[i_v] / yield_error_v3[i_v] : -10;
          bool v3_zero_error = !(yield_error_v3[i_v] > 0);

          if (yield_v3[i_v] <= 0)
          {
            yield_v3[i_v] = 0;
            yield_error_v3[i_v] = 0;
            chi2_ndf_for_q3_v3[i_cen][i_q][i_pt][i_v] = -10;
            sigma_v3[i_cen][i_q][i_pt][i_v] = -10;
          }
          if (sigma_v3[i_cen][i_q][i_pt][i_v] < 0)
          {
            const char *reason_v3 = v3_zero_error ? "zero_error" : "neg_yield";
            bad_bin_log << "v3 | " << i_cen << " | " << cen_name[i_cen]
            << " | " << i_q << " | " << i_pt << " | " << pt_name[i_pt]
            << " | " << i_v << " | " << reason_v3 << " | "
            << (v3_is_mid(i_v) ? "MID" : "TAIL") << "\n";
          }

          h_v3_hist[i_cen][i_q][i_pt]->SetBinContent(i_v + 1, yield_v3[i_v]);
          h_v3_hist[i_cen][i_q][i_pt]->SetBinError(i_v + 1, yield_error_v3[i_v]);

          c1->Clear();
          delete fitFcn_v3;
          delete fv3_sig;
          delete fv3_swap;
          delete fv3_bkg;
          delete fv3_kk;
          delete leg_v3;
          delete texCENT_v3;
          delete texPT_v3;
          delete texVbin_v3;
          delete texY_v3;
          delete texC_v3;

        } // ---- v3 loop ---
        delete f_q3;
        f_q3 = nullptr;

        yield_v2_graph[i_cen][i_q][i_pt] = new TGraphErrors(N_VBINS_V2, v2_x, yield_v2, v2_x_err, yield_error_v2);

        yield_v2_graph[i_cen][i_q][i_pt]->SetNameTitle(
            Form("v2_graph_%s_%s_q2bin%d", pt_name[i_pt], cen_name[i_cen], i_q),
            Form("v2_graph_%s_%s_q2bin%d", pt_name[i_pt], cen_name[i_cen], i_q));
        yield_v3_graph[i_cen][i_q][i_pt] = new TGraphErrors(N_VBINS_V3, v3_x, yield_v3, v3_x_err, yield_error_v3);

        yield_v3_graph[i_cen][i_q][i_pt]->SetNameTitle(
            Form("v3_graph_%s_%s_q3bin%d", pt_name[i_pt], cen_name[i_cen], i_q),
            Form("v3_graph_%s_%s_q3bin%d", pt_name[i_pt], cen_name[i_cen], i_q));
      } //-- NQBINS --

      // delete f_def;
      // f_def = nullptr;
    } // -- PTBIN ---
  } //-- CENTBIN --

  // int target_cen_group = (target1pct >= 0) ? (target1pct / 10) : -1;

  for (int i_cen = 0; i_cen < N_CENTBINS; i_cen++)
  {
    if (target_cent >= 0 && i_cen != target_cen_group)
      continue;

    for (int i_q = 0; i_q < N_QBINS; i_q++)
    {
      for (int i_pt = 0; i_pt < N_PTBINS; i_pt++)
      {
        // if (target_pt >= 0 && i_pt != target_pt) continue;
        pt_x[i_pt] = (pt_edges[i_pt] + pt_edges[i_pt + 1]) / 2;
        pt_x_error[i_pt] = (-1.0 * pt_edges[i_pt] + pt_edges[i_pt + 1]) / 2;

        mean_val_v2[i_pt] = h_v2_hist[i_cen][i_q][i_pt]->GetMean();
        //mean_error_v2[i_pt] = h_v2_hist[i_cen][i_q][i_pt]->GetMeanError();
        mean_error_v2[i_pt] = WeightedMeanError(N_VBINS_V2, v2_x_store[i_cen][i_pt], yield_v2, yield_error_v2, v2_x_err_store[i_cen][i_pt]);

        mean_val_v3[i_pt] = h_v3_hist[i_cen][i_q][i_pt]->GetMean();
        //mean_error_v3[i_pt] = h_v3_hist[i_cen][i_q][i_pt]->GetMeanError();
        mean_error_v3[i_pt] = WeightedMeanError(N_VBINS_V3, v3_x_store[i_cen][i_pt], yield_v3, yield_error_v3, v3_x_err_store[i_cen][i_pt]);

        chi2_ndf_for_q2_v2_graph[i_cen][i_q][i_pt] = new TGraph(N_VBINS_V2, v2_x_store[i_cen][i_pt], chi2_ndf_for_q2_v2[i_cen][i_q][i_pt]);

        chi2_ndf_for_q2_v2_graph[i_cen][i_q][i_pt]->SetNameTitle(
            Form("chi2ndf_graph_%s_%s_q2bin_%d", pt_name[i_pt], cen_name[i_cen], i_q),
            Form("chi2ndf_graph_%s_%s_q2bin_%d", pt_name[i_pt], cen_name[i_cen], i_q));

        // chi2_ndf_for_q2_v2_graph[i_cen][i_q][i_pt]->Write();
        chi2_ndf_for_q2_v2_graph[i_cen][i_q][i_pt]->GetXaxis()->SetTitle("v2");
        chi2_ndf_for_q2_v2_graph[i_cen][i_q][i_pt]->GetYaxis()->SetTitle("chi2/ndf");
        chi2_ndf_for_q2_v2_graph[i_cen][i_q][i_pt]->SetMarkerStyle(20); // <-- add
        chi2_ndf_for_q2_v2_graph[i_cen][i_q][i_pt]->SetMarkerSize(0.8); // <-- add
        chi2_ndf_for_q2_v2_graph[i_cen][i_q][i_pt]->SetDrawOption("ACP");
        // chi2_ndf_for_q2_v2_graph[i_cen][i_q][i_pt]->Draw("ACP");
        // c1->SaveAs(Form("prompt_sigma_chi2_plot/chi2_pt_%s_cen%s_q2bin_%d.pdf", pt_name[i_pt].c_str(), cen_name[i_cen].c_str(), i_q));

        chi2_ndf_for_q3_v3_graph[i_cen][i_q][i_pt] = new TGraph(N_VBINS_V3, v3_x_store[i_cen][i_pt], chi2_ndf_for_q3_v3[i_cen][i_q][i_pt]);

        chi2_ndf_for_q3_v3_graph[i_cen][i_q][i_pt]->SetNameTitle(
            Form("chi2ndf_graph_%s_%s_q3bin_%d", pt_name[i_pt], cen_name[i_cen], i_q),
            Form("chi2ndf_graph_%s_%s_q3bin_%d", pt_name[i_pt], cen_name[i_cen], i_q));

        // chi2_ndf_for_q3_v3_graph[i_cen][i_q][i_pt]->Write();
        chi2_ndf_for_q3_v3_graph[i_cen][i_q][i_pt]->GetXaxis()->SetTitle("v3");
        chi2_ndf_for_q3_v3_graph[i_cen][i_q][i_pt]->GetYaxis()->SetTitle("chi2/ndf");
        chi2_ndf_for_q3_v3_graph[i_cen][i_q][i_pt]->SetMarkerStyle(20);
        chi2_ndf_for_q3_v3_graph[i_cen][i_q][i_pt]->SetMarkerSize(0.8);
        chi2_ndf_for_q3_v3_graph[i_cen][i_q][i_pt]->SetDrawOption("ACP");
        // chi2_ndf_for_q3_v3_graph[i_cen][i_q][i_pt]->Draw("ACP");
        // c1->SaveAs(Form("prompt_sigma_chi2_plot/chi2_pt_%s_cen%s_q3bin_%d.pdf", pt_name[i_pt].c_str(), cen_name[i_cen].c_str(), i_q));

        sigma_v2_fit[i_cen][i_q][i_pt] = new TGraph(N_VBINS_V2, v2_x_store[i_cen][i_pt], sigma_v2[i_cen][i_q][i_pt]);
        sigma_v2_fit[i_cen][i_q][i_pt]->SetNameTitle(
            Form("sigma_v2_graph_%s_%s_q2bin_%d", pt_name[i_pt], cen_name[i_cen], i_q),
            Form("sigma_v2_graph_%s_%s_q2bin_%d", pt_name[i_pt], cen_name[i_cen], i_q));

        // sigma_v2_fit[i_cen][i_q][i_pt]->Write();
        sigma_v2_fit[i_cen][i_q][i_pt]->GetXaxis()->SetTitle("v2");
        sigma_v2_fit[i_cen][i_q][i_pt]->GetYaxis()->SetTitle("Significance");
        sigma_v2_fit[i_cen][i_q][i_pt]->SetMarkerStyle(20);
        sigma_v2_fit[i_cen][i_q][i_pt]->SetMarkerSize(0.8);
        sigma_v2_fit[i_cen][i_q][i_pt]->SetDrawOption("ACP");
        // sigma_v2_fit[i_cen][i_q][i_pt]->Draw("ACP");
        // c1->SaveAs(Form("prompt_sigma_chi2_plot/sigma_pt_%s_cen%s_q2bin_%d.pdf", pt_name[i_pt].c_str(), cen_name[i_cen].c_str(), i_q));

        sigma_v3_fit[i_cen][i_q][i_pt] = new TGraph(N_VBINS_V3, v3_x_store[i_cen][i_pt], sigma_v3[i_cen][i_q][i_pt]);

        sigma_v3_fit[i_cen][i_q][i_pt]->SetNameTitle(
            Form("sigma_v3_graph_%s_%s_q3bin_%d", pt_name[i_pt], cen_name[i_cen], i_q),
            Form("sigma_v3_graph_%s_%s_q3bin_%d", pt_name[i_pt], cen_name[i_cen], i_q));

        // sigma_v3_fit[i_cen][i_q][i_pt]->Write();
        sigma_v3_fit[i_cen][i_q][i_pt]->GetXaxis()->SetTitle("v3");
        sigma_v3_fit[i_cen][i_q][i_pt]->GetYaxis()->SetTitle("Significance");
        sigma_v3_fit[i_cen][i_q][i_pt]->SetMarkerStyle(20);
        sigma_v3_fit[i_cen][i_q][i_pt]->SetMarkerSize(0.8);
        sigma_v3_fit[i_cen][i_q][i_pt]->SetDrawOption("ACP");
        // sigma_v3_fit[i_cen][i_q][i_pt]->Draw("ACP");
        // c1->SaveAs(Form("prompt_sigma_chi2_plot/sigma_pt_%s_cen%s_q3bin_%d.pdf", pt_name[i_pt].c_str(), cen_name[i_cen].c_str(), i_q));
      }
      v2_cen[i_cen][i_q] = new TGraphErrors(N_PTBINS, pt_x, mean_val_v2, pt_x_error, mean_error_v2);
      // v2_cen[i_cen][i_q]->SetNameTitle(Form("v2_graph_%s_q2bin_%d", cen_name[i_cen].c_str(), i_q), Form("v2_graph_%s_q2bin_%d", cen_name[i_cen].c_str(), i_q));
      v2_cen[i_cen][i_q]->SetNameTitle(
          Form("v2_graph_%s_q2bin_%d", cen_name[i_cen], i_q),
          Form("v2_graph_%s_q2bin_%d", cen_name[i_cen], i_q));

      v3_cen[i_cen][i_q] = new TGraphErrors(N_PTBINS, pt_x, mean_val_v3, pt_x_error, mean_error_v3);
      // v3_cen[i_cen][i_q]->SetNameTitle(Form("v3_graph_%s_q3bin_%d", cen_name[i_cen].c_str(), i_q), Form("v3_graph_%s_q3bin_%d", cen_name[i_cen].c_str(), i_q));
      v3_cen[i_cen][i_q]->SetNameTitle(
          Form("v3_graph_%s_q3bin_%d", cen_name[i_cen], i_q),
          Form("v3_graph_%s_q3bin_%d", cen_name[i_cen], i_q)); // no .c_str()

    } // -- QBIN ---

    //-----------------------------
    // v2/v3 vs pT inclusive q bin
    //-----------------------------

    for (int i_pt = 0; i_pt < N_PTBINS; i_pt++)
    {
      // --- inverse-variance weighted ---
      double sum_v2 = 0, sum_v2_err2 = 0;
      double sum_v3 = 0, sum_v3_err2 = 0;
      // --- simple average ---
      double sum_v2_s = 0, sum_v2_err2_s = 0;
      double sum_v3_s = 0, sum_v3_err2_s = 0;
      int n_v2 = 0, n_v3 = 0;

      for (int i_q = 0; i_q < N_QBINS; i_q++)
      {
        if (v2_cen[i_cen][i_q])
        {
          double y, x;
          v2_cen[i_cen][i_q]->GetPoint(i_pt, x, y);
          double ey = v2_cen[i_cen][i_q]->GetErrorY(i_pt);
          if (ey > 0)
          {
            // inverse-variance
            double w = 1.0 / (ey * ey);
            sum_v2 += w * y;
            sum_v2_err2 += w;
            // simple
            sum_v2_s += y;
            sum_v2_err2_s += ey * ey;
            n_v2++;
          }
        }
        if (v3_cen[i_cen][i_q])
        {
          double y, x;
          v3_cen[i_cen][i_q]->GetPoint(i_pt, x, y);
          double ey = v3_cen[i_cen][i_q]->GetErrorY(i_pt);
          if (ey > 0)
          {
            // inverse-variance
            double w = 1.0 / (ey * ey);
            sum_v3 += w * y;
            sum_v3_err2 += w;
            // simple
            sum_v3_s += y;
            sum_v3_err2_s += ey * ey;
            n_v3++;
          }
        }
      } //----QBIN ----
      // inverse-variance weighted mean
      mean_v2_incl[i_pt] = (sum_v2_err2 > 0) ? sum_v2 / sum_v2_err2 : 0;
      mean_v2_incl_err[i_pt] = (sum_v2_err2 > 0) ? 1.0 / sqrt(sum_v2_err2) : 0;
      mean_v3_incl[i_pt] = (sum_v3_err2 > 0) ? sum_v3 / sum_v3_err2 : 0;
      mean_v3_incl_err[i_pt] = (sum_v3_err2 > 0) ? 1.0 / sqrt(sum_v3_err2) : 0;
      // simple average: mean = sum/N, error = sqrt(sum_err2)/N
      mean_v2_incl_simple[i_pt] = (n_v2 > 0) ? sum_v2_s / n_v2 : 0;
      mean_v2_incl_err_simple[i_pt] = (n_v2 > 0) ? sqrt(sum_v2_err2_s) / n_v2 : 0;
      mean_v3_incl_simple[i_pt] = (n_v3 > 0) ? sum_v3_s / n_v3 : 0;
      mean_v3_incl_err_simple[i_pt] = (n_v3 > 0) ? sqrt(sum_v3_err2_s) / n_v3 : 0;
    } // --- PTBIN ----

    // inverse-variance weighted graphs
    v2_cen_inclusive[i_cen] = new TGraphErrors(N_PTBINS, pt_x, mean_v2_incl, pt_x_error, mean_v2_incl_err);
    v2_cen_inclusive[i_cen]->SetNameTitle(
        Form("v2_graph_%s_qinclusive", cen_name[i_cen]),
        Form("v2_graph_%s_qinclusive", cen_name[i_cen]));
    v3_cen_inclusive[i_cen] = new TGraphErrors(N_PTBINS, pt_x, mean_v3_incl, pt_x_error, mean_v3_incl_err);
    v3_cen_inclusive[i_cen]->SetNameTitle(
        Form("v3_graph_%s_qinclusive", cen_name[i_cen]),
        Form("v3_graph_%s_qinclusive", cen_name[i_cen]));

    // simple average graphs
    v2_cen_inclusive_simple[i_cen] = new TGraphErrors(N_PTBINS, pt_x, mean_v2_incl_simple, pt_x_error, mean_v2_incl_err_simple);
    v2_cen_inclusive_simple[i_cen]->SetNameTitle(
        Form("v2_graph_%s_qinclusive_simpleavg", cen_name[i_cen]),
        Form("v2_graph_%s_qinclusive_simpleavg", cen_name[i_cen]));
    v3_cen_inclusive_simple[i_cen] = new TGraphErrors(N_PTBINS, pt_x, mean_v3_incl_simple, pt_x_error, mean_v3_incl_err_simple);
    v3_cen_inclusive_simple[i_cen]->SetNameTitle(
        Form("v3_graph_%s_qinclusive_simpleavg", cen_name[i_cen]),
        Form("v3_graph_%s_qinclusive_simpleavg", cen_name[i_cen]));

  } // --- CENT LOOP ---

  // write output ROOT file

  // build output name including optional cent/pt and SLURM ids

  outf->cd();

  for (int i_cen = 0; i_cen < N_CENTBINS; i_cen++)
  {
    if (target_cent >= 0 && i_cen != target_cen_group)
      continue;

    for (int i_q = 0; i_q < N_QBINS; i_q++)
    {
      for (int i_pt = 0; i_pt < N_PTBINS; i_pt++)
      {
        // if (target_pt >= 0 && i_pt != target_pt) continue;
        dir_yield_vn->cd();
        if (yield_v2_graph[i_cen][i_q][i_pt])
          yield_v2_graph[i_cen][i_q][i_pt]->Write();
        if (yield_v3_graph[i_cen][i_q][i_pt])
          yield_v3_graph[i_cen][i_q][i_pt]->Write();

        dir_hist_v2v3dist->cd();
        if (h_v2dist_grouped[i_cen][i_q][i_pt])
          h_v2dist_grouped[i_cen][i_q][i_pt]->Write();
        if (h_v3dist_grouped[i_cen][i_q][i_pt])
          h_v3dist_grouped[i_cen][i_q][i_pt]->Write();

        //-- SP v2/v3 (w/o avg) histograms
        dir_hist_sp->cd();
        if (h_v2_hist[i_cen][i_q][i_pt])
          h_v2_hist[i_cen][i_q][i_pt]->Write();
        if (h_v3_hist[i_cen][i_q][i_pt])
          h_v3_hist[i_cen][i_q][i_pt]->Write();

        // write quality / diagnostics graphs
        dir_diag_chi2_sigma->cd();
        if (chi2_ndf_for_q2_v2_graph[i_cen][i_q][i_pt])
          chi2_ndf_for_q2_v2_graph[i_cen][i_q][i_pt]->Write();
        if (chi2_ndf_for_q3_v3_graph[i_cen][i_q][i_pt])
          chi2_ndf_for_q3_v3_graph[i_cen][i_q][i_pt]->Write();
        if (sigma_v2_fit[i_cen][i_q][i_pt])
          sigma_v2_fit[i_cen][i_q][i_pt]->Write();
        if (sigma_v3_fit[i_cen][i_q][i_pt])
          sigma_v3_fit[i_cen][i_q][i_pt]->Write();
      }

      dir_summary->cd();
      //--- v2/v3 vs pt graphs for each q2/q3 bin [0..9] ---
      if (v2_cen[i_cen][i_q])
        v2_cen[i_cen][i_q]->Write();
      if (v3_cen[i_cen][i_q])
        v3_cen[i_cen][i_q]->Write();
    }
    //--- v2/v3 vs pt graphs for inclusive q2/q3 bin ---
    dir_summary->cd();
    if (v2_cen_inclusive[i_cen])
      v2_cen_inclusive[i_cen]->Write();
    if (v3_cen_inclusive[i_cen])
      v3_cen_inclusive[i_cen]->Write();
    if (v2_cen_inclusive_simple[i_cen])
      v2_cen_inclusive_simple[i_cen]->Write();
    if (v3_cen_inclusive_simple[i_cen])
      v3_cen_inclusive_simple[i_cen]->Write();
  }

  // outf->Write();
  outf->Close();

  delete tex;
  tex = nullptr;
  delete c1;
  c1 = nullptr;

  nt_Data_file->Close();
  inf_MC->Close();
  return 0;
}

int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    fit_mass_and_flow();
  }
  else if (argc == 2)
  {
    int idx = atoi(argv[1]);
    fit_mass_and_flow(idx);
  }
  else
  {
    std::cout << "Usage: ./fit_mass_and_flow [target_cent]" << std::endl;
    return 1;
  }
  return 0;
}
