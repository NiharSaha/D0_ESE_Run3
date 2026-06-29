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

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Quantiles/for_20Qbin/BDTHandler.h"
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Quantiles/for_20Qbin/BDTHandler.cc"
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Quantiles/for_20Qbin/quantile_12NUQbin_diffq2q3_cuts_2023_MB0to31.h"

using namespace std;

void flow_Analysis_latest(TString input_txt, TString output_path, int istart, int iend)
{

  BDTHandler bdtHandler;

  const int N_QBINS = 12;
  const int N_CENTBINS_1 = 50;

  const double MAX_PT_ANA = 100.0;
  const double MIN_PT_ANA = 1.0;
  const double MAX_Y_ANA = 2.4;

  // auto file_res = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/D0_Resolution_20Qbin_MB0to31_Jun16_FullStat/ROOT/Resolution_20Qbin_out_combined.root");
  auto file_res = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/D0_Resolution_12NUQbin_diffq2q3_MB0to31_Jun24_FullStat/ROOT/Resolution_12NUQbin_diffq2q3_out_combined.root");

  Double_t v2_den_Dy_plus[N_CENTBINS_1][N_QBINS];
  Double_t v2_den_Dy_minus[N_CENTBINS_1][N_QBINS];

  Double_t v3_den_Dy_plus[N_CENTBINS_1][N_QBINS];
  Double_t v3_den_Dy_minus[N_CENTBINS_1][N_QBINS];

  Float_t q2_hf_total, q3_hf_total;

  TH1D *hQ2Q2_HFmHFp[N_CENTBINS_1][N_QBINS];
  TH1D *hQ2Q2_HFmTrk[N_CENTBINS_1][N_QBINS];
  TH1D *hQ2Q2_HFpTrk[N_CENTBINS_1][N_QBINS];
  TH1D *hQ3Q3_HFmHFp[N_CENTBINS_1][N_QBINS];
  TH1D *hQ3Q3_HFmTrk[N_CENTBINS_1][N_QBINS];
  TH1D *hQ3Q3_HFpTrk[N_CENTBINS_1][N_QBINS];
  TString hname;

  for (Int_t i = 0; i < N_CENTBINS_1; i++)
  {
    for (Int_t j = 0; j < N_QBINS; j++)
    {
      hname = "Q2Q2_HFmHFp_Re_cen" + to_string(i) + "_" + to_string(i + 1) + "_q2bin_" + to_string(j);
      hQ2Q2_HFmHFp[i][j] = (TH1D *)file_res->Get(hname);
      hname = "Q2Q2_HFmTrk_Re_cen" + to_string(i) + "_" + to_string(i + 1) + "_q2bin_" + to_string(j);
      hQ2Q2_HFmTrk[i][j] = (TH1D *)file_res->Get(hname);
      hname = "Q2Q2_HFpTrk_Re_cen" + to_string(i) + "_" + to_string(i + 1) + "_q2bin_" + to_string(j);
      hQ2Q2_HFpTrk[i][j] = (TH1D *)file_res->Get(hname);

      hname = "Q3Q3_HFmHFp_Re_cen" + to_string(i) + "_" + to_string(i + 1) + "_q3bin_" + to_string(j);
      hQ3Q3_HFmHFp[i][j] = (TH1D *)file_res->Get(hname);
      hname = "Q3Q3_HFmTrk_Re_cen" + to_string(i) + "_" + to_string(i + 1) + "_q3bin_" + to_string(j);
      hQ3Q3_HFmTrk[i][j] = (TH1D *)file_res->Get(hname);
      hname = "Q3Q3_HFpTrk_Re_cen" + to_string(i) + "_" + to_string(i + 1) + "_q3bin_" + to_string(j);
      hQ3Q3_HFpTrk[i][j] = (TH1D *)file_res->Get(hname);

      v2_den_Dy_plus[i][j] = TMath::Sqrt((hQ2Q2_HFmHFp[i][j]->GetMean() * hQ2Q2_HFmTrk[i][j]->GetMean()) / (hQ2Q2_HFpTrk[i][j]->GetMean()));
      v2_den_Dy_minus[i][j] = TMath::Sqrt((hQ2Q2_HFmHFp[i][j]->GetMean() * hQ2Q2_HFpTrk[i][j]->GetMean()) / (hQ2Q2_HFmTrk[i][j]->GetMean()));

      v3_den_Dy_plus[i][j] = TMath::Sqrt((hQ3Q3_HFmHFp[i][j]->GetMean() * hQ3Q3_HFmTrk[i][j]->GetMean()) / (hQ3Q3_HFpTrk[i][j]->GetMean()));
      v3_den_Dy_minus[i][j] = TMath::Sqrt((hQ3Q3_HFmHFp[i][j]->GetMean() * hQ3Q3_HFpTrk[i][j]->GetMean()) / (hQ3Q3_HFmTrk[i][j]->GetMean()));
    }
  }
  file_res->Close();

  ifstream file_stream(input_txt.Data());
  TString outfile = TString::Format("%s/ROOT/flow_Analysis_%d_%d.root", output_path.Data(), istart, iend);
  TFile *fout = new TFile(outfile, "RECREATE");
  auto nt = new TNtuple("nt", "nt", "cent:q2_hf_total:q3_hf_total:pT:y:dca:mass:v2:v3");


  // --- Event-level counters per (1% cent bin, q2/q3 iq bin) ---
  Long64_t n_evts_q2[N_CENTBINS_1][N_QBINS] = {};
  Long64_t n_evts_q3[N_CENTBINS_1][N_QBINS] = {};

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

    std::cout << ">>> Processing ifile=" << ifile << " : " << filename << std::endl;

    TTree *tree = (TTree *)fin->Get("d0Analyzer/VCNtuple_D02kpi");
    TTree *t_eventinfoana = (TTree *)fin->Get("eventinfoana/EventInfoNtuple");
    tree->AddFriend(t_eventinfoana);

    // Standard Variable setup
    Int_t centrality;
    Int_t candSize;
    Float_t ephfmQ[3], ephfpQ[3], ephfmSumW[3], ephfpSumW[3], ephfpAngle[3], ephfmAngle[3], eptkAngle[2], eptkQ[2];
    vector<float> *pT = 0, *phi = 0, *mass = 0, *y = 0, *eta = 0, *mva = 0, *dca = 0;

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
    for (const auto &p : {"candSize", "centrality", "ephfpAngle", "ephfmAngle", "ephfpQ", "ephfmQ", "eptkAngle", "eptkQ", "ephfmSumW", "ephfpSumW", "y", "eta", "pT", "phi", "mass", "mva", "twoTrackDCA"})
      tree->SetBranchStatus(p, 1);

    Int_t n_entries = tree->GetEntries();

    for (Long64_t ii = 0; ii < n_entries; ii++)
    {
      tree->GetEntry(ii);

      if (ii % 100000 == 0)
        printf("Processing entry %d of %d : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);

      Int_t cent = centrality / 2;
      if (cent < 0 || cent >= N_CENTBINS_1)
        continue;

      TComplex aux_Q2_HFm(ephfmQ[1] * TMath::Cos(2.0 * ephfmAngle[1]), ephfmQ[1] * TMath::Sin(2.0 * ephfmAngle[1]), 0);
      TComplex aux_Q2_HFp(ephfpQ[1] * TMath::Cos(2.0 * ephfpAngle[1]), ephfpQ[1] * TMath::Sin(2.0 * ephfpAngle[1]), 0);
      TComplex aux_Q2_Trk(eptkQ[0] * TMath::Cos(2.0 * eptkAngle[0]), eptkQ[0] * TMath::Sin(2.0 * eptkAngle[0]), 0);

      TComplex aux_Q3_HFm(ephfmQ[2] * TMath::Cos(3.0 * ephfmAngle[2]), ephfmQ[2] * TMath::Sin(3.0 * ephfmAngle[2]), 0);
      TComplex aux_Q3_HFp(ephfpQ[2] * TMath::Cos(3.0 * ephfpAngle[2]), ephfpQ[2] * TMath::Sin(3.0 * ephfpAngle[2]), 0);
      TComplex aux_Q3_Trk(eptkQ[1] * TMath::Cos(3.0 * eptkAngle[1]), eptkQ[1] * TMath::Sin(3.0 * eptkAngle[1]), 0);

      float sumW2 = ephfpSumW[1] + ephfmSumW[1];
      float sumW3 = ephfpSumW[2] + ephfmSumW[2];

      if (sumW2 <= 0 || sumW3 <= 0)
        continue;
      q2_hf_total = (ephfpQ[1] + ephfmQ[1]) / sumW2;
      q3_hf_total = (ephfpQ[2] + ephfmQ[2]) / sumW3;

      int q2_idx = -1;
      int q3_idx = -1;
      int i_cen_match = -1;

      i_cen_match = cent;

      q2_idx = TMath::BinarySearch(N_QBINS + 1, q2_cuts[cent], (double)q2_hf_total);
      q3_idx = TMath::BinarySearch(N_QBINS + 1, q3_cuts[cent], (double)q3_hf_total);

      // BinarySearch returns -1 if value < first edge: clamp to valid range
      q2_idx = std::max(0, std::min(N_QBINS - 1, (int)q2_idx));
      q3_idx = std::max(0, std::min(N_QBINS - 1, (int)q3_idx));

      // --- Count unique events per (cent, q2_idx) and (cent, q3_idx) ---
      n_evts_q2[cent][q2_idx]++;
      n_evts_q3[cent][q3_idx]++;

      // 3. Pre-fetch Resolution Factors for this event
      Double_t res2_plus = v2_den_Dy_plus[i_cen_match][q2_idx];
      Double_t res2_minus = v2_den_Dy_minus[i_cen_match][q2_idx];
      Double_t res3_plus = v3_den_Dy_plus[i_cen_match][q3_idx];
      Double_t res3_minus = v3_den_Dy_minus[i_cen_match][q3_idx];

      for (int icand = 0; icand < candSize; ++icand)
      {

        float Y = y->at(icand);
        float Pt = pT->at(icand);
        float Mass = mass->at(icand);
        float Phi = phi->at(icand);
        float Mva = mva->at(icand);
        float Dca = dca->at(icand);

        if (Pt >= MAX_PT_ANA || Pt < MIN_PT_ANA || fabs(Y) >= MAX_Y_ANA)
          continue;

        // I need to use 2*cent because of format
        double bdt_cut = bdtHandler.getBDTCut(Y, 2 * cent, Pt);
        if (Mva <= bdt_cut)
          continue;

        TComplex aux_Q2_D0(TMath::Cos(2.0 * Phi), TMath::Sin(2.0 * Phi), 0);
        TComplex aux_Q3_D0(TMath::Cos(3.0 * Phi), TMath::Sin(3.0 * Phi), 0);

        Double_t v2_obs, v3_obs;
        Double_t res2_final, res3_final;

        if (Y > 0)
        {
          v2_obs = (aux_Q2_D0 * TComplex::Conjugate(aux_Q2_HFm)).Re();
          v3_obs = (aux_Q3_D0 * TComplex::Conjugate(aux_Q3_HFm)).Re();
          res2_final = res2_plus;
          res3_final = res3_plus;
        }
        else
        {
          v2_obs = (aux_Q2_D0 * TComplex::Conjugate(aux_Q2_HFp)).Re();
          v3_obs = (aux_Q3_D0 * TComplex::Conjugate(aux_Q3_HFp)).Re();
          res2_final = res2_minus;
          res3_final = res3_minus;
        }

        if (res2_final > 0 && res3_final > 0)
        {
          nt->Fill(cent, q2_hf_total, q3_hf_total, Pt, Y, Dca, Mass, v2_obs / res2_final, v3_obs / res3_final);
        }

      } // -- cand loop --
    } // -- Event loop ---

    fin->Close();
    delete fin;
    ifile++;
  } // -- file loop ---

  //============================
  // Write all histograms!!!
  //============================

  fout->cd();
  nt->Write();

  fout->cd();
  fout->Write(0, TObject::kOverwrite);
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

    flow_Analysis_latest(input_txt, output_path, istart, iend);
  }
  else
  {
    std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
    return 1;
  }
  return 0;
}
