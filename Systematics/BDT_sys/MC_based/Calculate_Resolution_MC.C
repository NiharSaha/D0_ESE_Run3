// =============================================================================
//  Calculate_Resolution_MC.C
//
//  Scalar-product 3-sub-event resolution ingredients from the prompt-D0 MC,
//  for the BDT systematic study.  Structured exactly like the data /
//  charged-particle resolution step (charge_vn/for_20Qbin/Calculate_Resolution_qbin.C):
//  a standalone job that runs over a slice of the file list and writes the
//  <Q_A . Q_B> histograms; the per-job files are then hadd-ed.
//
//  Differences vs the data version:
//    * q2/q3 INCLUSIVE   (no ESE quantile binning)
//    * centrality INCLUSIVE over the whole MC sample (all centralities)
//      (centrality is not physically meaningful in MC -> pool for statistics;
//       the same single resolution is used for every centrality class later)
//    * each event is weighted by  genWeight x nested-sample stitch weight
//      (USE_MC_WEIGHT in Analysis_bin_BDT.h; identical weight is used for the
//       flow numerator, so it cancels in the BDT systematic ratio)
//
//  Output histograms (one each, weighted):
//     Q2Q2_HFmHFp_Re  Q2Q2_HFmTrk_Re  Q2Q2_HFpTrk_Re
//     Q3Q3_HFmHFp_Re  Q3Q3_HFmTrk_Re  Q3Q3_HFpTrk_Re
//
//  Build / run:
//     g++ Calculate_Resolution_MC.C $(root-config --cflags --libs) -Wall -O2 -o <exe>
//     ./<exe>  <input_txt>  <output_path>  <istart>  <iend>
//  output:  <output_path>/ROOT/Resolution_MC_<istart>_<iend>.root
// =============================================================================

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TH1.h"
#include "TMath.h"
#include "TComplex.h"

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/BDT_sys/MC_based/Analysis_bin_BDT.h"

using namespace std;

void Calculate_Resolution_MC(TString input_txt, TString output_path, int istart, int iend)
{
  TH1::SetDefaultSumw2();
  // include under/overflow in GetMean(): the <Q.Q> distributions have long tails
  // and we must not bias the resolution by dropping the tail entries.
  TH1::StatOverflows(kTRUE);

  gSystem->Exec(TString::Format("mkdir -p %s/ROOT", output_path.Data()));
  TString outfile = TString::Format("%s/ROOT/Resolution_MC_%d_%d.root",
                                    output_path.Data(), istart, iend);

  // ~10 units / bin, range widened w.r.t. the data histograms to keep the tails
  // on-histogram (GetMean is exact regardless, thanks to StatOverflows above).
  TH1D *hQ2_HFmHFp = new TH1D("Q2Q2_HFmHFp_Re", "", 12000, -60000, 60000);
  TH1D *hQ2_HFmTrk = new TH1D("Q2Q2_HFmTrk_Re", "", 12000, -60000, 60000);
  TH1D *hQ2_HFpTrk = new TH1D("Q2Q2_HFpTrk_Re", "", 12000, -60000, 60000);
  TH1D *hQ3_HFmHFp = new TH1D("Q3Q3_HFmHFp_Re", "", 8000, -40000, 40000);
  TH1D *hQ3_HFmTrk = new TH1D("Q3Q3_HFmTrk_Re", "", 8000, -40000, 40000);
  TH1D *hQ3_HFpTrk = new TH1D("Q3Q3_HFpTrk_Re", "", 8000, -40000, 40000);
  TH1D *h_nevt = new TH1D("h_nevt", "weighted events used;;", 1, 0, 1);

  ifstream file_stream(input_txt.Data());
  string filename;
  int ifile = 0;

  while (file_stream >> filename)
  {
    if (ifile < istart) { ifile++; continue; }
    if (ifile >= iend) break;

    TFile *fin = TFile::Open(filename.c_str());
    if (!fin || fin->IsZombie())
    {
      cout << "Warning: skipping bad file: " << filename << endl;
      if (fin) { fin->Close(); delete fin; }
      ifile++;
      continue;
    }
    cout << ">>> Resolution pass: ifile=" << ifile << " : " << filename << endl;

    TTree *tree      = (TTree *)fin->Get("d0Analyzer/VCNtuple_D02kpi");
    TTree *tree_gen  = (TTree *)fin->Get("d0Analyzer/AllGensNtuple");
    TTree *t_evtinfo = (TTree *)fin->Get("eventinfoana/EventInfoNtuple");
    if (!tree || !t_evtinfo || (USE_MC_WEIGHT && !tree_gen))
    {
      cout << "Skipping: missing tree(s) in " << filename << endl;
      fin->Close(); delete fin; ifile++;
      continue;
    }
    if (USE_MC_WEIGHT)
      tree->AddFriend(tree_gen);
    tree->AddFriend(t_evtinfo);

    Int_t centrality = 0;
    Int_t N_genD0s = 0;
    vector<float> *genD0_pt = nullptr;
    Float_t genWeight = 1.f;
    Float_t ephfmQ[3], ephfpQ[3], ephfmSumW[3], ephfpSumW[3];
    Float_t ephfpAngle[3], ephfmAngle[3], eptkAngle[2], eptkQ[2];

    tree->SetBranchStatus("*", 0);
    for (const auto &p : {"centrality",
                          "ephfpAngle", "ephfmAngle", "ephfpQ", "ephfmQ",
                          "ephfpSumW", "ephfmSumW", "eptkAngle", "eptkQ"})
      tree->SetBranchStatus(p, 1);
    if (USE_MC_WEIGHT)
      for (const auto &p : {"N_genD0s", "genD0_pt", "genWeight"})
        tree->SetBranchStatus(p, 1);

    tree->SetBranchAddress("centrality", &centrality);
    tree->SetBranchAddress("ephfpAngle", ephfpAngle);
    tree->SetBranchAddress("ephfmAngle", ephfmAngle);
    tree->SetBranchAddress("ephfmQ", ephfmQ);
    tree->SetBranchAddress("ephfpQ", ephfpQ);
    tree->SetBranchAddress("ephfmSumW", ephfmSumW);
    tree->SetBranchAddress("ephfpSumW", ephfpSumW);
    tree->SetBranchAddress("eptkAngle", eptkAngle);
    tree->SetBranchAddress("eptkQ", eptkQ);
    if (USE_MC_WEIGHT)
    {
      tree->SetBranchAddress("N_genD0s", &N_genD0s);
      tree->SetBranchAddress("genD0_pt", &genD0_pt);
      tree->SetBranchAddress("genWeight", &genWeight);
    }

    Long64_t n_entries = tree->GetEntries();
    for (Long64_t ii = 0; ii < n_entries; ii++)
    {
      tree->GetEntry(ii);
      if (ii % 20000 == 0)
        printf("  entry %lld / %lld (%.1f%%)\n", ii, n_entries, 100. * ii / n_entries);

      if (centrality < 0 || centrality >= 200) // sanity only -- centrality inclusive
        continue;

      double w = 1.0;
      if (USE_MC_WEIGHT)
      {
        double p_gen_max = 0.0;
        for (int ig = 0; ig < N_genD0s; ++ig)
          p_gen_max = max(p_gen_max, (double)genD0_pt->at(ig));
        w = genWeight * getStitchWeight(p_gen_max);
        if (w <= 0.0)
          continue;
      }

      float sumW2 = ephfpSumW[1] + ephfmSumW[1];
      float sumW3 = ephfpSumW[2] + ephfmSumW[2];
      if (sumW2 <= 0 || sumW3 <= 0)
        continue;

      TComplex Q2_HFm(ephfmQ[1] * TMath::Cos(2. * ephfmAngle[1]), ephfmQ[1] * TMath::Sin(2. * ephfmAngle[1]));
      TComplex Q2_HFp(ephfpQ[1] * TMath::Cos(2. * ephfpAngle[1]), ephfpQ[1] * TMath::Sin(2. * ephfpAngle[1]));
      TComplex Q2_Trk(eptkQ[0] * TMath::Cos(2. * eptkAngle[0]), eptkQ[0] * TMath::Sin(2. * eptkAngle[0]));

      TComplex Q3_HFm(ephfmQ[2] * TMath::Cos(3. * ephfmAngle[2]), ephfmQ[2] * TMath::Sin(3. * ephfmAngle[2]));
      TComplex Q3_HFp(ephfpQ[2] * TMath::Cos(3. * ephfpAngle[2]), ephfpQ[2] * TMath::Sin(3. * ephfpAngle[2]));
      TComplex Q3_Trk(eptkQ[1] * TMath::Cos(3. * eptkAngle[1]), eptkQ[1] * TMath::Sin(3. * eptkAngle[1]));

      h_nevt->Fill(0.5, w);
      hQ2_HFmHFp->Fill((Q2_HFm * TComplex::Conjugate(Q2_HFp)).Re(), w);
      hQ2_HFmTrk->Fill((Q2_HFm * TComplex::Conjugate(Q2_Trk)).Re(), w);
      hQ2_HFpTrk->Fill((Q2_HFp * TComplex::Conjugate(Q2_Trk)).Re(), w);
      hQ3_HFmHFp->Fill((Q3_HFm * TComplex::Conjugate(Q3_HFp)).Re(), w);
      hQ3_HFmTrk->Fill((Q3_HFm * TComplex::Conjugate(Q3_Trk)).Re(), w);
      hQ3_HFpTrk->Fill((Q3_HFp * TComplex::Conjugate(Q3_Trk)).Re(), w);
    }

    fin->Close();
    delete fin;
    ifile++;
  }

  TFile *fout = new TFile(outfile, "RECREATE");
  fout->cd();
  hQ2_HFmHFp->Write();
  hQ2_HFmTrk->Write();
  hQ2_HFpTrk->Write();
  hQ3_HFmHFp->Write();
  hQ3_HFmTrk->Write();
  hQ3_HFpTrk->Write();
  h_nevt->Write();
  fout->Write(0, TObject::kOverwrite);
  fout->Close();
  cout << ">>> Resolution ingredients saved to: " << outfile << endl;
}

int main(int argc, char *argv[])
{
  if (argc != 5)
  {
    cout << "Usage: " << argv[0] << " input_txt output_path istart iend" << endl;
    return 1;
  }
  Calculate_Resolution_MC(argv[1], argv[2], std::stoi(argv[3]), std::stoi(argv[4]));
  return 0;
}
