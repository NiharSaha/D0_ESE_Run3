#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH3.h"
#include "TMath.h"

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/Prompt_frac/BDTHandler.h"
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/Prompt_frac/BDTHandler.cc"

using namespace std;
using namespace std::chrono;

void getDCA_fromData(TString input_txt, TString output_path, int istart, int iend)
{
  TH1::SetDefaultSumw2();
  auto start = high_resolution_clock::now();

  BDTHandler bdtHandler;

  const double MAX_PT_ANA = 100.0;
  const double MIN_PT_ANA = 2.0;
  const double MAX_Y_ANA = 1.0;

  const int N_CENTBIN = 5;
  Int_t min_centbin[N_CENTBIN] = {0, 10, 20, 30, 40};
  Int_t max_centbin[N_CENTBIN] = {10, 20, 30, 40, 50};
  Double_t centbinning[N_CENTBIN + 1] = {0, 10, 20, 30, 40, 50};
  const char *label_centbin[N_CENTBIN] = {"cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"};

  const int N_PTBIN = 9;
  Double_t ptbinning[N_PTBIN + 1] = {2, 3, 4, 5, 6, 8, 10, 15, 30, 100};
  const char *label_pTbin[N_PTBIN] = {"pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8", "pT8to10", "pT10to15", "pT15to30", "pT30to100"};

  // DCA binning: coarser tail for low pT (< 5 GeV), finer for higher pT.
  const Int_t dca_bins_hi = 18;
  Double_t DCA_hi[dca_bins_hi + 1] = {0,0.0005,0.0011,0.0014,0.0024,0.0029,0.0039,0.0045,0.0059,0.0067,0.0085,0.01179,0.016,0.0214,0.028,0.0366,0.0475,0.079,0.135};
  const char *label_dca_hi[dca_bins_hi] = {"dca_0to0p0005", "dca_0p0005to0p0011", "dca_0p0011to0p0014", "dca_0p0014to0p0024", "dca_0p0024to0p0029", "dca_0p0029to0p0039",
                                           "dca_0p0039to0p0045", "dca_0p0045to0p0059", "dca_0p0059to0p0067", "dca_0p0067to0p0085", "dca_0p0085to0p01179", "dca_0p01179to0p016",
                                           "dca_0p016to0p0214", "dca_0p0214to0p028", "dca_0p028to0p0366", "dca_0p0366to0p0475", "dca_0p0475to0p079", "dca_0p079to0p135"};

  const Int_t dca_bins_lo = 15;
  Double_t DCA_lo[dca_bins_lo + 1] = {0,0.0005,0.0011,0.0014,0.0024,0.0029,0.0039,0.0045,0.0059,0.0067,0.0085,0.01179,0.016,0.0214,0.0366,0.135};
  const char *label_dca_lo[dca_bins_lo] = {"dca_0to0p0005", "dca_0p0005to0p0011", "dca_0p0011to0p0014", "dca_0p0014to0p0024", "dca_0p0024to0p0029", "dca_0p0029to0p0039",
                                           "dca_0p0039to0p0045", "dca_0p0045to0p0059", "dca_0p0059to0p0067", "dca_0p0067to0p0085", "dca_0p0085to0p01179", "dca_0p01179to0p016",
                                           "dca_0p016to0p0214", "dca_0p0214to0p0366", "dca_0p0366to0p135"};

  const Int_t dca_bins_max = dca_bins_hi; // static array dimension

  // pT bins with lower edge < 5 GeV (pT2to3, pT3to4, pT4to5) use the low-pT binning
  auto useLowPtBins = [&](int ipt) { return ptbinning[ipt] < 5.0; };
  auto dcaNbins     = [&](int ipt) { return useLowPtBins(ipt) ? dca_bins_lo : dca_bins_hi; };
  auto dcaEdges     = [&](int ipt) -> Double_t * { return useLowPtBins(ipt) ? DCA_lo : DCA_hi; };
  auto dcaLabel     = [&](int ipt, int n) { return useLowPtBins(ipt) ? label_dca_lo[n] : label_dca_hi[n]; };

  // both binnings share the same first/last edge -> common acceptance window
  const Double_t DCA_EDGE_LO = DCA_hi[0];
  const Double_t DCA_EDGE_HI = DCA_hi[dca_bins_hi];

  TH1F *hcent[N_CENTBIN];
  TH1F *h_DCA[N_CENTBIN][N_PTBIN];
  TH1F *h_dmass[N_CENTBIN][N_PTBIN][dca_bins_max];
  TH1F *h_dmass_inclusive[N_CENTBIN][N_PTBIN];

  for (int icent = 0; icent < N_CENTBIN; icent++)
  {
    hcent[icent] = new TH1F(Form("hcent_%s", label_centbin[icent]), Form("hcent_%s", label_centbin[icent]), 100, 0.0, 100.0);
    for (int ipt = 0; ipt < N_PTBIN; ipt++)
    {
      h_DCA[icent][ipt] = new TH1F(Form("DCA_%s_%s", label_centbin[icent], label_pTbin[ipt]), Form("DCA_%s_%s", label_centbin[icent], label_pTbin[ipt]), dcaNbins(ipt), dcaEdges(ipt));

      h_dmass_inclusive[icent][ipt] = new TH1F(Form("dmass_inclusive_%s_%s", label_centbin[icent], label_pTbin[ipt]), Form("dmass_inclusive_%s_%s (All DCA)", label_centbin[icent], label_pTbin[ipt]), 48, 1.75, 1.99);

      for (int n = 0; n < dcaNbins(ipt); n++)
      {
        h_dmass[icent][ipt][n] = new TH1F(
            Form("dmass_%s_%s_%s", label_centbin[icent], label_pTbin[ipt], dcaLabel(ipt, n)),
            Form("dmass_%s_%s_%s", label_centbin[icent], label_pTbin[ipt], dcaLabel(ipt, n)),
            48, 1.75, 1.99);
      }
    }
  }

  // Used purely as a fast 3D bin-index lookup (cent, pT, DCA), same trick as the old file.
  TH3D *h_binning = new TH3D("centbins", "", N_CENTBIN, centbinning, N_PTBIN, ptbinning, dca_bins_hi, DCA_hi);

  ifstream file_stream(input_txt.Data());
  if (!file_stream.is_open())
  {
    cout << "Error: Could not open input list " << input_txt << endl;
    return;
  }

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
      cout << "Warning: Skipping bad file: " << filename << endl;
      if (fin)
      {
        fin->Close();
        delete fin;
      }
      ifile++;
      continue;
    }

    cout << ">>> Processing ifile=" << ifile << " : " << filename << endl;

    TTree *tree = (TTree *)fin->Get("d0Analyzer/VCNtuple_D02kpi");
    TTree *t_eventinfoana = (TTree *)fin->Get("eventinfoana/EventInfoNtuple");
    if (!tree || !t_eventinfoana)
    {
      cout << "Skipping: Tree(s) not found in " << filename << endl;
      fin->Close();
      delete fin;
      ifile++;
      continue;
    }
    tree->AddFriend(t_eventinfoana);

    Int_t candSize;
    Int_t centrality;
    vector<float> *pT = 0, *y = 0, *mass = 0, *mva = 0, *dca = 0;

    tree->SetBranchAddress("candSize", &candSize);
    tree->SetBranchAddress("centrality", &centrality);
    tree->SetBranchAddress("pT", &pT);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("mva", &mva);
    tree->SetBranchAddress("ip3d", &dca);

    tree->SetBranchStatus("*", 0);
    for (const auto &p : {"candSize", "centrality", "pT", "y", "mass", "mva", "ip3d"})
      tree->SetBranchStatus(p, 1);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t ievt = 0; ievt < nEntries; ievt++)
    {
      tree->GetEntry(ievt);

      Int_t cent = centrality / 2;
      if (cent < min_centbin[0] || cent >= max_centbin[N_CENTBIN - 1])
        continue;

      Int_t i_centbin = (h_binning->GetXaxis()->FindBin(cent)) - 1;
      if (i_centbin < 0 || i_centbin >= N_CENTBIN)
        continue;

      hcent[i_centbin]->Fill(cent);

      for (int icand = 0; icand < candSize; icand++)
      {
        float Pt = pT->at(icand);
        float Y = y->at(icand);
        float Mass = mass->at(icand);
        float Mva = mva->at(icand);
        float Dca = dca->at(icand);

        // Global acceptance, same as flow_Analysis_latest()
        if (Pt >= MAX_PT_ANA || Pt < MIN_PT_ANA || fabs(Y) >= MAX_Y_ANA)
          continue;

        // Local pT / DCA range for this histogram set
        if (Pt < ptbinning[0] || Pt >= ptbinning[N_PTBIN])
          continue;
        if (Dca < DCA_EDGE_LO || Dca >= DCA_EDGE_HI)
          continue;

        // BDT selection -- centrality convention needs *2 to match flow_Analysis_latest()
        double bdt_cut = bdtHandler.getBDTCut(Y, 2 * cent, Pt);
        if (Mva <= bdt_cut)
          continue;

        Int_t i_pTbin = (h_binning->GetYaxis()->FindBin(Pt)) - 1;
        if (i_pTbin < 0 || i_pTbin >= N_PTBIN) continue;

        // DCA slice index from this pT bin's own (possibly low-pT) binning
        Int_t i_dcabin = h_DCA[i_centbin][i_pTbin]->GetXaxis()->FindBin(Dca) - 1;
        if (i_dcabin < 0 || i_dcabin >= dcaNbins(i_pTbin)) continue;

        h_DCA[i_centbin][i_pTbin]->Fill(Dca);
        h_dmass_inclusive[i_centbin][i_pTbin]->Fill(Mass);
        h_dmass[i_centbin][i_pTbin][i_dcabin]->Fill(Mass);
      }
    }

    fin->Close();
    delete fin;
    ifile++;
  }
  file_stream.close();

  gSystem->mkdir(Form("%s/ROOT", output_path.Data()), kTRUE);
  TString outfilename = Form("Data_Mass_DCA_%d_%d", istart, iend);
  TString outfile = TString::Format("%s/ROOT/%s.root", output_path.Data(), outfilename.Data());
  TFile *fout = new TFile(outfile, "RECREATE");

  fout->cd();
  for (int icent = 0; icent < N_CENTBIN; icent++)
    hcent[icent]->Write();

  for (int icent = 0; icent < N_CENTBIN; icent++)
  {
    for (int ipt = 0; ipt < N_PTBIN; ipt++)
    {
      fout->mkdir(Form("out_Mass_DCA_%s_%s", label_centbin[icent], label_pTbin[ipt]));
      fout->cd(Form("out_Mass_DCA_%s_%s", label_centbin[icent], label_pTbin[ipt]));
      h_DCA[icent][ipt]->Write();
      h_dmass_inclusive[icent][ipt]->Write();
      for (int idca = 0; idca < dcaNbins(ipt); idca++)
        h_dmass[icent][ipt][idca]->Write();
    }
  }

  fout->Close();

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<minutes>(stop - start);
  cout << "Total time taken: " << duration.count() << " minutes" << endl;
  cout << "Successfully saved to " << outfile << endl;
}

int main(int argc, char **argv)
{
  // Now matches the istart/iend convention used by getDCA_fromMC and
  // flow_Analysis_latest, instead of the old single-index itxtoutFile scheme.
  if (argc == 5)
  {
    TString input_txt = argv[1];
    TString output_path = argv[2];
    int istart = std::stoi(argv[3]);
    int iend = std::stoi(argv[4]);

    getDCA_fromData(input_txt, output_path, istart, iend);
  }
  else
  {
    std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
    return 1;
  }
  return 0;
}
