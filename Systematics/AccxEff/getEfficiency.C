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

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/AccxEff/BDTHandler.h"
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/AccxEff/BDTHandler.cc"

using namespace std;
using namespace std::chrono;

const int N_CENTBIN = 5;
Int_t min_centbin[N_CENTBIN] = {0, 10, 20, 30, 40};
Int_t max_centbin[N_CENTBIN] = {10, 20, 30, 40, 50};
Double_t centbinning[N_CENTBIN + 1] = {0, 10, 20, 30, 40, 50};

const int N_CENTBIN_eff = 10;
Double_t centbinning_[N_CENTBIN_eff + 1] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50.};

const int N_PTBIN = 9;
Double_t min_pTbin[N_PTBIN] = {2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 30.0};
Double_t max_pTbin[N_PTBIN] = {3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 30.0, 100.0};
Double_t ptbinning[N_PTBIN + 1] = {2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 30.0, 100.0};

const int N_PTBIN_eff = 18;
Double_t ptbinning_[N_PTBIN_eff + 1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 12.5, 15.0, 20.0, 30.0, 50.0, 100.0};

const double MAX_PT_ANA = 100.0;
const double MIN_PT_ANA = 2.0;
const double MAX_Y_ANA = 1.0;

TH1D *h_ptgen_cent010 = new TH1D("ptgen_cent010", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent010 = new TH1D("ptreco_cent010", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent010_woBDT = new TH1D("ptreco_cent010_woBDT", ";pT", N_PTBIN_eff, ptbinning_);

TH1D *h_ptgen_cent1020 = new TH1D("ptgen_cent1020", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent1020 = new TH1D("ptreco_cent1020", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent1020_woBDT = new TH1D("ptreco_cent1020_woBDT", ";pT", N_PTBIN_eff, ptbinning_);

TH1D *h_ptgen_cent2030 = new TH1D("ptgen_cent2030", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent2030 = new TH1D("ptreco_cent2030", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent2030_woBDT = new TH1D("ptreco_cent2030_woBDT", ";pT", N_PTBIN_eff, ptbinning_);

TH1D *h_ptgen_cent3040 = new TH1D("ptgen_cent3040", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent3040 = new TH1D("ptreco_cent3040", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent3040_woBDT = new TH1D("ptreco_cent3040_woBDT", ";pT", N_PTBIN_eff, ptbinning_);

TH1D *h_ptgen_cent4050 = new TH1D("ptgen_cent4050", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent4050 = new TH1D("ptreco_cent4050", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent4050_woBDT = new TH1D("ptreco_cent4050_woBDT", ";pT", N_PTBIN_eff, ptbinning_);

TH2D *h_centptreco_woBDT = new TH2D("centptreco_woBDT", ";cent;pT", N_CENTBIN_eff, centbinning_, N_PTBIN_eff, ptbinning_);
TH2D *h_centptreco = new TH2D("centptreco", ";cent;pT", N_CENTBIN_eff, centbinning_, N_PTBIN_eff, ptbinning_);
TH2D *h_centptgen = new TH2D("centptgen", ";cent;pT", N_CENTBIN_eff, centbinning_, N_PTBIN_eff, ptbinning_);

void getEfficiency(TString input_txt, TString output_path, int istart, int iend)
{
  BDTHandler bdtHandler;

  TH1::StatOverflows(kTRUE);
  TH2::StatOverflows(kTRUE);
  TH3::StatOverflows(kTRUE);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  ifstream file_stream(input_txt.Data());
  TString outfile = TString::Format("%s/ROOT/Eff_%d_%d.root", output_path.Data(), istart, iend);
  TFile *fout = new TFile(outfile, "RECREATE");

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

    TTree *tree = (TTree *)fin->Get("d0Analyzer/VertexCompositeNtuple");
    TTree *tree_gen = (TTree *)fin->Get("d0Analyzer/AllGensNtuple");
    TTree *t_eventinfoana = (TTree *)fin->Get("eventinfoana/EventInfoNtuple");
    tree->AddFriend(tree_gen);
    tree->AddFriend(t_eventinfoana);

    const Int_t Max_cand = 10000;
    Int_t candSize, centrality, genD0_centrality, N_genD0s;
    Int_t matchGEN[Max_cand], isSwap[Max_cand];
    Float_t pT[Max_cand], mass[Max_cand], y[Max_cand], eta[Max_cand], mva[Max_cand], genD0_pt[Max_cand], genD0_y[Max_cand], genD0_eta[Max_cand];

    tree->SetBranchAddress("candSize", &candSize);
    tree->SetBranchAddress("centrality", &centrality);
    tree->SetBranchAddress("pT", &pT);
    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("mva", &mva);
    tree->SetBranchAddress("matchGEN", &matchGEN);
    tree->SetBranchAddress("isSwap", &isSwap);
    tree->SetBranchAddress("genD0_centrality", &genD0_centrality);
    tree->SetBranchAddress("N_genD0s", &N_genD0s);
    tree->SetBranchAddress("genD0_pt", &genD0_pt);
    tree->SetBranchAddress("genD0_y", &genD0_y);
    tree->SetBranchAddress("genD0_eta", &genD0_eta);

    tree->SetBranchStatus("*", 0);
    for (const auto &p : {"candSize", "centrality", "y", "eta", "pT", "mass", "mva", "matchGEN", "isSwap", "genD0_centrality", "N_genD0s", "genD0_pt", "genD0_y", "genD0_eta"})
      tree->SetBranchStatus(p, 1);

    Int_t nevent = tree->GetEntries();
    std::cout << "nevent=" << nevent << std::endl;

    for (int ievt = 0; ievt < nevent; ievt++)
    {
      tree->GetEntry(ievt);

      if (ievt % 1000 == 0)
        printf("current entry of event loop = %d out of %d : %.3f %%\n", ievt, nevent, (Double_t)ievt / nevent * 100);

      Int_t cent = centrality / 2;
      if (cent >= max_centbin[N_CENTBIN - 1]) continue;

      for (int icand = 0; icand < candSize; icand++)
      {

        if (pT[icand] >= MAX_PT_ANA || pT[icand] < MIN_PT_ANA || fabs(y[icand]) >= MAX_Y_ANA) continue;

        double bdt_cut = bdtHandler.getBDTCut(y[icand], 2 * cent, pT[icand]);

        cout<<"Candidate " << icand << ": pT=" << pT[icand] << ", y=" << y[icand] << ", cent=" << cent << ", bdt_cut=" << bdt_cut << ", mva=" << mva[icand] << endl;

        bool RECO_cut = (matchGEN[icand] == 1 && isSwap[icand] == 0 && mva[icand] > bdt_cut && pT[icand] >= MIN_PT_ANA && pT[icand] < MAX_PT_ANA && fabs(y[icand]) < MAX_Y_ANA);
        bool RECO_cut_woBDT = (matchGEN[icand] == 1 && isSwap[icand] == 0 && pT[icand] >= MIN_PT_ANA && pT[icand] < MAX_PT_ANA && fabs(y[icand]) < MAX_Y_ANA);

        if (RECO_cut)
        {
          if (cent >= 0 && cent < 10)
            h_ptreco_cent010->Fill(pT[icand]);
          if (cent >= 10 && cent < 20)
            h_ptreco_cent1020->Fill(pT[icand]);
          if (cent >= 20 && cent < 30)
            h_ptreco_cent2030->Fill(pT[icand]);
          if (cent >= 30 && cent < 40)
            h_ptreco_cent3040->Fill(pT[icand]);
          if (cent >= 40 && cent < 50)
            h_ptreco_cent4050->Fill(pT[icand]);
          if (cent >= 0 && cent < 50)
            h_centptreco->Fill(cent, pT[icand]);
        }

        if (RECO_cut_woBDT)
        {
          if (cent >= 0 && cent < 10)
            h_ptreco_cent010_woBDT->Fill(pT[icand]);
          if (cent >= 10 && cent < 20)
            h_ptreco_cent1020_woBDT->Fill(pT[icand]);
          if (cent >= 20 && cent < 30)
            h_ptreco_cent2030_woBDT->Fill(pT[icand]);
          if (cent >= 30 && cent < 40)
            h_ptreco_cent3040_woBDT->Fill(pT[icand]);
          if (cent >= 40 && cent < 50)
            h_ptreco_cent4050_woBDT->Fill(pT[icand]);
          if (cent >= 0 && cent < 50)
            h_centptreco_woBDT->Fill(cent, pT[icand]);
        }

      } // End of candidate loop

      // Gen loop
      for (int igen = 0; igen < N_genD0s; igen++)
      {

        bool GEN_cut = (genD0_pt[igen] > 0 && fabs(genD0_y[igen]) < MAX_Y_ANA && genD0_pt[igen] >= MIN_PT_ANA && genD0_pt[igen] < MAX_PT_ANA);

        if (GEN_cut)
        {
          if (cent >= 0 && cent < 10)
            h_ptgen_cent010->Fill(genD0_pt[igen]);
          if (cent >= 10 && cent < 20)
            h_ptgen_cent1020->Fill(genD0_pt[igen]);
          if (cent >= 20 && cent < 30)
            h_ptgen_cent2030->Fill(genD0_pt[igen]);
          if (cent >= 30 && cent < 40)
            h_ptgen_cent3040->Fill(genD0_pt[igen]);
          if (cent >= 40 && cent < 50)
            h_ptgen_cent4050->Fill(genD0_pt[igen]);
          if (cent >= 0 && cent < 50)
            h_centptgen->Fill(cent, genD0_pt[igen]);
        }
      } // End of gen loop
    } // End of event loop

    fin->Close();
    delete fin;
    ifile++;
  } // End of file loop

  //--------------------------------------
  // Write all histograms and efficiency
  //--------------------------------------
  fout->cd();

  h_centptreco->Write();
  h_centptreco_woBDT->Write();
  h_centptgen->Write();

  h_ptgen_cent010->Write();
  h_ptgen_cent1020->Write();
  h_ptgen_cent2030->Write();
  h_ptgen_cent3040->Write();
  h_ptgen_cent4050->Write();

  h_ptreco_cent010->Write();
  h_ptreco_cent1020->Write();
  h_ptreco_cent2030->Write();
  h_ptreco_cent3040->Write();
  h_ptreco_cent4050->Write();

  h_ptreco_cent010_woBDT->Write();
  h_ptreco_cent1020_woBDT->Write();
  h_ptreco_cent2030_woBDT->Write();
  h_ptreco_cent3040_woBDT->Write();
  h_ptreco_cent4050_woBDT->Write();

  fout->Close();

} // THE END

int main(int argc, char *argv[])
{
  if (argc == 5) // Expecting 4 arguments: input_txt, output_path, istart, iend

  {
    TString input_txt = argv[1];
    TString output_path = argv[2];
    int istart = std::stoi(argv[3]);
    int iend = std::stoi(argv[4]);

    getEfficiency(input_txt, output_path, istart, iend);
  }
  else
  {
    std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
    return 1;
  }
  return 0;
}
