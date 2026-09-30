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
#include <cstdio>
#include <unistd.h>

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

Int_t cent_case[N_CENTBIN_eff] = {5, 5, 15, 15, 25, 25, 35, 35, 45, 45};

const int N_PTBIN = 9;
Double_t min_pTbin[N_PTBIN] = {2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 30.0};
Double_t max_pTbin[N_PTBIN] = {3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 30.0, 100.0};
Double_t ptbinning[N_PTBIN + 1] = {2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 30.0, 100.0};

const int N_PTBIN_eff = 23;
Double_t ptbinning_[N_PTBIN_eff + 1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 70.0, 100.0};

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

TH1D *h_pthat = new TH1D("pthat", ";#hat{p}_{T} (GeV)", 300, 0., 300.);

TH1D *h_ptreco_cent[N_CENTBIN] = {h_ptreco_cent010, h_ptreco_cent1020, h_ptreco_cent2030, h_ptreco_cent3040, h_ptreco_cent4050};
TH1D *h_ptreco_cent_woBDT[N_CENTBIN] = {h_ptreco_cent010_woBDT, h_ptreco_cent1020_woBDT, h_ptreco_cent2030_woBDT, h_ptreco_cent3040_woBDT, h_ptreco_cent4050_woBDT};
TH1D *h_ptgen_cent[N_CENTBIN] = {h_ptgen_cent010, h_ptgen_cent1020, h_ptgen_cent2030, h_ptgen_cent3040, h_ptgen_cent4050};

// pT-weighted twins (totalWeight * getPtWeight(pt)) of every RECO histogram above, for the
// AccxEff x pT-weight nominal efficiency -- the unweighted (totalWeight-only) histograms above
// are kept as-is and become the "no pT-weight" systematic variation. Reco-side only, NOT gen:
// Efficiency = Reco/Gen is a ratio within each bin, so applying the same per-bin w(pT) to BOTH
// sides is a common factor that cancels almost exactly (up to bin-migration effects), making the
// weighted and unweighted efficiency curves come out nearly identical -- defeating the whole
// point. Applying it only to Reco gives Eff_w(bin) = w(bin) * Reco(bin)/Gen(bin) =
// w(bin) * Eff_unweighted(bin), i.e. literally AccxEff x pT-weight, which is what was actually
// intended (confirmed with the user, 2026-09-24). Both the nominal and systematic efficiency
// divide by the SAME unweighted Gen histograms above -- no separate gen _ptw twins needed.
// h_pthat has no per-fill pT (it's the event's own pthat, not a candidate pT) so no twin either.
TH2D *h_centptreco_ptw = new TH2D("centptreco_ptw", ";cent;pT", N_CENTBIN_eff, centbinning_, N_PTBIN_eff, ptbinning_);
TH2D *h_centptreco_woBDT_ptw = new TH2D("centptreco_woBDT_ptw", ";cent;pT", N_CENTBIN_eff, centbinning_, N_PTBIN_eff, ptbinning_);

TH1D *h_ptreco_cent010_ptw = new TH1D("ptreco_cent010_ptw", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent010_woBDT_ptw = new TH1D("ptreco_cent010_woBDT_ptw", ";pT", N_PTBIN_eff, ptbinning_);

TH1D *h_ptreco_cent1020_ptw = new TH1D("ptreco_cent1020_ptw", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent1020_woBDT_ptw = new TH1D("ptreco_cent1020_woBDT_ptw", ";pT", N_PTBIN_eff, ptbinning_);

TH1D *h_ptreco_cent2030_ptw = new TH1D("ptreco_cent2030_ptw", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent2030_woBDT_ptw = new TH1D("ptreco_cent2030_woBDT_ptw", ";pT", N_PTBIN_eff, ptbinning_);

TH1D *h_ptreco_cent3040_ptw = new TH1D("ptreco_cent3040_ptw", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent3040_woBDT_ptw = new TH1D("ptreco_cent3040_woBDT_ptw", ";pT", N_PTBIN_eff, ptbinning_);

TH1D *h_ptreco_cent4050_ptw = new TH1D("ptreco_cent4050_ptw", ";pT", N_PTBIN_eff, ptbinning_);
TH1D *h_ptreco_cent4050_woBDT_ptw = new TH1D("ptreco_cent4050_woBDT_ptw", ";pT", N_PTBIN_eff, ptbinning_);

TH1D *h_ptreco_cent_ptw[N_CENTBIN] = {h_ptreco_cent010_ptw, h_ptreco_cent1020_ptw, h_ptreco_cent2030_ptw, h_ptreco_cent3040_ptw, h_ptreco_cent4050_ptw};
TH1D *h_ptreco_cent_woBDT_ptw[N_CENTBIN] = {h_ptreco_cent010_woBDT_ptw, h_ptreco_cent1020_woBDT_ptw, h_ptreco_cent2030_woBDT_ptw, h_ptreco_cent3040_woBDT_ptw, h_ptreco_cent4050_woBDT_ptw};

// Cross-section-based stitching weights for the nested "Dpt > X GeV" generator-level
// D0-pT-filtered samples that make up inputFiles_Prompt_MC_Sept5_noDpt8.txt. Each sample requires
// at least one generated D0 above its own threshold, so the samples are *nested*
// (Dpt>80 subset of Dpt>60 subset of ... subset of Dpt>0), not disjoint pT-hat bins.
// Effective luminosity L_i = N_i / (xsec_i * filterEff_i); an event whose leading generated
// D0 has true pT = p is combined with weight = 1 / sum_{i : thr_i <= p} L_i.
//
// Dpt_xsec_mb[]/Dpt_filterEff[] updated 2026-09-24 with directly-measured GenXsecAnalyzer
// values (thr=10,20,40,60,80 measured 2026-09-23; thr=30 re-measured 2026-09-24 with a larger,
// 200k-event test batch -- the original 100k-event/1-survivor run was statistically unreliable,
// 100% relative uncertainty). See get_PtWeight.C's Dpt_xsecAfterFilter_pb[] for the equivalent
// single-number (xsec_mb*filterEff, in pb) form and full per-sample provenance notes. thr=0,1
// remain the original, unverified McM values -- deferred, same as in get_PtWeight.C.
//
// thr=8 dropped entirely (2026-09-24, one-off ablation): its GEN fragment was confirmed to be
// byte-identical to thr=1's (MomMinPt=1.0, not 8.0 -- a production mislabeling, full story in
// get_PtWeight.C), and its assumed xsec/filterEff (3.08e1 mb vs thr=1's 3.67e3 mb, a ~119x drop
// that doesn't match the two samples actually sharing the same real filter) was suspect. Removed
// from the arrays below (was index 2) rather than zeroed at runtime, to match its files also
// being removed from the input list -- keeps N_DPT_SAMPLES, the arrays, and the actual file list
// all consistent with each other. Previous (McM-only, as provided 2026-09-06, still including
// thr=8) kept commented out for comparison:
// const int N_DPT_SAMPLES = 9;
// Double_t Dpt_threshold[N_DPT_SAMPLES] = {0, 1, 8, 10, 20, 30, 40, 60, 80};
// Double_t Dpt_xsec_mb[N_DPT_SAMPLES]   = {4.63e3, 3.67e3, 3.08e1, 1.15e1,   7.1e-1,  1.15e-1, 3.5e-2,  4.57e-3, 1.11e-3};
// Double_t Dpt_filterEff[N_DPT_SAMPLES] = {7.43e-3, 5.9e-3, 5.8e-3, 3.9e-4,  3.9e-4,  6.3e-5,  2.3e-4,  1.42e-4, 1.11e-4};
// Double_t Dpt_Nevts[N_DPT_SAMPLES] = {1.0e6, 2.0e6, 1.5e5, 1.0e5, 5.0e4, 2.0e4, 2.0e4, 2.0e4, 2.0e4};
const int N_DPT_SAMPLES = 8;
Double_t Dpt_threshold[N_DPT_SAMPLES] = {0, 1, 10, 20, 30, 40, 60, 80};
Double_t Dpt_xsec_mb[N_DPT_SAMPLES]   = {4.63e3, 3.67e3, 29.36, 1.809,   1.810,   0.1539,  0.03214, 0.01002};
Double_t Dpt_filterEff[N_DPT_SAMPLES] = {7.43e-3, 5.9e-3, 3.400e-4, 3.800e-4, 5.500e-5, 1.200e-4, 1.700e-4, 7.000e-5};
Double_t Dpt_Nevts[N_DPT_SAMPLES] = {1.0e6, 2.0e6, 1.0e5, 5.0e4, 2.0e4, 2.0e4, 2.0e4, 2.0e4};
Double_t Dpt_Leff[N_DPT_SAMPLES];

double getStitchWeight(double p_gen)
{
  double sumL = 0;
  for (int i = 0; i < N_DPT_SAMPLES; i++)
    if (Dpt_threshold[i] <= p_gen)
      sumL += Dpt_Leff[i];
  return (sumL > 0) ? 1.0 / sumL : 0.0;
}

// pT-weighting factor (Data/MC ratio), derived in get_PtWeight_thr1only.C (2026-09-24): thr=1
// MC (pooled thr=1+thr=8, after get_MassSpectra_PtWeight.C's findDptSampleIndex() routing fix)
// vs. Data mass-fit yields, pol4 fit over [2,30) GeV, frozen at the [25,30) bin's own fit value
// beyond that (30-100 GeV all share that one constant, 1.284). Same fine pT binning as
// ptbinning_ above -- one value per bin, looked up by candidate/gen pT. Used to build the
// AccxEff x pT-weight nominal efficiency; AccxEff alone (unweighted, the existing
// totalWeight-only histograms below) is kept as the systematic variation.
//
// Previous (thr=0-only, from get_PtWeight_thr0only.C, superseded now that thr=8's fragment bug
// is fixed and thr=1 has substantially more pooled statistics) kept commented out for reference:
// Double_t ptWeightValue[N_PTBIN_eff] = {
//     1.2845, 1.1747, 1.0806, 1.0014, 0.93583, 0.88311, 0.84225, 0.81234, 0.7925, 0.78186,
//     0.77958, 0.78486, 0.80524, 0.85178, 0.97092, 1.1787, 1.4491, 1.5181, 1.1479, 1.1479,
//     1.1479, 1.1479, 1.1479};
//
// Previous (thr=1-derived, full fitted/frozen curve applied at all pT) kept commented out --
// superseded 2026-09-24 by capping the weight's reach to pT<10 GeV only (weight=1.0 above that),
// per the assumption that pT-weighting mainly matters at low pT:
// Double_t ptWeightValue[N_PTBIN_eff] = {
//     1.2646, 1.159, 1.0694, 0.99456, 0.93352, 0.88524, 0.84873, 0.82302, 0.80719, 0.80034,
//     0.80159, 0.81012, 0.83479, 0.88586, 1.0095, 1.2165, 1.476, 1.5481, 1.284, 1.284,
//     1.284, 1.284, 1.284};
Double_t ptWeightValue[N_PTBIN_eff] = {
    1.2646, 1.159, 1.0694, 0.99456, 0.93352, 0.88524, 0.84873, 0.82302, 0.80719, 0.80034,
    0.80159, 0.81012, 0.83479, 0.88586, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0};

double getPtWeight(double pt)
{
  for (int i = 0; i < N_PTBIN_eff; i++)
    if (pt >= ptbinning_[i] && pt < ptbinning_[i + 1])
      return ptWeightValue[i];
  return 1.0; // outside [2,100) GeV -- shouldn't occur given RECO_cut/GEN_cut already restrict pt
}

bool isNonPromptMother(int pdg)
    {
        int absPdg = std::abs(pdg);
        return ((absPdg / 100) % 10 == 5) || ((absPdg / 1000) % 10 == 5);
    }


void sanitizeBdtCutsFile(const TString &path)
{
  ifstream fin(path.Data());
  if (!fin.is_open())
    return;

  std::ostringstream cleaned;
  string line;
  while (getline(fin, line))
  {
    if (line.empty())
      continue;
    if (line[0] == '#')
    {
      cleaned << line << "\n";
      continue;
    }
    istringstream ss(line);
    int a, b, c;
    double d;
    if (!(ss >> a >> b >> c >> d))
      continue;
    if (a < 0 || a >= y_bins || b < 0 || b >= cent_bins || c < 0 || c >= pT_bins)
      continue;
    cleaned << line << "\n";
  }
  fin.close();

  TString tmpPath = TString::Format("%s.sanitize_tmp_%d", path.Data(), getpid());
  ofstream fout(tmpPath.Data());
  fout << cleaned.str();
  fout.close();

  std::rename(tmpPath.Data(), path.Data());
}

void getEfficiency(TString input_txt, TString output_path, int istart, int iend)
{
  sanitizeBdtCutsFile("/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/AccxEff/bdt_cuts.csv");

  BDTHandler bdtHandler;

  for (int i = 0; i < N_DPT_SAMPLES; i++)
    Dpt_Leff[i] = Dpt_Nevts[i] / (Dpt_xsec_mb[i] * Dpt_filterEff[i]);


  TH1::StatOverflows(kTRUE);
  TH2::StatOverflows(kTRUE);
  TH3::StatOverflows(kTRUE);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  ifstream file_stream(input_txt.Data());
  TString outfile = TString::Format("%s/ROOT/Eff_RecoGen_weight_%d_%d.root", output_path.Data(), istart, iend);
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

    TTree *tree = (TTree *)fin->Get("d0Analyzer/VCNtuple_D02kpi");
    TTree *tree_gen = (TTree *)fin->Get("d0Analyzer/AllGensNtuple");
    TTree *t_eventinfoana = (TTree *)fin->Get("eventinfoana/EventInfoNtuple");

    if (!tree || !tree_gen || !t_eventinfoana)
    {
      std::cout << "Skipping: Tree(s) missing in " << filename << std::endl;
      fin->Close();
      delete fin;
      ifile++;
      continue;
    }

    tree->AddFriend(tree_gen);
    tree->AddFriend(t_eventinfoana);

    Int_t candSize = 0;
    Int_t centrality = 0;
    Int_t N_genD0s = 0;
    Float_t pthat = 0;
    Float_t genWeight = 1;

    std::vector<Bool_t> *matchGEN = nullptr;
    std::vector<Bool_t> *isSwap = nullptr;
    std::vector<Float_t> *pT = nullptr;
    std::vector<Float_t> *mass = nullptr;
    std::vector<Float_t> *y = nullptr;
    std::vector<Float_t> *mva = nullptr;
    std::vector<bool> *Dgen_isPrompt = nullptr; // per-candidate, gen-matched: true if this gen D0 descends from a b hadron

    std::vector<Float_t> *genD0_pt = nullptr;
    std::vector<Float_t> *genD0_y = nullptr;
    std::vector<Int_t> *genD0_MotherPdgID = nullptr;

    tree->SetBranchStatus("*", 0);
    for (const auto &p : {"candSize", "centrality", "y", "pT", "mass", "mva", "ip3d",
                           "matchGEN", "isSwap", "Dgen_isPrompt",
                           "N_genD0s", "genD0_pt", "genD0_y", "genD0_MotherPdgID",
                           "pthat", "genWeight"})
    {
      tree->SetBranchStatus(p, 1);
    }

    tree->SetBranchAddress("centrality", &centrality);
    tree->SetBranchAddress("candSize", &candSize);
    tree->SetBranchAddress("N_genD0s", &N_genD0s);
    tree->SetBranchAddress("pT", &pT);
    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("mva", &mva);
    tree->SetBranchAddress("matchGEN", &matchGEN);
    tree->SetBranchAddress("isSwap", &isSwap);
    tree->SetBranchAddress("Dgen_isPrompt", &Dgen_isPrompt);
    tree->SetBranchAddress("genD0_pt", &genD0_pt);
    tree->SetBranchAddress("genD0_y", &genD0_y);
    tree->SetBranchAddress("genD0_MotherPdgID", &genD0_MotherPdgID);
    tree->SetBranchAddress("pthat", &pthat);
    tree->SetBranchAddress("genWeight", &genWeight);

    Int_t nevent = tree->GetEntries();
    std::cout << "nevent=" << nevent << std::endl;

    for (int ievt = 0; ievt < nevent; ievt++)
    {
      tree->GetEntry(ievt);

      if (ievt % 1000 == 0)
        printf("current entry of event loop = %d out of %d : %.3f %%\n", ievt, nevent, (Double_t)ievt / nevent * 100);

      if (centrality < 0 || centrality >= 200)
        continue;

      // Leading generated D0 pT in the event decides which of the nested Dpt-threshold
      // samples this event was eligible to appear in; the same event-level weight is
      // reused below for both the reco (numerator) and gen (denominator) fills.
      double p_gen_max = 0;
      for (int igen = 0; igen < N_genD0s; igen++)
        p_gen_max = std::max(p_gen_max, (double)genD0_pt->at(igen));

      double totalWeight = genWeight * getStitchWeight(p_gen_max);

      h_pthat->Fill(pthat, totalWeight);

      for (int icand = 0; icand < candSize; icand++)
      {

        if (matchGEN->at(icand) != 1 || isSwap->at(icand) != 0)
          continue;

        float Pt = pT->at(icand);
        float Y = y->at(icand);
        float Mva = mva->at(icand);
        

        // Direct per-candidate ancestry flag from the gen-matched branch (full ancestry,
        // not just the immediate mother's PDG code)
        if (!Dgen_isPrompt->at(icand))
          continue;

        bool RECO_cut = (Pt >= MIN_PT_ANA && Pt < MAX_PT_ANA && fabs(Y) < MAX_Y_ANA);
        bool RECO_cut_woBDT = RECO_cut;

        // AccxEff x pT-weight nominal: same totalWeight, additionally scaled by the pT-weight
        // factor for this candidate's own Pt. AccxEff alone (totalWeight, filled above/below
        // unchanged) remains the "no pT-weight" systematic variation.
        double totalWeight_ptw = totalWeight * getPtWeight(Pt);

        if (RECO_cut)
        {
          for (int ibin = 0; ibin < N_CENTBIN_eff; ibin++)
          {
            if (Mva <= bdtHandler.getBDTCut(Y, 2 * cent_case[ibin], Pt))
              continue;
            double cent_center = 0.5 * (centbinning_[ibin] + centbinning_[ibin + 1]);
            h_centptreco->Fill(cent_center, Pt, totalWeight);
            h_centptreco_ptw->Fill(cent_center, Pt, totalWeight_ptw);
          }

          for (int iwide = 0; iwide < N_CENTBIN; iwide++)
          {
            if (Mva <= bdtHandler.getBDTCut(Y, 2 * cent_case[2 * iwide], Pt))
              continue;
            h_ptreco_cent[iwide]->Fill(Pt, totalWeight);
            h_ptreco_cent_ptw[iwide]->Fill(Pt, totalWeight_ptw);
          }
        }

        if (RECO_cut_woBDT)
        {
          for (int ibin = 0; ibin < N_CENTBIN_eff; ibin++)
          {
            double cent_center = 0.5 * (centbinning_[ibin] + centbinning_[ibin + 1]);
            h_centptreco_woBDT->Fill(cent_center, Pt, totalWeight);
            h_centptreco_woBDT_ptw->Fill(cent_center, Pt, totalWeight_ptw);
          }

          for (int iwide = 0; iwide < N_CENTBIN; iwide++)
          {
            h_ptreco_cent_woBDT[iwide]->Fill(Pt, totalWeight);
            h_ptreco_cent_woBDT_ptw[iwide]->Fill(Pt, totalWeight_ptw);
          }
        }

      } // End of candidate loop

      // Gen loop (denominator): same prompt/non-prompt definition as the numerator,
      // using the immediate-mother PDG code since AllGensNtuple has no full-ancestry flag
      for (int igen = 0; igen < N_genD0s; igen++)
      {
        int motherId = std::abs(genD0_MotherPdgID->at(igen));
        if (isNonPromptMother(motherId))
          continue;

        bool GEN_cut = (fabs(genD0_y->at(igen)) < MAX_Y_ANA && genD0_pt->at(igen) >= MIN_PT_ANA && genD0_pt->at(igen) < MAX_PT_ANA);

        if (GEN_cut)
        {
          // No pT-weighting here -- see the declaration comment above the _ptw reco histograms
          // for why: it would cancel against the reco-side weighting in the efficiency ratio.
          // Both the nominal (w-corrected) and systematic efficiency divide by these same,
          // single set of gen histograms.
          for (int ibin = 0; ibin < N_CENTBIN_eff; ibin++)
          {
            double cent_center = 0.5 * (centbinning_[ibin] + centbinning_[ibin + 1]);
            h_centptgen->Fill(cent_center, genD0_pt->at(igen), totalWeight);
          }

          for (int iwide = 0; iwide < N_CENTBIN; iwide++)
          {
            h_ptgen_cent[iwide]->Fill(genD0_pt->at(igen), totalWeight);
          }
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
  h_pthat->Write();

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

  // pT-weighted twins (AccxEff x pT-weight nominal), reco-side only -- see the declarations
  // above for why (h_pthat and the gen histograms have no _ptw twin).
  h_centptreco_ptw->Write();
  h_centptreco_woBDT_ptw->Write();

  h_ptreco_cent010_ptw->Write();
  h_ptreco_cent1020_ptw->Write();
  h_ptreco_cent2030_ptw->Write();
  h_ptreco_cent3040_ptw->Write();
  h_ptreco_cent4050_ptw->Write();

  h_ptreco_cent010_woBDT_ptw->Write();
  h_ptreco_cent1020_woBDT_ptw->Write();
  h_ptreco_cent2030_woBDT_ptw->Write();
  h_ptreco_cent3040_woBDT_ptw->Write();
  h_ptreco_cent4050_woBDT_ptw->Write();

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
