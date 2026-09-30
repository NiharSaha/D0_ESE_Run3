// Step 1/2 of the pT-weighting derivation for AccxEff.
//
// Compiled batch executable (same pattern as getEfficiency.C / getDCA_fromData.C):
// reads ONE chunk (istart..iend) of either the data or the (prompt) MC TTree file
// list, applies the optimised BDT cuts (same BDTHandler), and produces:
//   - mode "data": one mass histogram per fine pT bin (ptbinning_ below), restricted
//     to the same 0-50% centrality / |y|<1 / 2<pT<100 acceptance used elsewhere in
//     this directory -> DataMass_<istart>_<iend>.root. Fit in get_PtWeight_fromMass.C
//     (macro 2) after hadd'ing all chunks into one file.
//   - mode "mc": one weighted reco pT spectrum PER Dpt-threshold sample (9 histograms,
//     one per nested "Dpt > X GeV" sample, weighted by genWeight only -- no cross-sample
//     stitching at this stage), plus a parallel set of 9 unweighted (weight=1) histograms
//     with the same selection for comparing against genWeight's effect on the shape ->
//     MCSpectrum_<istart>_<iend>.root. No mass fit is needed for MC since truth-matching
//     already removes combinatorial background. Which of the 9 histograms a file's events
//     fill is decided per FILE (from the "Dpt{N}" tag in its path), not per event:
//     combining the 9 samples into one smooth spectrum (properly weighting each by its own
//     per-bin statistical power, since a blanket eligibility-sum approach lets a sparse
//     high-threshold sample distort a bin the instant it crosses its own threshold,
//     regardless of whether it has any real statistics there) is done downstream in
//     get_PtWeight_fromMass.C, after hadd'ing all chunks together.
//
// Data uses the candidate's own event centrality (branch "centrality", raw CMS 0-200
// units) to look up the right BDT working point, restricted to 0-50%. MC does not: its
// per-event centrality has no physical meaning (there is no real underlying-event
// multiplicity behind it) and MC statistics are too limited to bin by it anyway, so MC
// candidates are taken with no centrality cut at all and, exactly as in getEfficiency.C,
// the BDT cut is evaluated at 5 representative centrality points (the 5/15/25/35/45%
// midpoints of the 0-50% decades) rather than at the candidate's own centrality -- a
// candidate is filled once for each of those 5 points it passes.
//
// Usage: ./exe <data|mc> input_txt output_path istart iend
// (compile with: g++ get_MassSpectra_PtWeight.C $(root-config --cflags --libs) -Wall -O2 -o exe)

#include <cstdlib>
#include <cctype>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TH1.h"
#include "TMath.h"

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/AccxEff/BDTHandler.h"
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/AccxEff/BDTHandler.cc"

using namespace std;

// --- analysis acceptance, same as getEfficiency.C / getDCA_fromData.C ---
const double MIN_PT_ANA = 2.0;
const double MAX_PT_ANA = 100.0;
const double MAX_Y_ANA = 1.0;
const Int_t MIN_CENT = 0;  // percent
const Int_t MAX_CENT = 50; // percent

// MC representative-centrality points for the BDT cut (see note above), same convention
// as getEfficiency.C's cent_case wide-bin midpoints: 5/15/25/35/45% = midpoints of the
// five 0-50% decades.
const int N_CENT_REPR = 5;
Int_t centRepr[N_CENT_REPR] = {5, 15, 25, 35, 45};

// --- fine pT binning used for the pT-weight derivation (finer than the AccxEff pT bins) ---
const Int_t N_PTBIN_eff = 23;
Double_t ptbinning_[N_PTBIN_eff + 1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 70.0, 100.0};
const char *ptLabel_[N_PTBIN_eff] = {
    "pT2to2p5", "pT2p5to3", "pT3to3p5", "pT3p5to4", "pT4to4p5", "pT4p5to5",
    "pT5to5p5", "pT5p5to6", "pT6to6p5", "pT6p5to7", "pT7to7p5", "pT7p5to8",
    "pT8to9", "pT9to10", "pT10to12p5", "pT12p5to15", "pT15to20", "pT20to25",
    "pT25to30", "pT30to40", "pT40to50", "pT50to70", "pT70to100"};

const Double_t MASS_LO = 1.75;
const Double_t MASS_HI = 1.99;
const Int_t MASS_NBIN = 48;

// --- the nested "Dpt > X GeV" MC sample thresholds (same 9 samples as getEfficiency.C /
// inputFiles_Prompt_MC_Sept5.txt). Cross-section/filterEff/Nevts and the resulting Leff-based
// combination now live in get_PtWeight_fromMass.C, where the 9 per-sample histograms this file
// writes get combined -- this file only needs the threshold values themselves, to identify
// which sample a given input file belongs to. ---
const int N_DPT_SAMPLES = 9;
Double_t Dpt_threshold[N_DPT_SAMPLES] = {0, 1, 8, 10, 20, 30, 40, 60, 80};

// Which of the 9 Dpt-threshold samples does this file belong to? Determined from the "Dpt{N}"
// tag in its path (e.g. ".../PromptD0ToKPi_Dpt8_Sept5/...") -- confirmed this session that the
// input file list groups files in contiguous per-sample blocks, but NOT sorted by threshold
// value, and a single chunk job's file range can span a sample boundary, so this has to be
// resolved per file, not assumed constant for a whole chunk.
int findDptSampleIndex(const string &filename)
{
  size_t pos = filename.find("Dpt");
  if (pos == string::npos)
    return -1;
  pos += 3;
  size_t end = pos;
  while (end < filename.size() && isdigit((unsigned char)filename[end]))
    end++;
  if (end == pos)
    return -1;
  int thr = std::stoi(filename.substr(pos, end - pos));

  // Dpt8's GEN fragment was confirmed (2026-09-23) to be byte-identical to Dpt1's (MomMinPt=1.0,
  // not 8.0) -- a production mislabeling, not a real higher-pT filter (full story in
  // get_PtWeight.C). Route its files directly into thr=1's slot here, at the source: this pools
  // its raw statistics with thr=1's own (more precision for the same effective sample) AND fixes
  // a second-order bug this caused -- fillMCPtSpectrum()'s sub-threshold cut below uses
  // Dpt_threshold[sampleIdx], so as long as Dpt8 files kept their own index, that cut was Pt<8,
  // wrongly discarding their perfectly legitimate pT=1-8 GeV entries (confirmed: the Sept24_v1
  // MC output has zero thr=8 entries below 8 GeV). Relabeling to thr=1 here makes the cut Pt<1
  // instead, matching what the sample's real filter actually justifies.
  if (thr == 8)
    thr = 1;

  for (int i = 0; i < N_DPT_SAMPLES; i++)
    if ((int)Dpt_threshold[i] == thr)
      return i;
  return -1;
}

bool isNonPromptMother(int pdg)
{
  int absPdg = std::abs(pdg);
  return ((absPdg / 100) % 10 == 5) || ((absPdg / 1000) % 10 == 5);
}

// Same guard as getEfficiency.C: BDTHandler::loadCuts() has no bounds checking, so
// scrub the CSV before the BDTHandler ctor parses it.
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

Int_t findFineBin(Double_t pt)
{
  for (int i = 0; i < N_PTBIN_eff; i++)
    if (pt >= ptbinning_[i] && pt < ptbinning_[i + 1])
      return i;
  return -1;
}

void fillDataMassHistograms(TString data_txt, Long64_t istart, Long64_t iend, BDTHandler &bdtHandler, TH1F *h_dmass[N_PTBIN_eff])
{
  ifstream file_stream(data_txt.Data());
  if (!file_stream.is_open())
  {
    cout << "Error: could not open data file list " << data_txt << endl;
    return;
  }

  string filename;
  Long64_t ifile = 0;
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
      cout << "Warning: Skipping bad data file: " << filename << endl;
      if (fin)
      {
        fin->Close();
        delete fin;
      }
      ifile++;
      continue;
    }

    cout << ">>> [DATA] ifile=" << ifile << " : " << filename << endl;

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

    Int_t candSize = 0;
    Int_t centrality = 0;
    vector<Float_t> *pT = nullptr;
    vector<Float_t> *y = nullptr;
    vector<Float_t> *mass = nullptr;
    vector<Float_t> *mva = nullptr;

    tree->SetBranchStatus("*", 0);
    for (const auto &p : {"candSize", "centrality", "pT", "y", "mass", "mva"})
      tree->SetBranchStatus(p, 1);

    tree->SetBranchAddress("candSize", &candSize);
    tree->SetBranchAddress("centrality", &centrality);
    tree->SetBranchAddress("pT", &pT);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("mva", &mva);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t ievt = 0; ievt < nEntries; ievt++)
    {
      tree->GetEntry(ievt);

      Int_t cent = centrality / 2;
      if (cent < MIN_CENT || cent >= MAX_CENT)
        continue;

      for (int icand = 0; icand < candSize; icand++)
      {
        float Pt = pT->at(icand);
        float Y = y->at(icand);
        float Mass = mass->at(icand);
        float Mva = mva->at(icand);

        if (Pt < MIN_PT_ANA || Pt >= MAX_PT_ANA || fabs(Y) >= MAX_Y_ANA)
          continue;

        double bdt_cut = bdtHandler.getBDTCut(Y, 2 * cent, Pt);
        if (Mva <= bdt_cut)
          continue;

        Int_t ipt = findFineBin(Pt);
        if (ipt < 0)
          continue;

        h_dmass[ipt]->Fill(Mass);
      }
    }

    fin->Close();
    delete fin;
    ifile++;
  }
  file_stream.close();
}

void fillMCPtSpectrum(TString mc_txt, Long64_t istart, Long64_t iend, BDTHandler &bdtHandler,
                       TH1D *h_pt_mc[N_DPT_SAMPLES], TH1D *h_pt_mc_noGenWeight[N_DPT_SAMPLES])
{
  ifstream file_stream(mc_txt.Data());
  if (!file_stream.is_open())
  {
    cout << "Error: could not open MC file list " << mc_txt << endl;
    return;
  }

  string filename;
  Long64_t ifile = 0;
  while (file_stream >> filename)
  {
    if (ifile < istart)
    {
      ifile++;
      continue;
    }
    if (ifile >= iend)
      break;

    const int sampleIdx = findDptSampleIndex(filename);
    if (sampleIdx < 0)
    {
      cout << "Warning: could not identify Dpt sample for " << filename << ", skipping" << endl;
      ifile++;
      continue;
    }

    TFile *fin = TFile::Open(filename.c_str());
    if (!fin || fin->IsZombie())
    {
      cout << "Warning: Skipping bad MC file: " << filename << endl;
      if (fin)
      {
        fin->Close();
        delete fin;
      }
      ifile++;
      continue;
    }

    cout << ">>> [MC] ifile=" << ifile << " (Dpt" << (int)Dpt_threshold[sampleIdx] << ") : " << filename << endl;

    TTree *tree = (TTree *)fin->Get("d0Analyzer/VCNtuple_D02kpi");
    TTree *t_eventinfoana = (TTree *)fin->Get("eventinfoana/EventInfoNtuple");

    if (!tree || !t_eventinfoana)
    {
      cout << "Skipping: Tree(s) missing in " << filename << endl;
      fin->Close();
      delete fin;
      ifile++;
      continue;
    }

    tree->AddFriend(t_eventinfoana);

    Int_t candSize = 0;
    Float_t genWeight = 1;

    vector<Bool_t> *matchGEN = nullptr;
    vector<Bool_t> *isSwap = nullptr;
    vector<Float_t> *pT = nullptr;
    vector<Float_t> *y = nullptr;
    vector<Float_t> *mva = nullptr;
    vector<bool> *Dgen_isPrompt = nullptr;

    tree->SetBranchStatus("*", 0);
    for (const auto &p : {"candSize", "y", "pT", "mva",
                           "matchGEN", "isSwap", "Dgen_isPrompt", "genWeight"})
      tree->SetBranchStatus(p, 1);

    tree->SetBranchAddress("candSize", &candSize);
    tree->SetBranchAddress("pT", &pT);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("mva", &mva);
    tree->SetBranchAddress("matchGEN", &matchGEN);
    tree->SetBranchAddress("isSwap", &isSwap);
    tree->SetBranchAddress("Dgen_isPrompt", &Dgen_isPrompt);
    tree->SetBranchAddress("genWeight", &genWeight);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t ievt = 0; ievt < nEntries; ievt++)
    {
      tree->GetEntry(ievt);

      // No centrality cut for MC (see note above the mode dispatch): every event is used.
      // No stitching weight here either: this event unambiguously belongs to the sample
      // identified by the file it came from (sampleIdx, fixed per file above), so it's just
      // genWeight -- cross-sample combination happens downstream, in get_PtWeight_fromMass.C.
      for (int icand = 0; icand < candSize; icand++)
      {
        if (matchGEN->at(icand) != 1 || isSwap->at(icand) != 0)
          continue;
        if (!Dgen_isPrompt->at(icand))
          continue;

        float Pt = pT->at(icand);
        float Y = y->at(icand);
        float Mva = mva->at(icand);

        if (Pt < MIN_PT_ANA || Pt >= MAX_PT_ANA || fabs(Y) >= MAX_Y_ANA)
          continue;

        // This sample's own generator filter only guarantees SOME D0 in the event has
        // pT >= Dpt_threshold[sampleIdx] -- not that THIS candidate does. Below its own
        // threshold, a candidate is very likely an unrelated companion D0 from the same (rare,
        // filter-selected) event topology, not an unbiased sample of the marginal pT spectrum
        // (confirmed this session by inspecting real Dpt80 entries at pT=2-3 GeV). Drop those
        // here at the source instead of relying on combineDptSamples()'s downstream eligibility
        // gate to discard them -- same effect, but the saved per-sample spectra are then
        // directly usable/inspectable without a separate combination step.
        if (Pt < Dpt_threshold[sampleIdx])
          continue;

        // Same methodology as getEfficiency.C: evaluate the BDT cut at each of the 5
        // representative centrality points rather than the candidate's own (meaningless,
        // and here uncut) centrality, filling once per point passed.
        for (int irep = 0; irep < N_CENT_REPR; irep++)
        {
          double bdt_cut = bdtHandler.getBDTCut(Y, 2 * centRepr[irep], Pt);
          if (Mva <= bdt_cut)
            continue;

          h_pt_mc[sampleIdx]->Fill(Pt, genWeight);
          h_pt_mc_noGenWeight[sampleIdx]->Fill(Pt); // same selection, weight=1, for comparison
        }
      }
    }

    fin->Close();
    delete fin;
    ifile++;
  }
  file_stream.close();
}

// One chunk of the data file list -> DataMass_<istart>_<iend>.root
// (per-fine-pT-bin mass histograms only; fit in macro 2 after hadd'ing chunks).
void runData(TString input_txt, TString output_path, int istart, int iend, BDTHandler &bdtHandler)
{
  TH1F *h_dmass[N_PTBIN_eff];
  for (int ipt = 0; ipt < N_PTBIN_eff; ipt++)
  {
    h_dmass[ipt] = new TH1F(Form("dmass_%s", ptLabel_[ipt]), Form("dmass_%s", ptLabel_[ipt]), MASS_NBIN, MASS_LO, MASS_HI);
    h_dmass[ipt]->SetDirectory(nullptr);
  }

  fillDataMassHistograms(input_txt, istart, iend, bdtHandler, h_dmass);

  TString outfile = TString::Format("%s/ROOT/DataMass_%d_%d.root", output_path.Data(), istart, iend);
  TFile *fout = new TFile(outfile, "RECREATE");
  fout->cd();
  for (int ipt = 0; ipt < N_PTBIN_eff; ipt++)
    h_dmass[ipt]->Write();
  fout->Write();
  fout->Close();

  cout << ">>> [DATA] Saved " << outfile << endl;
}

// One chunk of the MC file list -> MCSpectrum_<istart>_<iend>.root
// (one weighted reco pT spectrum PER Dpt-threshold sample, ready to hadd -- no fit needed for
// MC, and no cross-sample combination here either; see get_PtWeight_fromMass.C for that).
void runMC(TString input_txt, TString output_path, int istart, int iend, BDTHandler &bdtHandler)
{
  TH1D *h_pt_mc[N_DPT_SAMPLES];
  TH1D *h_pt_mc_noGenWeight[N_DPT_SAMPLES]; // same selection, weight=1 -- for comparing against genWeight's effect
  for (int i = 0; i < N_DPT_SAMPLES; i++)
  {
    h_pt_mc[i] = new TH1D(Form("h_pt_mc_raw_thr%d", (int)Dpt_threshold[i]), ";p_{T} (GeV);weighted reco counts / GeV", N_PTBIN_eff, ptbinning_);
    h_pt_mc[i]->SetDirectory(nullptr);
    h_pt_mc_noGenWeight[i] = new TH1D(Form("h_pt_mc_raw_thr%d_noGenWeight", (int)Dpt_threshold[i]), ";p_{T} (GeV);reco counts / GeV (no genWeight)", N_PTBIN_eff, ptbinning_);
    h_pt_mc_noGenWeight[i]->SetDirectory(nullptr);
  }

  fillMCPtSpectrum(input_txt, istart, iend, bdtHandler, h_pt_mc, h_pt_mc_noGenWeight);

  // ptbinning_ is variable-width (0.5 GeV at low pT, widening to 30 GeV at high pT), so raw
  // (un-normalized) bin content visually "bumps" at every width change even though the true
  // differential rate is falling smoothly -- purely an artifact of aggregating more phase space
  // into a wider bin, confirmed this session by comparing raw content to content/width directly.
  // Divide by bin width here (safe to do per-chunk: this is a fixed, chunk-independent factor
  // per bin, so it commutes with hadd's later bin-by-bin addition across chunks) so the saved
  // spectra are directly plottable without looking distorted, and macro 2 can use them as-is.
  for (int i = 0; i < N_DPT_SAMPLES; i++)
  {
    h_pt_mc[i]->Scale(1.0, "width");
    h_pt_mc_noGenWeight[i]->Scale(1.0, "width");
  }

  TString outfile = TString::Format("%s/ROOT/MCSpectrum_%d_%d.root", output_path.Data(), istart, iend);
  TFile *fout = new TFile(outfile, "RECREATE");
  fout->cd();
  for (int i = 0; i < N_DPT_SAMPLES; i++)
  {
    h_pt_mc[i]->Write();
    h_pt_mc_noGenWeight[i]->Write();
  }
  fout->Write();
  fout->Close();

  cout << ">>> [MC] Saved " << outfile << endl;
}

int main(int argc, char *argv[])
{
  if (argc != 6)
  {
    cout << "Usage: ./exe <data|mc> input_txt output_path istart iend" << endl;
    return 1;
  }

  TString mode = argv[1];
  TString input_txt = argv[2];
  TString output_path = argv[3];
  int istart = std::stoi(argv[4]);
  int iend = std::stoi(argv[5]);

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2();

  sanitizeBdtCutsFile("/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/AccxEff/bdt_cuts.csv");
  BDTHandler bdtHandler;

  if (mode == "data")
    runData(input_txt, output_path, istart, iend, bdtHandler);
  else if (mode == "mc")
    runMC(input_txt, output_path, istart, iend, bdtHandler);
  else
  {
    cout << "Unknown mode '" << mode << "'. Use 'data' or 'mc'." << endl;
    return 1;
  }

  return 0;
}
