// =============================================================================
//  scan_significance_vs_BDT.C
//
//  Data-driven justification for the size of a BDT-cut systematic variation.
//
//  Reads the TMVA training/test trees in BDT_training/TMVA.root
//  (dataset/TestTree + dataset/TrainTree), which carry, per candidate:
//     classID (0=signal, 1=background), the evaluated BDT score,
//     and npT / ncent / ny -- the SAME kinematic variables BDTHandler::get_bins
//     slices the analysis cut table on.
//  -> lets us rebuild a per-(y,cent,pT) proxy figure-of-merit-vs-cut curve in
//     the exact BDTHandler binning (2 y-bins x 4 cent-bins x 13 pT-bins).
//
//  Proxy FoM = S/sqrt(S+B) (unweighted training-sample counts; weight branch
//  is 1 for every event in this file).
//
//  IMPORTANT geometry finding (confirmed on real curves before committing to
//  this method): the analysis's nominal BDT cut (bdt_cuts.csv) sits well past
//  this proxy FoM's own maximum, on its monotonically-declining tail -- most
//  likely because the actual working point was tuned with a mass-fit-based
//  significance = S/S_error, which is a different functional form that peaks
//  elsewhere. Consequence: walking *tighter* from nominal, the proxy FoM keeps
//  falling (well-posed -- a 1%-drop point exists nearby). Walking *looser*
//  from nominal, the proxy FoM only *rises* back toward its own maximum (not
//  well-posed -- there is no nearby 1%-drop point in that direction). So:
//     delta_tight = cut_tight - nominal   (found from the well-posed tight-side
//                                          1% search)
//     cut_loose   = nominal - delta_tight (same shift, mirrored)
//  i.e. the window MAGNITUDE is derived from the one direction where the local
//  curve behavior is meaningful, then applied symmetrically.
//
//  Coverage caveat (found by inspection, confirmed below): the training only
//  covers 1 <= pT < 20 GeV (BDTHandler pT-bins 0-8). pT-bins 9-12 (20-900 GeV)
//  have zero training statistics -- bdt_cuts.csv already uses hand-set
//  placeholder cuts there. For those bins delta_tight from the last covered
//  bin (pT-bin 8, 15-20 GeV, same y/cent) is applied (absolute BDT-score
//  shift) to the nominal cut. A fractional/relative shift would be meaningless
//  on a sentinel value like 1.0 (reject-all) or -1.0 (accept-all), so sentinel
//  nominal cuts (|v|>=1), in any pT bin, are left unshifted and flagged --
//  there is no cut to perturb.
//
//  Outputs (into this directory):
//     significance_scan.root      -- per-bin FoM-vs-cut TGraph + raw hS/hB
//     significance_plots/*.png    -- FoM-vs-cut, pT-bins 0-8 overlaid, per (y,cent)
//     bdt_cuts_loose.csv          -- same format as bdt_cuts.csv
//     bdt_cuts_tight.csv          -- same format as bdt_cuts.csv
//     a printed summary table (also useful piped to a log file)
//
//  Nothing outside Systematics/BDT_sys/ is read or written; BDTHandler.{h,cc}
//  and bdt_cuts.csv (already in this directory) are reused UNMODIFIED.
//
//  Build / run:
//     g++ scan_significance_vs_BDT.C $(root-config --cflags --libs) -Wall -O2 -o scan_significance_vs_BDT.exe
//     ./scan_significance_vs_BDT.exe
// =============================================================================

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/BDT_sys/MC_based/BDTHandler.h"
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/BDT_sys/MC_based/BDTHandler.cc"

using namespace std;

// ---------------------------------------------------------------------------
//  Bin edges -- must match BDTHandler::get_bins() (private, so duplicated here)
// ---------------------------------------------------------------------------
static const int NY = 2, NCENT = 4, NPT = 13;

static const double pt_lo[NPT] = {1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 40, 60, 100};
static const double pt_hi[NPT] = {2, 3, 4, 5, 6, 8, 10, 15, 20, 40, 60, 100, 900};
static const char *pt_label[NPT] = {"1-2", "2-3", "3-4", "4-5", "5-6", "6-8", "8-10",
                                    "10-15", "15-20", "20-40", "40-60", "60-100", "100-900"};
static const int PT_LAST_COVERED = 8; // 15-20 GeV: last bin with training statistics

static const int cent_lo[NCENT] = {0, 20, 60, 100}; // 2 x cent% units
static const int cent_hi[NCENT] = {20, 60, 100, 180};
static const char *cent_label[NCENT] = {"0-10%", "10-30%", "30-50%", "50-90%"};

static const char *y_label[NY] = {"|y|<1", "1<=|y|<3"};

static int y_bin_of(float y)
{
  float ay = fabs(y);
  if (ay < 1) return 0;
  if (ay < 3) return 1;
  return -1;
}
static int cent_bin_of(int c2)
{
  for (int i = 0; i < NCENT; ++i)
    if (c2 >= cent_lo[i] && c2 < cent_hi[i]) return i;
  return -1;
}
static int pt_bin_of(float pt)
{
  for (int i = 0; i < NPT; ++i)
    if (pt >= pt_lo[i] && pt < pt_hi[i]) return i;
  return -1;
}
static bool isSentinel(double v) { return fabs(v) >= 1.0 - 1e-6; }

// ---------------------------------------------------------------------------
//  BDT-score histogram binning for the scan
// ---------------------------------------------------------------------------
static const int NSCAN = 150;
static const double SCAN_LO = -0.65, SCAN_HI = 0.90;

// minimum-statistics guard (applies uniformly -- also catches thin y=1/low-pT
// corners inside the nominally-covered pT<20 GeV range, not just pT>=20 GeV)
static const double MIN_SIG = 100, MIN_BKG = 200;
static const double FOM_DROP = 0.01; // 1% relative drop sizes the window

// ---------------------------------------------------------------------------
//  proxy FoM curve for one (y,cent,pT) bin
// ---------------------------------------------------------------------------
struct ScanCurve
{
  bool derived = false; // passed the minimum-statistics guard
  double nSig = 0, nBkg = 0;
  int nb = 0;
  vector<double> cutAt, fom; // valid index range 1..nb
};

static ScanCurve buildCurve(TH1D *hS, TH1D *hB)
{
  ScanCurve r;
  int nb = hS->GetNbinsX();
  r.nb = nb;
  r.nSig = hS->Integral(0, nb + 1);
  r.nBkg = hB->Integral(0, nb + 1);
  if (r.nSig < MIN_SIG || r.nBkg < MIN_BKG)
    return r; // derived stays false

  vector<double> Stail(nb + 2, 0.0), Btail(nb + 2, 0.0);
  Stail[nb + 1] = hS->GetBinContent(nb + 1); // overflow
  Btail[nb + 1] = hB->GetBinContent(nb + 1);
  for (int i = nb; i >= 1; --i)
  {
    Stail[i] = Stail[i + 1] + hS->GetBinContent(i);
    Btail[i] = Btail[i + 1] + hB->GetBinContent(i);
  }
  r.fom.assign(nb + 1, 0.0);
  r.cutAt.assign(nb + 1, 0.0);
  for (int i = 1; i <= nb; ++i)
  {
    double S = Stail[i], B = Btail[i];
    r.fom[i] = (S + B > 0) ? S / sqrt(S + B) : 0.0;
    r.cutAt[i] = hS->GetXaxis()->GetBinLowEdge(i); // "accept if score > cutAt[i]"
  }
  r.derived = true;
  return r;
}

// scan's own maximum (diagnostic only -- shows how far the proxy FoM's peak
// sits from the analysis's actual nominal cut; not used to define the window)
static void scanPeak(const ScanCurve &c, double &cut_opt, double &fom_max)
{
  cut_opt = 0; fom_max = -1;
  for (int i = 1; i <= c.nb; ++i)
    if (c.fom[i] > fom_max) { fom_max = c.fom[i]; cut_opt = c.cutAt[i]; }
}

struct Window
{
  bool ok = false;          // false -> nominal is a sentinel, or curve not derived
  bool tightBoundary = false;
  double delta_tight = 0;   // >=0, from the well-posed (tight) direction
  double cut_loose = 0, cut_tight = 0;
  double fom_nom = 0;
};

// Anchor at the ACTUAL nominal cut. Only the tight direction is well-posed
// here (the proxy FoM's own maximum sits looser than nominal, so walking
// looser from nominal only climbs back toward that maximum -- see header
// note). delta_tight is measured on the well-posed side and mirrored.
static Window windowAtNominal(const ScanCurve &c, double nominal)
{
  Window w;
  if (isSentinel(nominal) || !c.derived) return w;

  int i_nom = 1; double best = 1e18;
  for (int i = 1; i <= c.nb; ++i)
  {
    double d = fabs(c.cutAt[i] - nominal);
    if (d < best) { best = d; i_nom = i; }
  }
  w.fom_nom = c.fom[i_nom];
  if (w.fom_nom <= 0) return w;

  double thr = (1.0 - FOM_DROP) * w.fom_nom;
  int iT = i_nom;
  while (iT < c.nb && c.fom[iT] > thr) iT++;
  w.tightBoundary = (iT == c.nb && c.fom[iT] > thr);

  w.delta_tight = fabs(c.cutAt[iT] - nominal);
  w.cut_tight = nominal + w.delta_tight;
  w.cut_loose = nominal - w.delta_tight;
  w.ok = true;
  return w;
}

static void fillFromTree(TTree *t, TH1D *hS[NY][NCENT][NPT], TH1D *hB[NY][NCENT][NPT])
{
  if (!t) return;
  Int_t classID = -1;
  Float_t bdt = 0, npt = 0, ncent = 0, ny = 0, weight = 1;
  t->SetBranchAddress("classID", &classID);
  t->SetBranchAddress("BDT", &bdt);
  t->SetBranchAddress("npT", &npt);
  t->SetBranchAddress("ncent", &ncent);
  t->SetBranchAddress("ny", &ny);
  t->SetBranchAddress("weight", &weight);

  Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; ++i)
  {
    t->GetEntry(i);
    int yb = y_bin_of(ny);
    int cb = cent_bin_of((int)ncent);
    int pb = pt_bin_of(npt);
    if (yb < 0 || cb < 0 || pb < 0) continue;
    if (classID == 0) hS[yb][cb][pb]->Fill(bdt, weight);
    else if (classID == 1) hB[yb][cb][pb]->Fill(bdt, weight);
  }
}

static void writeCsv(const TString &path, double cut[NY][NCENT][NPT])
{
  ofstream f(path.Data());
  f << "# this file contains all the bdt_cuts\n";
  f << "# generated by scan_significance_vs_BDT.C from BDT_training/TMVA.root\n";
  f << "# proxy significance = S/sqrt(S+B) on TestTree+TrainTree (unweighted counts);\n";
  f << "# window anchored at the nominal cut, sized by a 1% relative FoM drop on the\n";
  f << "# well-posed (tighter) side and mirrored -- see printed table for per-bin status.\n";
  f << "# pT>=20 GeV bins extrapolated (absolute score shift) from the 15-20 GeV bin.\n";
  f << "# y_bins: 0-(0_1) 1-(1_3)\n";
  f << "# cent_bins: 0-(0_10) 1-(10_30) 2-(30_50) 3-(50_90)\n";
  f << "# pt_bins 0-(1_2) 1-(2_3) 2-(3_4) 3-(4_5) 4-(5_6) 5-(6_8) 6-(8_10) 7-(10_15) 8-(15_20) 9-(20-40) 10-(40-60) 11-(60-100) 12-(100-900)\n";
  f << "#\n";
  f << "#y_bin, cent_bin, pt_bin, cut_value\n";
  for (int y = 0; y < NY; ++y)
    for (int c = 0; c < NCENT; ++c)
      for (int p = 0; p < NPT; ++p)
        f << y << " " << c << " " << p << " " << cut[y][c][p] << "\n";
  f.close();
  cout << ">>> wrote " << path << endl;
}

void scan_significance_vs_BDT()
{
  gStyle->SetOptStat(0);
  const TString outdir = "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/BDT_sys";
  const TString tmvaFile = outdir + "/BDT_training/TMVA.root";
  gSystem->Exec("mkdir -p " + outdir + "/significance_plots");

  TFile *fin = TFile::Open(tmvaFile);
  if (!fin || fin->IsZombie()) { cerr << "FATAL: cannot open " << tmvaFile << endl; return; }
  TTree *testT = (TTree *)fin->Get("dataset/TestTree");
  TTree *trainT = (TTree *)fin->Get("dataset/TrainTree");
  if (!testT || !trainT) { cerr << "FATAL: TestTree/TrainTree not found in " << tmvaFile << endl; return; }

  BDTHandler bdtHandler; // reused unmodified -> reads the LOCAL bdt_cuts.csv (nominal)

  // -------------------------------------------------------------------------
  //  book + fill BDT-score histograms per (y,cent,pT) bin
  // -------------------------------------------------------------------------
  TH1D *hS[NY][NCENT][NPT], *hB[NY][NCENT][NPT];
  for (int y = 0; y < NY; ++y)
    for (int c = 0; c < NCENT; ++c)
      for (int p = 0; p < NPT; ++p)
      {
        hS[y][c][p] = new TH1D(Form("hS_y%d_c%d_p%d", y, c, p), "", NSCAN, SCAN_LO, SCAN_HI);
        hB[y][c][p] = new TH1D(Form("hB_y%d_c%d_p%d", y, c, p), "", NSCAN, SCAN_LO, SCAN_HI);
        hS[y][c][p]->SetDirectory(nullptr);
        hB[y][c][p]->SetDirectory(nullptr);
      }
  fillFromTree(testT, hS, hB);
  fillFromTree(trainT, hS, hB);
  fin->Close();

  // -------------------------------------------------------------------------
  //  scan each bin, build nominal / loose / tight tables
  // -------------------------------------------------------------------------
  ScanCurve curve[NY][NCENT][NPT];
  double nominal[NY][NCENT][NPT], loose[NY][NCENT][NPT], tight[NY][NCENT][NPT];
  double cut_opt[NY][NCENT][NPT] = {}, fom_max[NY][NCENT][NPT] = {};
  double delta_tight[NY][NCENT][NPT] = {};
  string status[NY][NCENT][NPT];

  for (int y = 0; y < NY; ++y)
    for (int c = 0; c < NCENT; ++c)
      for (int p = 0; p < NPT; ++p)
      {
        float y_rep = (y == 0) ? 0.5f : 2.0f;
        float cent_rep = 0.5f * (cent_lo[c] + cent_hi[c]);
        float pt_rep = pt_lo[p] + 0.5f;
        nominal[y][c][p] = bdtHandler.getBDTCut(y_rep, (int)cent_rep, pt_rep);

        if (p <= PT_LAST_COVERED)
        {
          curve[y][c][p] = buildCurve(hS[y][c][p], hB[y][c][p]);
          if (!curve[y][c][p].derived)
          {
            status[y][c][p] = "insufficient-stats";
            loose[y][c][p] = tight[y][c][p] = nominal[y][c][p];
            continue;
          }
          scanPeak(curve[y][c][p], cut_opt[y][c][p], fom_max[y][c][p]);
          if (isSentinel(nominal[y][c][p]))
          {
            status[y][c][p] = "sentinel-nominal";
            loose[y][c][p] = tight[y][c][p] = nominal[y][c][p];
            continue;
          }
          Window w = windowAtNominal(curve[y][c][p], nominal[y][c][p]);
          if (!w.ok)
          {
            status[y][c][p] = "window-failed";
            loose[y][c][p] = tight[y][c][p] = nominal[y][c][p];
            continue;
          }
          delta_tight[y][c][p] = w.delta_tight;
          loose[y][c][p] = w.cut_loose;
          tight[y][c][p] = w.cut_tight;
          status[y][c][p] = w.tightBoundary ? "derived (tight-boundary)" : "derived";
        }
        else
        {
          // pT >= 20 GeV: no training statistics; extrapolate delta_tight from pT-bin 8
          const string &srcStatus = status[y][c][PT_LAST_COVERED];
          bool srcOk = (srcStatus == "derived" || srcStatus == "derived (tight-boundary)");
          if (isSentinel(nominal[y][c][p]))
          {
            status[y][c][p] = "no-coverage-sentinel-kept";
            loose[y][c][p] = tight[y][c][p] = nominal[y][c][p];
          }
          else if (srcOk)
          {
            double d = delta_tight[y][c][PT_LAST_COVERED];
            delta_tight[y][c][p] = d;
            loose[y][c][p] = nominal[y][c][p] - d;
            tight[y][c][p] = nominal[y][c][p] + d;
            status[y][c][p] = "no-coverage-extrapolated";
          }
          else
          {
            status[y][c][p] = "no-coverage-no-source";
            loose[y][c][p] = tight[y][c][p] = nominal[y][c][p];
          }
        }
      }

  // -------------------------------------------------------------------------
  //  outputs: root file (graphs), csvs, printed table
  // -------------------------------------------------------------------------
  TFile *fout = new TFile(outdir + "/significance_scan.root", "RECREATE");
  fout->cd();
  for (int y = 0; y < NY; ++y)
    for (int c = 0; c < NCENT; ++c)
      for (int p = 0; p <= PT_LAST_COVERED; ++p)
      {
        hS[y][c][p]->Write();
        hB[y][c][p]->Write();
        if (curve[y][c][p].derived)
        {
          TGraph *g = new TGraph(curve[y][c][p].nb, &curve[y][c][p].cutAt[1], &curve[y][c][p].fom[1]);
          g->SetName(Form("fom_y%d_c%d_p%d", y, c, p));
          g->SetTitle(Form("FoM vs BDT cut: y=%s cent=%s pT=%s GeV;BDT cut;S/#sqrt{S+B}",
                           y_label[y], cent_label[c], pt_label[p]));
          g->Write();
        }
      }

  cout << fixed << setprecision(4);
  cout << "\n===================================================================================================================\n";
  cout << left << setw(10) << "y" << setw(10) << "cent" << setw(10) << "pT"
       << setw(10) << "nSig" << setw(10) << "nBkg"
       << setw(12) << "nominal" << setw(12) << "scan_peak"
       << setw(12) << "loose" << setw(12) << "tight" << setw(10) << "delta"
       << setw(26) << "status" << "\n";
  cout << "-------------------------------------------------------------------------------------------------------------------\n";
  int nStatus[8] = {0};
  const char *statusNames[8] = {"derived", "derived (tight-boundary)", "insufficient-stats",
                                "sentinel-nominal", "no-coverage-extrapolated",
                                "no-coverage-sentinel-kept", "no-coverage-no-source", "window-failed"};
  for (int y = 0; y < NY; ++y)
    for (int c = 0; c < NCENT; ++c)
      for (int p = 0; p < NPT; ++p)
      {
        const auto &st = status[y][c][p];
        for (int k = 0; k < 8; ++k) if (st == statusNames[k]) nStatus[k]++;

        cout << left << setw(10) << y_label[y] << setw(10) << cent_label[c] << setw(10) << pt_label[p]
             << setw(10) << Form("%.0f", (p <= PT_LAST_COVERED) ? curve[y][c][p].nSig : -1.0)
             << setw(10) << Form("%.0f", (p <= PT_LAST_COVERED) ? curve[y][c][p].nBkg : -1.0)
             << setw(12) << nominal[y][c][p]
             << setw(12) << ((p <= PT_LAST_COVERED && curve[y][c][p].derived) ? cut_opt[y][c][p] : 0.0)
             << setw(12) << loose[y][c][p]
             << setw(12) << tight[y][c][p]
             << setw(10) << delta_tight[y][c][p]
             << setw(26) << st << "\n";
      }
  cout << "===================================================================================================================\n";
  cout << "bins: ";
  for (int k = 0; k < 8; ++k) cout << statusNames[k] << "=" << nStatus[k] << "  ";
  cout << " (total=" << NY * NCENT * NPT << ")\n\n";

  writeCsv(outdir + "/bdt_cuts_loose.csv", loose);
  writeCsv(outdir + "/bdt_cuts_tight.csv", tight);

  fout->Write(0, TObject::kOverwrite);
  fout->Close();

  // -------------------------------------------------------------------------
  //  plots: per (y,cent), overlay the 9 covered pT-bin FoM-vs-cut curves
  // -------------------------------------------------------------------------
  int colors[9] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
                   kOrange + 1, kCyan + 2, kGray + 2, kViolet + 1};
  const char *centTag[NCENT] = {"cent0to10", "cent10to30", "cent30to50", "cent50to90"};
  for (int y = 0; y < NY; ++y)
    for (int c = 0; c < NCENT; ++c)
    {
      TCanvas *cv = new TCanvas(Form("c_y%d_c%d", y, c), "", 800, 600);
      TLegend *leg = new TLegend(0.65, 0.15, 0.89, 0.55);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.028);

      // common y-axis range across all pT curves in this (y,cent) panel, so no
      // curve gets clipped by whichever one happens to be drawn first
      double yMax = 0;
      for (int p = 0; p <= PT_LAST_COVERED; ++p)
        if (curve[y][c][p].derived)
          for (int i = 1; i <= curve[y][c][p].nb; ++i)
            yMax = max(yMax, curve[y][c][p].fom[i]);

      bool first = true;
      for (int p = 0; p <= PT_LAST_COVERED; ++p)
      {
        if (!curve[y][c][p].derived) continue;
        TGraph *g = new TGraph(curve[y][c][p].nb, &curve[y][c][p].cutAt[1], &curve[y][c][p].fom[1]);
        g->SetLineColor(colors[p]);
        g->SetLineWidth(2);
        g->SetTitle(Form("y=%s  cent=%s;BDT cut;S/#sqrt{S+B}", y_label[y], cent_label[c]));
        if (first) { g->SetMinimum(0); g->SetMaximum(1.1 * yMax); }
        g->Draw(first ? "AL" : "L SAME");
        leg->AddEntry(g, Form("pT %s GeV", pt_label[p]), "l");
        first = false;
      }
      if (!first) { leg->Draw(); cv->SaveAs(Form("%s/significance_plots/fom_y%d_%s.png", outdir.Data(), y, centTag[c])); }
      delete cv;
    }

  cout << ">>> significance_scan.root, bdt_cuts_loose.csv, bdt_cuts_tight.csv, "
       << "significance_plots/*.png written to " << outdir << endl;
}

int main()
{
  scan_significance_vs_BDT();
  return 0;
}
