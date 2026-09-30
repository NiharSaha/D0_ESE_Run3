// Standalone pT-weight derivation using ONLY the thr=1 MC sample -- no cross-sample combination,
// no reliance on the other Dpt-threshold samples' xsec/filterEff calibration (where several bugs
// were found this session: thr=8's fragment mislabeling, thr=30's low-statistics filterEff, and
// a still-unresolved Nevts inconsistency for thr>=10). Supersedes the earlier thr=0-only version
// of this macro (get_PtWeight_thr0only.C): thr=8's GEN fragment was confirmed byte-identical to
// thr=1's (MomMinPt=1.0, not 8.0), so get_MassSpectra_PtWeight.C now routes thr=8's files
// directly into thr=1's histogram at the source (2026-09-24) -- thr=1 here is really a pooled
// thr=1+thr=8 sample, with substantially more raw statistics than thr=0 alone ever had (e.g.
// 15081 vs. 5430 raw entries at [2,2.5) GeV, and now reaches all the way to [70,100) with a few
// entries, which thr=0 never did). thr=1 tracks Data well out to ~30-40 GeV, but still runs out
// of raw statistics beyond that. For that region, fit the Data/MC ratio over the statistically
// trustworthy range and extrapolate it to get a reasonable ratio factor at higher pT. A plain
// power law can't do this: the ratio falls to a minimum around 6-8 GeV and rises again, which a
// monotonic function structurally cannot follow (confirmed: chi2/ndf ~200 with a power law) --
// see the pol4 fit below instead.
//
// Does not modify get_MassSpectra_PtWeight.C / get_PtWeight.C beyond reading their existing
// outputs (the already-fitted Data yields from DataPtSpectrum_out.root, and the thr=1 raw MC
// spectrum from the hadd'ed MCSpectrum_out_combined.root).

#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <Math/MinimizerOptions.h>
#include <iostream>

using namespace std;

// pol4 evaluated at min(pT, par[5]) -- i.e. the fitted pol4 shape for pT < par[5], frozen at
// its value at par[5] beyond that. par[5] is set to fitMax after the fit below, so this gives
// exactly "use the fitted curve in the trusted range, hold the last value constant beyond it"
// as requested, instead of letting the raw polynomial run away past where it's constrained.
Double_t pol4Capped(Double_t *x, Double_t *par)
{
  Double_t ptEval = (x[0] < par[5]) ? x[0] : par[5];
  return par[0] + par[1] * ptEval + par[2] * ptEval * ptEval + par[3] * ptEval * ptEval * ptEval + par[4] * ptEval * ptEval * ptEval * ptEval;
}

// Minimum raw (un-weighted) MC entries in a bin to trust it for the ratio fit -- below
// this, Poisson noise dominates and the bin is left to the extrapolated fit instead. 20 raw
// entries corresponds to ~22% relative statistical uncertainty on that bin alone.
const Double_t MIN_RAW_ENTRIES_TRUSTED = 20.0;

void get_PtWeight_thr1only(
    TString data_infile = "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/AccxEff/DataPtSpectrum_out.root",
    // Sept24_v1, not v2: v1 predates the source-level thr=8->thr=1 routing fix, so its thr=1
    // histogram is pure Dpt1 statistics only (283974 entries), with none of Dpt8's folded in --
    // per request, Dpt8 (low raw stats, 787 entries even after the fix) is deliberately excluded
    // here rather than pooled in, unlike getEfficiency.C's ptWeightValue[] derivation.
    TString mc_infile = "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/PtWeight_Sept24_v1/ROOT/MCSpectrum_out_combined.root",
    TString outfile_name = "PtWeight_thr1only_out.root")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);

  TFile *fData = TFile::Open(data_infile);
  TFile *fMC = TFile::Open(mc_infile);
  if (!fData || fData->IsZombie() || !fMC || fMC->IsZombie())
  {
    cerr << "Error opening input files (" << data_infile << ", " << mc_infile << ")" << endl;
    return;
  }

  TH1D *h_pt_data_raw = (TH1D *)fData->Get("h_pt_data_raw");
  TH1D *h_pt_mc_thr1 = (TH1D *)fMC->Get("h_pt_mc_raw_thr1");             // already density (content/GeV)
  TH1D *h_pt_mc_thr1_noW = (TH1D *)fMC->Get("h_pt_mc_raw_thr1_noGenWeight");
  if (!h_pt_data_raw || !h_pt_mc_thr1 || !h_pt_mc_thr1_noW)
  {
    cerr << "Missing expected histograms in input files" << endl;
    return;
  }

  TH1D *h_pt_data = (TH1D *)h_pt_data_raw->Clone("h_pt_data");
  TH1D *h_pt_mc = (TH1D *)h_pt_mc_thr1->Clone("h_pt_mc");

  h_pt_data->Scale(1.0 / h_pt_data->Integral("width"));
  h_pt_mc->Scale(1.0 / h_pt_mc->Integral("width"));

  const Int_t nBins = h_pt_data->GetNbinsX();

  TH1D *h_ratio = (TH1D *)h_pt_data->Clone("h_ratio_data_over_mc");
  h_ratio->Divide(h_pt_mc);

  // Determine the trustworthy fit range from thr=1's own raw (un-weighted) entry counts, but
  // cap it at 30 GeV regardless -- kept from the thr=0-only version's setup (where the [30,40)
  // bin's stray statistics pulled the fit into a steep, fragile rise right at its own edge) for
  // now. thr=1's pooled thr=1+thr=8 statistics are substantially better in this region (88 raw
  // entries at [30,40), 65 at [40,50), vs. thr=0's 26/10) -- worth revisiting whether this cap
  // should move to 50 GeV, but left untouched here pending that separate decision.
  Double_t fitMax = h_pt_data->GetXaxis()->GetXmin();
  for (Int_t b = 1; b <= nBins; b++)
  {
    Double_t rawEntries = h_pt_mc_thr1_noW->GetBinContent(b) * h_pt_mc_thr1_noW->GetXaxis()->GetBinWidth(b);
    Double_t binHi = h_pt_data->GetXaxis()->GetBinUpEdge(b);
    if (rawEntries >= MIN_RAW_ENTRIES_TRUSTED && binHi <= 30.0)
      fitMax = binHi;
  }
  const Double_t fitMin = h_pt_data->GetXaxis()->GetXmin();
  cout << ">>> Trustworthy fit range (>= " << MIN_RAW_ENTRIES_TRUSTED << " raw thr=1 entries/bin, "
       << "capped at 30 GeV): ["
       << fitMin << ", " << fitMax << ") GeV" << endl;

  // A plain power law (p0*pT^p1) is strictly monotonic and cannot describe this ratio, which
  // falls to a minimum around 6-8 GeV and then rises again (confirmed directly in the fit
  // residuals: chi2/ndf ~200). Use a low-order polynomial instead -- pol2 is the lowest order
  // that can capture a single dip-then-rise, and staying low-order matters here specifically
  // because we're extrapolating a long way past the fit range (40 -> 100 GeV): a higher-order
  // polynomial (e.g. pol4, as used for the full-range fit in get_PtWeight_fromMass.C, where it
  // never needs to extrapolate beyond real MC statistics) risks curving unpredictably once
  // evaluated far outside where it was constrained.
  // pol2 (one turning point) wasn't enough: the measured ratio has TWO turning points in the
  // trusted range -- falls to a minimum around 7 GeV, rises to a local maximum around 15-20
  // GeV, then falls again out to ~30 GeV. pol3 (confirmed: chi2/ndf ~11) captures that shape
  // well in-range, but its negative leading coefficient sends it unphysically negative almost
  // immediately past the fit range. Using pol4 instead, per request -- standard chi2 fit (NOT
  // "L"/log-likelihood): h_ratio's bin content is Data/MC, a ratio of two densities with
  // properly Gaussian-propagated errors, not a raw Poisson event count -- "L" assumes the
  // latter and fails outright the moment the fit function dips non-positive during minimization.
  // NOTE: a higher-order polynomial fits the trusted range more tightly but its extrapolation
  // beyond that range is correspondingly LESS constrained/predictable, not more -- see the
  // per-bin printout below for how it actually behaves past fitMax.
  // "N": don't auto-attach this (uncapped) function to h_ratio's draw list -- the capped
  // f_ptWeight below is what actually gets drawn/used, and leaving the raw fit attached made
  // it show up too (a stray descending line past fitMax, separate from the capped plateau).
  TF1 *f_ratioFit = new TF1("f_ratioFit", "pol4", fitMin, 100.0);
  TFitResultPtr fitRes = h_ratio->Fit(f_ratioFit, "N Q S R", "", fitMin, fitMax);

  cout << "\n>>> Polynomial fit (Data/MC = p0 + p1*pT + p2*pT^2 + p3*pT^3 + p4*pT^4) over [" << fitMin << ", " << fitMax
       << "), status = " << fitRes->Status() << endl;
  cout << "    p0 = " << f_ratioFit->GetParameter(0) << " +- " << f_ratioFit->GetParError(0) << endl;
  cout << "    p1 = " << f_ratioFit->GetParameter(1) << " +- " << f_ratioFit->GetParError(1) << endl;
  cout << "    p2 = " << f_ratioFit->GetParameter(2) << " +- " << f_ratioFit->GetParError(2) << endl;
  cout << "    p3 = " << f_ratioFit->GetParameter(3) << " +- " << f_ratioFit->GetParError(3) << endl;
  cout << "    p4 = " << f_ratioFit->GetParameter(4) << " +- " << f_ratioFit->GetParError(4) << endl;
  cout << "    chi2/ndf = " << f_ratioFit->GetChisquare() << " / " << f_ratioFit->GetNDF() << endl;

  // The pT-weighting factor to actually use: pol4 in the trusted range, frozen beyond fitMax at
  // the pol4 fit's value at the LAST TRUSTED BIN'S OWN CENTER (not at the fitMax edge itself --
  // the curve is often still visibly moving right at that edge, so evaluating exactly there can
  // undershoot/overshoot the last trusted bin's own characteristic value; freezing at that bin's
  // center instead uses the value the fit already reports for that bin, per request).
  Double_t lastTrustedBinCenter = fitMax;
  for (Int_t b = 1; b <= nBins; b++)
  {
    if (h_pt_data->GetXaxis()->GetBinUpEdge(b) == fitMax)
      lastTrustedBinCenter = 0.5 * (h_pt_data->GetXaxis()->GetBinLowEdge(b) + h_pt_data->GetXaxis()->GetBinUpEdge(b));
  }

  TF1 *f_ptWeight = new TF1("f_ptWeight", pol4Capped, fitMin, 100.0, 6);
  for (Int_t p = 0; p < 5; p++)
    f_ptWeight->SetParameter(p, f_ratioFit->GetParameter(p));
  f_ptWeight->FixParameter(5, lastTrustedBinCenter);
  const Double_t frozenValue = f_ptWeight->Eval(lastTrustedBinCenter);
  cout << "\n>>> pT weight frozen beyond fitMax=" << fitMax << " GeV, using the pol4 fit's value at "
       << "the last trusted bin's center (" << lastTrustedBinCenter << " GeV) = " << frozenValue << endl;

  TH1D *h_ptWeight = (TH1D *)h_pt_data->Clone("h_ptWeight");
  h_ptWeight->Reset();
  h_ptWeight->GetYaxis()->SetTitle("pT weight (Data/MC, thr=1 + pol4, frozen beyond fitMax)");
  for (Int_t b = 1; b <= nBins; b++)
  {
    Double_t center = 0.5 * (h_ptWeight->GetXaxis()->GetBinLowEdge(b) + h_ptWeight->GetXaxis()->GetBinUpEdge(b));
    h_ptWeight->SetBinContent(b, f_ptWeight->Eval(center));
  }

  cout << "\n>>> pT WEIGHTING FACTOR per bin (this is what to use downstream):\n";
  cout << "    [pT_lo, pT_hi)      pt_weight    source\n";
  for (Int_t b = 1; b <= nBins; b++)
  {
    Double_t lo = h_pt_data->GetXaxis()->GetBinLowEdge(b);
    Double_t hi = h_pt_data->GetXaxis()->GetBinUpEdge(b);
    Double_t rawEntries = h_pt_mc_thr1_noW->GetBinContent(b) * h_pt_mc_thr1_noW->GetXaxis()->GetBinWidth(b);
    Double_t ratioVal = h_pt_mc->GetBinContent(b) > 0 ? h_ratio->GetBinContent(b) : -1;
    Double_t weightVal = h_ptWeight->GetBinContent(b);
    printf("[%6.1f,%7.1f)  %12.5g   %s (raw_thr1_entries=%.0f, measured_ratio=%.4g)\n",
           lo, hi, weightVal,
           (hi <= fitMax) ? "pol4 fit" : "FROZEN at fitMax",
           rawEntries, ratioVal);
  }

  // --- overlay canvas ---
  TFile *fout = new TFile(outfile_name, "RECREATE");

  TCanvas *c1 = new TCanvas("cOverlay", "cOverlay", 800, 800);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();

  h_pt_data->SetTitle("");
  h_pt_data->SetMarkerStyle(20);
  h_pt_data->SetMarkerColor(kBlack);
  h_pt_data->SetLineColor(kBlack);
  h_pt_data->GetYaxis()->SetTitle("Normalized dN/dp_{T}");
  h_pt_data->GetXaxis()->SetLabelSize(0);
  h_pt_data->Draw("E1");

  h_pt_mc->SetLineColor(kRed + 1);
  h_pt_mc->SetLineWidth(2);
  h_pt_mc->SetMarkerStyle(0);
  h_pt_mc->Draw("HIST SAME");

  TLegend *leg = new TLegend(0.5, 0.65, 0.88, 0.85);
  leg->AddEntry(h_pt_data, "Data", "lep");
  leg->AddEntry(h_pt_mc, "MC (thr=1 only)", "l");
  leg->SetBorderSize(0);
  leg->Draw();

  c1->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();
  pad2->cd();

  h_ratio->SetTitle("");
  h_ratio->GetYaxis()->SetTitle("Data / MC");
  h_ratio->GetYaxis()->SetRangeUser(0, 1.2 * frozenValue);
  h_ratio->GetXaxis()->SetTitle("p_{T} (GeV)");
  h_ratio->SetMarkerStyle(20);
  h_ratio->Draw("E1");

  f_ptWeight->SetLineColor(kBlue + 1);
  f_ptWeight->SetLineStyle(2);
  f_ptWeight->SetLineWidth(2);
  f_ptWeight->SetNpx(1000);
  f_ptWeight->Draw("SAME");

  TLine *line1 = new TLine(h_ratio->GetXaxis()->GetXmin(), 1, h_ratio->GetXaxis()->GetXmax(), 1);
  line1->SetLineStyle(3);
  line1->Draw();
  TLine *lineFitMax = new TLine(fitMax, 0, fitMax, 1.2 * frozenValue);
  lineFitMax->SetLineStyle(2);
  lineFitMax->SetLineColor(kGray + 2);
  lineFitMax->Draw();

  c1->Write();
  c1->Print("PtWeight_thr1only_Overlay.png");

  h_pt_data->Write();
  h_pt_mc->Write();
  h_ratio->Write();
  h_ptWeight->Write();
  f_ratioFit->Write("f_ratioFit_pol4_raw");
  f_ptWeight->Write("f_ptWeight_capped");
  fitRes->Write("fitResult_ratioFit");

  fout->Write();
  fout->Close();

  cout << "\n>>> Done. Full output (histograms, fit function, fit result, canvas) saved to "
       << outfile_name << " -- open with a TBrowser or TFile::Open() to inspect directly.\n";
  cout << ">>> Overlay PNG also saved to PtWeight_thr1only_Overlay.png\n";
}
