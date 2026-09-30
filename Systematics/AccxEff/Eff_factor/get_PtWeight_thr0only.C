// Standalone pT-weight derivation using ONLY the thr=0 (unfiltered, "Dpt>0") MC sample --
// no cross-sample combination, no reliance on the other 8 Dpt-threshold samples' xsec/filterEff
// calibration (where several bugs were found this session: thr=8's fragment mislabeling,
// thr=30's low-statistics filterEff, and a still-unresolved Nevts inconsistency for thr>=10).
// thr=0 alone tracks Data well out to ~30-40 GeV (confirmed directly, ratio ~0.75-1.6), but has
// essentially no raw statistics beyond that (10 raw entries at [40,50), 0 beyond). For that
// region, fit the Data/MC ratio over the statistically trustworthy range and extrapolate it to
// get a reasonable ratio factor at higher pT. A plain power law can't do this: the ratio falls
// to a minimum around 6-8 GeV and rises again, which a monotonic function structurally cannot
// follow (confirmed: chi2/ndf ~200 with a power law) -- see the pol2 fit below instead.
//
// Does not modify or depend on get_MassSpectra_PtWeight.C / get_PtWeight_fromMass.C beyond
// reading their existing outputs (the already-fitted Data yields from PtWeight_out.root, and
// the thr=0 raw MC spectrum from the hadd'ed MCSpectrum_out_combined.root).

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
#include <TMath.h>
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

void get_PtWeight_thr0only(
    TString data_infile = "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/AccxEff/Eff_factor/DataPtSpectrum_out.root",
    TString mc_infile = "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/PtWeight_Sept24_v1/ROOT/MCSpectrum_out_combined.root",
    TString outfile_name = "PtWeight_thr0only_out.root")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);

  TFile *fData = TFile::Open(data_infile);
  TFile *fMC = TFile::Open(mc_infile);
  if (!fData || fData->IsZombie() || !fMC || fMC->IsZombie())
  {
    cerr << "Error opening input files (" << data_infile << ", " << mc_infile << ")" << endl;
    return;
  }

  TH1D *h_pt_data_raw = (TH1D *)fData->Get("h_pt_data_raw");
  TH1D *h_pt_mc_thr0 = (TH1D *)fMC->Get("h_pt_mc_raw_thr0");             // already density (content/GeV)
  TH1D *h_pt_mc_thr0_noW = (TH1D *)fMC->Get("h_pt_mc_raw_thr0_noGenWeight");
  if (!h_pt_data_raw || !h_pt_mc_thr0 || !h_pt_mc_thr0_noW)
  {
    cerr << "Missing expected histograms in input files" << endl;
    return;
  }

  TH1D *h_pt_data = (TH1D *)h_pt_data_raw->Clone("h_pt_data");
  TH1D *h_pt_mc = (TH1D *)h_pt_mc_thr0->Clone("h_pt_mc");

  h_pt_data->Scale(1.0 / h_pt_data->Integral("width"));
  h_pt_mc->Scale(1.0 / h_pt_mc->Integral("width"));

  const Int_t nBins = h_pt_data->GetNbinsX();

  TH1D *h_ratio = (TH1D *)h_pt_data->Clone("h_ratio_data_over_mc");
  h_ratio->Divide(h_pt_mc);

  // Determine the trustworthy fit range from thr=0's own raw (un-weighted) entry counts, but
  // cap it at 30 GeV regardless: the [30,40) bin technically clears MIN_RAW_ENTRIES_TRUSTED (26
  // raw entries) but is still noisy enough (measured ratio 3.35 vs. a smooth ~1.2-1.4 trend in
  // the neighboring bins) that including it pulled the pol4 fit into a steep, fragile rise right
  // at its own edge -- freezing there gave a plateau (5.16) driven almost entirely by that one
  // bin. Per request, use [25,30) instead as the last trusted bin (101 raw entries, ~10%
  // relative uncertainty) and freeze from 30 GeV onward at the fit's value there.
  Double_t fitMax = h_pt_data->GetXaxis()->GetXmin();
  for (Int_t b = 1; b <= nBins; b++)
  {
    Double_t rawEntries = h_pt_mc_thr0_noW->GetBinContent(b) * h_pt_mc_thr0_noW->GetXaxis()->GetBinWidth(b);
    Double_t binHi = h_pt_data->GetXaxis()->GetBinUpEdge(b);
    if (rawEntries >= MIN_RAW_ENTRIES_TRUSTED && binHi <= 30.0)
      fitMax = binHi;
  }
  const Double_t fitMin = h_pt_data->GetXaxis()->GetXmin();
  cout << ">>> Trustworthy fit range (>= " << MIN_RAW_ENTRIES_TRUSTED << " raw thr=0 entries/bin, "
       << "capped at 30 GeV to exclude the noisy [30,40) bin): ["
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
  // the curve is still visibly descending right at that edge, so evaluating exactly there
  // (0.828) undershoot the [25,30) bin's own characteristic value; freezing at that bin's
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
  h_ptWeight->GetYaxis()->SetTitle("pT weight (Data/MC, thr=0 + pol4, frozen beyond fitMax)");
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
    Double_t rawEntries = h_pt_mc_thr0_noW->GetBinContent(b) * h_pt_mc_thr0_noW->GetXaxis()->GetBinWidth(b);
    Double_t ratioVal = h_pt_mc->GetBinContent(b) > 0 ? h_ratio->GetBinContent(b) : -1;
    Double_t weightVal = h_ptWeight->GetBinContent(b);
    printf("[%6.1f,%7.1f)  %12.5g   %s (raw_thr0_entries=%.0f, measured_ratio=%.4g)\n",
           lo, hi, weightVal,
           (hi <= fitMax) ? "pol4 fit" : "FROZEN at fitMax",
           rawEntries, ratioVal);
  }

  // --- overlay canvas (display only; fit/weights above still use the full pT range) ---
  TFile *fout = new TFile(outfile_name, "RECREATE");

  // Only the plotted x-range is restricted; normalization and fit are unchanged.
  const Double_t plotXmin = h_pt_data->GetXaxis()->GetXmin();
  const Double_t plotXmax = 10.0;

  // Draw clones so the histograms written to the output file keep their full axis range.
  TH1D *h_data_plot = (TH1D *)h_pt_data->Clone("h_data_plot");
  TH1D *h_mc_plot = (TH1D *)h_pt_mc->Clone("h_mc_plot");
  TH1D *h_ratio_plot = (TH1D *)h_ratio->Clone("h_ratio_plot");
  for (TH1D *h : {h_data_plot, h_mc_plot, h_ratio_plot})
  {
    h->SetStats(0);
    h->SetTitle("");
    h->GetXaxis()->SetRangeUser(plotXmin, plotXmax - 1e-6);
  }

  // y-ranges from the bins actually shown
  Double_t yMax = 0, yMin = 1e30, rMax = 0, rMin = 1e30;
  for (Int_t b = 1; b <= nBins; b++)
  {
    if (h_pt_data->GetXaxis()->GetBinUpEdge(b) > plotXmax + 1e-6)
      continue;
    for (TH1D *h : {h_pt_data, h_pt_mc})
    {
      if (h->GetBinContent(b) <= 0)
        continue;
      yMax = TMath::Max(yMax, h->GetBinContent(b) + h->GetBinError(b));
      yMin = TMath::Min(yMin, h->GetBinContent(b) - h->GetBinError(b));
    }
    if (h_pt_mc->GetBinContent(b) > 0)
    {
      rMax = TMath::Max(rMax, h_ratio->GetBinContent(b) + h_ratio->GetBinError(b));
      rMin = TMath::Min(rMin, h_ratio->GetBinContent(b) - h_ratio->GetBinError(b));
    }
  }

  TCanvas *c1 = new TCanvas("cOverlay", "cOverlay", 800, 800);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetLeftMargin(0.14);
  pad1->SetRightMargin(0.04);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.02);
  pad1->SetTicks(1, 1);
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();

  h_data_plot->SetMarkerStyle(20);
  h_data_plot->SetMarkerSize(1.1);
  h_data_plot->SetMarkerColor(kBlack);
  h_data_plot->SetLineColor(kBlack);
  h_data_plot->GetYaxis()->SetTitle("Normalized dN/dp_{T} (GeV/c)^{-1}");
  h_data_plot->GetYaxis()->SetTitleSize(0.05);
  h_data_plot->GetYaxis()->SetTitleOffset(1.3);
  h_data_plot->GetYaxis()->SetLabelSize(0.045);
  h_data_plot->GetXaxis()->SetLabelSize(0);
  h_data_plot->GetXaxis()->SetTitleSize(0);
  h_data_plot->SetMinimum(0.5 * yMin);
  h_data_plot->SetMaximum(3.0 * yMax);
  h_data_plot->Draw("E1");

  h_mc_plot->SetLineColor(kRed + 1);
  h_mc_plot->SetLineWidth(2);
  h_mc_plot->SetMarkerStyle(0);
  h_mc_plot->Draw("HIST SAME");
  h_data_plot->Draw("E1 SAME"); // keep data points on top

  TLegend *leg = new TLegend(0.60, 0.70, 0.92, 0.88);
  leg->AddEntry(h_data_plot, "Data", "lep");
  leg->AddEntry(h_mc_plot, "MC", "l");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->Draw();

  c1->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
  pad2->SetLeftMargin(0.14);
  pad2->SetRightMargin(0.04);
  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.35);
  pad2->SetTicks(1, 1);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  // pad2 is 0.3/0.7 the height of pad1 -> scale text sizes so they match on the canvas
  const Double_t sizeScale = 0.7 / 0.3;
  h_ratio_plot->SetMarkerStyle(20);
  h_ratio_plot->SetMarkerSize(1.1);
  h_ratio_plot->SetMarkerColor(kBlack);
  h_ratio_plot->SetLineColor(kBlack);
  h_ratio_plot->GetYaxis()->SetTitle("Data / MC");
  h_ratio_plot->GetYaxis()->SetTitleSize(0.05 * sizeScale);
  h_ratio_plot->GetYaxis()->SetTitleOffset(1.3 / sizeScale);
  h_ratio_plot->GetYaxis()->SetLabelSize(0.045 * sizeScale);
  h_ratio_plot->GetYaxis()->SetNdivisions(505);
  h_ratio_plot->GetYaxis()->CenterTitle();
  h_ratio_plot->GetYaxis()->SetRangeUser(TMath::Max(0.0, rMin - 0.1), rMax + 0.1);
  h_ratio_plot->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_ratio_plot->GetXaxis()->SetTitleSize(0.05 * sizeScale);
  h_ratio_plot->GetXaxis()->SetTitleOffset(1.0);
  h_ratio_plot->GetXaxis()->SetLabelSize(0.045 * sizeScale);
  h_ratio_plot->GetXaxis()->SetTickLength(0.03 * sizeScale);
  h_ratio_plot->Draw("E1");

  f_ptWeight->SetLineColor(kBlue + 1);
  f_ptWeight->SetLineStyle(2);
  f_ptWeight->SetLineWidth(2);
  f_ptWeight->SetNpx(1000);
  f_ptWeight->SetRange(plotXmin, plotXmax);
  f_ptWeight->Draw("SAME");

  TLine *line1 = new TLine(plotXmin, 1, plotXmax, 1);
  line1->SetLineStyle(3);
  line1->Draw();
  h_ratio_plot->Draw("E1 SAME");

  c1->Write();
  c1->Print("PtWeight_thr0only_Overlay.png");

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
  cout << ">>> Overlay PNG also saved to PtWeight_thr0only_Overlay.png\n";
}
