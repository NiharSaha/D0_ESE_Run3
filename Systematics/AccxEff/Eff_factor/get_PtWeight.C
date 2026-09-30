// Step 3 (ratio + fit) of the pT-weighting derivation for AccxEff.
//
// Reads: (a) h_pt_data_raw from get_DataPtSpectrum_fromMass.C's output (step 2 -- the Data
// mass-fit pT spectrum, split out into its own macro since it's the slow part and only needs
// rerunning when the Data mass histograms actually change), and (b) the per-Dpt-sample MC pT
// histograms hadd'ed from get_MassSpectra_PtWeight.C (step 1). Combines the 9 MC samples,
// builds the Data/MC ratio, and fits it with a polynomial; that fit is the pT-weighting factor
// to feed into the AccxEff calculation.

#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <iostream>
#include <Math/MinimizerOptions.h>
#include <cmath>

using namespace std;

// --- fine pT binning, must match get_MassSpectra_PtWeight.C / get_DataPtSpectrum_fromMass.C ---
const Int_t N_PTBIN_eff = 23;
Double_t ptbinning_[N_PTBIN_eff + 1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 70.0, 100.0};

// --- the 9 nested "Dpt > X GeV" MC samples (same set as getEfficiency.C / Analysis_bin_BDT.h),
// combined here per pT bin instead of get_MassSpectra_PtWeight.C's old global eligibility-sum
// getStitchWeight() scheme. That scheme let a sparse high-threshold sample's full luminosity
// distort a bin the instant it crossed its own threshold, regardless of whether it actually had
// any statistical power there -- confirmed this session as the cause of the sharp, unphysical
// break in the MC spectrum after 8 GeV. get_MassSpectra_PtWeight.C now writes one histogram per
// sample (h_pt_mc_raw_thr<N>, unstitched, weight=genWeight only); combineDptSamples() below
// combines them here via per-bin inverse-variance weighting instead. ---
const int N_DPT_SAMPLES = 9;
Double_t Dpt_threshold[N_DPT_SAMPLES] = {0, 1, 8, 10, 20, 30, 40, 60, 80};

// Dpt_xsecAfterFilter_pb[i]: the full generator-level cross section AFTER all gen filters
// (parton filter + pTHatMin, where applicable, + D0 MomMinPt/rapidity), in pb, taken directly
// from GenXsecAnalyzer's "After filter: final cross section" line for that sample's own
// gen-fragment. This replaces the previous Dpt_xsec_mb[]/Dpt_filterEff[] two-number split,
// which is where both the thr=30 (typo'd filterEff, since fixed by this re-measurement) and
// thr=8 (fragment byte-identical to thr=1, confirmed and handled by the merge below) bugs were
// hiding -- a single directly-measured number per sample removes that whole failure mode.
//
// thr=0,1: NOT yet re-measured via GenXsecAnalyzer (still the original McM-derived xsec_mb *
// filterEff, converted to pb) -- deferred since the pT<8 GeV region already matches Data.
// thr=8: merged into thr=1 (see the Dpt8 fragment bug, confirmed 2026-09-23); entry unused
// (h_pt_mc_sample[2] is zeroed below before combineDptSamples() runs).
// thr=10,20,40,60,80: measured via GenXsecAnalyzer (user, 2026-09-23).
// thr=30: re-measured with a larger test batch (user, 2026-09-24): 200000 events, 11 survivors,
// filterEff = 5.500e-05 +- 1.658e-05 (~30% relative uncertainty) -> after-filter sigma =
// 9.957e4 +- 3.002e4 pb. Supersedes the original 100000-event/1-survivor run (100% relative
// uncertainty, 1.809e4 pb), which was flagged unreliable and is why this was redone.
Double_t Dpt_xsecAfterFilter_pb[N_DPT_SAMPLES] = {
    4.63e3 * 7.43e-3 * 1e9, // thr=0  (unchanged, unverified)
    3.67e3 * 5.9e-3 * 1e9,  // thr=1  (unchanged, unverified)
    0.0,                    // thr=8  (merged into thr=1, unused)
    9.984e6,                // thr=10 (GenXsecAnalyzer)
    6.873e5,                // thr=20 (GenXsecAnalyzer)
    9.957e4,                // thr=30 (GenXsecAnalyzer, 200k-event rerun, ~30% rel. unc.)
    1.847e4,                // thr=40 (GenXsecAnalyzer)
    5.463e3,                // thr=60 (GenXsecAnalyzer)
    7.015e2,                // thr=80 (GenXsecAnalyzer)
};

// Nevts_i: actual raw generated events in the production ntuples this analysis processes --
// NOT the small test-batch sizes (50000/100000) used above just to measure filterEff/xsec via
// GenXsecAnalyzer. Unchanged from before this round of fixes; flag if these don't match the
// real production sample sizes.
Double_t Dpt_Nevts[N_DPT_SAMPLES] = {1.0e6, 2.0e6, 1.5e5, 1.0e5, 5.0e4, 2.0e4, 2.0e4, 2.0e4, 2.0e4};

// Combine the 9 per-sample MC histograms into one, per bin, via inverse-variance weighting --
// but ONLY among samples that are "eligible" for that bin (Dpt_threshold[i] <= bin's low edge).
//
// Why the eligibility gate is required, not optional (found by inspecting real Sept23_v5 MC
// output): fillMCPtSpectrum() fills every matchGEN/prompt/non-swap RECONSTRUCTED candidate in an
// event, not just the one leading D0 whose generated pT satisfied the sample's own Dpt-threshold
// filter. So a Dpt80 file (leading generated D0 pT>80) still legitimately contains occasional
// OTHER, unrelated reconstructed D0 candidates at pT=2-3 GeV (e.g. a second charm hadron in the
// same rare hard-scatter event). Those low-pT entries are real, but they estimate a *different,
// much rarer, conditional* rate ("D0 in this bin AND a companion D0 with pT>80 elsewhere"), not
// the marginal inclusive rate that Dpt0/Dpt1 measure. Naive inverse-variance weighting across ALL
// 9 samples cannot tell these apart: because Dpt80's Leff is ~5.5e6x larger than Dpt0's, even a
// handful of these stray/unrelated entries get a deceptively tiny *absolute* error after /Leff,
// giving them 15-20 orders of magnitude more weight than a well-measured low-threshold sample --
// confirmed directly against real v5 output, where a 5-entry Dpt80 stray in bin [2,2.5) got
// weight ~5e21 vs Dpt0's 5430-entry weight of ~1.6e5. Increasing a minimum-entry-count is not a
// fix either: the dominant weight ratio is set by (Leff_i/Leff_j)^2, ~10-20 orders of magnitude
// between the sparsest and densest samples, so no realistic entry count closes that gap. Only
// restricting each sample to bins its own filter can actually speak for (threshold <= bin's pT)
// removes the contamination at its source, while inverse-variance weighting *among the eligible
// samples* still gives smooth, non-discontinuous transitions at each threshold crossing (unlike
// the original getStitchWeight()'s blunt eligibility-luminosity-sum, which jumped discretely).
TH1D *combineDptSamples(TH1D *h_pt_mc_sample[N_DPT_SAMPLES])
{
  Double_t Dpt_Leff[N_DPT_SAMPLES];
  for (int i = 0; i < N_DPT_SAMPLES; i++)
    Dpt_Leff[i] = Dpt_Nevts[i] / Dpt_xsecAfterFilter_pb[i];

  TH1D *h_combined = (TH1D *)h_pt_mc_sample[0]->Clone("h_pt_mc_raw_combined");
  h_combined->Reset();

  for (int ipt = 0; ipt < N_PTBIN_eff; ipt++)
  {
    const Double_t pt_lo = h_pt_mc_sample[0]->GetXaxis()->GetBinLowEdge(ipt + 1);

    Double_t sumWInv2 = 0.0;   // sum of 1/sigma_i^2
    Double_t sumRateWInv2 = 0.0; // sum of rate_i/sigma_i^2
    for (int i = 0; i < N_DPT_SAMPLES; i++)
    {
      if (Dpt_threshold[i] > pt_lo)
        continue; // sample's own generator filter can't speak for this pT region -- excluded

      const Double_t content = h_pt_mc_sample[i]->GetBinContent(ipt + 1);
      const Double_t error = h_pt_mc_sample[i]->GetBinError(ipt + 1);
      if (content <= 0.0 || error <= 0.0)
        continue; // no (or no usable) statistical power from this sample in this bin

      const Double_t rate = content / Dpt_Leff[i];
      const Double_t rateErr = error / Dpt_Leff[i];
      const Double_t wInv2 = 1.0 / (rateErr * rateErr);

      sumWInv2 += wInv2;
      sumRateWInv2 += rate * wInv2;
    }

    if (sumWInv2 <= 0.0)
      continue; // no eligible sample had any statistical power here; bin stays at 0

    h_combined->SetBinContent(ipt + 1, sumRateWInv2 / sumWInv2);
    h_combined->SetBinError(ipt + 1, std::sqrt(1.0 / sumWInv2));
  }

  return h_combined;
}

void get_PtWeight(TString data_pt_infile = "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/AccxEff/DataPtSpectrum_out.root",
                            TString mc_infile = "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/PtWeight_Sept24_v2/ROOT/MCSpectrum_out_combined.root",
                            TString outfile_name = "PtWeight_out.root",
                            Int_t poly_degree = 4)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TFile *fin = TFile::Open(data_pt_infile);
  TFile *finMC = TFile::Open(mc_infile);
  if (!fin || fin->IsZombie() || !finMC || finMC->IsZombie())
  {
    cerr << "Error opening input files (" << data_pt_infile << ", " << mc_infile << ")" << endl;
    return;
  }

  TFile *fout = new TFile(outfile_name, "RECREATE");

  TH1D *h_pt_data_raw = (TH1D *)fin->Get("h_pt_data_raw");
  if (!h_pt_data_raw)
  {
    cerr << "Missing h_pt_data_raw in " << data_pt_infile
         << " -- run get_DataPtSpectrum_fromMass.C first." << endl;
    return;
  }
  h_pt_data_raw = (TH1D *)h_pt_data_raw->Clone("h_pt_data_raw");

  TH1D *h_pt_data = (TH1D *)h_pt_data_raw->Clone("h_pt_data");
  TH1D *h_pt_mc = new TH1D("h_pt_mc", ";p_{T} (GeV);dN/dp_{T}", N_PTBIN_eff, ptbinning_);

  TH1D *h_pt_mc_sample[N_DPT_SAMPLES];
  for (int i = 0; i < N_DPT_SAMPLES; i++)
  {
    h_pt_mc_sample[i] = (TH1D *)finMC->Get(Form("h_pt_mc_raw_thr%d", (int)Dpt_threshold[i]));
    if (!h_pt_mc_sample[i])
    {
      cerr << "Missing h_pt_mc_raw_thr" << (int)Dpt_threshold[i] << " in " << mc_infile << endl;
      return;
    }
  }
  // Dpt8's GEN fragment was confirmed (2026-09-23) to be byte-identical to Dpt1's (MomMinPt=1.0,
  // not 8.0) -- a production mislabeling, not a genuine higher-pT filter. Its recorded xsec/
  // filterEff were booked assuming the (nonexistent) pT>8 cut, so they cannot be used to
  // normalize its raw events. Treat Dpt8's statistics as additional Dpt1-pool events instead:
  // merge its histogram into thr=1's, fold its Nevts into Dpt1's Nevts, and zero it out so it
  // can no longer contribute as an independent (wrongly-normalized) sample downstream.
  h_pt_mc_sample[1]->Add(h_pt_mc_sample[2]);
  Dpt_Nevts[1] += Dpt_Nevts[2];
  h_pt_mc_sample[2]->Reset();

  TH1D *h_pt_mc_raw = combineDptSamples(h_pt_mc_sample);

  // h_pt_mc_raw is built from get_MassSpectra_PtWeight.C's per-sample histograms, which are
  // already divided by bin width at the source (Scale(1.0, "width") in runMC()) -- do NOT
  // divide by bin width again here, or the MC spectrum ends up double-width-normalized.
  for (int ipt = 0; ipt < N_PTBIN_eff; ipt++)
  {
    h_pt_mc->SetBinContent(ipt + 1, h_pt_mc_raw->GetBinContent(ipt + 1));
    h_pt_mc->SetBinError(ipt + 1, h_pt_mc_raw->GetBinError(ipt + 1));
  }

  // Keep the raw (un-normalized) yield spectra too -- normalization below is only
  // for the ratio/weight fit and the presentation overlay.
  TH1D *h_pt_mc_raw_out = (TH1D *)h_pt_mc->Clone("h_pt_mc_raw");

  // --- normalize both spectra to equal (unit) area before the ratio/weight fit ---
  const Double_t data_integral = h_pt_data->Integral("width");
  const Double_t mc_integral = h_pt_mc->Integral("width");
  if (data_integral > 0)
    h_pt_data->Scale(1.0 / data_integral);
  if (mc_integral > 0)
    h_pt_mc->Scale(1.0 / mc_integral);

  // --- ratio + polynomial fit ---
  // The old MC_TRUSTED_PT_MAX range restriction (only fitting pT<8 GeV, since above that the
  // global eligibility-sum stitching was unphysical) is no longer needed now that
  // combineDptSamples() properly weights each Dpt sample by its own per-bin statistical power --
  // the full range should be trustworthy.
  TH1D *h_ratio = (TH1D *)h_pt_data->Clone("h_ratio_data_over_mc");
  h_ratio->SetTitle(";p_{T} (GeV);Data / MC");
  h_ratio->Divide(h_pt_mc);

  TF1 *ptWeightFit = new TF1("ptWeightFit", Form("pol%d", poly_degree), ptbinning_[0], ptbinning_[N_PTBIN_eff]);
  TFitResultPtr ratioFitRes = h_ratio->Fit(ptWeightFit, "S Q", "", ptbinning_[0], ptbinning_[N_PTBIN_eff]);
  ptWeightFit->SetLineColor(kRed);
  ptWeightFit->SetLineWidth(2);

  cout << "\n>>> pT-weight polynomial fit (pol" << poly_degree << "), status = "
       << (ratioFitRes.Get() ? ratioFitRes->Status() : -1) << endl;
  for (int p = 0; p <= poly_degree; p++)
    cout << "    par[" << p << "] = " << ptWeightFit->GetParameter(p)
         << " +- " << ptWeightFit->GetParError(p) << endl;

  // --- overlay (Data vs MC, normalized) + ratio canvas, for presentation ---
  TCanvas *cOverlay = new TCanvas("cOverlay_DataMC_Pt", "Data vs MC p_{T} spectra", 800, 800);
  TPad *padTop = new TPad("padTop", "padTop", 0, 0.3, 1, 1.0);
  TPad *padBot = new TPad("padBot", "padBot", 0, 0.0, 1, 0.3);
  padTop->SetBottomMargin(0.02);
  padTop->SetLogy();
  padTop->SetLogx();
  padBot->SetTopMargin(0.02);
  padBot->SetBottomMargin(0.35);
  padBot->SetLogx();
  padTop->Draw();
  padBot->Draw();

  padTop->cd();
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

  TLegend *leg = new TLegend(0.65, 0.72, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_pt_data, "Data", "lep");
  leg->AddEntry(h_pt_mc, "MC", "l");
  leg->Draw();

  padBot->cd();
  h_ratio->SetTitle("");
  h_ratio->GetYaxis()->SetTitle("Data / MC");
  h_ratio->GetXaxis()->SetTitle("p_{T} (GeV)");
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->GetYaxis()->SetTitleSize(0.1);
  h_ratio->GetYaxis()->SetTitleOffset(0.45);
  h_ratio->GetYaxis()->SetLabelSize(0.08);
  h_ratio->GetXaxis()->SetTitleSize(0.1);
  h_ratio->GetXaxis()->SetLabelSize(0.08);
  h_ratio->SetMarkerStyle(20);
  h_ratio->Draw("E1");
  ptWeightFit->Draw("SAME");

  TLine *unity = new TLine(ptbinning_[0], 1.0, ptbinning_[N_PTBIN_eff], 1.0);
  unity->SetLineStyle(2);
  unity->SetLineColor(kGray + 2);
  unity->Draw("SAME");

  cOverlay->cd();

  fout->cd();
  h_pt_data->Write();
  h_pt_mc->Write();
  h_pt_data_raw->Write();
  h_pt_mc_raw_out->Write();
  h_ratio->Write();
  ptWeightFit->Write();
  cOverlay->Write();

  TString overlayPng = outfile_name;
  overlayPng.ReplaceAll(".root", "_DataMC_Overlay.png");
  cOverlay->SaveAs(overlayPng);

  fout->Write();
  fout->Close();

  cout << "\n>>> Done! pT-weight output saved to " << outfile_name << endl;
}
