#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TDirectory.h>
#include <TPaveText.h>
#include <TMath.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <TParameter.h>
#include <Math/PdfFuncMathCore.h>
#include <Math/MinimizerOptions.h>
#include <cmath>

using namespace std;

// Pivot for the polynomial background. Expressing the background as
// b0 + b1*(m - MASS_PIVOT) + b2*(m - MASS_PIVOT)^2 (instead of powers of m)
// de-correlates b0/b1/b2 over the narrow fit window and makes b0 the actual
// background level (counts/bin) at the centre of the window.
const Double_t MASS_PIVOT = 1.87;

Double_t TotalFuncModel(Double_t *x, Double_t *par)
{
  const double xx = x[0];

  const double signalYield = par[0];
  const double signalMean = par[1];
  const double sigma1 = par[2];
  const double sigma2 = par[3];
  const double narrowFrac = par[4];
  const double signalSplit = par[5];
  const double widthScale = par[6];
  const double sigma3 = par[7];

  const double b0 = par[8];
  const double b1 = par[9];
  const double b2 = par[10];

  const double tailFrac = par[11]; // Crystal Ball tail size relative to signalYield
  const double tail1Scale = par[12];
  const double tail2Scale = par[13];

  const double sigma1_eff = sigma1 * (1.0 + widthScale);
  const double sigma2_eff = sigma2 * (1.0 + widthScale);
  const double sigma3_eff = sigma3 * (1.0 + widthScale);

  const double gaus1 = TMath::Gaus(xx, signalMean, sigma1_eff, kTRUE) / (TMath::Sqrt(2.0 * TMath::Pi()) * sigma1_eff);
  const double gaus2 = TMath::Gaus(xx, signalMean, sigma2_eff, kTRUE) / (TMath::Sqrt(2.0 * TMath::Pi()) * sigma2_eff);
  const double gaus3 = TMath::Gaus(xx, signalMean, sigma3_eff, kTRUE) / (TMath::Sqrt(2.0 * TMath::Pi()) * sigma3_eff);

  const double signalPart = signalYield * ((signalSplit * (narrowFrac * gaus1 + (1.0 - narrowFrac) * gaus2)) + (1.0 - signalSplit) * gaus3);

  const double cb1 = ROOT::Math::crystalball_function(xx, 2.2, 17, 0.0267 * (1.0 + widthScale), 1.96 * (1.0 + tail1Scale));
  const double cb2 = ROOT::Math::crystalball_function(xx, 0.34, 5, 0.0146 * (1.0 + widthScale), 1.7734 * (1.0 + tail2Scale));

  const double tb = xx - MASS_PIVOT;
  const double background = b0 + b1 * tb + b2 * tb * tb;

  // The Crystal Ball tails describe the K-pi swapped D0 and therefore scale
  // with the true D0 yield. Tying them to signalYield lets the tail shrink
  // automatically in the low-S/B (high-DCA) bins instead of staying frozen at
  // the inclusive value.
  return signalPart + background + signalYield * tailFrac * (cb1 + 4.0 * cb2);
}

Double_t doubleGauss_norm(Double_t *x, Double_t *par)
{
  const double xx = x[0];
  const double signalYield = par[0];
  const double mean = par[1];
  const double sigma1 = par[2];
  const double sigma2 = par[3];
  const double narrowFrac = par[4];
  const double widthScale = par[5];

  const double sigma1_eff = sigma1 * (1.0 + widthScale);
  const double sigma2_eff = sigma2 * (1.0 + widthScale);

  const double gaus1 = TMath::Gaus(xx, mean, sigma1_eff, kTRUE) / (TMath::Sqrt(2.0 * TMath::Pi()) * sigma1_eff);
  const double gaus2 = TMath::Gaus(xx, mean, sigma2_eff, kTRUE) / (TMath::Sqrt(2.0 * TMath::Pi()) * sigma2_eff);

  return signalYield * (narrowFrac * gaus1 + (1.0 - narrowFrac) * gaus2);
}

const Double_t D0_MASS = 1.86484;
const Double_t MASS_FIT_MIN = 1.76;
const Double_t MASS_FIT_MAX = 1.98;

// --- Per-DCA-bin background strategy ------------------------------------------
// The background treatment is chosen from a sideband estimate of the D0 signal
// significance  S / sqrt(S+B)  in a fixed +-SIG_WINDOW mass region (identical
// calculation for every bin), NOT from a hardcoded bin index. S and B are the
// sideband-subtracted signal and the background under the peak.
//   Tier 0 : significance >= BKG_TIER0_MIN_SIGNIF  (and signalFrac not tiny)
//            -> quadratic background (b0, b1, b2 free)
//   Tier 1 : significance >= BKG_TIER1_MIN_SIGNIF
//            -> linear background (b0, b1 free, b2 = 0)
//   Tier 2 : significance <  BKG_TIER1_MIN_SIGNIF
//            -> linear background with the SLOPE fixed from the sideband points
//               (b0 and yield free, b1 from the two sidebands, b2 = 0)
//   Tier 3 : escalation-only fallback
//            -> background fully fixed from the sidebands, only the yield floats
// Rationale: below a few-sigma sideband significance the polynomial background
// shape and the yield are not simultaneously constrainable from that bin alone,
// so the shape is taken from the sideband regions of the same histogram. In
// practice this selects the last 2-3 DCA bins but adapts per (cent,pt). A bin
// also escalates one tier whenever the fit-quality gate below fails.
const Double_t SIG_WINDOW = 0.030;             // +-30 MeV signal region (~2 sigma)
const Double_t BKG_TIER0_MIN_SIGNIF = 25.0;
const Double_t BKG_TIER1_MIN_SIGNIF = 8.0;
const Double_t BKG_TIER0_MIN_SIGFRAC = 0.01;   // extra guard on the free quadratic
// Only genuine convergence failures escalate the tier. With O(10^6) counts per
// bin a fixed-shape fit legitimately reaches chi2/ndf ~ 10, so the ceiling here
// is deliberately loose and only catches pathological (chi2/ndf ~ 10^3) fits.
const Double_t FIT_CHI2NDF_CEILING = 60.0;

void get_Data_DCA_fromMass()
{
  TFile *file = new TFile("Data_Mass_DCA_out_combined_Sept9.root");
  TFile *fileMC = TFile::Open("D0_MCtemplate_out_combined.root");
  const char *date = "Sept9";
  TFile *f_out = new TFile(Form("HIST_DATA_DCA_%s.root", date), "recreate");
  bool Massplot = false;

  Int_t pt_min = 0;
  Int_t pt_max = 9;

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  // Give the minimiser more headroom for the low-S/B DCA bins.
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(20000);

  const int N_CENTBIN = 5;
  const char *label_centbin[N_CENTBIN] = {"cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"};

  const Int_t PTBIN = 9;
  const char *pt_label[PTBIN] = {"pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8", "pT8to10", "pT10to15", "pT15to30", "pT30to100"};
  const char *mc_pt_label[PTBIN] = {"pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8", "pT8to10", "pT10to15", "pT15to20", "pT20to40"};

  // Two DCA binnings, matching the mass-distribution producer: pT bins whose
  // lower edge is < 5 GeV (pT2to3, pT3to4, pT4to5) use a coarser high-DCA tail
  // (dca_bins_lo) so the low-statistics high-DCA bins fit better; every other pT
  // bin uses the finer dca_bins_hi. Both binnings share the same first/last edge
  // and identical edges up to 0.0214 cm.
  const Int_t dca_bins_hi = 18;
  Double_t DCA_hi[dca_bins_hi + 1] = {0, 0.0005, 0.0011, 0.0014, 0.0024, 0.0029, 0.0039, 0.0045, 0.0059, 0.0067, 0.0085, 0.01179, 0.016, 0.0214, 0.028, 0.0366, 0.0475, 0.079, 0.135};
  const char *label_dca_hi[dca_bins_hi] = {"dca_0to0p0005", "dca_0p0005to0p0011", "dca_0p0011to0p0014", "dca_0p0014to0p0024", "dca_0p0024to0p0029", "dca_0p0029to0p0039",
                                           "dca_0p0039to0p0045", "dca_0p0045to0p0059", "dca_0p0059to0p0067", "dca_0p0067to0p0085", "dca_0p0085to0p01179", "dca_0p01179to0p016",
                                           "dca_0p016to0p0214", "dca_0p0214to0p028", "dca_0p028to0p0366", "dca_0p0366to0p0475", "dca_0p0475to0p079", "dca_0p079to0p135"};

  const Int_t dca_bins_lo = 15;
  Double_t DCA_lo[dca_bins_lo + 1] = {0, 0.0005, 0.0011, 0.0014, 0.0024, 0.0029, 0.0039, 0.0045, 0.0059, 0.0067, 0.0085, 0.01179, 0.016, 0.0214, 0.0366, 0.135};
  const char *label_dca_lo[dca_bins_lo] = {"dca_0to0p0005", "dca_0p0005to0p0011", "dca_0p0011to0p0014", "dca_0p0014to0p0024", "dca_0p0024to0p0029", "dca_0p0029to0p0039",
                                           "dca_0p0039to0p0045", "dca_0p0045to0p0059", "dca_0p0059to0p0067", "dca_0p0067to0p0085", "dca_0p0085to0p01179", "dca_0p01179to0p016",
                                           "dca_0p016to0p0214", "dca_0p0214to0p0366", "dca_0p0366to0p135"};

  const Int_t dca_bins_max = dca_bins_hi; // static-array dimension (largest of the two)

  // numeric pT lower edges, in the same order as pt_label above
  const Double_t pt_low_edge[PTBIN] = {2, 3, 4, 5, 6, 8, 10, 15, 30};
  auto useLowPtBins = [&](int ipt) { return pt_low_edge[ipt] < 5.0; };
  auto dcaNbins = [&](int ipt) { return useLowPtBins(ipt) ? dca_bins_lo : dca_bins_hi; };
  auto dcaEdges = [&](int ipt) -> Double_t * { return useLowPtBins(ipt) ? DCA_lo : DCA_hi; };
  auto dcaLabel = [&](int ipt, int n) { return useLowPtBins(ipt) ? label_dca_lo[n] : label_dca_hi[n]; };

  TCanvas *canvas[N_CENTBIN][PTBIN];
  TDirectory *dir[N_CENTBIN][PTBIN];
  TH1F *h_dmass[N_CENTBIN][PTBIN][dca_bins_max];
  TH1D *hdca[N_CENTBIN][PTBIN];

  TCanvas *c[N_CENTBIN][PTBIN][dca_bins_max];
  TF1 *total_func[N_CENTBIN][PTBIN][dca_bins_max];
  TF1 *doubGauss[N_CENTBIN][PTBIN][dca_bins_max];

  Double_t BinWidth[N_CENTBIN][PTBIN][dca_bins_max];
  Double_t Yield[N_CENTBIN][PTBIN][dca_bins_max];
  Double_t Yield_error[N_CENTBIN][PTBIN][dca_bins_max];

  for (int icent = 0; icent < N_CENTBIN; icent++)
  {
    for (int ipt = 0; ipt < PTBIN; ipt++)
    {
      hdca[icent][ipt] = new TH1D(Form("hdca_%s_%s", label_centbin[icent], pt_label[ipt]),
                                  Form("hdca_%s_%s", label_centbin[icent], pt_label[ipt]), dcaNbins(ipt), dcaEdges(ipt));
      hdca[icent][ipt]->SetDirectory(nullptr); // avoid auto-attaching to f_out; written explicitly below
    }
  }

  

  for (int icent = 0; icent < N_CENTBIN; icent++)
  {
    for (int ipt = pt_min; ipt < pt_max; ipt++)
    {
      dir[icent][ipt] = (TDirectory *)file->Get(Form("out_Mass_DCA_%s_%s", label_centbin[icent], pt_label[ipt]));
      if (!dir[icent][ipt])
        continue;

      // canvas[icent][ipt] = new TCanvas(Form("c%d%d", icent, ipt), Form("c%d%d", icent, ipt), 600, 600);

      f_out->cd();
      // f_out->mkdir(Form("MassFits_%s_%s", label_centbin[icent], pt_label[ipt]));
      TDirectory *massDir = f_out->mkdir(Form("MassFits_%s_%s", label_centbin[icent], pt_label[ipt]));
      //---------

      Double_t signalFrac = 0.5;
      Double_t signalMean = D0_MASS;
      Double_t signalSigma1 = 0.008;
      Double_t signalSigma2 = 0.018;
      Double_t signalScale = 0.0;

      TH1F *h_mass_inclusive = (TH1F *)dir[icent][ipt]->Get(Form("dmass_inclusive_%s_%s", label_centbin[icent], pt_label[ipt]));

      if (!h_mass_inclusive)
      {
        cout << "Missing inclusive histogram for "
             << label_centbin[icent] << " " << pt_label[ipt] << endl;
        continue;
      }

      TF1 *inclusiveFit = new TF1(Form("inclusiveFit_%d_%d", icent, ipt), TotalFuncModel, MASS_FIT_MIN, MASS_FIT_MAX, 14);

      TH1F *h_mc_match_signal = (TH1F *)fileMC->Get(Form("hMass_Signal_%s", mc_pt_label[ipt]));
      TH1F *h_mc_match_all = (TH1F *)fileMC->Get(Form("hMass_Signal_Plus_Swap_%s", mc_pt_label[ipt]));
      if (!h_mc_match_signal || !h_mc_match_all)
      {
        cerr << "Missing MC templates for " << mc_pt_label[ipt] << endl;
        delete inclusiveFit;
        continue;
      }

      Double_t inclusiveMax = h_mass_inclusive->GetMaximum();
      const Double_t inclusiveIntegral = h_mass_inclusive->Integral() * h_mass_inclusive->GetBinWidth(1);

      // Double_t inclusiveParameters[14] = {inclusiveIntegral, D0_MASS, 0.008, 0.018, 0.7, 0.5, 0.0, 0.012, inclusiveMax * 0.5, 0.0, 0.0, 20.0, 0.0, 0.0};
      // inclusiveFit->SetParameters(inclusiveParameters);

      inclusiveFit->SetParameter(0, 100.);    // overall normalization
      inclusiveFit->SetParameter(1, D0_MASS); // D0 mass mean
      inclusiveFit->SetParameter(2, 0.03);    // wide Gaussian sigma
      inclusiveFit->SetParameter(3, 0.005);   // narrow Gaussian sigma
      inclusiveFit->SetParameter(4, 0.1);     // narrow Gaussian fraction
      inclusiveFit->FixParameter(5, 1.0);     // signal fraction
      inclusiveFit->FixParameter(6, 0.0);     // global width smearing
      inclusiveFit->FixParameter(7, 0.1);     // swap Gaussian sigma
      inclusiveFit->FixParameter(8, 0.0);     // background constant
      inclusiveFit->FixParameter(9, 0.0);     // background linear term
      inclusiveFit->FixParameter(10, 0.0);    // background quadratic term
      inclusiveFit->FixParameter(11, 0.0);    // Crystal Ball tail fraction (rel. to signalYield)
      inclusiveFit->FixParameter(12, 0.0);    // Crystal Ball 1 mean shift
      inclusiveFit->FixParameter(13, 0.0);    // Crystal Ball 2 mean shift

      if (ipt < 5)
      {
        inclusiveFit->SetParLimits(2, 0.01, 0.5);
        inclusiveFit->SetParLimits(3, 0.001, 0.25);
      }
      else
      {
        inclusiveFit->SetParLimits(2, 0.005, 0.15);
        inclusiveFit->SetParLimits(3, 0.001, 0.08);
      }
      inclusiveFit->SetParLimits(4, 0.0, 1.0);
      inclusiveFit->SetParLimits(5, 0.0, 1.0);

      // Fit the gen-matched signal MC, exactly as in the flow analysis.
      inclusiveFit->FixParameter(1, D0_MASS);
      h_mc_match_signal->Fit(inclusiveFit, "q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      h_mc_match_signal->Fit(inclusiveFit, "q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      inclusiveFit->ReleaseParameter(1);
      h_mc_match_signal->Fit(inclusiveFit, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      h_mc_match_signal->Fit(inclusiveFit, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      h_mc_match_signal->Fit(inclusiveFit, "L q ", "", MASS_FIT_MIN, MASS_FIT_MAX);

      inclusiveFit->FixParameter(1, inclusiveFit->GetParameter(1));
      inclusiveFit->FixParameter(2, inclusiveFit->GetParameter(2));
      inclusiveFit->FixParameter(3, inclusiveFit->GetParameter(3));
      inclusiveFit->FixParameter(4, inclusiveFit->GetParameter(4));

      // Fit signal-plus-swap MC to determine the swap fraction and width.
      inclusiveFit->ReleaseParameter(5);
      inclusiveFit->ReleaseParameter(7);
      inclusiveFit->SetParameter(7, 0.1);
      for (int iteration = 0; iteration < 5; iteration++)
        h_mc_match_all->Fit(inclusiveFit, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);

      inclusiveFit->FixParameter(5, inclusiveFit->GetParameter(5));
      inclusiveFit->FixParameter(7, inclusiveFit->GetParameter(7));
      inclusiveFit->FixParameter(1, inclusiveFit->GetParameter(1));
      inclusiveFit->FixParameter(2, inclusiveFit->GetParameter(2));
      inclusiveFit->FixParameter(3, inclusiveFit->GetParameter(3));
      inclusiveFit->FixParameter(4, inclusiveFit->GetParameter(4));
      inclusiveFit->FixParameter(6, 0.0);

      // Release polynomial background parameters for the inclusive data fit.
      inclusiveFit->ReleaseParameter(8);
      inclusiveFit->ReleaseParameter(9);
      inclusiveFit->ReleaseParameter(10);

      // First fit with mean and width smearing fixed, then release both.
      inclusiveFit->FixParameter(1, D0_MASS);
      inclusiveFit->FixParameter(6, 0.0);
      h_mass_inclusive->Fit(inclusiveFit, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      h_mass_inclusive->Fit(inclusiveFit, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      inclusiveFit->ReleaseParameter(1);
      inclusiveFit->ReleaseParameter(6);
      h_mass_inclusive->Fit(inclusiveFit, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      h_mass_inclusive->Fit(inclusiveFit, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      h_mass_inclusive->Fit(inclusiveFit, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      inclusiveFit->FixParameter(6, inclusiveFit->GetParameter(6));
      inclusiveFit->FixParameter(1, inclusiveFit->GetParameter(1));

      // Final inclusive fit: release the mean, Crystal Ball terms, and refit.
      h_mass_inclusive->Fit(inclusiveFit, "q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      h_mass_inclusive->Fit(inclusiveFit, "q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      inclusiveFit->ReleaseParameter(1);
      h_mass_inclusive->Fit(inclusiveFit, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      h_mass_inclusive->Fit(inclusiveFit, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      h_mass_inclusive->Fit(inclusiveFit, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
      inclusiveFit->ReleaseParameter(11);
      inclusiveFit->SetParameter(11, 0.5); // par[11] is now a tail fraction of signalYield
      inclusiveFit->SetParLimits(11, 0.0, 20.0);
      h_mass_inclusive->Fit(inclusiveFit, "L q ", "", MASS_FIT_MIN, MASS_FIT_MAX);
      h_mass_inclusive->Fit(inclusiveFit, "L q ", "", MASS_FIT_MIN, MASS_FIT_MAX);
      inclusiveFit->SetParLimits(12, -0.03, 0.03);
      inclusiveFit->SetParLimits(13, -0.03, 0.03);
      inclusiveFit->ReleaseParameter(12);
      inclusiveFit->ReleaseParameter(13);
      h_mass_inclusive->Fit(inclusiveFit, "L q ", "", MASS_FIT_MIN, MASS_FIT_MAX);
      h_mass_inclusive->Fit(inclusiveFit, "L q ", "", MASS_FIT_MIN, MASS_FIT_MAX);
      h_mass_inclusive->Fit(inclusiveFit, "L q ", "", MASS_FIT_MIN, MASS_FIT_MAX);

      TFitResultPtr inclusiveResult = h_mass_inclusive->Fit(inclusiveFit, "L q S", "", MASS_FIT_MIN, MASS_FIT_MAX);

      if (inclusiveResult.Get() && inclusiveResult->Status() == 0)
      {
        signalFrac = inclusiveFit->GetParameter(4);
        signalMean = inclusiveFit->GetParameter(1);
        signalSigma1 = inclusiveFit->GetParameter(2);
        signalSigma2 = inclusiveFit->GetParameter(3);
      }
      else
      {
        cout << "Inclusive fit failed for "
             << label_centbin[icent] << " " << pt_label[ipt]
             << ". Using default signal parameters." << endl;
      }

      TH1F *h_inclusive_fit = (TH1F *)h_mass_inclusive->Clone(Form("dmass_inclusive_fit_%s_%s", label_centbin[icent], pt_label[ipt]));
      h_inclusive_fit->SetTitle(Form("Inclusive mass fit: %s %s", label_centbin[icent], pt_label[ipt]));
      h_inclusive_fit->GetXaxis()->SetRangeUser(MASS_FIT_MIN, MASS_FIT_MAX);
      h_inclusive_fit->SetDirectory(nullptr);

      TF1 *savedInclusiveFunc = (TF1 *)inclusiveFit->Clone(Form("inclusive_total_func_%s_%s", label_centbin[icent], pt_label[ipt]));
      savedInclusiveFunc->SetLineWidth(1);
      savedInclusiveFunc->SetLineColor(kRed);

      h_inclusive_fit->GetListOfFunctions()->Add(savedInclusiveFunc);

      const Int_t inclusiveNDF = inclusiveFit->GetNDF();
      const Double_t inclusiveChi2NDF = (inclusiveNDF > 0) ? inclusiveFit->GetChisquare() / inclusiveNDF : -1.0;
      const Double_t inclusiveYield = inclusiveFit->GetParameter(0) / h_mass_inclusive->GetBinWidth(1);   // Normalize yield to per bin width
      const Double_t inclusiveYieldErr = inclusiveFit->GetParError(0) / h_mass_inclusive->GetBinWidth(1); // Normalize yield error to per bin width

      TPaveText *inclusivePaveText = new TPaveText(0.15, 0.70, 0.48, 0.88, "brNDC");
      inclusivePaveText->SetBorderSize(0);
      inclusivePaveText->SetFillColor(0);
      inclusivePaveText->SetFillStyle(0);
      inclusivePaveText->SetTextFont(42);
      inclusivePaveText->SetTextSize(0.04);
      inclusivePaveText->SetTextAlign(12);
      const Double_t inclusiveSignificance = (inclusiveYieldErr > 0) ? inclusiveYield / inclusiveYieldErr : 0.0;
      inclusivePaveText->AddText(Form("Yield = %.1f #pm %.1f", inclusiveYield, inclusiveYieldErr));
      inclusivePaveText->AddText(Form("#chi^{2}/NDF = %.2f", inclusiveChi2NDF));
      inclusivePaveText->AddText(Form("Significance = %.2f", inclusiveSignificance));
      h_inclusive_fit->GetListOfFunctions()->Add(inclusivePaveText);

      /*TList *inclusiveInfo = new TList();
      inclusiveInfo->SetName("InclusiveFitInfo");
      inclusiveInfo->SetOwner(kTRUE);
      inclusiveInfo->Add(new TParameter<Double_t>("Yield", inclusiveFit->GetParameter(0)));
      inclusiveInfo->Add(new TParameter<Double_t>("Mean", inclusiveFit->GetParameter(1)));
      inclusiveInfo->Add(new TParameter<Double_t>("Sigma1", inclusiveFit->GetParameter(2)));
      inclusiveInfo->Add(new TParameter<Double_t>("Sigma2", inclusiveFit->GetParameter(3)));
      inclusiveInfo->Add(new TParameter<Double_t>("SignalFraction", inclusiveFit->GetParameter(4)));
      inclusiveInfo->Add(new TParameter<Double_t>("SignalSplit", inclusiveFit->GetParameter(5)));
      inclusiveInfo->Add(new TParameter<Double_t>("WidthScale", inclusiveFit->GetParameter(6)));
      inclusiveInfo->Add(new TParameter<Double_t>("Sigma3", inclusiveFit->GetParameter(7)));
      inclusiveInfo->Add(new TParameter<Double_t>("BackgroundConstant", inclusiveFit->GetParameter(8)));
      inclusiveInfo->Add(new TParameter<Double_t>("Chi2NDF", inclusiveChi2NDF));
      inclusiveInfo->Add(new TParameter<Int_t>("FitStatus", inclusiveResult.Get() ? inclusiveResult->Status() : -1));
      h_inclusive_fit->GetListOfFunctions()->Add(inclusiveInfo);*/

      massDir->cd();
      h_inclusive_fit->Write("", TObject::kOverwrite);
      delete h_inclusive_fit;

      //-----

      for (int idca = 0; idca < dcaNbins(ipt); idca++)
      {
        const char *dcaLab = dcaLabel(ipt, idca);
        if (Massplot)
          c[icent][ipt][idca] = new TCanvas(Form("c%d%d%d", icent, ipt, idca), "", 600, 600);

        h_dmass[icent][ipt][idca] = (TH1F *)dir[icent][ipt]->Get(Form("dmass_%s_%s_%s", label_centbin[icent], pt_label[ipt], dcaLab));
        if (!h_dmass[icent][ipt][idca])
          continue;

        total_func[icent][ipt][idca] = new TF1(Form("total_func_%d%d%d", icent, ipt, idca), TotalFuncModel, MASS_FIT_MIN, MASS_FIT_MAX, 14);
        doubGauss[icent][ipt][idca] = new TF1(Form("doubGauss_%d%d%d", icent, ipt, idca), doubleGauss_norm, MASS_FIT_MIN, MASS_FIT_MAX, 6);

        if (Massplot)
          c[icent][ipt][idca]->cd();

        TF1 *f = total_func[icent][ipt][idca];
        TH1F *h = h_dmass[icent][ipt][idca];

        // --- histogram quantities over the fit window ---
        const Int_t fitBinLo = h->FindBin(MASS_FIT_MIN + 1e-6);
        const Int_t fitBinHi = h->FindBin(MASS_FIT_MAX - 1e-6);
        const Double_t bin_w = h->GetBinWidth(1);
        const Double_t h_integral = h->Integral(fitBinLo, fitBinHi);
        const Double_t nEntries = h->GetEntries();
        const Double_t nIntegral = h_integral;

        // --- sideband estimate of the D0 signal significance (same recipe for every bin) ---
        const Double_t fitMean = inclusiveFit->GetParameter(1);
        const Int_t peakLowBin = h->FindBin(fitMean - SIG_WINDOW);
        const Int_t peakHighBin = h->FindBin(fitMean + SIG_WINDOW);
        const Int_t leftLo = fitBinLo;
        const Int_t leftHi = peakLowBin - 1;
        const Int_t rightLo = peakHighBin + 1;
        const Int_t rightHi = fitBinHi;

        const Int_t sidebandBins = std::max(1, (leftHi - leftLo + 1) + (rightHi - rightLo + 1));
        const Double_t sidebandIntegral = h->Integral(leftLo, leftHi) + h->Integral(rightLo, rightHi);
        const Double_t sidebandDensity = std::max(0.0, sidebandIntegral / sidebandBins); // counts/bin

        const Int_t peakBins = std::max(1, peakHighBin - peakLowBin + 1);
        const Double_t peakTotal = std::max(1.0, h->Integral(peakLowBin, peakHighBin));
        const Double_t peakSignal = std::max(0.0, peakTotal - sidebandDensity * peakBins);

        const Double_t signalFrac = peakSignal / peakTotal;
        const Double_t significance = peakSignal / std::sqrt(peakTotal);

        // local background slope from the two sidebands (seed for b1)
        const Double_t leftCenter = 0.5 * (h->GetBinLowEdge(leftLo) + h->GetBinLowEdge(leftHi + 1));
        const Double_t rightCenter = 0.5 * (h->GetBinLowEdge(rightLo) + h->GetBinLowEdge(rightHi + 1));
        const Double_t leftDensity = (leftHi >= leftLo) ? h->Integral(leftLo, leftHi) / (leftHi - leftLo + 1) : sidebandDensity;
        const Double_t rightDensity = (rightHi >= rightLo) ? h->Integral(rightLo, rightHi) / (rightHi - rightLo + 1) : sidebandDensity;
        const Double_t initial_b1 = (rightCenter > leftCenter) ? (rightDensity - leftDensity) / (rightCenter - leftCenter) : 0.0;

        const Double_t signalInitial = std::max(1.0, std::min(h_integral, peakSignal));

        // --- signal shape held fixed from the inclusive/MC fit (params 1-7, 11-13) ---
        Double_t totalParameters[14] = {
            signalInitial,
            inclusiveFit->GetParameter(1),
            inclusiveFit->GetParameter(2),
            inclusiveFit->GetParameter(3),
            inclusiveFit->GetParameter(4),
            inclusiveFit->GetParameter(5),
            inclusiveFit->GetParameter(6),
            inclusiveFit->GetParameter(7),
            std::max(1.0, sidebandDensity),
            initial_b1,
            0.0,
            inclusiveFit->GetParameter(11),
            inclusiveFit->GetParameter(12),
            inclusiveFit->GetParameter(13)};
        f->SetParameters(totalParameters);
        for (int p = 1; p <= 7; ++p)
          f->FixParameter(p, inclusiveFit->GetParameter(p));
        f->FixParameter(11, inclusiveFit->GetParameter(11)); // tail fraction scales with the fitted yield
        f->FixParameter(12, inclusiveFit->GetParameter(12));
        f->FixParameter(13, inclusiveFit->GetParameter(13));

        // --- background tier from the sideband significance (documented criterion) ---
        //   0: quadratic (b0,b1,b2 free)  1: linear (b0,b1 free)
        //   2: linear, slope fixed from sidebands (b0 free)  3: background fully fixed (escalation only)
        int bkgTier;
        if (significance >= BKG_TIER0_MIN_SIGNIF && signalFrac >= BKG_TIER0_MIN_SIGFRAC)
          bkgTier = 0;
        else if (significance >= BKG_TIER1_MIN_SIGNIF)
          bkgTier = 1;
        else
          bkgTier = 2;

        const Double_t b0hi = 5.0 * std::max(1.0, sidebandDensity);
        const Double_t b1span = 5.0 * std::max(std::fabs(initial_b1), 10.0 * std::max(1.0, sidebandDensity));
        const Double_t b2span = b1span / 0.02;

        auto applyTier = [&](int tier)
        {
          f->SetParameter(0, signalInitial);
          f->SetParLimits(0, 0.0, 1.e9);

          if (tier <= 2) // b0 (background level at the pivot) floats
          {
            f->SetParameter(8, std::max(1.0, sidebandDensity));
            f->ReleaseParameter(8);
            f->SetParLimits(8, 0.0, b0hi);
          }
          else // tier 3: background fully fixed from the sidebands
          {
            f->FixParameter(8, std::max(1.0, sidebandDensity));
          }

          if (tier <= 1) // slope floats
          {
            f->SetParameter(9, initial_b1);
            f->ReleaseParameter(9);
            f->SetParLimits(9, -b1span, b1span);
          }
          else // slope fixed to the two-sideband estimate
          {
            f->FixParameter(9, initial_b1);
          }

          if (tier == 0) // curvature floats
          {
            f->SetParameter(10, 0.0);
            f->ReleaseParameter(10);
            f->SetParLimits(10, -b2span, b2span);
          }
          else
          {
            f->FixParameter(10, 0.0);
          }
        };

        auto fitPass = [&]() -> TFitResultPtr
        {
          h->Fit(f, "L I Q", "", MASS_FIT_MIN, MASS_FIT_MAX);
          h->Fit(f, "L I Q", "", MASS_FIT_MIN, MASS_FIT_MAX);
          return h->Fit(f, "L I Q S", "", MASS_FIT_MIN, MASS_FIT_MAX);
        };

        auto isGood = [&](TFitResultPtr r) -> bool
        {
          if (!r.Get() || !r->IsValid())
            return false;
          const Int_t nd = f->GetNDF();
          if (nd <= 0)
            return false;
          const Double_t y = f->GetParameter(0);
          const Double_t ye = f->GetParError(0);
          if (!std::isfinite(y) || !std::isfinite(ye) || ye <= 0.0)
            return false;
          if (y > 1.0 && ye < 1e-4 * y) // yield pinned at a limit
            return false;
          if (f->GetChisquare() / nd > FIT_CHI2NDF_CEILING)
            return false;
          return true;
        };

        if (bkgTier > 0)
          cout << "BKG TIER " << bkgTier << ": " << label_centbin[icent] << " " << pt_label[ipt] << " " << dcaLab
               << " | entries = " << nEntries << ", signalFrac = " << signalFrac
               << ", significance = " << significance << endl;

        int usedTier = bkgTier;
        applyTier(bkgTier);
        TFitResultPtr fitRes = fitPass();
        for (int t = bkgTier + 1; t <= 3 && !isGood(fitRes); ++t)
        {
          usedTier = t;
          applyTier(t);
          fitRes = fitPass();
        }
        const bool qualityPass = isGood(fitRes);

        total_func[icent][ipt][idca]->SetNpx(2000);
        total_func[icent][ipt][idca]->SetLineWidth(1);
        total_func[icent][ipt][idca]->SetLineColor(kRed);

        if (Massplot)
        {
          h_dmass[icent][ipt][idca]->SetTitle(Form("hdmass_%s_%s_%s", label_centbin[icent], pt_label[ipt], dcaLab));
          h_dmass[icent][ipt][idca]->SetMinimum(0);
          h_dmass[icent][ipt][idca]->Draw("E1");
        }

        // doubleGauss_norm par 5 is widthScale, which is par 6 in TotalFuncModel.
        for (int p = 0; p < 5; p++)
        {
          doubGauss[icent][ipt][idca]->SetParameter(p, total_func[icent][ipt][idca]->GetParameter(p));
        }
        doubGauss[icent][ipt][idca]->SetParameter(5, total_func[icent][ipt][idca]->GetParameter(6));
        doubGauss[icent][ipt][idca]->SetNpx(2000);
        doubGauss[icent][ipt][idca]->SetLineWidth(2);
        doubGauss[icent][ipt][idca]->SetLineColor(kGreen);
        if (Massplot)
          doubGauss[icent][ipt][idca]->Draw("SAME");

        const Int_t fitStatus = fitRes.Get() ? fitRes->Status() : -1;

        if (!qualityPass)
        {
          cout << "FIT NOT CLEAN: "
               << label_centbin[icent] << " "
               << pt_label[ipt] << " "
               << dcaLab
               << " | status = " << fitStatus
               << ", usedTier = " << usedTier
               << ", signalFrac = " << signalFrac
               << ", entries = " << nEntries
               << ", integral = " << nIntegral << endl;
        }

        Yield[icent][ipt][idca] = total_func[icent][ipt][idca]->GetParameter(0) / bin_w;
        Yield_error[icent][ipt][idca] = total_func[icent][ipt][idca]->GetParError(0) / bin_w;

        if (Yield[icent][ipt][idca] < 0)
        {
          Yield[icent][ipt][idca] = 0.0;
          Yield_error[icent][ipt][idca] = 0.0;
        }

        // CRITICAL FIX: Divide yield by the DCA bin width so the plot represents dN/dDCA
        Double_t dca_bin_width = dcaEdges(ipt)[idca + 1] - dcaEdges(ipt)[idca];
        hdca[icent][ipt]->SetBinContent(idca + 1, Yield[icent][ipt][idca] / dca_bin_width);
        hdca[icent][ipt]->SetBinError(idca + 1, Yield_error[icent][ipt][idca] / dca_bin_width);

        /*if (Massplot)
        {
          f_out->cd(Form("MassFits_%s_%s", label_centbin[icent], pt_label[ipt]));
          c[icent][ipt][idca]->Write();
        }*/

        TH1F *h_fit = (TH1F *)h_dmass[icent][ipt][idca]->Clone(Form("dmass_fit_%s_%s_%s", label_centbin[icent], pt_label[ipt], dcaLab));
        h_fit->GetXaxis()->SetRangeUser(MASS_FIT_MIN, MASS_FIT_MAX);
        h_fit->SetTitle(Form("Mass fit: %s %s %s", label_centbin[icent], pt_label[ipt], dcaLab));
        h_fit->SetDirectory(nullptr);

        TF1 *savedTotalFunc = (TF1 *)total_func[icent][ipt][idca]->Clone(Form("total_func_%s", dcaLab));

        h_fit->GetListOfFunctions()->Add(savedTotalFunc);

        const Int_t ndf = total_func[icent][ipt][idca]->GetNDF();
        const Double_t chi2ndf = (ndf > 0) ? total_func[icent][ipt][idca]->GetChisquare() / ndf : -1.0;

        TPaveText *paveText = new TPaveText(0.15, 0.70, 0.48, 0.88, "brNDC");
        paveText->SetBorderSize(0); // No border
        paveText->SetFillColor(0);  // Transparent background
        paveText->SetFillStyle(0);
        paveText->SetTextFont(42);
        paveText->SetTextSize(0.04);
        paveText->SetTextAlign(12); // Left-aligned

        const Double_t fitSignificance = (Yield_error[icent][ipt][idca] > 0) ? Yield[icent][ipt][idca] / Yield_error[icent][ipt][idca] : 0.0;
        paveText->AddText(Form("Yield = %.1f #pm %.1f", Yield[icent][ipt][idca], Yield_error[icent][ipt][idca]));
        paveText->AddText(Form("#chi^{2}/NDF = %.2f", chi2ndf));
        paveText->AddText(Form("Significance = %.2f", fitSignificance));

        h_fit->GetListOfFunctions()->Add(paveText);

        TList *fitInfo = new TList;
        fitInfo->SetName("FitInfo");
        fitInfo->SetOwner(kTRUE);
        fitInfo->Add(new TParameter<Double_t>("Yield", Yield[icent][ipt][idca]));
        fitInfo->Add(new TParameter<Double_t>("YieldError", Yield_error[icent][ipt][idca]));
        fitInfo->Add(new TParameter<Double_t>("Chi2NDF", chi2ndf));
        fitInfo->Add(new TParameter<Int_t>("FitStatus", fitRes.Get() ? fitRes->Status() : -1));
        fitInfo->Add(new TParameter<Int_t>("BkgTier", usedTier));
        fitInfo->Add(new TParameter<Int_t>("QualityPass", qualityPass ? 1 : 0));
        fitInfo->Add(new TParameter<Double_t>("SignalFraction", signalFrac));
        fitInfo->Add(new TParameter<Double_t>("Significance", significance));
        h_fit->GetListOfFunctions()->Add(fitInfo);

        massDir->cd();
        h_fit->Write("", TObject::kOverwrite);
        delete h_fit;

      } //---DCA loop

    } //----pt loop---
  } //---cent loop ---

  f_out->cd();
  for (int icent = 0; icent < N_CENTBIN; icent++)
  {
    for (int i_pt = pt_min; i_pt < pt_max; i_pt++)
    {
      hdca[icent][i_pt]->Write();
    }
  }
  f_out->Write();
  cout << "\n>>> Done! All data DCA distributions saved to " << f_out->GetName() << endl;
  f_out->Close();
}
