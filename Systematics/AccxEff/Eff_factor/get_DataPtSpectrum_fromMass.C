#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TPaveText.h>
#include <TParameter.h>
#include <TList.h>
#include <TDirectory.h>
#include <TMath.h>
#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <Math/PdfFuncMathCore.h>
#include <Math/MinimizerOptions.h>
#include <cmath>

using namespace std;

// --- mass model, identical to get_Data_DCA_fromMass.C / get_PtWeight_fromMass.C ---
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

  const double tailFrac = par[11];
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

  return signalPart + background + signalYield * tailFrac * (cb1 + 4.0 * cb2);
}

const Double_t D0_MASS = 1.86484;
const Double_t MASS_FIT_MIN = 1.76;
const Double_t MASS_FIT_MAX = 1.98;

// --- fine pT binning, must match get_MassSpectra_PtWeight.C ---
const Int_t N_PTBIN_eff = 23;
Double_t ptbinning_[N_PTBIN_eff + 1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 70.0, 100.0};
const char *ptLabel_[N_PTBIN_eff] = {
    "pT2to2p5", "pT2p5to3", "pT3to3p5", "pT3p5to4", "pT4to4p5", "pT4p5to5",
    "pT5to5p5", "pT5p5to6", "pT6to6p5", "pT6p5to7", "pT7to7p5", "pT7p5to8",
    "pT8to9", "pT9to10", "pT10to12p5", "pT12p5to15", "pT15to20", "pT20to25",
    "pT25to30", "pT30to40", "pT40to50", "pT50to70", "pT70to100"};

// --- MC signal-template pT ranges available in D0_MCtemplate_out_combined.root ---
const Int_t N_MC_TEMPLATE = 10;
const char *mcTemplateLabel[N_MC_TEMPLATE] = {"pT1to2", "pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8", "pT8to10", "pT10to15", "pT15to20", "pT20to40"};
Double_t mcTemplateLow[N_MC_TEMPLATE] = {1, 2, 3, 4, 5, 6, 8, 10, 15, 20};
Double_t mcTemplateHigh[N_MC_TEMPLATE] = {2, 3, 4, 5, 6, 8, 10, 15, 20, 1e4}; // last bucket open-ended: covers our bins up to 100 GeV

const char *mcTemplateForFineBin(int ipt)
{
  double center = 0.5 * (ptbinning_[ipt] + ptbinning_[ipt + 1]);
  for (int i = 0; i < N_MC_TEMPLATE; i++)
    if (center >= mcTemplateLow[i] && center < mcTemplateHigh[i])
      return mcTemplateLabel[i];
  return mcTemplateLabel[N_MC_TEMPLATE - 1];
}

// Signal shape derived once per unique MC template (several fine pT bins share the
// same coarse template -- see mcTemplateForFineBin above).
struct SignalShape
{
  Double_t par[14];
};

// Estimate a linear background (level + slope, evaluated about MASS_PIVOT) from
// the mass sidebands, away from the D0 peak. The per-bin fit below used to start
// the background from a blanket guess (0.3*peak, zero slope); that is a poor seed
// wherever the combinatorial background is steep, which is what was causing the
// low-pT mass fits to converge badly. Seeding from the actual sideband shape fixes
// that without touching the fit model itself.
void estimateSidebandBackground(TH1F *h, Double_t &b0est, Double_t &b1est)
{
  const Double_t peakHalfWidth = 0.05; // GeV, mass window excluded around the D0 peak
  const Int_t loLo = h->FindBin(MASS_FIT_MIN + 1e-6);
  const Int_t loHi = h->FindBin(D0_MASS - peakHalfWidth);
  const Int_t hiLo = h->FindBin(D0_MASS + peakHalfWidth);
  const Int_t hiHi = h->FindBin(MASS_FIT_MAX - 1e-6);

  Double_t sumLo = 0.0, nLo = 0.0, sumHi = 0.0, nHi = 0.0;
  for (int b = loLo; b <= loHi; b++)
  {
    sumLo += h->GetBinContent(b);
    nLo += 1.0;
  }
  for (int b = hiLo; b <= hiHi; b++)
  {
    sumHi += h->GetBinContent(b);
    nHi += 1.0;
  }

  if (nLo < 1.0 || nHi < 1.0)
  {
    b0est = std::max(1.0, h->GetMaximum() * 0.3);
    b1est = 0.0;
    return;
  }

  const Double_t mLo = 0.5 * (MASS_FIT_MIN + (D0_MASS - peakHalfWidth));
  const Double_t mHi = 0.5 * ((D0_MASS + peakHalfWidth) + MASS_FIT_MAX);
  const Double_t bkgLo = sumLo / nLo;
  const Double_t bkgHi = sumHi / nHi;

  b1est = (bkgHi - bkgLo) / (mHi - mLo);
  b0est = bkgLo + b1est * (MASS_PIVOT - mLo);
  if (b0est < 0.0)
    b0est = std::max(1.0, h->GetMaximum() * 0.1);
}

SignalShape deriveSignalShapeFromMC(TFile *fileMC, const char *label)
{
  SignalShape shape;

  TH1F *h_mc_match_signal = (TH1F *)fileMC->Get(Form("hMass_Signal_%s", label));
  TH1F *h_mc_match_all = (TH1F *)fileMC->Get(Form("hMass_Signal_Plus_Swap_%s", label));
  if (!h_mc_match_signal || !h_mc_match_all)
  {
    cerr << "Missing MC templates for " << label << endl;
    for (int p = 0; p < 14; p++)
      shape.par[p] = 0.0;
    return shape;
  }

  TF1 *f = new TF1(Form("shapeFit_%s", label), TotalFuncModel, MASS_FIT_MIN, MASS_FIT_MAX, 14);

  f->SetParameter(0, 100.);
  f->SetParameter(1, D0_MASS);
  f->SetParameter(2, 0.03);
  f->SetParameter(3, 0.005);
  f->SetParameter(4, 0.1);
  f->FixParameter(5, 1.0);
  f->FixParameter(6, 0.0);
  f->FixParameter(7, 0.1);
  f->FixParameter(8, 0.0);
  f->FixParameter(9, 0.0);
  f->FixParameter(10, 0.0);
  f->FixParameter(11, 0.0);
  f->FixParameter(12, 0.0);
  f->FixParameter(13, 0.0);

  // Bounds tightened to the physical D0 mass-resolution range (~0.01-0.03 GeV, as seen
  // consistently across the other MC templates). The old 0.005-0.5 / 0.001-0.25 ranges
  // were loose enough that the shared high-pT template ("pT20to40", used for every fine
  // pT bin from 20 up to 100 GeV) converged to a spurious sigma1 ~ 0.2 GeV -- about 10x
  // too wide -- which then got fixed into every one of those data-bin fits and was the
  // root cause of the high-pT fits diverging.
  f->SetParLimits(2, 0.005, 0.08);
  f->SetParLimits(3, 0.001, 0.06);
  f->SetParLimits(4, 0.0, 1.0);
  f->SetParLimits(5, 0.0, 1.0);

  f->FixParameter(1, D0_MASS);
  h_mc_match_signal->Fit(f, "q", "", MASS_FIT_MIN, MASS_FIT_MAX);
  h_mc_match_signal->Fit(f, "q", "", MASS_FIT_MIN, MASS_FIT_MAX);
  f->ReleaseParameter(1);
  h_mc_match_signal->Fit(f, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
  h_mc_match_signal->Fit(f, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
  h_mc_match_signal->Fit(f, "L q ", "", MASS_FIT_MIN, MASS_FIT_MAX);

  f->FixParameter(1, f->GetParameter(1));
  f->FixParameter(2, f->GetParameter(2));
  f->FixParameter(3, f->GetParameter(3));
  f->FixParameter(4, f->GetParameter(4));

  f->ReleaseParameter(5);
  f->ReleaseParameter(7);
  f->SetParameter(7, 0.1);
  f->SetParLimits(7, 0.02, 0.3); // sigma3 (wide-Gaussian component) had no bound at all before
  for (int iteration = 0; iteration < 5; iteration++)
    h_mc_match_all->Fit(f, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);

  f->FixParameter(5, f->GetParameter(5));
  f->FixParameter(7, f->GetParameter(7));

  for (int p = 0; p < 14; p++)
    shape.par[p] = f->GetParameter(p);

  delete f;
  return shape;
}

void get_DataPtSpectrum_fromMass(
    TString data_infile = "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/PtWeight_Sept22_v1/ROOT/DataMass_out_combined.root",
    TString mc_template_file = "D0_MCtemplate_out_combined.root",
    TString outfile_name = "DataPtSpectrum_out.root")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(20000);
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2); // more robust minimization for the low-stat / low-pT bins

  TFile *fin = TFile::Open(data_infile);
  TFile *fileMC = TFile::Open(mc_template_file);
  if (!fin || fin->IsZombie() || !fileMC || fileMC->IsZombie())
  {
    cerr << "Error opening input files (" << data_infile << ", " << mc_template_file << ")" << endl;
    return;
  }

  TFile *fout = new TFile(outfile_name, "RECREATE");
  TDirectory *massDir = fout->mkdir("MassFits");

  TH1D *h_pt_data_raw = new TH1D("h_pt_data_raw", ";p_{T} (GeV);dN/dp_{T} (raw yield)", N_PTBIN_eff, ptbinning_);

  map<string, SignalShape> shapeCache;

  for (int ipt = 0; ipt < N_PTBIN_eff; ipt++)
  {
    TH1F *h = (TH1F *)fin->Get(Form("dmass_%s", ptLabel_[ipt]));
    if (!h)
    {
      cout << "Missing dmass_" << ptLabel_[ipt] << endl;
      continue;
    }

    const char *templateLabel = mcTemplateForFineBin(ipt);
    if (shapeCache.find(templateLabel) == shapeCache.end())
      shapeCache[templateLabel] = deriveSignalShapeFromMC(fileMC, templateLabel);
    const SignalShape &shape = shapeCache[templateLabel];

    TF1 *f = new TF1(Form("total_func_%s", ptLabel_[ipt]), TotalFuncModel, MASS_FIT_MIN, MASS_FIT_MAX, 14);

    const Int_t fitBinLo = h->FindBin(MASS_FIT_MIN + 1e-6);
    const Int_t fitBinHi = h->FindBin(MASS_FIT_MAX - 1e-6);
    const Double_t bin_w = h->GetBinWidth(1);
    const Double_t h_integral = h->Integral(fitBinLo, fitBinHi);
    const Double_t inclusiveMax = h->GetMaximum();

    Double_t b0est, b1est;
    estimateSidebandBackground(h, b0est, b1est);

    Double_t initPars[14] = {
        h_integral * bin_w, shape.par[1], shape.par[2], shape.par[3], shape.par[4],
        shape.par[5], shape.par[6], shape.par[7],
        b0est, b1est, 0.0, 0.0, 0.0, 0.0};
    f->SetParameters(initPars);

    // Signal shape fixed from MC (mean, widths, fractions, swap sigma).
    for (int p = 1; p <= 7; p++)
      f->FixParameter(p, shape.par[p]);
    // Crystal-ball ("swap") peak positions/widths start fixed at their default (0 =
    // hardcoded universal shape); optionally released for a final tail-only refinement
    // further below, once the rest of the fit has already converged well.
    f->FixParameter(12, 0.0);
    f->FixParameter(13, 0.0);

    // Free quadratic background + free yield (this is the per-pT-bin analogue of the
    // "inclusive" mass fit in get_Data_DCA_fromMass.C -- there is no DCA loop here,
    // so every bin gets the same free-background treatment as that macro's tier 0).
    f->SetParameter(0, std::max(1.0, h_integral * bin_w));
    f->SetParLimits(0, 0.0, 1.0e9);
    f->SetParameter(8, std::max(1.0, b0est));
    f->SetParLimits(8, 0.0, 5.0 * std::max(1.0, inclusiveMax));
    // b1/b2 limits scale with the bin's own content (like b0's above) instead of the old
    // flat +-1e6/+-1e7: the low-pT bins carry O(10^8) weighted entries (PbPb combinatorics),
    // and a flat +-1e6 cap on the background slope was being hit exactly, starving the
    // background of the freedom it needed and driving chi2/ndf into the thousands. The
    // fitted b1/b2 in those bins only ever land within ~1x the bin's own peak height, so
    // +-15x gives a wide safety margin without leaving so much slack that a low-statistics
    // bin (e.g. pT50to70, ~44k entries) can wander off into an unphysical, wildly overfit
    // background -- which is what happened when this was left at +-50x/+-500x.
    const Double_t bkgScale = std::max(1.0, inclusiveMax);
    f->SetParameter(9, b1est);
    f->SetParLimits(9, -15.0 * bkgScale, 15.0 * bkgScale);
    f->SetParameter(10, 0.0);
    f->SetParLimits(10, -15.0 * bkgScale, 15.0 * bkgScale);
    f->FixParameter(11, 0.0); // released below, once the background/yield/mean have settled

    h->Fit(f, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
    h->Fit(f, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
    // Release the tail/mean terms for the final pass, same as the reference macro.
    f->ReleaseParameter(1);
    // Keep the peak anchored physically -- previously unbounded, which let the low-stat
    // high-pT bins (sharing the pathological wide template above) run away: one bin's
    // fitted mean diverged to -1.1e6 GeV.
    f->SetParLimits(1, D0_MASS - 0.03, D0_MASS + 0.03);
    // Also let the overall width scale float here. It was fixed at 0 (i.e. use MC's own
    // width completely unmodified) everywhere -- through the MC-shape derivation AND every
    // per-bin data fit -- even though the model carries this parameter specifically to absorb
    // a resolution difference between MC and data. With it stuck at 0, every pT bin's fit was
    // systematically undershooting the peak bin (checked: +0.4% to +8% too low, consistently
    // in the same direction across every bin with real statistics) -- the signature of a
    // fixed width that's too wide relative to the real data peak, propped up with extra yield
    // in the wings instead of actually matching the sharper true peak.
    f->ReleaseParameter(6);
    f->SetParLimits(6, -0.3, 0.3);
    h->Fit(f, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
    h->Fit(f, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
    f->ReleaseParameter(11);
    f->SetParameter(11, 0.5);
    f->SetParLimits(11, 0.0, 20.0);
    h->Fit(f, "L q", "", MASS_FIT_MIN, MASS_FIT_MAX);
    TFitResultPtr fitRes = h->Fit(f, "L q S", "", MASS_FIT_MIN, MASS_FIT_MAX);

    // Final, optional refinement: let the crystal-ball ("swap") peak positions/widths
    // nudge slightly to match this bin's own residual tail structure -- they were
    // otherwise always sitting at their hardcoded universal values. Kept tight (+-15%)
    // since this is a secondary shape correction, not the primary peak, and only kept
    // if it actually improves the fit (clean convergence + lower chi2) -- otherwise the
    // already-good fit above is restored, so this step can only help, never regress it.
    if (fitRes.Get() && fitRes->Status() == 0)
    {
      Double_t parsBeforeTail[14];
      for (int p = 0; p < 14; p++)
        parsBeforeTail[p] = f->GetParameter(p);
      const Double_t chi2Before = f->GetChisquare();

      f->ReleaseParameter(12);
      f->SetParLimits(12, -0.15, 0.15);
      f->ReleaseParameter(13);
      f->SetParLimits(13, -0.15, 0.15);
      TFitResultPtr fitResTail = h->Fit(f, "L q S", "", MASS_FIT_MIN, MASS_FIT_MAX);

      if (fitResTail.Get() && fitResTail->Status() == 0 && f->GetChisquare() <= chi2Before)
      {
        fitRes = fitResTail;
      }
      else
      {
        f->SetParameters(parsBeforeTail);
        f->FixParameter(12, parsBeforeTail[12]);
        f->FixParameter(13, parsBeforeTail[13]);
        fitRes = h->Fit(f, "L q S", "", MASS_FIT_MIN, MASS_FIT_MAX);
      }
    }

    // Low-stat bins occasionally don't converge cleanly on the first pass. Retry from a
    // freshly reset (perturbed) seed each time -- NOT from wherever the failed attempt
    // left b1/b2/yield -- since a non-converged fit can leave the background walked off
    // into an unphysical region, and re-fitting from there tends to compound the problem
    // rather than fix it (this is what was pushing yield/b2 to absurd values in the
    // low-statistics high-pT bins).
    for (int nRetries = 0; (!fitRes.Get() || fitRes->Status() != 0) && nRetries < 3; nRetries++)
    {
      const Double_t perturb = 1.0 + 0.2 * (nRetries + 1);
      f->SetParameter(0, std::max(1.0, h_integral * bin_w));
      f->SetParameter(8, std::max(1.0, b0est) * perturb);
      f->SetParameter(9, b1est);
      f->SetParameter(10, 0.0);
      fitRes = h->Fit(f, "L q S", "", MASS_FIT_MIN, MASS_FIT_MAX);
    }

    f->SetNpx(2000);
    f->SetLineWidth(1);
    f->SetLineColor(kRed);

    Double_t yield = f->GetParameter(0) / bin_w;
    Double_t yieldErr = f->GetParError(0) / bin_w;
    if (yield < 0)
    {
      // No signal above background in this pT bin: treat like get_Data_DCA_fromMass.C
      // does for a negative-yield fit and zero both the value and its error.
      yield = 0.0;
      yieldErr = 0.0;
    }
    const Int_t ndf = f->GetNDF();
    const Double_t chi2ndf = (ndf > 0) ? f->GetChisquare() / ndf : -1.0;
    const Int_t fitStatus = fitRes.Get() ? fitRes->Status() : -1;

    if (!fitRes.Get() || fitRes->Status() != 0)
      cout << "FIT NOT CLEAN: " << ptLabel_[ipt] << " | status = " << fitStatus << endl;

    const Double_t ptBinW = ptbinning_[ipt + 1] - ptbinning_[ipt];
    h_pt_data_raw->SetBinContent(ipt + 1, yield / ptBinW);
    h_pt_data_raw->SetBinError(ipt + 1, yieldErr / ptBinW);

    // Save the mass fit for visual inspection.
    TH1F *h_fit = (TH1F *)h->Clone(Form("dmass_fit_%s", ptLabel_[ipt]));
    h_fit->GetXaxis()->SetRangeUser(MASS_FIT_MIN, MASS_FIT_MAX);
    h_fit->SetTitle(Form("Mass fit: %s (MC template %s)", ptLabel_[ipt], templateLabel));
    h_fit->SetDirectory(nullptr);

    TF1 *savedFunc = (TF1 *)f->Clone(Form("total_func_saved_%s", ptLabel_[ipt]));
    h_fit->GetListOfFunctions()->Add(savedFunc);

    TPaveText *pave = new TPaveText(0.15, 0.70, 0.48, 0.88, "brNDC");
    pave->SetBorderSize(0);
    pave->SetFillColor(0);
    pave->SetFillStyle(0);
    pave->SetTextFont(42);
    pave->SetTextSize(0.04);
    pave->SetTextAlign(12);
    const Double_t significance = (yieldErr > 0) ? yield / yieldErr : 0.0;
    pave->AddText(Form("Yield = %.1f #pm %.1f", yield, yieldErr));
    pave->AddText(Form("#chi^{2}/NDF = %.2f", chi2ndf));
    pave->AddText(Form("Significance = %.2f", significance));
    h_fit->GetListOfFunctions()->Add(pave);

    TList *fitInfo = new TList;
    fitInfo->SetName("FitInfo");
    fitInfo->SetOwner(kTRUE);
    fitInfo->Add(new TParameter<Double_t>("Yield", yield));
    fitInfo->Add(new TParameter<Double_t>("YieldError", yieldErr));
    fitInfo->Add(new TParameter<Double_t>("Chi2NDF", chi2ndf));
    fitInfo->Add(new TParameter<Int_t>("FitStatus", fitStatus));
    h_fit->GetListOfFunctions()->Add(fitInfo);

    massDir->cd();
    h_fit->Write("", TObject::kOverwrite);
    delete h_fit;
  }

  fout->cd();
  h_pt_data_raw->Write();
  fout->Write();
  fout->Close();

  cout << "\n>>> Done! Data pT spectrum (h_pt_data_raw) saved to " << outfile_name
       << " -- feed this into get_PtWeight_fromMass.C's data_pt_infile.\n";
}
