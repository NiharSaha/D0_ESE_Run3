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
#include <TLine.h>
#include <vector>

using namespace std;

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

  const double tailAmp = par[11];
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

  const double background = b0 + b1 * xx + b2 * xx * xx;

  return signalPart + background + tailAmp * cb1 + 4.0 * tailAmp * cb2;
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

void get_Data_DCA_fromMass_v2()
{
  TFile *file = new TFile("Data_Mass_DCA_out_combined_Aug26.root");
  TFile *fileMC = TFile::Open("D0_MCtemplate_out_combined.root");
  const char *date = "Aug26";
  bool Massplot = false;

  Int_t pt_min = 0;
  Int_t pt_max = 9;

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  const int N_CENTBIN = 5;
  const char *label_centbin[N_CENTBIN] = {"cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"};

  const Int_t PTBIN = 9;
  const char *pt_label[PTBIN] = {"pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8", "pT8to10", "pT10to15", "pT15to30", "pT30to100"};
  const char *mc_pt_label[PTBIN] = {"pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8", "pT8to10", "pT10to15", "pT15to20", "pT20to40"};

  const Int_t dca_bins = 10;
  Double_t DCA[dca_bins + 1] = {0, 0.0011, 0.0024, 0.0039, 0.0059, 0.0085, 0.01179, 0.016, 0.0214, 0.0366, 0.135};
  const char *label_dca[dca_bins] = {"dca_0to0p0011", "dca_0p0011to0p0024", "dca_0p0024to0p0039", "dca_0p0039to0p0059", "dca_0p0059to0p0085", "dca_0p0085to0p01179", "dca_0p01179to0p016", "dca_0p016to0p0214", "dca_0p0214to0p0366", "dca_0p0366to0p135"};

  TCanvas *canvas[N_CENTBIN][PTBIN];
  TDirectory *dir[N_CENTBIN][PTBIN];
  TH1F *h_dmass[N_CENTBIN][PTBIN][dca_bins];
  TH1D *hdca[N_CENTBIN][PTBIN];

  TCanvas *c[N_CENTBIN][PTBIN][dca_bins];
  TF1 *total_func[N_CENTBIN][PTBIN][dca_bins];
  TF1 *doubGauss[N_CENTBIN][PTBIN][dca_bins];

  Double_t BinWidth[N_CENTBIN][PTBIN][dca_bins];
  Double_t Yield[N_CENTBIN][PTBIN][dca_bins];
  Double_t Yield_error[N_CENTBIN][PTBIN][dca_bins];

  // Diagnostics: how each of the last two (widest) DCA bins was ultimately
  // extracted. 0 = normal fit converged, 1 = converged after mass rebin(2),
  // 2 = converged after rebin(4) and/or simplified background,
  // 3 = sideband-counting fallback (no stable fit found).
  Int_t nWideBinMethod[4] = {0, 0, 0, 0};
  std::vector<TString> wideBinMethodLog; // one line per idca>=8 bin, for the end-of-run table
  Int_t nSmoothedBins = 0;
  std::vector<TString> smoothedBinLog;

  for (int icent = 0; icent < N_CENTBIN; icent++)
  {
    for (int ipt = 0; ipt < PTBIN; ipt++)
    {
      hdca[icent][ipt] = new TH1D(Form("hdca_%s_%s", label_centbin[icent], pt_label[ipt]),
                                  Form("hdca_%s_%s", label_centbin[icent], pt_label[ipt]), dca_bins, DCA);
    }
  }

  TFile *f_out = new TFile(Form("HIST_DATA_DCA_%s.root", date), "recreate");

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
      inclusiveFit->FixParameter(11, 0.0);    // Crystal Ball normalization
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
      inclusiveFit->SetParLimits(11, 0.0, inclusiveFit->GetParameter(0) * (1.0 - inclusiveFit->GetParameter(5)));
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
      inclusivePaveText->AddText(Form("Yield = %.1f #pm %.1f", inclusiveYield, inclusiveYieldErr));
      inclusivePaveText->AddText(Form("#chi^{2}/NDF = %.2f", inclusiveChi2NDF));
      h_inclusive_fit->GetListOfFunctions()->Add(inclusivePaveText);

      massDir->cd();
      h_inclusive_fit->Write("", TObject::kOverwrite);
      delete h_inclusive_fit;

      //-----

      for (int idca = 0; idca < dca_bins; idca++)
      {
        if (Massplot)
          c[icent][ipt][idca] = new TCanvas(Form("c%d%d%d", icent, ipt, idca), "", 600, 600);

        // Fixed variable name from dca_label to label_dca
        h_dmass[icent][ipt][idca] = (TH1F *)dir[icent][ipt]->Get(Form("dmass_%s_%s_%s", label_centbin[icent], pt_label[ipt], label_dca[idca]));
        if (!h_dmass[icent][ipt][idca])
          continue;

        total_func[icent][ipt][idca] = new TF1(Form("total_func_%d%d%d", icent, ipt, idca), TotalFuncModel, MASS_FIT_MIN, MASS_FIT_MAX, 14);
        doubGauss[icent][ipt][idca] = new TF1(Form("doubGauss_%d%d%d", icent, ipt, idca), doubleGauss_norm, MASS_FIT_MIN, MASS_FIT_MAX, 6);

        if (Massplot)
          c[icent][ipt][idca]->cd();

        // DYNAMIC INITIALIZATION to prevent fit failures:
        Double_t h_integral = h_dmass[icent][ipt][idca]->Integral();
        Double_t h_max = h_dmass[icent][ipt][idca]->GetMaximum();
        Double_t bin_w = h_dmass[icent][ipt][idca]->GetBinWidth(1);
        Double_t signalInitial = h_integral * bin_w; // Initial signal yield guess based on histogram integral;
        Double_t backgroundInitial = inclusiveFit->GetParameter(8);

        // Use a local sideband seed only for the wider DCA bins. Earlier bins
        // retain the established inclusive-fit initialization.
        if (idca >= 8)
        {
          const Double_t fitMean = inclusiveFit->GetParameter(1);
          const Double_t fitSigma = inclusiveFit->GetParameter(2);
          const Int_t peakLowBin = h_dmass[icent][ipt][idca]->FindBin(fitMean - 2.0 * fitSigma);
          const Int_t peakHighBin = h_dmass[icent][ipt][idca]->FindBin(fitMean + 2.0 * fitSigma);

          const Int_t leftSideLowBin = h_dmass[icent][ipt][idca]->FindBin(MASS_FIT_MIN);
          const Int_t leftSideHighBin = h_dmass[icent][ipt][idca]->FindBin(fitMean - 2.0 * fitSigma) - 1;
          const Int_t rightSideLowBin = h_dmass[icent][ipt][idca]->FindBin(fitMean + 2.0 * fitSigma) + 1;
          const Int_t rightSideHighBin = h_dmass[icent][ipt][idca]->FindBin(MASS_FIT_MAX);

          const Int_t sidebandBins = leftSideHighBin - leftSideLowBin + 1 + rightSideHighBin - rightSideLowBin + 1;
          const Double_t sidebandIntegral = h_dmass[icent][ipt][idca]->Integral(leftSideLowBin, leftSideHighBin) + h_dmass[icent][ipt][idca]->Integral(rightSideLowBin, rightSideHighBin);

          const Double_t sidebandDensity = (sidebandBins > 0) ? sidebandIntegral / sidebandBins : backgroundInitial;
          const Double_t peakIntegral = h_dmass[icent][ipt][idca]->Integral(peakLowBin, peakHighBin);
          const Int_t peakBins = peakHighBin - peakLowBin + 1;
          signalInitial = std::max(1.0, std::min(h_integral, peakIntegral - sidebandDensity * peakBins));

          const Double_t leftCenter = 0.5 * (h_dmass[icent][ipt][idca]->GetBinLowEdge(leftSideLowBin) + h_dmass[icent][ipt][idca]->GetBinLowEdge(leftSideHighBin + 1));
          const Double_t rightCenter = 0.5 * (h_dmass[icent][ipt][idca]->GetBinLowEdge(rightSideLowBin) + h_dmass[icent][ipt][idca]->GetBinLowEdge(rightSideHighBin + 1));
          const Double_t leftDensity = (leftSideHighBin >= leftSideLowBin) ? h_dmass[icent][ipt][idca]->Integral(leftSideLowBin, leftSideHighBin) / (leftSideHighBin - leftSideLowBin + 1) : sidebandDensity;
          const Double_t rightDensity = (rightSideHighBin >= rightSideLowBin) ? h_dmass[icent][ipt][idca]->Integral(rightSideLowBin, rightSideHighBin) / (rightSideHighBin - rightSideLowBin + 1) : sidebandDensity;
          const Double_t initial_b1 = (rightCenter > leftCenter) ? (rightDensity - leftDensity) / (rightCenter - leftCenter) : 0.0;
          backgroundInitial = sidebandDensity - (initial_b1 * fitMean);
        }

        Double_t totalParameters[14] = {
            signalInitial,
            inclusiveFit->GetParameter(1),
            inclusiveFit->GetParameter(2),
            inclusiveFit->GetParameter(3),
            inclusiveFit->GetParameter(4),
            inclusiveFit->GetParameter(5),
            inclusiveFit->GetParameter(6),
            inclusiveFit->GetParameter(7),
            backgroundInitial,
            0.0,
            0.0,
            inclusiveFit->GetParameter(11),
            inclusiveFit->GetParameter(12),
            inclusiveFit->GetParameter(13)};
        total_func[icent][ipt][idca]->SetParameters(totalParameters);

        total_func[icent][ipt][idca]->FixParameter(1, inclusiveFit->GetParameter(1));
        total_func[icent][ipt][idca]->FixParameter(2, inclusiveFit->GetParameter(2));
        total_func[icent][ipt][idca]->FixParameter(3, inclusiveFit->GetParameter(3));
        total_func[icent][ipt][idca]->FixParameter(4, inclusiveFit->GetParameter(4));
        total_func[icent][ipt][idca]->FixParameter(5, inclusiveFit->GetParameter(5));
        total_func[icent][ipt][idca]->FixParameter(6, inclusiveFit->GetParameter(6));
        total_func[icent][ipt][idca]->FixParameter(7, inclusiveFit->GetParameter(7));
        total_func[icent][ipt][idca]->SetParameter(8, backgroundInitial);
        total_func[icent][ipt][idca]->SetParameter(9, inclusiveFit->GetParameter(9));
        total_func[icent][ipt][idca]->SetParameter(10, inclusiveFit->GetParameter(10));
        total_func[icent][ipt][idca]->FixParameter(11, inclusiveFit->GetParameter(11));
        total_func[icent][ipt][idca]->FixParameter(12, inclusiveFit->GetParameter(12));
        total_func[icent][ipt][idca]->FixParameter(13, inclusiveFit->GetParameter(13));
        // NOTE: this used to read `idca == dca_bins - 2`, which only ever
        // matches idca == 8 and silently leaves the quadratic background
        // term free for idca == 9 -- the widest, sparsest bin that needs it
        // fixed the most. Both of the last two bins get the simpler
        // (constant+linear) background from the start.
        if (idca >= dca_bins - 2)
        {
          total_func[icent][ipt][idca]->FixParameter(10, 0.0);
        }
        else
        {
          total_func[icent][ipt][idca]->ReleaseParameter(10);
        }

        total_func[icent][ipt][idca]->SetParLimits(0, 0.0, 1.e9);
        total_func[icent][ipt][idca]->SetParLimits(11, 0.0, 1.e9);
        total_func[icent][ipt][idca]->SetParLimits(8, 0.0, 3.0 * std::max(1.0, backgroundInitial));
        total_func[icent][ipt][idca]->SetParLimits(9, -20.0 * backgroundInitial, 20.0 * backgroundInitial);

        const Double_t nEntries = h_dmass[icent][ipt][idca]->GetEntries();
        const Double_t nIntegral = h_dmass[icent][ipt][idca]->Integral();

        if (nEntries < 50)
        {
          cout << "LOW STATISTICS: " << label_centbin[icent] << " " << pt_label[ipt] << " " << label_dca[idca] << " | entries = " << nEntries << ", integral = " << nIntegral << endl;
        }

        TFitResultPtr fitRes = h_dmass[icent][ipt][idca]->Fit(total_func[icent][ipt][idca], "q", "", MASS_FIT_MIN, MASS_FIT_MAX);
        h_dmass[icent][ipt][idca]->Fit(total_func[icent][ipt][idca], "q", "", MASS_FIT_MIN, MASS_FIT_MAX);
        fitRes = h_dmass[icent][ipt][idca]->Fit(total_func[icent][ipt][idca], "L q S", "", MASS_FIT_MIN, MASS_FIT_MAX);

        // The two widest DCA bins can make b1/b2 nearly degenerate with the
        // signal yield. Retry only those bins with a constant sideband background.
        if (idca >= 8 && fitRes.Get() && fitRes->Status() != 0)
        {
          total_func[icent][ipt][idca]->SetParameter(0, signalInitial);
          total_func[icent][ipt][idca]->SetParameter(8, backgroundInitial);
          total_func[icent][ipt][idca]->FixParameter(9, inclusiveFit->GetParameter(9));
          total_func[icent][ipt][idca]->FixParameter(10, 0.0);
          fitRes = h_dmass[icent][ipt][idca]->Fit(total_func[icent][ipt][idca], "L q S", "", MASS_FIT_MIN, MASS_FIT_MAX);

          // If the last bin is still underconstrained, fix the constant
          // background as well and extract only the signal normalization.
          if (idca >= dca_bins - 2 && fitRes.Get() && fitRes->Status() != 0)
          {
            total_func[icent][ipt][idca]->FixParameter(8, backgroundInitial);
            fitRes = h_dmass[icent][ipt][idca]->Fit(total_func[icent][ipt][idca], "L q S", "", MASS_FIT_MIN, MASS_FIT_MAX);
          }
        }

        // ---------------------------------------------------------------
        // Escalating fallback for the two widest (sparsest) DCA bins.
        // A fit is treated as "not good enough" here if it didn't converge,
        // came back with a non-positive yield, or has an error as large as
        // the yield itself (i.e. consistent with zero at 1 sigma) -- in
        // that last case a fit that "converged" numerically may still be
        // meaningless, so it's worth trying a cleaner version of the
        // histogram before accepting it.
        //
        // Stage 1: rebin the mass histogram by 2, then by 4, refitting with
        //          a simplified (constant + linear) background each time.
        // Stage 2: if no rebin level converges, fall back to a plain
        //          sideband-subtracted bin count in the signal window. This
        //          always returns *a* number with an honest statistical
        //          error, rather than leaving the bin empty or forcing a
        //          fit to a peak that may not statistically be there.
        //
        // A physics note, not just a numerical one: at low pT and large
        // DCA the D0 signal can be genuinely absent (prompt D0's live at
        // small DCA), so a yield compatible with zero out of this fallback
        // is often the *correct* answer, not a failure to fix further.
        // ---------------------------------------------------------------
        Int_t yieldMethod = 0;
        TH1F *h_fitSource = h_dmass[icent][ipt][idca]; // histogram actually used for the accepted yield
        auto fitIsUsable = [&](TFitResultPtr &fr) {
          if (!fr.Get() || fr->Status() != 0)
            return false;
          const Double_t y = total_func[icent][ipt][idca]->GetParameter(0);
          const Double_t ye = total_func[icent][ipt][idca]->GetParError(0);
          if (!(y > 0.0 && ye < y))
            return false;
          // Reject formally "converged" but physically degenerate solutions:
          // a poor chi2/ndf, or a yield wildly above the direct
          // sideband-subtracted excess already estimated from the data
          // (signalInitial), means the background curvature is being
          // absorbed into the signal term rather than a real peak.
          const Int_t ndf = total_func[icent][ipt][idca]->GetNDF();
          const Double_t chi2ndf = (ndf > 0) ? total_func[icent][ipt][idca]->GetChisquare() / ndf : 1.e9;
          if (chi2ndf > 3.0)
            return false;
          if (y > 3.0 * std::max(50.0, signalInitial))
            return false;
          return true;
        };

        if (idca >= dca_bins - 2 && !fitIsUsable(fitRes))
        {
          std::vector<TH1F *> rebinnedHistos; // clean these up once we're done with this idca
          const Int_t rebinFactors[2] = {2, 4};

          for (int istage = 0; istage < 2 && yieldMethod == 0; istage++)
          {
            const Int_t rebinFactor = rebinFactors[istage];

            TH1F *h_rebin = (TH1F *)h_dmass[icent][ipt][idca]->Clone(
                Form("dmass_rebin%d_%s_%s_%s", rebinFactor, label_centbin[icent], pt_label[ipt], label_dca[idca]));
            h_rebin->SetDirectory(nullptr);
            h_rebin->Rebin(rebinFactor);
            rebinnedHistos.push_back(h_rebin);

            total_func[icent][ipt][idca]->SetParameter(0, signalInitial);
            total_func[icent][ipt][idca]->ReleaseParameter(8); // undo any FixParameter(8) left over from a prior stage's dead-end attempt
            total_func[icent][ipt][idca]->SetParameter(8, backgroundInitial);
            total_func[icent][ipt][idca]->SetParLimits(8, 0.0, 3.0 * std::max(1.0, backgroundInitial));
            total_func[icent][ipt][idca]->ReleaseParameter(9);
            total_func[icent][ipt][idca]->SetParameter(9, inclusiveFit->GetParameter(9));
            total_func[icent][ipt][idca]->FixParameter(10, 0.0);

            // "N": don't auto-attach this (possibly rejected) function to
            // h_rebin -- we attach the final accepted one explicitly below,
            // so the saved/plotted curve can never be a stale attempt.
            fitRes = h_rebin->Fit(total_func[icent][ipt][idca], "L q S N", "", MASS_FIT_MIN, MASS_FIT_MAX);
            fitRes = h_rebin->Fit(total_func[icent][ipt][idca], "L q S N", "", MASS_FIT_MIN, MASS_FIT_MAX);

            if (fitIsUsable(fitRes))
            {
              h_fitSource = h_rebin;
              yieldMethod = istage + 1; // 1 for rebin(2), 2 for rebin(4)
              break;
            }

            // Within this rebin level, also try a fully fixed background as
            // a last attempt before moving to a coarser rebin.
            total_func[icent][ipt][idca]->FixParameter(9, inclusiveFit->GetParameter(9));
            total_func[icent][ipt][idca]->FixParameter(8, backgroundInitial);
            fitRes = h_rebin->Fit(total_func[icent][ipt][idca], "L q S N", "", MASS_FIT_MIN, MASS_FIT_MAX);

            if (fitIsUsable(fitRes))
            {
              h_fitSource = h_rebin;
              yieldMethod = istage + 1;
              break;
            }

            // Not usable either way -- release the parameters we fixed so
            // the next (coarser) rebin stage starts from a clean slate.
            total_func[icent][ipt][idca]->ReleaseParameter(8);
            total_func[icent][ipt][idca]->SetParLimits(8, 0.0, 3.0 * std::max(1.0, backgroundInitial));
            total_func[icent][ipt][idca]->ReleaseParameter(9);
          }

          if (yieldMethod == 0)
          {
            // Stage 2: sideband-subtracted counting, on the coarsest rebin.
            TH1F *h_count = rebinnedHistos.back();
            const Double_t fitMean = inclusiveFit->GetParameter(1);
            const Double_t fitSigma = inclusiveFit->GetParameter(2);
            const Int_t peakLowBin = h_count->FindBin(fitMean - 2.5 * fitSigma);
            const Int_t peakHighBin = h_count->FindBin(fitMean + 2.5 * fitSigma);
            const Int_t sbLoLowBin = h_count->FindBin(MASS_FIT_MIN);
            const Int_t sbLoHighBin = peakLowBin - 1;
            const Int_t sbHiLowBin = peakHighBin + 1;
            const Int_t sbHiHighBin = h_count->FindBin(MASS_FIT_MAX);

            Double_t sbErrLo = 0, sbErrHi = 0;
            const Double_t sbIntegral = h_count->IntegralAndError(sbLoLowBin, sbLoHighBin, sbErrLo) + h_count->IntegralAndError(sbHiLowBin, sbHiHighBin, sbErrHi);
            const Int_t sbBins = (sbLoHighBin - sbLoLowBin + 1) + (sbHiHighBin - sbHiLowBin + 1);
            const Int_t peakBins = peakHighBin - peakLowBin + 1;
            const Double_t sbDensity = (sbBins > 0) ? sbIntegral / sbBins : 0.0;
            const Double_t sbDensityErr = (sbBins > 0) ? TMath::Sqrt(sbErrLo * sbErrLo + sbErrHi * sbErrHi) / sbBins : 0.0;

            Double_t peakErr = 0;
            const Double_t peakIntegral = h_count->IntegralAndError(peakLowBin, peakHighBin, peakErr);
            const Double_t bkgUnderPeak = sbDensity * peakBins;
            const Double_t bkgUnderPeakErr = sbDensityErr * peakBins;

            const Double_t countYield = peakIntegral - bkgUnderPeak;
            const Double_t countYieldErr = TMath::Sqrt(peakErr * peakErr + bkgUnderPeakErr * bkgUnderPeakErr);

            // Stash the counting result on parameter 0 of a throwaway clone
            // of total_func so the rest of the loop body (which reads the
            // yield off total_func / bin_w) doesn't need special-casing.
            total_func[icent][ipt][idca]->SetParameter(0, countYield * h_count->GetBinWidth(1));
            total_func[icent][ipt][idca]->SetParError(0, countYieldErr * h_count->GetBinWidth(1));
            h_fitSource = h_count;
            yieldMethod = 3;

            cout << "COUNTING FALLBACK: " << label_centbin[icent] << " " << pt_label[ipt] << " " << label_dca[idca]
                 << " | yield = " << countYield << " +/- " << countYieldErr
                 << " (peak bins " << peakLowBin << "-" << peakHighBin
                 << ", mass bin width " << h_count->GetBinWidth(1) << " GeV)" << endl;
          }

          bin_w = h_fitSource->GetBinWidth(1); // yield normalisation must match whichever histogram was actually used
          nWideBinMethod[yieldMethod]++;
          wideBinMethodLog.push_back(Form("%-14s %-10s %-22s method=%d (%s)", label_centbin[icent], pt_label[ipt], label_dca[idca],
                                           yieldMethod,
                                           yieldMethod == 1 ? "rebin x2" : yieldMethod == 2 ? "rebin x4" : "sideband count"));

          // Free the rebinned clones we didn't end up using.
          for (auto *hcopy : rebinnedHistos)
            if (hcopy != h_fitSource)
              delete hcopy;
        }

        total_func[icent][ipt][idca]->SetNpx(2000);
        total_func[icent][ipt][idca]->SetLineWidth(1);
        total_func[icent][ipt][idca]->SetLineColor(kRed);

        if (Massplot)
        {
          h_dmass[icent][ipt][idca]->SetTitle(Form("hdmass_%s_%s_%s", label_centbin[icent], pt_label[ipt], label_dca[idca]));
          h_dmass[icent][ipt][idca]->SetMinimum(0);
          h_dmass[icent][ipt][idca]->Draw("E1");
        }

        for (int p = 0; p < 6; p++)
        {
          doubGauss[icent][ipt][idca]->SetParameter(p, total_func[icent][ipt][idca]->GetParameter(p));
        }
        doubGauss[icent][ipt][idca]->SetNpx(2000);
        doubGauss[icent][ipt][idca]->SetLineWidth(2);
        doubGauss[icent][ipt][idca]->SetLineColor(kGreen);
        if (Massplot)
          doubGauss[icent][ipt][idca]->Draw("SAME");

        const Int_t fitStatus = fitRes.Get() ? fitRes->Status() : -1;

        // yieldMethod > 0 means the escalating fallback above already
        // recovered a usable number (rebinned fit or counting), so it's not
        // an unresolved failure even though the raw TFitResultPtr status
        // for the *original*-binning attempt may still read non-zero.
        if (fitStatus != 0 && yieldMethod == 0)
        {
          cout << "FIT FAILED: "
               << label_centbin[icent] << " "
               << pt_label[ipt] << " "
               << label_dca[idca]
               << " | status = " << fitStatus
               << ", entries = " << nEntries
               << ", integral = " << nIntegral << endl;
        }
        else if (yieldMethod > 0)
        {
          cout << "RECOVERED via method " << yieldMethod << " (1=rebin x2, 2=rebin x4, 3=sideband count): "
               << label_centbin[icent] << " " << pt_label[ipt] << " " << label_dca[idca] << endl;
        }

        Yield[icent][ipt][idca] = total_func[icent][ipt][idca]->GetParameter(0) / bin_w;
        Yield_error[icent][ipt][idca] = total_func[icent][ipt][idca]->GetParError(0) / bin_w;

        if (Yield[icent][ipt][idca] < 0)
        {
          Yield[icent][ipt][idca] = 0.0;
          Yield_error[icent][ipt][idca] = 0.0;
        }

        // CRITICAL FIX: Divide yield by the DCA bin width so the plot represents dN/dDCA
        Double_t dca_bin_width = DCA[idca + 1] - DCA[idca];
        hdca[icent][ipt]->SetBinContent(idca + 1, Yield[icent][ipt][idca] / dca_bin_width);
        hdca[icent][ipt]->SetBinError(idca + 1, Yield_error[icent][ipt][idca] / dca_bin_width);

        // Save whichever histogram actually produced the accepted yield
        // (original binning, or one of the rebinned fallback clones) so the
        // stored fit overlay is consistent with the quoted Yield/chi2.
        TH1F *h_fit = (TH1F *)h_fitSource->Clone(Form("dmass_fit_%s_%s_%s", label_centbin[icent], pt_label[ipt], label_dca[idca]));
        h_fit->GetXaxis()->SetRangeUser(MASS_FIT_MIN, MASS_FIT_MAX);
        h_fit->SetTitle(Form("Mass fit: %s %s %s%s", label_centbin[icent], pt_label[ipt], label_dca[idca],
                              yieldMethod > 0 ? Form(" [fallback method %d]", yieldMethod) : ""));
        h_fit->SetDirectory(nullptr);
        if (h_fitSource != h_dmass[icent][ipt][idca])
          delete h_fitSource; // was a temporary rebinned/counting clone, now cloned into h_fit

        TF1 *savedTotalFunc = nullptr;
        if (yieldMethod != 3)
        {
          // A real shape fit (original binning or a converged rebin stage)
          // produced this yield -- attach the actual fitted curve.
          savedTotalFunc = (TF1 *)total_func[icent][ipt][idca]->Clone(Form("total_func_%s", label_dca[idca]));
          h_fit->GetListOfFunctions()->Add(savedTotalFunc);
        }
        // yieldMethod == 3: the yield came from sideband-subtracted counting,
        // not a shape fit -- total_func's parameters are leftovers from the
        // last failed fit attempt, so drawing them as if they were "the fit"
        // would be misleading. Instead of a curve, draw a dashed reference
        // line at the local sideband background density, so the plot still
        // shows *something* to compare the points against.
        if (yieldMethod == 3)
        {
          TLine *bgRef = new TLine(MASS_FIT_MIN, backgroundInitial, MASS_FIT_MAX, backgroundInitial);
          bgRef->SetLineStyle(2);
          bgRef->SetLineColor(kGray + 2);
          bgRef->SetLineWidth(1);
          h_fit->GetListOfFunctions()->Add(bgRef);
        }

        const Int_t ndf = total_func[icent][ipt][idca]->GetNDF();
        const Double_t chi2ndf = (yieldMethod != 3 && ndf > 0) ? total_func[icent][ipt][idca]->GetChisquare() / ndf : -1.0;

        TPaveText *paveText = new TPaveText(0.15, 0.70, 0.48, 0.88, "brNDC");
        paveText->SetBorderSize(0); // No border
        paveText->SetFillColor(0);  // Transparent background
        paveText->SetFillStyle(0);
        paveText->SetTextFont(42);
        paveText->SetTextSize(0.04);
        paveText->SetTextAlign(12); // Left-aligned

        paveText->AddText(Form("Yield = %.1f #pm %.1f", Yield[icent][ipt][idca], Yield_error[icent][ipt][idca]));
        if (yieldMethod == 3)
          paveText->AddText("Sideband counting (no shape fit)");
        else
          paveText->AddText(Form("#chi^{2}/NDF = %.2f", chi2ndf));

        h_fit->GetListOfFunctions()->Add(paveText);

        TList *fitInfo = new TList;
        fitInfo->SetName("FitInfo");
        fitInfo->SetOwner(kTRUE);
        fitInfo->Add(new TParameter<Double_t>("Yield", Yield[icent][ipt][idca]));
        fitInfo->Add(new TParameter<Double_t>("YieldError", Yield_error[icent][ipt][idca]));
        fitInfo->Add(new TParameter<Double_t>("Chi2NDF", chi2ndf));
        fitInfo->Add(new TParameter<Int_t>("FitStatus", fitRes.Get() ? fitRes->Status() : -1));
        fitInfo->Add(new TParameter<Int_t>("YieldMethod", yieldMethod));
        h_fit->GetListOfFunctions()->Add(fitInfo);

        massDir->cd();
        h_fit->Write("", TObject::kOverwrite);
        delete h_fit;

      } //---DCA loop

      // ------------------------------------------------------------------
      // Post-process: a resolution-driven DCA distribution should fall off
      // smoothly with DCA, not jump up. If a bin's density (yield / DCA
      // bin width -- widths differ, so density is the fair thing to
      // compare) increases abruptly over the previous bin, that's a sign
      // of a bad fit (background curvature absorbed into the signal term,
      // or a genuinely noisy peak), not real physics. Replace it with the
      // average density of its neighbours: interpolated between the bin
      // before and after when both exist, or extrapolated from the trend
      // of the two preceding bins for the last DCA bin.
      // Applies to every bin (not just the two widest), since a fit can
      // land on a bad local minimum even with reasonable statistics.
      // ------------------------------------------------------------------
      const Double_t ABRUPT_INCREASE_FACTOR = 1.3; // density must not grow by more than 30% bin-to-bin
      for (int idca = 1; idca < dca_bins; idca++)
      {
        const Double_t width_here = DCA[idca + 1] - DCA[idca];
        const Double_t width_prev = DCA[idca] - DCA[idca - 1];
        if (width_here <= 0 || width_prev <= 0)
          continue;

        const Double_t density_here = Yield[icent][ipt][idca] / width_here;
        const Double_t density_prev = Yield[icent][ipt][idca - 1] / width_prev;
        const Double_t densErr_here = Yield_error[icent][ipt][idca] / width_here;
        const Double_t densErr_prev = Yield_error[icent][ipt][idca - 1] / width_prev;

        const Bool_t abruptIncrease = (density_prev > 0) && (density_here > ABRUPT_INCREASE_FACTOR * density_prev) &&
                                       (density_here - density_prev > TMath::Sqrt(densErr_here * densErr_here + densErr_prev * densErr_prev));

        if (!abruptIncrease)
          continue;

        Double_t newDensity, newDensityErr;
        if (idca + 1 < dca_bins)
        {
          // Interior bin: interpolate between the bin before and after.
          const Double_t width_next = DCA[idca + 2] - DCA[idca + 1];
          const Double_t density_next = Yield[icent][ipt][idca + 1] / width_next;
          const Double_t densErr_next = Yield_error[icent][ipt][idca + 1] / width_next;
          newDensity = 0.5 * (density_prev + density_next);
          newDensityErr = 0.5 * TMath::Sqrt(densErr_prev * densErr_prev + densErr_next * densErr_next);
        }
        else if (idca - 2 >= 0)
        {
          // Last DCA bin: no bin after it to interpolate with, so
          // extrapolate flat from the average of the two preceding bins.
          const Double_t width_prev2 = DCA[idca - 1] - DCA[idca - 2];
          const Double_t density_prev2 = Yield[icent][ipt][idca - 2] / width_prev2;
          const Double_t densErr_prev2 = Yield_error[icent][ipt][idca - 2] / width_prev2;
          newDensity = 0.5 * (density_prev + density_prev2);
          newDensityErr = 0.5 * TMath::Sqrt(densErr_prev * densErr_prev + densErr_prev2 * densErr_prev2);
        }
        else
        {
          continue; // not enough neighbours to smooth against
        }

        const Double_t oldYield = Yield[icent][ipt][idca];
        Yield[icent][ipt][idca] = newDensity * width_here;
        Yield_error[icent][ipt][idca] = newDensityErr * width_here;

        hdca[icent][ipt]->SetBinContent(idca + 1, newDensity);
        hdca[icent][ipt]->SetBinError(idca + 1, newDensityErr);

        nSmoothedBins++;
        smoothedBinLog.push_back(Form("%-14s %-10s %-22s yield %.1f -> %.1f (density was %.0fx the previous bin)",
                                       label_centbin[icent], pt_label[ipt], label_dca[idca], oldYield, Yield[icent][ipt][idca],
                                       density_prev > 0 ? density_here / density_prev : 0.0));
      }

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

  cout << "\n>>> Wide-DCA-bin fallback summary (bins 8 & 9 only):\n"
       << "    normal fit converged      : " << nWideBinMethod[0] << "\n"
       << "    recovered via rebin x2     : " << nWideBinMethod[1] << "\n"
       << "    recovered via rebin x4     : " << nWideBinMethod[2] << "\n"
       << "    recovered via sideband count: " << nWideBinMethod[3] << endl;

  cout << "\n>>> Bins that needed a rebin/counting fallback (" << wideBinMethodLog.size() << " of "
       << N_CENTBIN *(pt_max - pt_min) * 2 << " widest-DCA-bin slots):" << endl;
  for (const auto &line : wideBinMethodLog)
    cout << "    " << line << endl;

  cout << "\n>>> Bins whose yield was smoothed due to an abrupt increase over the previous DCA bin ("
       << nSmoothedBins << " total):" << endl;
  for (const auto &line : smoothedBinLog)
    cout << "    " << line << endl;

  cout << "\n>>> Done! All data DCA distributions saved to " << f_out->GetName() << endl;
  f_out->Close();
}