#include <iostream>
#include <limits>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TStyle.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <THStack.h>
#include <TROOT.h>
#include <TDirectory.h>
#include <TString.h>
#include <TSystem.h>
#include <TParameter.h>

using namespace std;

//===============================================
// Input/output files
//===============================================
TFile *promptMC_DCA = TFile::Open("Hist_DCA_Prompt_out_combined_Sept6.root");
TFile *nonpromptMC_DCA = TFile::Open("Hist_DCA_NonPrompt_out_combined_Sept6.root");
TFile *data_DCA = TFile::Open("HIST_DATA_DCA_Aug26.root");

TFile *f_out = new TFile("Output_Chi2_vs_SF_AllBins_Sept6.root", "RECREATE");

const char *kPlotDir = "plots_Chi2_vs_SF_v2";

const double kDcaCut = 0.01;  // cm
const double kPullMax = 5.0;
//===============================================

const int bins = 10;
Double_t DCA[bins + 1] = {0, 0.0011, 0.0024, 0.0039, 0.0059, 0.0085, 0.01179, 0.016, 0.0214, 0.0366, 0.135};

// Global pointers for the fit function
TH1D *h_prompt_MC = nullptr;
TH1D *h_nonprompt_MC = nullptr;

Double_t ftotal(Double_t *x, Double_t *par)
{
	Int_t bin = h_nonprompt_MC->GetXaxis()->FindBin(x[0]);
	if (bin < 1) bin = 1;
	if (bin > h_nonprompt_MC->GetNbinsX()) bin = h_nonprompt_MC->GetNbinsX();

	Double_t prompt_norm = h_prompt_MC->Integral("width");
	Double_t nonprompt_norm = h_nonprompt_MC->Integral("width");

	Double_t sr = 0.0;
	Double_t br = 0.0;

	if (prompt_norm > 0.0)
		sr = par[1] * par[0] * (h_prompt_MC->GetBinContent(bin) / prompt_norm);

	if (nonprompt_norm > 0.0)
		br = par[1] * (1.0 - par[0]) * (h_nonprompt_MC->GetBinContent(bin) / nonprompt_norm);

	return sr + br;
}


bool RequireFile(TFile *f, const char *tag)
{
	if (!f || f->IsZombie())
	{
		cerr << "ERROR: could not open " << tag << " -- check the filename." << endl;
		return false;
	}
	return true;
}

TF1 *fitHist(TH1D *h, const char *name, bool quiet, bool &ok)
{
	double minRange = h->GetXaxis()->GetXmin();
	double maxRange = h->GetXaxis()->GetXmax();
	double yield = h->Integral("width");

	TF1 *fitFunc = new TF1(name, ftotal, minRange, maxRange, 2);
	fitFunc->SetLineColor(kRed);
	fitFunc->SetLineWidth(3);
	fitFunc->SetNpx(2000);

	fitFunc->SetParameter(0, 0.8);
	fitFunc->SetParameter(1, yield);
	fitFunc->SetParLimits(0, 0.0, 1.0);
	fitFunc->SetParLimits(1, 0.0, 2.0 * yield);

	TString opt = quiet ? "SRIQN" : "SERIN";
	TFitResultPtr res = h->Fit(fitFunc, opt.Data(), "", minRange, maxRange);

	const Int_t ndf = fitFunc->GetNDF();
	ok = res.Get() && res->Status() == 0 && ndf > 0;
	return fitFunc;
}

void BinWidthNormalization(TH1 *hist)
{
	Int_t numBins = hist->GetNbinsX();
	for (Int_t i = 1; i <= numBins; ++i)
	{
		hist->SetBinContent(i, hist->GetBinContent(i) / hist->GetBinWidth(i));
		hist->SetBinError(i, hist->GetBinError(i) / hist->GetBinWidth(i));
	}
}


int getBestSFIndex(int size, const Double_t *chi2, const bool *valid)
{
	int minIndex = -1;
	Double_t minValue = std::numeric_limits<Double_t>::infinity();
	for (int i = 0; i < size; ++i)
	{
		if (valid[i] && chi2[i] < minValue)
		{
			minValue = chi2[i];
			minIndex = i;
		}
	}
	return minIndex;
}


bool AtParLimit(TF1 *f, int ipar, double tol = 1e-4)
{
	double lo = 0.0, hi = 0.0;
	f->GetParLimits(ipar, lo, hi);
	if (lo == hi)
		return false; // parameter isn't limited
	const double val = f->GetParameter(ipar);
	return (val - lo < tol) || (hi - val < tol);
}


double PositiveFloor(TH1 *h, double fallback)
{
	double floorVal = std::numeric_limits<double>::infinity();
	for (int b = 1; b <= h->GetNbinsX(); ++b)
	{
		double v = h->GetBinContent(b);
		if (v > 0)
			floorVal = TMath::Min(floorVal, v);
	}
	return std::isfinite(floorVal) ? floorVal : fallback;
}


double FractionBelow(TH1 *h, double xCut)
{
	double total = 0.0, below = 0.0;
	for (int b = 1; b <= h->GetNbinsX(); ++b)
	{
		const double lo = h->GetBinLowEdge(b);
		const double w = h->GetBinWidth(b);
		const double area = h->GetBinContent(b) * w; // density * width = yield
		total += area;
		if (lo + w <= xCut)
			below += area;
		else if (lo < xCut)
			below += h->GetBinContent(b) * (xCut - lo);
	}
	return (total > 0.0) ? below / total : 0.0;
}

// Prompt fraction among candidates with DCA < xCut, derived from the whole-range
// fitted prompt fraction (fWhole) and the shapes of the unit-normalised prompt /
// non-prompt templates:
//   f(<xCut) = fWhole*P(<xCut) / [ fWhole*P(<xCut) + (1-fWhole)*NP(<xCut) ]
// The overall normalisation par[1] cancels. Only the fWhole uncertainty is
// propagated (df/dfWhole = P*NP / denom^2); template statistical uncertainty on
// P/NP is neglected, consistent with how the fit itself treats the templates.
double PromptFractionBelow(double fWhole, TH1 *hPrompt, TH1 *hNonprompt,
                           double xCut, double fWholeErr, double &fCutErr)
{
	const double pB = FractionBelow(hPrompt, xCut);
	const double npB = FractionBelow(hNonprompt, xCut);
	const double denom = fWhole * pB + (1.0 - fWhole) * npB;
	if (denom <= 0.0)
	{
		fCutErr = 0.0;
		return 0.0;
	}
	fCutErr = std::fabs(pB * npB / (denom * denom)) * fWholeErr;
	return (fWhole * pB) / denom;
}

TCanvas *DrawTemplateFit(const char *cname,
                         TH1D *h_data,
						 TF1 *fit_func,
                         TH1D *h_prompt,
						 TH1D *h_nonprompt,
                         Double_t sf_val,
						 Double_t chi2ndf,
                         const char *label_cent,
						 const char *label_pt,
                         int n_dca_bins,
						 Double_t *dca_edges,
                         Double_t cmp_sf = -1.0,
						 Double_t cmp_chi2ndf = -1.0)
{
    TCanvas *c = new TCanvas(cname, cname, 0, 53, 900, 800);
    c->cd();

    TPad *pad1 = new TPad(Form("pad1_%s", cname), "", 0, 0.25, 1, 1);
    pad1->SetTopMargin(0.08); pad1->SetBottomMargin(0.0); pad1->SetLogy(1);
    pad1->Draw(); pad1->cd();

    TH1D *p_sc  = (TH1D *)h_prompt->Clone(Form("p_sc_%s", cname));
    TH1D *np_sc = (TH1D *)h_nonprompt->Clone(Form("np_sc_%s", cname));
    p_sc ->Scale((fit_func->GetParameter(0)* fit_func->GetParameter(1)) / p_sc ->Integral("width"));
    np_sc->Scale(((1.0 - fit_func->GetParameter(0)) * fit_func->GetParameter(1)) / np_sc->Integral("width"));

    TH1D *hd = (TH1D *)h_data->Clone(Form("hd_%s", cname));
    hd->SetStats(0); hd->SetMarkerStyle(20); hd->SetMarkerSize(1.3); hd->SetLineColor(1);
    hd->GetXaxis()->SetLabelSize(0); hd->GetXaxis()->SetTitleSize(0);
    hd->GetYaxis()->SetTitle("dN/dDCA (cm^{-1})"); hd->GetYaxis()->CenterTitle(true);
    hd->GetYaxis()->SetLabelFont(42); hd->GetYaxis()->SetTitleFont(132);
    hd->GetYaxis()->SetTitleSize(0.06); hd->GetYaxis()->SetTitleOffset(0.9);

    // v2: size the log-scale range to what's actually on the canvas (data,
    // and the stacked prompt+non-prompt templates) instead of a fixed
    // "*100 / floor at 1.0" window that could clip real points off-scale.
    double stackMax = 0.0;
    for (int b = 1; b <= p_sc->GetNbinsX(); ++b)
        stackMax = TMath::Max(stackMax, p_sc->GetBinContent(b) + np_sc->GetBinContent(b));
    double plotMax = TMath::Max(hd->GetMaximum(), stackMax);
    double plotMin = PositiveFloor(hd, 1e-3);
    hd->SetMaximum(plotMax * 50.0);
    hd->SetMinimum(plotMin * 0.3);
    hd->Draw("E1");

    p_sc->SetStats(0);  p_sc->SetFillColor(42);  p_sc->SetFillStyle(1001);  p_sc->SetLineColor(1); p_sc->SetLineWidth(1); p_sc->SetMarkerStyle(0);
    np_sc->SetStats(0); np_sc->SetFillColor(51); np_sc->SetFillStyle(1001); np_sc->SetLineColor(1); np_sc->SetLineWidth(1); np_sc->SetMarkerStyle(0);
    THStack *hs = new THStack(Form("hs_%s", cname), "");
    hs->Add(np_sc); hs->Add(p_sc); hs->Draw("HIST SAME");
    fit_func->SetLineColor(kRed); fit_func->SetLineWidth(3); fit_func->Draw("L SAME");
    hd->Draw("E1 SAME");

    TLatex *tcms = new TLatex(0.1, 0.93, "CMS #it{Preliminary}");
    tcms->SetNDC(); tcms->SetTextFont(42); tcms->SetTextSize(0.05); tcms->Draw();
    TLatex *tlumi = new TLatex(0.9, 0.93, "PbPb 5.36 TeV");
    tlumi->SetNDC(); tlumi->SetTextAlign(31); tlumi->SetTextFont(42); tlumi->SetTextSize(0.05); tlumi->Draw();

    double fCutErr = 0.0;
    const double fCut = PromptFractionBelow(fit_func->GetParameter(0), h_prompt, h_nonprompt,
                                            kDcaCut, fit_func->GetParError(0), fCutErr);

    TLatex *ti = new TLatex();
    ti->SetNDC(); ti->SetTextAlign(12); ti->SetTextFont(42); ti->SetTextSize(0.038);
    ti->DrawLatex(0.15, 0.86, Form("Cent: %s%%, p_{T}: %s GeV/c", label_cent, label_pt));
    ti->DrawLatex(0.15, 0.80, Form("f_{prompt}(whole) = %.3f #pm %.3f", fit_func->GetParameter(0), fit_func->GetParError(0)));
    ti->DrawLatex(0.15, 0.74, Form("f_{prompt}(DCA < %.4g) = %.3f #pm %.3f", kDcaCut, fCut, fCutErr));
    if (cmp_sf > 0)
        ti->DrawLatex(0.15, 0.68, Form("#chi^{2}/NDF = %.2f  (best SF=%.2f: %.2f)", chi2ndf, cmp_sf, cmp_chi2ndf));
    else
        ti->DrawLatex(0.15, 0.68, Form("#chi^{2}/NDF = %.2f  (SF = %.2f)", chi2ndf, sf_val));

    TLegend *leg = new TLegend(0.66, 0.66, 0.90, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    TLegendEntry *e;
    e = leg->AddEntry(hd, "Data", "pe"); e->SetLineColor(kBlack); e->SetMarkerColor(kBlack); e->SetMarkerStyle(20); e->SetMarkerSize(1.3); e->SetTextFont(42);
    e = leg->AddEntry(p_sc,  "MC Prompt D^{0}",    "f"); e->SetFillColor(42); e->SetFillStyle(1001); e->SetLineColor(1); e->SetTextFont(42);
    e = leg->AddEntry(np_sc, "MC NonPrompt D^{0}", "f"); e->SetFillColor(51); e->SetFillStyle(1001); e->SetLineColor(1); e->SetTextFont(42);
    TString fit_label = (cmp_sf > 0) ? Form("Total Fit (SF=%.2f)", sf_val) : "Total Fit";
    e = leg->AddEntry(fit_func, fit_label.Data(), "l"); e->SetLineColor(kRed); e->SetLineWidth(3); e->SetTextFont(42);
    leg->Draw();

    c->cd();
    TPad *pad2 = new TPad(Form("pad2_%s", cname), "", 0, 0, 1, 0.25);
    pad2->SetTopMargin(0.0); pad2->SetBottomMargin(0.35); pad2->SetGridy(); pad2->Draw(); pad2->cd();

    // Pull = (data - model) / sqrt(data_err^2 + prompt_err^2)
    TH1D *pred = (TH1D *)p_sc->Clone(Form("pred_%s", cname));
    pred->Add(np_sc);

    TH1D *pull = new TH1D(Form("pull_%s", cname), "", n_dca_bins, dca_edges);
    double maxAbsPull = 0.0;
    for (int k = 1; k <= pull->GetNbinsX(); k++) {
        double dc = h_data->GetBinContent(k), de = h_data->GetBinError(k);
        double mc = pred->GetBinContent(k);
        double pe = p_sc->GetBinError(k); // prompt-template error only 
        double sig = std::sqrt(de * de + pe * pe);
        if (sig <= 0.0) continue;
        double pv = (dc - mc) / sig;
        pull->SetBinContent(k, pv);
        pull->SetBinError(k, 0.0);
        maxAbsPull = TMath::Max(maxAbsPull, std::fabs(pv));
    }
    double pullAxis = TMath::Max(kPullMax, 1.15 * maxAbsPull);

    pull->SetStats(0); pull->SetMarkerStyle(20); pull->SetMarkerSize(1.2);
    pull->SetMarkerColor(1); pull->SetLineColor(1);
    pull->GetXaxis()->SetTitle("DCA (cm)"); pull->GetXaxis()->CenterTitle(true);
    pull->GetXaxis()->SetLabelFont(132); pull->GetXaxis()->SetTitleFont(132);
    pull->GetXaxis()->SetLabelSize(0.12); pull->GetXaxis()->SetTitleSize(0.15); pull->GetXaxis()->SetTitleOffset(1.0);
    pull->GetYaxis()->SetTitle("Pull"); pull->GetYaxis()->CenterTitle(true); pull->GetYaxis()->SetNdivisions(505);
    pull->GetYaxis()->SetLabelFont(132); pull->GetYaxis()->SetTitleFont(132);
    pull->GetYaxis()->SetLabelSize(0.1); pull->GetYaxis()->SetTitleSize(0.13); pull->GetYaxis()->SetTitleOffset(0.4);
    pull->GetYaxis()->SetRangeUser(-pullAxis, pullAxis);
    pull->Draw("P");

    TLine *zero = new TLine(pull->GetXaxis()->GetXmin(), 0.0, pull->GetXaxis()->GetXmax(), 0.0);
    zero->SetLineColor(kRed); zero->SetLineStyle(2); zero->Draw();
    for (double s : {-3.0, 3.0}) {
        if (std::fabs(s) < pullAxis) {
            TLine *g = new TLine(pull->GetXaxis()->GetXmin(), s, pull->GetXaxis()->GetXmax(), s);
            g->SetLineColor(kGray + 1); g->SetLineStyle(3); g->Draw();
        }
    }
    return c;
}

void Chi2_vs_SF_v2()
{
	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	TH1::SetDefaultSumw2(true);

	// v2: fail loudly and immediately on a bad filename instead of
	// segfaulting deep inside the (cent,pt) loop.
	if (!RequireFile(promptMC_DCA, "prompt MC DCA file") ||
	    !RequireFile(nonpromptMC_DCA, "non-prompt MC DCA file") ||
	    !RequireFile(data_DCA, "data DCA file") ||
	    !RequireFile(f_out, "output ROOT file"))
		return;

	gSystem->mkdir(kPlotDir, kTRUE);

	bool binWidthNorm = false; // Because the input DCA is already bin normalised!

	const Int_t N_CENTBIN = 5;
	const Int_t N_PTBIN = 9;

	const int cent_min = 0;
	const int cent_max = N_CENTBIN;
	const int ptbin_min = 0;
	const int ptbin_max = N_PTBIN;

	const char *cent[N_CENTBIN] = {"cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"};
	const char *pt[N_PTBIN] = {"pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8", "pT8to10", "pT10to15", "pT15to30", "pT30to100"};

	
	const Double_t PT_EDGES[N_PTBIN + 1] = {2, 3, 4, 5, 6, 8, 10, 15, 30, 100};

	const Int_t numSF = 101;
	Double_t ScaleFactor[numSF] = {
		0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59,
		0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69,
		0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79,
		0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89,
		0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
		1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09,
		1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19,
		1.20, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29,
		1.30, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39,
		1.40, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49,
		1.50};

	const char *sf[numSF] = {
		"x0p5", "x0p51", "x0p52", "x0p53", "x0p54", "x0p55", "x0p56", "x0p57", "x0p58", "x0p59",
		"x0p6", "x0p61", "x0p62", "x0p63", "x0p64", "x0p65", "x0p66", "x0p67", "x0p68", "x0p69",
		"x0p7", "x0p71", "x0p72", "x0p73", "x0p74", "x0p75", "x0p76", "x0p77", "x0p78", "x0p79",
		"x0p8", "x0p81", "x0p82", "x0p83", "x0p84", "x0p85", "x0p86", "x0p87", "x0p88", "x0p89",
		"x0p9", "x0p91", "x0p92", "x0p93", "x0p94", "x0p95", "x0p96", "x0p97", "x0p98", "x0p99",
		"x1p0", "x1p01", "x1p02", "x1p03", "x1p04", "x1p05", "x1p06", "x1p07", "x1p08", "x1p09",
		"x1p1", "x1p11", "x1p12", "x1p13", "x1p14", "x1p15", "x1p16", "x1p17", "x1p18", "x1p19",
		"x1p2", "x1p21", "x1p22", "x1p23", "x1p24", "x1p25", "x1p26", "x1p27", "x1p28", "x1p29",
		"x1p3", "x1p31", "x1p32", "x1p33", "x1p34", "x1p35", "x1p36", "x1p37", "x1p38", "x1p39",
		"x1p4", "x1p41", "x1p42", "x1p43", "x1p44", "x1p45", "x1p46", "x1p47", "x1p48", "x1p49",
		"x1p5"};
	const int sf_unity_idx = 50; // ScaleFactor[50] == 1.0

	int nBinsProcessed = 0, nBinsSkippedNoData = 0, nBinsSkippedNoValidSF = 0, nBinsNoUnityComparison = 0;
	int nBoundaryPegged = 0;
	std::vector<TString> boundaryLog; // bins whose best-SF prompt fraction sits at 0 or 1

	
	Double_t bestSF_result[N_CENTBIN][N_PTBIN];
	bool bestSF_valid[N_CENTBIN][N_PTBIN];
	bool bestSF_pegged[N_CENTBIN][N_PTBIN];
	// prompt fraction over the whole DCA range vs. restricted to DCA < kDcaCut
	Double_t pfWhole[N_CENTBIN][N_PTBIN], pfWholeErr[N_CENTBIN][N_PTBIN];
	Double_t pfCut[N_CENTBIN][N_PTBIN], pfCutErr[N_CENTBIN][N_PTBIN];
	for (int i = 0; i < N_CENTBIN; i++)
		for (int j = 0; j < N_PTBIN; j++)
			bestSF_valid[i][j] = bestSF_pegged[i][j] = false;

	for (int icent = cent_min; icent < cent_max; icent++)
	{
		for (int ipt = ptbin_min; ipt < ptbin_max; ipt++)
		{

			cout << "\n==============================================" << endl;
			cout << "Processing Bin: " << cent[icent] << " | " << pt[ipt] << endl;

			TH1D *h_data = (TH1D *)data_DCA->Get(Form("hdca_%s_%s", cent[icent], pt[ipt]));
			if (!h_data)
			{
				cout << "Skipping missing data histogram: hdca_" << cent[icent] << "_" << pt[ipt] << endl;
				nBinsSkippedNoData++;
				continue;
			}
			h_data->SetNameTitle(Form("data_dca_%d_%d", icent, ipt), "");
			if (binWidthNorm)
				BinWidthNormalization(h_data);
			h_data->GetXaxis()->SetRangeUser(0.0, DCA[bins]);

			Double_t chi2s_current[numSF];
			bool sfValid[numSF];
			TH1D *h_prompt_arr[numSF];
			TH1D *h_nonprompt_arr[numSF];

			// ---- Stage 1: fast, quiet scan over all scale factors ----
			for (int isf = 0; isf < numSF; isf++)
			{
				h_prompt_arr[isf] = (TH1D *)promptMC_DCA->Get(Form("hdca_%s_%s_%s", cent[icent], pt[ipt], sf[isf]));
				h_nonprompt_arr[isf] = (TH1D *)nonpromptMC_DCA->Get(Form("hdca_%s_%s_%s", cent[icent], pt[ipt], sf[isf]));

				if (!h_prompt_arr[isf] || !h_nonprompt_arr[isf])
				{
					chi2s_current[isf] = std::numeric_limits<Double_t>::infinity();
					sfValid[isf] = false;
					continue;
				}

				if (binWidthNorm)
				{
					BinWidthNormalization(h_prompt_arr[isf]);
					BinWidthNormalization(h_nonprompt_arr[isf]);
				}

				Double_t prompt_norm = h_prompt_arr[isf]->Integral("width");
				Double_t nonprompt_norm = h_nonprompt_arr[isf]->Integral("width");

				if (prompt_norm > 0.0)
					h_prompt_arr[isf]->Scale(1.0 / prompt_norm);

				if (nonprompt_norm > 0.0)
					h_nonprompt_arr[isf]->Scale(1.0 / nonprompt_norm);

				h_prompt_arr[isf]->GetXaxis()->SetRangeUser(0.0, DCA[bins]);
				h_nonprompt_arr[isf]->GetXaxis()->SetRangeUser(0.0, DCA[bins]);

				h_prompt_MC = h_prompt_arr[isf];
				h_nonprompt_MC = h_nonprompt_arr[isf];

				bool ok = false;
				TF1 *scanFit = fitHist(h_data, Form("fitFuncScan_%d_%d_%d", icent, ipt, isf), /*quiet=*/true, ok);
				chi2s_current[isf] = ok ? scanFit->GetChisquare() / scanFit->GetNDF()
				                        : std::numeric_limits<Double_t>::infinity();
				sfValid[isf] = ok;
				delete scanFit; // only the chi2/ndf from the scan is needed; the winner gets refit below
			}

			int best_isf = getBestSFIndex(numSF, chi2s_current, sfValid);

			if (best_isf < 0)
			{
				cout << "No valid scale factor found for " << cent[icent] << " " << pt[ipt] << endl;
				nBinsSkippedNoValidSF++;
				continue;
			}

			Double_t best_SF = ScaleFactor[best_isf];

			// ---- Stage 2: precise refit (Minos errors, verbose) of just the
			// winning SF -- this is the fit whose parameters/errors are
			// actually reported and plotted, so it's worth the extra cost. ----
			h_prompt_MC = h_prompt_arr[best_isf];
			h_nonprompt_MC = h_nonprompt_arr[best_isf];
			bool bestOk = false;
			TF1 *best_fit = fitHist(h_data, Form("fitFunc_best_%d_%d", icent, ipt), /*quiet=*/false, bestOk);
			if (bestOk)
				chi2s_current[best_isf] = best_fit->GetChisquare() / best_fit->GetNDF();

			// prompt fraction: whole range, and restricted to DCA < kDcaCut
			double pfCutErrBin = 0.0;
			const double pfCutBin = PromptFractionBelow(best_fit->GetParameter(0),
			                                            h_prompt_arr[best_isf], h_nonprompt_arr[best_isf],
			                                            kDcaCut, best_fit->GetParError(0), pfCutErrBin);
			pfWhole[icent][ipt] = best_fit->GetParameter(0);
			pfWholeErr[icent][ipt] = best_fit->GetParError(0);
			pfCut[icent][ipt] = pfCutBin;
			pfCutErr[icent][ipt] = pfCutErrBin;

			cout << ">>> Optimal Scale Factor: " << best_SF << " (chi2/NDF = " << chi2s_current[best_isf] << ")" << endl;
			cout << ">>> Prompt Fraction (whole)        : " << best_fit->GetParameter(0) << " +/- " << best_fit->GetParError(0) << endl;
			cout << ">>> Prompt Fraction (DCA < " << kDcaCut << ") : " << pfCutBin << " +/- " << pfCutErrBin << endl;

			const bool pegged = AtParLimit(best_fit, 0);
			bestSF_result[icent][ipt] = best_SF;
			bestSF_valid[icent][ipt] = true;
			bestSF_pegged[icent][ipt] = pegged;

			if (pegged)
			{
				nBoundaryPegged++;
				boundaryLog.push_back(Form("%-14s %-10s  fraction=%.3f (best SF=%.2f, chi2/NDF=%.1f) -- no template/SF combination fits this bin well",
				                            cent[icent], pt[ipt], best_fit->GetParameter(0), best_SF, chi2s_current[best_isf]));
				cout << "WARNING: prompt fraction pegged at its physical boundary for "
				     << cent[icent] << " " << pt[ipt] << " -- treat this bin's result with caution." << endl;
			}

			//=========================================
			// Draw template fit with best chi2/ndf
			//=========================================
			TString s_cent = cent[icent];
            s_cent.ReplaceAll("cent", "");
            s_cent.ReplaceAll("to", "-");
            TString s_pt = pt[ipt];
            s_pt.ReplaceAll("pT", "");
            s_pt.ReplaceAll("to", "-");

            TCanvas *c_best = DrawTemplateFit(
                Form("c_best_%d_%d", icent, ipt),
                h_data, best_fit,
                h_prompt_arr[best_isf], h_nonprompt_arr[best_isf],
                best_SF, chi2s_current[best_isf],
                s_cent.Data(), s_pt.Data(),
                bins, DCA);

			// ==========================================
            // Unscaled (SF=1.0) comparison canvas -- also refit precisely
            // (or reuse best_fit if SF=1.0 happens to be the winner) since
            // its prompt fraction and error are quoted on the plot too.
            // ==========================================
            TCanvas *c_unscaled = nullptr;
            TF1 *unity_fit = nullptr;
            if (sfValid[sf_unity_idx])
            {
                if (sf_unity_idx == best_isf)
                {
                    unity_fit = best_fit;
                }
                else
                {
                    h_prompt_MC = h_prompt_arr[sf_unity_idx];
                    h_nonprompt_MC = h_nonprompt_arr[sf_unity_idx];
                    bool unityOk = false;
                    unity_fit = fitHist(h_data, Form("fitFunc_unity_%d_%d", icent, ipt), /*quiet=*/false, unityOk);
                    if (unityOk)
                        chi2s_current[sf_unity_idx] = unity_fit->GetChisquare() / unity_fit->GetNDF();
                }
            }
            if (unity_fit)
                c_unscaled = DrawTemplateFit(
                    Form("c_unscaled_%d_%d", icent, ipt),
                    h_data, unity_fit,
                    h_prompt_arr[sf_unity_idx], h_nonprompt_arr[sf_unity_idx],
                    1.0, chi2s_current[sf_unity_idx],
                    s_cent.Data(), s_pt.Data(),
                    bins, DCA,
                    best_SF, chi2s_current[best_isf]);
            else
                nBinsNoUnityComparison++;


			// ==========================================
			// Draw the Chi2 vs SF
			// ==========================================
			TCanvas *c_sf = new TCanvas(Form("c_chi2_%d_%d", icent, ipt), Form("Chi2 vs SF %s %s", cent[icent], pt[ipt]), 900, 700);
			c_sf->SetLeftMargin(0.12);
			c_sf->SetRightMargin(0.05);
			c_sf->SetTopMargin(0.1);
			c_sf->SetBottomMargin(0.12);
			c_sf->cd();

			// v2: only plot SF points with a usable fit. Feeding the
			// invalid-SF sentinel value into the TGraph (9999 previously,
			// or an outright infinity now) would swamp the y-axis and
			// squash the actual chi2/ndf curve into a flat line.
			std::vector<Double_t> validSFvals, validChi2vals;
			for (int isf = 0; isf < numSF; isf++)
			{
				if (sfValid[isf])
				{
					validSFvals.push_back(ScaleFactor[isf]);
					validChi2vals.push_back(chi2s_current[isf]);
				}
			}

			TGraph *graph = new TGraph((int)validSFvals.size(), validSFvals.data(), validChi2vals.data());
			graph->SetName(Form("graph_chi2_%s_%s", cent[icent], pt[ipt]));
			graph->SetTitle("; Scale Factor; #chi^{2}/NDF");
			graph->SetMarkerStyle(20);
			graph->SetMarkerSize(0.9);
			graph->SetMarkerColor(kBlue + 1);
			graph->SetLineColor(kBlue + 1);
			graph->SetLineWidth(1);
			graph->GetXaxis()->SetLabelFont(42);
			graph->GetXaxis()->SetTitleFont(42);
			graph->GetXaxis()->SetTitleSize(0.05);
			graph->GetYaxis()->SetLabelFont(42);
			graph->GetYaxis()->SetTitleFont(42);
			graph->GetYaxis()->SetTitleSize(0.05);
			graph->GetYaxis()->SetTitleOffset(1.1);
			graph->Draw("AP");

			// Mark minimum point
			TGraph *gmin = new TGraph(1, &best_SF, &chi2s_current[best_isf]);
			gmin->SetMarkerStyle(29);
			gmin->SetMarkerColor(kRed);
			gmin->SetMarkerSize(2.5);
			gmin->Draw("P SAME");

			TLatex *tex_cms2 = new TLatex(0.12, 0.93, "CMS #it{Preliminary}");
			tex_cms2->SetNDC();
			tex_cms2->SetTextFont(42);
			tex_cms2->SetTextSize(0.05);
			tex_cms2->Draw();

			TLatex *tex_lumi2 = new TLatex(0.95, 0.93, "PbPb 5.36 TeV");
			tex_lumi2->SetNDC();
			tex_lumi2->SetTextAlign(31);
			tex_lumi2->SetTextFont(42);
			tex_lumi2->SetTextSize(0.05);
			tex_lumi2->Draw();

			TLatex *tex_info2 = new TLatex();
			tex_info2->SetNDC();
			tex_info2->SetTextFont(42);
			tex_info2->SetTextSize(0.045);
			tex_info2->SetTextAlign(12);
			tex_info2->DrawLatex(0.15, 0.85, Form("Cent: %s%%, p_{T}: %s GeV/c", s_cent.Data(), s_pt.Data()));
			tex_info2->DrawLatex(0.15, 0.79, Form("Best SF = %.2f,  #chi^{2}/NDF = %.2f", best_SF, chi2s_current[best_isf]));

			// =========================================================
			// Draw Shape Overlay Plot (Without Fit, Area Normalized)
			// =========================================================
			TCanvas *c_overlay = new TCanvas(Form("c_overlay_%d_%d", icent, ipt), Form("Overlay Shape %s %s", cent[icent], pt[ipt]), 900, 800);
			c_overlay->SetLeftMargin(0.12);
			c_overlay->SetRightMargin(0.05);
			c_overlay->SetTopMargin(0.1);
			c_overlay->SetBottomMargin(0.12);
			c_overlay->SetLogy();
			c_overlay->cd();

			TH1D *h_data_shape = (TH1D *)h_data->Clone(Form("data_shape_%d_%d", icent, ipt));
			TH1D *p_MC_shape = (TH1D *)h_prompt_arr[best_isf]->Clone(Form("prompt_shape_%d_%d", icent, ipt));
			TH1D *np_MC_shape = (TH1D *)h_nonprompt_arr[best_isf]->Clone(Form("nonprompt_shape_%d_%d", icent, ipt));

			if (h_data_shape->Integral("width") > 0)
				h_data_shape->Scale(1.0 / h_data_shape->Integral("width"));
			if (p_MC_shape->Integral("width") > 0)
				p_MC_shape->Scale(1.0 / p_MC_shape->Integral("width"));
			if (np_MC_shape->Integral("width") > 0)
				np_MC_shape->Scale(1.0 / np_MC_shape->Integral("width"));

			h_data_shape->SetStats(0); // EXPLICITLY TURN OFF STATS
			h_data_shape->GetYaxis()->SetTitle("Normalised entries");
			h_data_shape->GetXaxis()->SetTitle("DCA (cm)");
			h_data_shape->GetXaxis()->SetLabelSize(0.04);
			h_data_shape->GetXaxis()->SetTitleSize(0.05);
			h_data_shape->GetYaxis()->SetLabelSize(0.04);
			h_data_shape->GetYaxis()->SetTitleSize(0.05);
			h_data_shape->GetYaxis()->SetTitleOffset(1.2);
			h_data_shape->SetMarkerStyle(20);
			h_data_shape->SetMarkerColor(kBlack);
			h_data_shape->SetLineColor(kBlack);

			p_MC_shape->SetStats(0);
			p_MC_shape->SetLineColor(kBlue);
			p_MC_shape->SetLineWidth(2);
			p_MC_shape->SetMarkerStyle(0);
			p_MC_shape->SetFillStyle(0);

			np_MC_shape->SetStats(0);
			np_MC_shape->SetLineColor(kRed);
			np_MC_shape->SetLineWidth(2);
			np_MC_shape->SetMarkerStyle(0);
			np_MC_shape->SetFillStyle(0);

			double MaxValShape = TMath::Max(h_data_shape->GetMaximum(), TMath::Max(p_MC_shape->GetMaximum(), np_MC_shape->GetMaximum()));
			double MinValShape = TMath::Min(PositiveFloor(h_data_shape, 1e-3),
			                                 TMath::Min(PositiveFloor(p_MC_shape, 1e-3), PositiveFloor(np_MC_shape, 1e-3)));
			h_data_shape->SetMaximum(MaxValShape * 5.0);
			h_data_shape->SetMinimum(MinValShape * 0.3);
			h_data_shape->GetXaxis()->SetLabelFont(42);
			h_data_shape->GetXaxis()->SetTitleFont(42);
			h_data_shape->GetYaxis()->SetLabelFont(42);
			h_data_shape->GetYaxis()->SetTitleFont(42);

			h_data_shape->Draw("E1");
			p_MC_shape->Draw("HIST SAME");
			np_MC_shape->Draw("HIST SAME");

			TLatex *tex_cms_ov = new TLatex(0.12, 0.93, "CMS #it{Preliminary}");
			tex_cms_ov->SetNDC();
			tex_cms_ov->SetTextFont(42);
			tex_cms_ov->SetTextSize(0.05);
			tex_cms_ov->Draw();

			TLatex *tex_lumi_ov = new TLatex(0.95, 0.93, "PbPb 5.36 TeV");
			tex_lumi_ov->SetNDC();
			tex_lumi_ov->SetTextAlign(31);
			tex_lumi_ov->SetTextFont(42);
			tex_lumi_ov->SetTextSize(0.05);
			tex_lumi_ov->Draw();

			TLatex *tex_info_ov = new TLatex();
			tex_info_ov->SetNDC();
			tex_info_ov->SetTextFont(42);
			tex_info_ov->SetTextSize(0.045);
			tex_info_ov->SetTextAlign(12);
			tex_info_ov->DrawLatex(0.15, 0.85, Form("Cent: %s%%", s_cent.Data()));
			tex_info_ov->DrawLatex(0.15, 0.79, Form("p_{T}: %s GeV/c", s_pt.Data()));

			TLegend *leg_ov = new TLegend(0.55, 0.65, 0.88, 0.88);
			leg_ov->SetBorderSize(0);
			leg_ov->SetFillStyle(0);
			leg_ov->SetTextFont(42);
			leg_ov->SetTextSize(0.04);
			leg_ov->AddEntry(h_data_shape, "Data", "pe");
			leg_ov->AddEntry(p_MC_shape, "Prompt MC Shape", "l");
			leg_ov->AddEntry(np_MC_shape, "Non-Prompt MC Shape", "l");
			leg_ov->Draw();

			// ==========================================
			// WRITE TO FILE (ROOT file, plus flat PNGs for direct reuse)
			// ==========================================
			f_out->cd();
			TDirectory *dir = f_out->mkdir(Form("%s_%s", cent[icent], pt[ipt]));
			if (!dir)
			{
				cout << "WARNING: could not create output directory for "
				     << cent[icent] << " " << pt[ipt] << " -- skipping write." << endl;
				continue;
			}
			dir->cd();
            c_overlay->Write("Canvas_Overlay_Shape");
			c_sf->Write("Canvas_Chi2_vs_SF");
			if (c_unscaled)
				c_unscaled->Write("Canvas_Unscaled_Fit");
			c_best->Write("Canvas_BestFit");

			// scalar results for this (cent,pt) bin
			TParameter<double>("DCAcut", kDcaCut).Write();
			TParameter<double>("BestSF", best_SF).Write();
			TParameter<double>("Chi2NDF", chi2s_current[best_isf]).Write();
			TParameter<double>("PromptFraction_Whole", pfWhole[icent][ipt]).Write();
			TParameter<double>("PromptFraction_Whole_Err", pfWholeErr[icent][ipt]).Write();
			TParameter<double>("PromptFraction_DCAcut", pfCut[icent][ipt]).Write();
			TParameter<double>("PromptFraction_DCAcut_Err", pfCutErr[icent][ipt]).Write();

			TString pngTag = Form("%s/%s_%s", kPlotDir, cent[icent], pt[ipt]);
			c_best->SaveAs(pngTag + "_BestFit.png");
			if (c_unscaled)
				c_unscaled->SaveAs(pngTag + "_UnscaledFit.png");
			c_sf->SaveAs(pngTag + "_Chi2vsSF.png");
			c_overlay->SaveAs(pngTag + "_OverlayShape.png");

			nBinsProcessed++;
		}
	}

	// ==================================================================
	// Summary plot: best scale factor vs pT, one series per centrality bin.
	// Lets the SF trend/consistency across bins be checked at a glance
	// instead of digging through all 45 individual chi2-vs-SF canvases.
	// Points whose prompt fraction hit its [0,1] boundary are flagged with
	// an open black marker overlay -- their "best SF" comes from a fit that
	// could not describe the data well at any SF, so it's a number you have
	// but shouldn't lean on the same way as the others.
	// ==================================================================
	{
		Double_t ptCenter[N_PTBIN], ptHalfWidth[N_PTBIN];
		for (int j = 0; j < N_PTBIN; j++)
		{
			ptCenter[j] = 0.5 * (PT_EDGES[j] + PT_EDGES[j + 1]);
			ptHalfWidth[j] = 0.5 * (PT_EDGES[j + 1] - PT_EDGES[j]);
		}

		TCanvas *c_bestSF = new TCanvas("c_bestSF_vs_pT", "Best SF vs pT", 900, 750);
		c_bestSF->SetLeftMargin(0.12);
		c_bestSF->SetRightMargin(0.05);
		c_bestSF->SetTopMargin(0.1);
		c_bestSF->SetBottomMargin(0.12);
		c_bestSF->SetLogx();
		c_bestSF->cd();

		const Int_t sfColor[N_CENTBIN] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 2};
		const Int_t sfMarker[N_CENTBIN] = {20, 21, 22, 23, 33};
		const char *centLegendLabel[N_CENTBIN] = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%"};

		TMultiGraph *mg = new TMultiGraph();
		mg->SetName("mg_bestSF_vs_pT");

		TLegend *legSF = new TLegend(0.47, 0.60, 0.90, 0.88);
		legSF->SetBorderSize(0);
		legSF->SetFillStyle(0);
		legSF->SetTextFont(42);
		legSF->SetTextSize(0.032);

		std::vector<Double_t> peggedX, peggedY;

		for (int icent = 0; icent < N_CENTBIN; icent++)
		{
			std::vector<Double_t> xv, xe, yv, ye;
			for (int ipt = 0; ipt < N_PTBIN; ipt++)
			{
				if (!bestSF_valid[icent][ipt])
					continue;
				xv.push_back(ptCenter[ipt]);
				xe.push_back(ptHalfWidth[ipt]);
				yv.push_back(bestSF_result[icent][ipt]);
				ye.push_back(0.0);
				if (bestSF_pegged[icent][ipt])
				{
					peggedX.push_back(ptCenter[ipt]);
					peggedY.push_back(bestSF_result[icent][ipt]);
				}
			}
			if (xv.empty())
				continue;

			TGraphErrors *g = new TGraphErrors((int)xv.size(), xv.data(), yv.data(), xe.data(), ye.data());
			g->SetName(Form("gr_bestSF_%s", cent[icent]));
			g->SetMarkerStyle(sfMarker[icent]);
			g->SetMarkerColor(sfColor[icent]);
			g->SetLineColor(sfColor[icent]);
			g->SetMarkerSize(1.4);
			g->SetLineWidth(2);
			mg->Add(g, "PL");
			legSF->AddEntry(g, centLegendLabel[icent], "pl");
		}

		mg->SetTitle("; p_{T} (GeV/c); Best Scale Factor (SF)");
		mg->Draw("A");
		mg->GetXaxis()->SetLimits(PT_EDGES[0] * 0.8, PT_EDGES[N_PTBIN] * 1.2);
		mg->SetMinimum(0.45);
		mg->SetMaximum(1.55);
		mg->GetXaxis()->SetTitleFont(42); mg->GetXaxis()->SetLabelFont(42); mg->GetXaxis()->SetTitleSize(0.05);
		mg->GetYaxis()->SetTitleFont(42); mg->GetYaxis()->SetLabelFont(42); mg->GetYaxis()->SetTitleSize(0.05);
		mg->GetYaxis()->SetTitleOffset(1.1);
		mg->GetXaxis()->SetMoreLogLabels();
		gPad->Modified();
		gPad->Update();

		TLine *sfRef = new TLine(mg->GetXaxis()->GetXmin(), 1.0, mg->GetXaxis()->GetXmax(), 1.0);
		sfRef->SetLineStyle(2);
		sfRef->SetLineColor(kGray + 2);
		sfRef->Draw();

		if (!peggedX.empty())
		{
			TGraph *gFlag = new TGraph((int)peggedX.size(), peggedX.data(), peggedY.data());
			gFlag->SetMarkerStyle(24);
			gFlag->SetMarkerColor(kBlack);
			gFlag->SetMarkerSize(2.0);
			gFlag->Draw("P SAME");
			legSF->AddEntry(gFlag, "Fraction pegged at boundary", "p");
		}

		legSF->Draw();

		TLatex *tex_cms3 = new TLatex(0.12, 0.93, "CMS #it{Preliminary}");
		tex_cms3->SetNDC(); tex_cms3->SetTextFont(42); tex_cms3->SetTextSize(0.05); tex_cms3->Draw();
		TLatex *tex_lumi3 = new TLatex(0.95, 0.93, "PbPb 5.36 TeV");
		tex_lumi3->SetNDC(); tex_lumi3->SetTextAlign(31); tex_lumi3->SetTextFont(42); tex_lumi3->SetTextSize(0.05); tex_lumi3->Draw();

		f_out->cd();
		c_bestSF->Write("Canvas_BestSF_vs_pT");
		c_bestSF->SaveAs(Form("%s/Summary_BestSF_vs_pT.png", kPlotDir));
	}

	f_out->cd();
	f_out->Close();
	cout << "\n>>> Done! " << nBinsProcessed << " bins processed." << endl;
	cout << ">>> Skipped (no data histogram): " << nBinsSkippedNoData << endl;
	cout << ">>> Skipped (no valid scale factor / fit converged): " << nBinsSkippedNoValidSF << endl;
	cout << ">>> Bins without an SF=1.0 comparison plot: " << nBinsNoUnityComparison << endl;
	cout << "\n>>> Bins where the prompt fraction hit its [0,1] boundary (" << nBoundaryPegged << " total) -- "
	     << "these are not trustworthy prompt-fraction measurements as they stand:" << endl;
	for (const auto &line : boundaryLog)
		cout << "    " << line << endl;

	cout << "\n>>> Prompt fraction: whole DCA range vs. DCA < " << kDcaCut << " cm" << endl;
	cout << "    " << Form("%-12s %-10s %18s %18s %8s", "centrality", "pT",
	                       "f(whole)", "f(DCA<cut)", "ratio") << endl;
	for (int i = 0; i < N_CENTBIN; i++)
		for (int j = 0; j < N_PTBIN; j++)
		{
			if (!bestSF_valid[i][j])
				continue;
			const double ratio = (pfWhole[i][j] > 0.0) ? pfCut[i][j] / pfWhole[i][j] : 0.0;
			cout << "    " << Form("%-12s %-10s   %6.3f +/- %5.3f   %6.3f +/- %5.3f   %5.2f",
			                       cent[i], pt[j],
			                       pfWhole[i][j], pfWholeErr[i][j],
			                       pfCut[i][j], pfCutErr[i][j], ratio)
			     << endl;
		}

	cout << "\n>>> ROOT output :"<<f_out->GetName()<< endl;
	cout << ">>> PNG plots:   " << kPlotDir << "/" << endl;
}
