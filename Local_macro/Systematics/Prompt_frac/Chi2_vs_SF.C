#include <iostream>
#include <fstream>
#include <limits>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TROOT.h>
#include <TDirectory.h>
#include <TString.h>

using namespace std;

//===============================================
// Input output files
//===============================================
TFile *promptMC_DCA = TFile::Open("Hist_DCA_Prompt_out_combined_Sept10.root");
TFile *nonpromptMC_DCA = TFile::Open("Hist_DCA_NonPrompt_out_combined_Sept10.root");
TFile *data_DCA = TFile::Open("HIST_DATA_DCA_Sept10.root");

TFile *f_out = new TFile("Output_Chi2_vs_SF_AllBins_Sept10.root", "RECREATE");

const double dca_cut = 0.0085; // cm
// Per-bin f_prompt,1 (DCA < dca_cut), f_prompt,2 (whole DCA range), and the
// v_n unfolding factor (1-f_prompt,1)/(f_prompt,1-f_prompt,2), one text file
// row per (centrality, pT) bin.
const char *factorTxtName = "PromptFraction_vn_factor.txt";
// PLACEHOLDER for (v_n,1 - v_n,2), which isn't measured yet. Stands in just so
// the "factor" column can be turned into an actual v_n shift / systematic
// uncertainty estimate; replace with the real per-bin (v_n,1 - v_n,2) once
// available (at that point this should become an array, not one constant).
const double dvn_placeholder = 0.01;
//===============================================

// Two DCA binnings, matching the mass-fit producer (get_Data_DCA_fromMass.C):
// pT bins with lower edge < 5 GeV (pT2to3, pT3to4, pT4to5) use the coarser
// dca_bins_lo tail; every other pT bin uses dca_bins_hi. Both share every edge
// up to 0.0214 cm and the same first/last edge, so DCA_lo is an exact subset of
// DCA_hi (an 18-bin template can be Rebin()-ed onto DCA_lo losslessly).
const Int_t dca_bins_hi = 18;
Double_t DCA_hi[dca_bins_hi + 1] = {0,0.0005,0.0011,0.0014,0.0024,0.0029,0.0039,0.0045,0.0059,0.0067,0.0085,0.01179,0.016,0.0214,0.028,0.0366,0.0475,0.079,0.135};

const Int_t dca_bins_lo = 15;
Double_t DCA_lo[dca_bins_lo + 1] = {0,0.0005,0.0011,0.0014,0.0024,0.0029,0.0039,0.0045,0.0059,0.0067,0.0085,0.01179,0.016,0.0214,0.0366,0.135};

const Double_t DCA_MAX = 0.135; // common upper edge of both binnings

// pT lower edges, in the same order as pt[] in Chi2_vs_SF()
const Double_t PT_LOW_EDGE[9] = {2, 3, 4, 5, 6, 8, 10, 15, 30};
inline bool useLowPtBins(int ipt) { return PT_LOW_EDGE[ipt] < 5.0; }
inline int        dcaNbins(int ipt) { return useLowPtBins(ipt) ? dca_bins_lo : dca_bins_hi; }
inline Double_t  *dcaEdges(int ipt) { return useLowPtBins(ipt) ? DCA_lo : DCA_hi; }

// Rebin an 18-bin MC template onto this pT bin's DCA binning if needed (a no-op
// once the template producer writes the two binnings directly).
TH1D *MatchDcaBinning(TH1D *h, int ipt)
{
	if (!h || h->GetNbinsX() == dcaNbins(ipt))
		return h;
	return (TH1D *)h->Rebin(dcaNbins(ipt), Form("%s_rb", h->GetName()), dcaEdges(ipt));
}

// Global pointers for the fit function
TH1D *h_prompt_MC = nullptr;
TH1D *h_nonprompt_MC = nullptr;

Double_t ftotal(Double_t *x, Double_t *par)
{
	Int_t bin = h_nonprompt_MC->GetXaxis()->FindBin(x[0]);
	if (bin < 1) bin = 1;
	if (bin > h_nonprompt_MC->GetNbinsX()) bin = h_nonprompt_MC->GetNbinsX();

	//Double_t prompt_norm = h_prompt_MC->Integral("width");
	//Double_t nonprompt_norm = h_nonprompt_MC->Integral("width");

	Double_t prompt_norm = h_prompt_MC->Integral();
	Double_t nonprompt_norm = h_nonprompt_MC->Integral();
	

	Double_t sr = 0.0;
	Double_t br = 0.0;

	if (prompt_norm > 0.0)
		sr = par[1] * par[0] * (h_prompt_MC->GetBinContent(bin) / prompt_norm);

	if (nonprompt_norm > 0.0)
		br = par[1] * (1.0 - par[0]) * (h_nonprompt_MC->GetBinContent(bin) / nonprompt_norm);

	return sr + br;
}

TF1 *fitHist(TH1D *h, const char *name)
{
	double minRange = h->GetXaxis()->GetXmin();
	double maxRange = h->GetXaxis()->GetXmax();
	//double yield = h->Integral("width");
	double yield = h->Integral();

	TF1 *fitFunc = new TF1(name, ftotal, minRange, maxRange, 2);
	fitFunc->SetLineColor(kRed);
	fitFunc->SetLineWidth(3);
	fitFunc->SetNpx(2000);

	fitFunc->SetParameter(0, 0.8);
	fitFunc->SetParameter(1, yield);
	fitFunc->SetParLimits(0, 0.0, 1.0); // prompt fraction is physically bounded to [0,1]
	fitFunc->SetParLimits(1, 0.0, 2.0 * yield);

	//h->Fit(fitFunc, "SERIN", "", minRange, maxRange); 
	h->Fit(fitFunc, "E R 0 N", "", minRange, 0.08); 
	h->Fit(fitFunc, "E R 0 N", "", minRange, 0.08);
	h->Fit(fitFunc, "E R b N", "", minRange, 0.08);
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



// To find prompt fraction below DCA cut.
double FractionBelow(TH1D *h, double xCut)
{
	double total = 0.0, below = 0.0;
	for (int b = 1; b <= h->GetNbinsX(); ++b) {
		double lo = h->GetBinLowEdge(b), w = h->GetBinWidth(b);
		double area = h->GetBinContent(b) * w;
		total += area;
		if (lo + w <= xCut) below += area;
		else if (lo < xCut) below += h->GetBinContent(b) * (xCut - lo);
	}
	return (total > 0.0) ? below / total : 0.0;
}

//To find best SF based on chi2/ndf
int getBestSFIndex(int size, Double_t *Chi2_array)
{
	double minValue = 999999.0;
	int minIndex = 0;
	for (int i = 0; i < size; ++i)
	{
		if (Chi2_array[i] < minValue && Chi2_array[i] > 0)
		{
			minValue = Chi2_array[i];
			minIndex = i;
		}
	}
	return minIndex;
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
	p_sc->Sumw2();
	np_sc->Sumw2();
    //p_sc ->Scale((fit_func->GetParameter(0)* fit_func->GetParameter(1)) / p_sc ->Integral("width"));
    //np_sc->Scale(((1.0 - fit_func->GetParameter(0)) * fit_func->GetParameter(1)) / np_sc->Integral("width"));

	p_sc ->Scale((fit_func->GetParameter(0)* fit_func->GetParameter(1)) / p_sc ->Integral());
    np_sc->Scale(((1.0 - fit_func->GetParameter(0)) * fit_func->GetParameter(1)) / np_sc->Integral());

    // Y-axis floor: half the smallest nonzero non-prompt bin, so the non-prompt (and
    // therefore also prompt, since it never sits lower) distribution is never clipped,
    // whatever scale this particular bin happens to sit at.
    double npMin = 1e300;
    for (int b = 1; b <= np_sc->GetNbinsX(); ++b) {
        double v = np_sc->GetBinContent(b);
        if (v > 0.0 && v < npMin) npMin = v;
    }
    double yMin = (npMin < 1e300) ? npMin * 0.5 : 1e-5;

    TH1D *hd = (TH1D *)h_data->Clone(Form("hd_%s", cname));
    hd->SetStats(0); hd->SetMarkerStyle(20); hd->SetMarkerSize(1.3); hd->SetLineColor(1);
    hd->GetXaxis()->SetLabelSize(0); hd->GetXaxis()->SetTitleSize(0);
    hd->GetYaxis()->SetTitle("dN/dDCA (cm^{-1})"); hd->GetYaxis()->CenterTitle(true);
    hd->GetYaxis()->SetLabelFont(42); hd->GetYaxis()->SetTitleFont(132);
    hd->GetYaxis()->SetTitleSize(0.06); hd->GetYaxis()->SetTitleOffset(0.9);

    // Headroom above the peak has to scale with how many decades the (auto-computed)
    // floor sits below it, not a fixed multiplier -- otherwise a bin whose floor lands
    // many decades down (as yMin now can, since it adapts per bin) leaves the peak
    // sitting close to the top of the pad, colliding with the text block. Keeping the
    // peak at roughly 60% of the pad's log-height leaves the text block (top ~35%) clear
    // regardless of how many decades the floor spans.
    double peakVal = hd->GetMaximum();
    double decadesBelow = (yMin > 0.0 && peakVal > yMin) ? TMath::Log10(peakVal / yMin) : 8.0;
    double decadesAbove = decadesBelow * 0.6;
    hd->SetMaximum(peakVal * TMath::Power(10.0, decadesAbove));
    hd->SetMinimum(yMin);
    hd->Draw("E1");

    p_sc->SetStats(0);  p_sc->SetFillColor(42);  p_sc->SetFillStyle(1001);  p_sc->SetLineColor(1); p_sc->SetLineWidth(1); p_sc->SetMarkerStyle(0);
    np_sc->SetStats(0); np_sc->SetFillColor(51); np_sc->SetFillStyle(1001); np_sc->SetLineColor(1); np_sc->SetLineWidth(1); np_sc->SetMarkerStyle(0);
    THStack *hs = new THStack(Form("hs_%s", cname), "");
    hs->Add(np_sc); hs->Add(p_sc); hs->Draw("HIST SAME");

    fit_func->SetLineColor(kRed); fit_func->SetLineWidth(3); fit_func->Draw("L SAME");
    hd->Draw("E1 SAME");

    TLatex *tcms = new TLatex(0.1, 0.93, "#bf{CMS} #it{Preliminary}");
    tcms->SetNDC(); tcms->SetTextFont(42); tcms->SetTextSize(0.05); tcms->Draw();
    TLatex *tlumi = new TLatex(0.9, 0.93, "PbPb 5.36 TeV");
    tlumi->SetNDC(); tlumi->SetTextAlign(31); tlumi->SetTextFont(42); tlumi->SetTextSize(0.05); tlumi->Draw();

    TLatex *ti = new TLatex();
    ti->SetNDC(); ti->SetTextAlign(12); ti->SetTextFont(42); ti->SetTextSize(0.045);
    ti->DrawLatex(0.15, 0.85, Form("Cent: %s%%, p_{T}: %s GeV/c", label_cent, label_pt));
    if (cmp_sf > 0)
        ti->DrawLatex(0.15, 0.79, Form("#chi^{2}/NDF = %.2f  (best SF=%.2f)", chi2ndf, cmp_sf));
    else
        ti->DrawLatex(0.15, 0.79, Form("#chi^{2}/NDF = %.2f  (SF = %.2f)", chi2ndf, sf_val));
    ti->DrawLatex(0.15, 0.73, Form("f_{prompt} (all) = %.3f #pm %.3f", fit_func->GetParameter(0), fit_func->GetParError(0)));

    double pBelow = FractionBelow(h_prompt, dca_cut);
    double npBelow = FractionBelow(h_nonprompt, dca_cut);
    double fWhole = fit_func->GetParameter(0);
    double denom = fWhole * pBelow + (1.0 - fWhole) * npBelow;
    double fCut = (denom > 0) ? fWhole * pBelow / denom : 0.0;
    
    double fCutErr = (denom > 0) ? (pBelow * npBelow / (denom * denom)) * fit_func->GetParError(0) : 0.0;
    ti->DrawLatex(0.15, 0.67, Form("f_{prompt}(DCA<%.3g) = %.3f #pm %.3f", dca_cut, fCut, fCutErr));

    TLegend *leg = new TLegend(0.65, 0.65, 0.88, 0.88);
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
    pad2->SetTopMargin(0.0); pad2->SetBottomMargin(0.35); pad2->Draw(); pad2->cd();
    TH1D *pull = new TH1D(Form("pull_%s", cname), "", n_dca_bins, dca_edges);
    for (int k = 1; k <= pull->GetNbinsX(); k++) {
        double fc = fit_func->Eval(pull->GetXaxis()->GetBinCenter(k));
        double dc = h_data->GetBinContent(k), de = h_data->GetBinError(k);
        double pe = p_sc->GetBinError(k);
        double sig = TMath::Sqrt(de * de + pe * pe);
        if (sig == 0) continue;
        pull->SetBinContent(k, (dc - fc) / sig); pull->SetBinError(k, 0.0);
    }
    pull->SetStats(0); pull->SetMarkerStyle(20); pull->SetMarkerSize(0.8); pull->SetLineColor(1);
    pull->GetXaxis()->SetTitle("DCA (cm)"); pull->GetXaxis()->CenterTitle(true);
    pull->GetXaxis()->SetLabelFont(132); pull->GetXaxis()->SetTitleFont(132);
    pull->GetXaxis()->SetLabelSize(0.12); pull->GetXaxis()->SetTitleSize(0.15); pull->GetXaxis()->SetTitleOffset(1.0);
    pull->GetYaxis()->SetTitle("Pull"); pull->GetYaxis()->CenterTitle(true); pull->GetYaxis()->SetNdivisions(505);
    pull->GetYaxis()->SetLabelFont(132); pull->GetYaxis()->SetTitleFont(132);
    pull->GetYaxis()->SetLabelSize(0.1); pull->GetYaxis()->SetTitleSize(0.13); pull->GetYaxis()->SetTitleOffset(0.4);
    pull->GetYaxis()->SetRangeUser(-5.0, 5.0); pull->Draw("P");
    TLine *ref = new TLine(pull->GetXaxis()->GetXmin(), 0.0, pull->GetXaxis()->GetXmax(), 0.0);
    ref->SetLineColor(kRed); ref->SetLineStyle(2); ref->Draw();
    return c;
}





void Chi2_vs_SF()
{
	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	TH1::SetDefaultSumw2(true);

	bool binWidthNorm = false; // Because the input DCA is already bin normalised!

	// TEMP DEBUG (low-pT SF investigation, after widening the MC DCA SF sweep to 0.0-1.5):
	// saves 4 extra per-bin template-fit canvases (best SF +/-1 step, and fixed SF=0.0/0.5).
	// Set to false, or delete this line plus the two blocks in the main loop marked
	// "BEGIN/END TEMP DEBUG BLOCK", to remove once the low-pT fit issue is understood.
	bool saveExtraSFDebugPlots = true;

	const Int_t N_CENTBIN = 5;
	const Int_t N_PTBIN = 9;

	const int cent_min = 0;
	const int cent_max = N_CENTBIN;
	const int ptbin_min = 0;
	const int ptbin_max = N_PTBIN;

	const char *cent[N_CENTBIN] = {"cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"};
	const char *pt[N_PTBIN] = {"pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8", "pT8to10", "pT10to15", "pT15to30", "pT30to100"};
	// pT bin edges, used to place each point at the bin center with an x-error bar spanning the bin range
	// in the prompt-fraction-vs-pT summary graphs.
	const Double_t ptLow[N_PTBIN]  = {2, 3, 4, 5, 6, 8, 10, 15, 30};
	const Double_t ptHigh[N_PTBIN] = {3, 4, 5, 6, 8, 10, 15, 30, 100};

	// Scale factors 0.00 .. 1.50 in steps of 0.01, generated in a loop so there are
	// never missing points. numSF = 151, SF = 1.0 at index 100.
	const Int_t numSF = 151;
	Double_t ScaleFactor[numSF];
	TString sf[numSF];
	for (int i = 0; i < numSF; i++)
	{
		ScaleFactor[i] = i / 100.0;           // i is the value in hundredths (0..150)
		int whole = i / 100;
		int frac = i % 100;
		sf[i] = (frac % 10 == 0) ? Form("x%dp%d", whole, frac / 10)   // 0.50 -> x0p5
		                         : Form("x%dp%02d", whole, frac);      // 0.51 -> x0p51
	}

	// v_n unfolding factor, one row per (cent,pT) bin: f_prompt,1 = f_prompt(DCA <
	// dca_cut), f_prompt,2 = f_prompt(whole DCA range), factor = (1-f1)/(f1-f2).
	ofstream factorTxt(factorTxtName);
	factorTxt << "# f_prompt,1 = prompt fraction for DCA < " << dca_cut << " cm\n";
	factorTxt << "# f_prompt,2 = prompt fraction over the whole DCA range\n";
	factorTxt << "# factor     = (1 - f_prompt,1) / (f_prompt,1 - f_prompt,2)\n";
	factorTxt << "# vn_shift   = factor * (v_n,1 - v_n,2), with (v_n,1 - v_n,2) NOT YET\n";
	factorTxt << "#             MEASURED -- placeholder value used below: " << dvn_placeholder << "\n";
	factorTxt << "#             Replace dvn_placeholder with the real (v_n,1 - v_n,2) and\n";
	factorTxt << "#             recompute once it's available.\n";
	factorTxt << Form("%-12s %-10s %10s %10s %10s %10s %14s %16s\n",
	                   "centrality", "pT", "f_prompt,1", "err", "f_prompt,2", "err", "factor",
	                   Form("vn_shift(dvn=%.3g)", dvn_placeholder));

	// Accumulate prompt fraction (full range and DCA<cut) vs pT for every centrality bin, so all

	Double_t ptCenter_all[N_CENTBIN][N_PTBIN], ptHalfW_all[N_CENTBIN][N_PTBIN];
	Double_t fFull_all[N_CENTBIN][N_PTBIN], fFullErr_all[N_CENTBIN][N_PTBIN];
	Double_t fCut_all[N_CENTBIN][N_PTBIN], fCutErr_all[N_CENTBIN][N_PTBIN];
	Int_t nPtValid_all[N_CENTBIN];

	for (int icent = cent_min; icent < cent_max; icent++)
	{
		nPtValid_all[icent] = 0;

		for (int ipt = ptbin_min; ipt < ptbin_max; ipt++)
		{

			cout << "\n==============================================" << endl;
			cout << "Processing Bin: " << cent[icent] << " | " << pt[ipt] << endl;

			// DCA binning for this pT bin (coarse tail for pT < 5 GeV)
			const int      nb  = dcaNbins(ipt);
			Double_t *const edg = dcaEdges(ipt);

			TH1D *h_data = (TH1D *)data_DCA->Get(Form("hdca_%s_%s", cent[icent], pt[ipt]));
			if (!h_data)
			{
				cout << "Skipping missing data histogram: hdca_" << cent[icent] << "_" << pt[ipt] << endl;
				continue;
			}
			h_data->SetNameTitle(Form("data_dca_%d_%d", icent, ipt), "");
			if (binWidthNorm)
				BinWidthNormalization(h_data);
			h_data->GetXaxis()->SetRangeUser(0.0, DCA_MAX);

			Double_t chi2s_current[numSF];
			TF1 *fits[numSF];
			TH1D *h_prompt_arr[numSF];
			TH1D *h_nonprompt_arr[numSF];

			for (int isf = 0; isf < numSF; isf++)
			{
				h_prompt_arr[isf] = (TH1D *)promptMC_DCA->Get(Form("hdca_%s_%s_%s", cent[icent], pt[ipt], sf[isf].Data()));
				h_nonprompt_arr[isf] = (TH1D *)nonpromptMC_DCA->Get(Form("hdca_%s_%s_%s", cent[icent], pt[ipt], sf[isf].Data()));

				if (!h_prompt_arr[isf] || !h_nonprompt_arr[isf])
				{
					chi2s_current[isf] = 9999.0;
					fits[isf] = nullptr;
					continue;
				}

				// bring 18-bin templates onto this pT bin's DCA binning
				h_prompt_arr[isf]    = MatchDcaBinning(h_prompt_arr[isf], ipt);
				h_nonprompt_arr[isf] = MatchDcaBinning(h_nonprompt_arr[isf], ipt);

				if (binWidthNorm)
				{
					BinWidthNormalization(h_prompt_arr[isf]);
					BinWidthNormalization(h_nonprompt_arr[isf]);
				}

				// Double_t prompt_norm = h_prompt_arr[isf]->Integral("width");
				// Double_t nonprompt_norm = h_nonprompt_arr[isf]->Integral("width");

				Double_t prompt_norm = h_prompt_arr[isf]->Integral();
				Double_t nonprompt_norm = h_nonprompt_arr[isf]->Integral();

				if (prompt_norm > 0.0)
					h_prompt_arr[isf]->Scale(1.0 / prompt_norm);

				if (nonprompt_norm > 0.0)
					h_nonprompt_arr[isf]->Scale(1.0 / nonprompt_norm);

				h_prompt_arr[isf]->GetXaxis()->SetRangeUser(0.0, DCA_MAX);
				h_nonprompt_arr[isf]->GetXaxis()->SetRangeUser(0.0, DCA_MAX);

				h_prompt_MC = h_prompt_arr[isf];
				h_nonprompt_MC = h_nonprompt_arr[isf];

				fits[isf] = fitHist(h_data, Form("fitFunc_%d_%d_%d", icent, ipt, isf));
				chi2s_current[isf] = fits[isf]->GetChisquare() / fits[isf]->GetNDF();
			}

			int best_isf = getBestSFIndex(numSF, chi2s_current);

			if (best_isf < 0 || fits[best_isf] == nullptr || chi2s_current[best_isf] >= 9999.0)
			{
				cout << "No valid scale factor found for " << cent[icent] << " " << pt[ipt] << endl;
				continue;
			}

			Double_t best_SF = ScaleFactor[best_isf];
			TF1 *best_fit = fits[best_isf];

			h_prompt_MC = h_prompt_arr[best_isf];
			h_nonprompt_MC = h_nonprompt_arr[best_isf];

			cout << ">>> Optimal Scale Factor: " << best_SF << " (chi2/NDF = " << chi2s_current[best_isf] << ")" << endl;
			cout << ">>> Prompt Fraction: " << best_fit->GetParameter(0) << " +/- " << best_fit->GetParError(0) << endl;

			// Same f_prompt(DCA<cut) calculation as in DrawTemplateFit, kept here too so it
			double pBelow_val = FractionBelow(h_prompt_arr[best_isf], dca_cut);
			double npBelow_val = FractionBelow(h_nonprompt_arr[best_isf], dca_cut);
			double fWhole_val = best_fit->GetParameter(0);
			double denom_val = fWhole_val * pBelow_val + (1.0 - fWhole_val) * npBelow_val;
			double fCut_val = (denom_val > 0) ? fWhole_val * pBelow_val / denom_val : 0.0;
			double fWholeErr_val = best_fit->GetParError(0);
			// d(fCut)/d(fWhole) = pBelow*npBelow / denom^2 -- propagates the fit's own
			double fCutErr_val = (denom_val > 0) ? (pBelow_val * npBelow_val / (denom_val * denom_val)) * fWholeErr_val : 0.0;

			// v_n unfolding factor: f_prompt,1 = fCut_val (DCA < dca_cut), f_prompt,2 =
			// fWhole_val (whole DCA range). fCut_val and fWhole_val come from the same
			// fit (fCut_val is a deterministic reparametrisation of fWhole_val via the
			// template shapes), so no separate error is propagated for the factor itself.
			double vnFactor = (fCut_val - fWhole_val != 0.0)
			                       ? (1.0 - fCut_val) / (fCut_val - fWhole_val)
			                       : std::numeric_limits<double>::quiet_NaN();
			double vnShift = vnFactor * dvn_placeholder; // PLACEHOLDER (v_n,1-v_n,2) -- see header
			factorTxt << Form("%-12s %-10s %10.4f %10.4f %10.4f %10.4f %14.4f %16.4f\n",
			                   cent[icent], pt[ipt], fCut_val, fCutErr_val, fWhole_val, fWholeErr_val,
			                   vnFactor, vnShift);

			Int_t &nPtValid = nPtValid_all[icent];
			ptCenter_all[icent][nPtValid] = 0.5 * (ptLow[ipt] + ptHigh[ipt]);
			ptHalfW_all[icent][nPtValid] = 0.5 * (ptHigh[ipt] - ptLow[ipt]);
			fFull_all[icent][nPtValid] = fWhole_val;
			fFullErr_all[icent][nPtValid] = fWholeErr_val;
			fCut_all[icent][nPtValid] = fCut_val;
			fCutErr_all[icent][nPtValid] = fCutErr_val;
			nPtValid++;

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
                nb, edg);

			// ==========================================
            // Unscaled (SF=1.0) comparison canvas
            // ==========================================
            const int sf_unity_idx = 100; // ScaleFactor[100] == 1.0
            TCanvas *c_unscaled = nullptr;
            if (fits[sf_unity_idx] && h_prompt_arr[sf_unity_idx] && h_nonprompt_arr[sf_unity_idx])
                c_unscaled = DrawTemplateFit(
                    Form("c_unscaled_%d_%d", icent, ipt),
                    h_data, fits[sf_unity_idx],
                    h_prompt_arr[sf_unity_idx], h_nonprompt_arr[sf_unity_idx],
                    1.0, chi2s_current[sf_unity_idx],
                    s_cent.Data(), s_pt.Data(),
                    nb, edg,
                    best_SF, chi2s_current[best_isf]);

			// ===================== BEGIN TEMP DEBUG BLOCK =====================
			// Extra SF template-fit canvases for the low-pT investigation. To remove:
			// delete this block, the matching WRITE-section block below (same marker),
			// and the `saveExtraSFDebugPlots` declaration near the top of this function.
			TCanvas *c_sf_bestMinus1 = nullptr;
			TCanvas *c_sf_bestPlus1 = nullptr;
			TCanvas *c_sf_zero = nullptr;
			TCanvas *c_sf_half = nullptr;
			if (saveExtraSFDebugPlots)
			{
				const int sf_zero_idx = 0;  // ScaleFactor[0]  == 0.0
				const int sf_half_idx = 50; // ScaleFactor[50] == 0.5

				if (best_isf - 1 >= 0 && fits[best_isf - 1] && h_prompt_arr[best_isf - 1] && h_nonprompt_arr[best_isf - 1])
					c_sf_bestMinus1 = DrawTemplateFit(
						Form("c_sf_bestMinus1_%d_%d", icent, ipt),
						h_data, fits[best_isf - 1],
						h_prompt_arr[best_isf - 1], h_nonprompt_arr[best_isf - 1],
						ScaleFactor[best_isf - 1], chi2s_current[best_isf - 1],
						s_cent.Data(), s_pt.Data(),
						nb, edg,
						best_SF, chi2s_current[best_isf]);
				else
					cout << "DEBUG: SF(best-1) canvas skipped for " << cent[icent] << " " << pt[ipt] << " (best_isf=" << best_isf << ")" << endl;

				if (best_isf + 1 < numSF && fits[best_isf + 1] && h_prompt_arr[best_isf + 1] && h_nonprompt_arr[best_isf + 1])
					c_sf_bestPlus1 = DrawTemplateFit(
						Form("c_sf_bestPlus1_%d_%d", icent, ipt),
						h_data, fits[best_isf + 1],
						h_prompt_arr[best_isf + 1], h_nonprompt_arr[best_isf + 1],
						ScaleFactor[best_isf + 1], chi2s_current[best_isf + 1],
						s_cent.Data(), s_pt.Data(),
						nb, edg,
						best_SF, chi2s_current[best_isf]);
				else
					cout << "DEBUG: SF(best+1) canvas skipped for " << cent[icent] << " " << pt[ipt] << " (best_isf=" << best_isf << ")" << endl;

				if (fits[sf_zero_idx] && h_prompt_arr[sf_zero_idx] && h_nonprompt_arr[sf_zero_idx])
					c_sf_zero = DrawTemplateFit(
						Form("c_sf_zero_%d_%d", icent, ipt),
						h_data, fits[sf_zero_idx],
						h_prompt_arr[sf_zero_idx], h_nonprompt_arr[sf_zero_idx],
						ScaleFactor[sf_zero_idx], chi2s_current[sf_zero_idx],
						s_cent.Data(), s_pt.Data(),
						nb, edg,
						best_SF, chi2s_current[best_isf]);
				else
					cout << "DEBUG: SF=0.0 canvas skipped for " << cent[icent] << " " << pt[ipt]
					     << " (fits=" << (void*)fits[sf_zero_idx] << ", prompt=" << (void*)h_prompt_arr[sf_zero_idx]
					     << ", nonprompt=" << (void*)h_nonprompt_arr[sf_zero_idx] << ")" << endl;

				if (fits[sf_half_idx] && h_prompt_arr[sf_half_idx] && h_nonprompt_arr[sf_half_idx])
					c_sf_half = DrawTemplateFit(
						Form("c_sf_half_%d_%d", icent, ipt),
						h_data, fits[sf_half_idx],
						h_prompt_arr[sf_half_idx], h_nonprompt_arr[sf_half_idx],
						ScaleFactor[sf_half_idx], chi2s_current[sf_half_idx],
						s_cent.Data(), s_pt.Data(),
						nb, edg,
						best_SF, chi2s_current[best_isf]);
			}
			// ====================== END TEMP DEBUG BLOCK =======================

			// ==========================================
			// Draw the Chi2 vs SF
			// ==========================================
			TCanvas *c_sf = new TCanvas(Form("c_chi2_%d_%d", icent, ipt), Form("Chi2 vs SF %s %s", cent[icent], pt[ipt]), 900, 700);
			c_sf->SetLeftMargin(0.12);
			c_sf->SetRightMargin(0.05);
			c_sf->SetTopMargin(0.1);
			c_sf->SetBottomMargin(0.12);
			c_sf->cd();

			TGraph *graph = new TGraph(numSF, ScaleFactor, chi2s_current);
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

			TLatex *tex_cms2 = new TLatex(0.12, 0.93, "#bf{CMS} #it{Preliminary}");
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

			/*if (h_data_shape->Integral("width") > 0)
				h_data_shape->Scale(1.0 / h_data_shape->Integral("width"));
			if (p_MC_shape->Integral("width") > 0)
				p_MC_shape->Scale(1.0 / p_MC_shape->Integral("width"));
			if (np_MC_shape->Integral("width") > 0)
				np_MC_shape->Scale(1.0 / np_MC_shape->Integral("width"));
				*/

			if (h_data_shape->Integral() > 0)
				h_data_shape->Scale(1.0 / h_data_shape->Integral());
			if (p_MC_shape->Integral() > 0)
				p_MC_shape->Scale(1.0 / p_MC_shape->Integral());
			if (np_MC_shape->Integral() > 0)
				np_MC_shape->Scale(1.0 / np_MC_shape->Integral());

			h_data_shape->SetStats(0); // EXPLICITLY TURN OFF STATS
			//h_data_shape->SetTitle(Form("Shape Comparison %s %s", cent[icent], pt[ipt]));
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

			// Y-axis floor: half the smallest nonzero non-prompt-shape bin, same approach
			// as Canvas_BestFit, so the non-prompt (and prompt) shape is never clipped.
			double ovMin = 1e300;
			for (int b = 1; b <= np_MC_shape->GetNbinsX(); ++b) {
				double v = np_MC_shape->GetBinContent(b);
				if (v > 0.0 && v < ovMin) ovMin = v;
			}
			double ovYMin = (ovMin < 1e300) ? ovMin * 0.5 : 1e-5;

			double MaxValShape = TMath::Max(h_data_shape->GetMaximum(), TMath::Max(p_MC_shape->GetMaximum(), np_MC_shape->GetMaximum()));
			h_data_shape->SetMaximum(MaxValShape * 50.0);
			h_data_shape->SetMinimum(ovYMin);
			h_data_shape->GetXaxis()->SetLabelFont(42);
			h_data_shape->GetXaxis()->SetTitleFont(42);
			h_data_shape->GetYaxis()->SetLabelFont(42);
			h_data_shape->GetYaxis()->SetTitleFont(42);

			h_data_shape->Draw("E1");
			p_MC_shape->Draw("HIST SAME");
			np_MC_shape->Draw("HIST SAME");

			TLatex *tex_cms_ov = new TLatex(0.12, 0.93, "#bf{CMS} #it{Preliminary}");
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
			// Gen-level DCA overlay (SF=0.0, i.e. no resolution smearing applied --
			// stands in for the gen/truth-level DCA shape). Prompt/non-prompt only,
			// no data equivalent at gen level.
			// ==========================================
			TCanvas *c_overlay_gen = new TCanvas(Form("c_overlay_gen_%d_%d", icent, ipt), Form("Overlay Shape Gen %s %s", cent[icent], pt[ipt]), 900, 800);
			c_overlay_gen->SetLeftMargin(0.12);
			c_overlay_gen->SetRightMargin(0.05);
			c_overlay_gen->SetTopMargin(0.1);
			c_overlay_gen->SetBottomMargin(0.12);
			c_overlay_gen->SetLogy();
			c_overlay_gen->cd();

			const int sf_zero_idx_gen = 0; // ScaleFactor[0] == 0.0
			if (h_prompt_arr[sf_zero_idx_gen] && h_nonprompt_arr[sf_zero_idx_gen])
			{
				TH1D *p_gen_shape = (TH1D *)h_prompt_arr[sf_zero_idx_gen]->Clone(Form("prompt_gen_shape_%d_%d", icent, ipt));
				TH1D *np_gen_shape = (TH1D *)h_nonprompt_arr[sf_zero_idx_gen]->Clone(Form("nonprompt_gen_shape_%d_%d", icent, ipt));

				if (p_gen_shape->Integral() > 0)
					p_gen_shape->Scale(1.0 / p_gen_shape->Integral());
				if (np_gen_shape->Integral() > 0)
					np_gen_shape->Scale(1.0 / np_gen_shape->Integral());

				p_gen_shape->SetStats(0);
				p_gen_shape->SetTitle("");
				p_gen_shape->GetYaxis()->SetTitle("Normalised entries");
				p_gen_shape->GetXaxis()->SetTitle("DCA (cm)");
				p_gen_shape->GetXaxis()->SetLabelSize(0.04);
				p_gen_shape->GetXaxis()->SetTitleSize(0.05);
				p_gen_shape->GetYaxis()->SetLabelSize(0.04);
				p_gen_shape->GetYaxis()->SetTitleSize(0.05);
				p_gen_shape->GetYaxis()->SetTitleOffset(1.2);
				p_gen_shape->GetXaxis()->SetLabelFont(42);
				p_gen_shape->GetXaxis()->SetTitleFont(42);
				p_gen_shape->GetYaxis()->SetLabelFont(42);
				p_gen_shape->GetYaxis()->SetTitleFont(42);
				p_gen_shape->SetLineColor(kBlue);
				p_gen_shape->SetLineWidth(2);
				p_gen_shape->SetMarkerStyle(0);
				p_gen_shape->SetFillStyle(0);

				np_gen_shape->SetStats(0);
				np_gen_shape->SetLineColor(kRed);
				np_gen_shape->SetLineWidth(2);
				np_gen_shape->SetMarkerStyle(0);
				np_gen_shape->SetFillStyle(0);

				double genMin = 1e300;
				for (int b = 1; b <= np_gen_shape->GetNbinsX(); ++b) {
					double v = np_gen_shape->GetBinContent(b);
					if (v > 0.0 && v < genMin) genMin = v;
				}
				double genYMin = (genMin < 1e300) ? genMin * 0.5 : 1e-5;
				double MaxValGen = TMath::Max(p_gen_shape->GetMaximum(), np_gen_shape->GetMaximum());
				p_gen_shape->SetMaximum(MaxValGen * 50.0);
				p_gen_shape->SetMinimum(genYMin);

				p_gen_shape->Draw("HIST");
				np_gen_shape->Draw("HIST SAME");

				TLatex *tex_cms_gen = new TLatex(0.12, 0.93, "#bf{CMS} #it{Preliminary}");
				tex_cms_gen->SetNDC(); tex_cms_gen->SetTextFont(42); tex_cms_gen->SetTextSize(0.05); tex_cms_gen->Draw();

				TLatex *tex_lumi_gen = new TLatex(0.95, 0.93, "PbPb 5.36 TeV");
				tex_lumi_gen->SetNDC(); tex_lumi_gen->SetTextAlign(31); tex_lumi_gen->SetTextFont(42); tex_lumi_gen->SetTextSize(0.05); tex_lumi_gen->Draw();

				TLatex *tex_info_gen = new TLatex();
				tex_info_gen->SetNDC(); tex_info_gen->SetTextFont(42); tex_info_gen->SetTextSize(0.045); tex_info_gen->SetTextAlign(12);
				tex_info_gen->DrawLatex(0.15, 0.85, Form("Cent: %s%%", s_cent.Data()));
				tex_info_gen->DrawLatex(0.15, 0.79, Form("p_{T}: %s GeV/c", s_pt.Data()));
				tex_info_gen->DrawLatex(0.15, 0.73, "Gen-level DCA (SF = 0.0)");

				TLegend *leg_gen = new TLegend(0.55, 0.70, 0.88, 0.88);
				leg_gen->SetBorderSize(0);
				leg_gen->SetFillStyle(0);
				leg_gen->SetTextFont(42);
				leg_gen->SetTextSize(0.04);
				leg_gen->AddEntry(p_gen_shape, "Prompt Gen Shape", "l");
				leg_gen->AddEntry(np_gen_shape, "Non-Prompt Gen Shape", "l");
				leg_gen->Draw();
			}

			// ==========================================
			// WRITE TO FILE
			// ==========================================
			f_out->cd();
			TDirectory *dir = f_out->mkdir(Form("%s_%s", cent[icent], pt[ipt]));
			dir->cd();
            c_overlay->Write("Canvas_Overlay_Shape");
			c_overlay_gen->Write("Canvas_Overlay_Shape_Gen");
			c_sf->Write("Canvas_Chi2_vs_SF");
			c_unscaled->Write("Canvas_Unscaled_Fit");
			c_best->Write("Canvas_BestFit");

			// ===================== BEGIN TEMP DEBUG BLOCK =====================
			if (saveExtraSFDebugPlots)
			{
				if (c_sf_bestMinus1) c_sf_bestMinus1->Write("Canvas_SF_BestMinus1");
				if (c_sf_bestPlus1) c_sf_bestPlus1->Write("Canvas_SF_BestPlus1");
				if (c_sf_zero) c_sf_zero->Write("Canvas_SF_0p0");
				if (c_sf_half) c_sf_half->Write("Canvas_SF_0p5");
			}
			// ====================== END TEMP DEBUG BLOCK =======================
			
			
			
		}

	}

	// ==========================================
	// Prompt fraction (Full DCA vs DCA<cut) vs pT 
	// ==========================================
	for (int icent = cent_min; icent < cent_max; icent++)
	{
		Int_t nPtValid = nPtValid_all[icent];

		TCanvas *c_fpt = new TCanvas(Form("c_fpt_%s", cent[icent]), Form("Prompt Fraction vs pT %s", cent[icent]), 900, 700);
		c_fpt->SetLeftMargin(0.12);
		c_fpt->SetRightMargin(0.05);
		c_fpt->SetTopMargin(0.1);
		c_fpt->SetBottomMargin(0.12);
		c_fpt->cd();

		TGraphErrors *g_full = new TGraphErrors(nPtValid, ptCenter_all[icent], fFull_all[icent], ptHalfW_all[icent], fFullErr_all[icent]);
		g_full->SetName(Form("g_fullDCA_%s", cent[icent]));
		g_full->SetTitle("; p_{T} (GeV/c); Prompt fraction");
		g_full->SetMarkerStyle(20);
		g_full->SetMarkerSize(1.2);
		g_full->SetMarkerColor(kBlue + 1);
		g_full->SetLineColor(kBlue + 1);
		g_full->GetXaxis()->SetLimits(0.0, ptHigh[N_PTBIN - 1]);
		g_full->GetXaxis()->SetLabelFont(42);
		g_full->GetXaxis()->SetTitleFont(42);
		g_full->GetXaxis()->SetTitleSize(0.05);
		g_full->GetYaxis()->SetLabelFont(42);
		g_full->GetYaxis()->SetTitleFont(42);
		g_full->GetYaxis()->SetTitleSize(0.05);
		g_full->GetYaxis()->SetTitleOffset(1.1);
		g_full->GetYaxis()->SetRangeUser(0.0, 1.2);
		g_full->Draw("AP");

		TGraphErrors *g_cut = new TGraphErrors(nPtValid, ptCenter_all[icent], fCut_all[icent], ptHalfW_all[icent], fCutErr_all[icent]);
		g_cut->SetName(Form("g_dcaCutFrac_%s", cent[icent]));
		g_cut->SetMarkerStyle(21);
		g_cut->SetMarkerSize(1.2);
		g_cut->SetMarkerColor(kRed + 1);
		g_cut->SetLineColor(kRed + 1);
		g_cut->Draw("P SAME");

		TLine *refLine_fpt = new TLine(ptLow[0], 1.0, ptHigh[N_PTBIN - 1], 1.0);
		refLine_fpt->SetLineColor(kBlack);
		refLine_fpt->SetLineStyle(3); // dotted
		refLine_fpt->Draw();

		TLatex *tex_cms3 = new TLatex(0.12, 0.93, "#bf{CMS} #it{Preliminary}");
		tex_cms3->SetNDC();
		tex_cms3->SetTextFont(42);
		tex_cms3->SetTextSize(0.05);
		tex_cms3->Draw();

		TLatex *tex_lumi3 = new TLatex(0.95, 0.93, "PbPb 5.36 TeV");
		tex_lumi3->SetNDC();
		tex_lumi3->SetTextAlign(31);
		tex_lumi3->SetTextFont(42);
		tex_lumi3->SetTextSize(0.05);
		tex_lumi3->Draw();

		TString s_cent_fpt = cent[icent];
		s_cent_fpt.ReplaceAll("cent", "");
		s_cent_fpt.ReplaceAll("to", "-");
		TLatex *tex_info3 = new TLatex();
		tex_info3->SetNDC();
		tex_info3->SetTextFont(42);
		tex_info3->SetTextSize(0.045);
		tex_info3->SetTextAlign(12);
		tex_info3->DrawLatex(0.15, 0.85, Form("Cent: %s%%", s_cent_fpt.Data()));

		TLegend *leg_fpt = new TLegend(0.55, 0.2, 0.88, 0.35);
		leg_fpt->SetBorderSize(0);
		leg_fpt->SetFillStyle(0);
		leg_fpt->SetTextFont(42);
		leg_fpt->SetTextSize(0.04);
		leg_fpt->AddEntry(g_full, "Full DCA range", "pe");
		leg_fpt->AddEntry(g_cut, Form("DCA < %.3g cm", dca_cut), "pe");
		leg_fpt->Draw();

		f_out->cd();
		c_fpt->Write(Form("Canvas_PromptFraction_vs_pT_%s", cent[icent]));
	}

	f_out->cd();
	f_out->Close();
	factorTxt.close();
	cout << "\n>>> Done! All output saved to "<<f_out->GetName() << endl;
	cout << ">>> v_n factor table written to " << factorTxtName << endl;
}
