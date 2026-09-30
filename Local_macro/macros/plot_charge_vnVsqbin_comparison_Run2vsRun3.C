// plot_vn_comparison_Run2vsRun3.C
//
// Comparison of charged-particle vn vs qn-bin:
//   Run 3  (PbPb 5.36 TeV, 0.5 < pT < 3 GeV/c) — flow_Analysis_chg_out_combined_May11.root
//   Run 3  (PbPb 5.36 TeV, 1.0 < pT < 3 GeV/c) — same file
//   Run 2  (PbPb 5.02 TeV, 1.0 < pT < 3 GeV/c) — ch_v2_vs_q2.root  (v2 only)
//
// Centrality bins: 0-10, 10-20, 20-30, 30-40, 40-50 %
//
// Pass includeRun2 = false to skip Run2 entirely (Run3-only plots) — the Run2
// input files then aren't even opened, and output filenames say "Run3only"
// instead of "Run2vsRun3" so both variants can coexist.
//
// Output: plot_v2_vsq2_Run2vsRun3_May11.pdf  (or ..._Run3only_May11.pdf)
//         plot_v3_vsq3_Run2vsRun3_May11.pdf  (or ..._Run3only_May11.pdf)
//
// Usage:
//   root -l -q plot_charge_vnVsqbin_comparison_Run2vsRun3.C+
//   root -l -q 'plot_charge_vnVsqbin_comparison_Run2vsRun3.C+(false)'   // Run3 only

#include "TFile.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TString.h"

// ── global style ─────────────────────────────────────────────────────────────
void setStyle_comp()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
}

// ── Convert TProfile → TGraphErrors using actual bin centres ────────────────
TGraphErrors* profileToGraph_comp(TProfile* p, const char* name)
{
    int n = 0;
    for (int i = 1; i <= p->GetNbinsX(); i++)
        if (p->GetBinEntries(i) > 0) n++;

    TGraphErrors* g = new TGraphErrors(n);
    g->SetName(name);
    int pt = 0;
    for (int i = 1; i <= p->GetNbinsX(); i++) {
        if (p->GetBinEntries(i) <= 0) continue;
        g->SetPoint    (pt, p->GetBinCenter(i) - 0.5, p->GetBinContent(i));
        g->SetPointError(pt, 0.0,                         p->GetBinError(i));
        pt++;
    }
    return g;
}

// ── Convert TH1D → TGraphErrors, bin centres shifted to 0,1,...,9 index ──────
// Run2 v2 histogram binning is [-0.5,9.5] -> centres already at 0,1,...,9 (no shift).
// Run2 v3 histogram binning is [0,10]     -> centres at 0.5,1.5,...,9.5 (needs -0.5 shift).
TGraphErrors* histToGraph_comp(TH1D* h, const char* name, bool shiftHalf = false)
{
    int n = 0;
    for (int i = 1; i <= h->GetNbinsX(); i++)
        if (h->GetBinContent(i) != 0.0 || h->GetBinError(i) != 0.0) n++;

    TGraphErrors* g = new TGraphErrors(n);
    g->SetName(name);
    int pt = 0;
    for (int i = 1; i <= h->GetNbinsX(); i++) {
        if (h->GetBinContent(i) == 0.0 && h->GetBinError(i) == 0.0) continue;
        double x = shiftHalf ? h->GetBinCenter(i) - 0.5 : h->GetBinCenter(i);
        g->SetPoint    (pt, x, h->GetBinContent(i));
        g->SetPointError(pt, 0.0,               h->GetBinError(i));
        pt++;
    }
    return g;
}

// ── Draw one comparison pad ───────────────────────────────────────────────────
void drawCompPad(TVirtualPad* pad,
                 TGraphErrors* gRun3_05to3,
                 TGraphErrors* gRun3_1to3,
                 TGraphErrors* gRun2,          // may be nullptr (v3 has no Run2)
                 int           ivn,
                 const char*   centLabel,
                 double        yMin,
                 double        yMax,
                 bool          legendTopRight = false)
{
    pad->SetLeftMargin  (0.16);
    pad->SetRightMargin (0.04);
    pad->SetBottomMargin(0.14);
    pad->SetTopMargin   (0.10);
    pad->cd();

    // Axis frame: Run3 spans 0,1,...,11 (12 q-bins); Run2 (10 q-bins) sits at 0,1,...,9
    TH2F* hf = new TH2F(Form("hf_%s", gRun3_05to3->GetName()), "",
                         12, -0.5, 11.5, 100, yMin, yMax);
    hf->GetXaxis()->SetTitle(Form("q_{%d} bin", ivn));
    hf->GetYaxis()->SetTitle(Form("v_{%d}(h^{#pm})", ivn));
    hf->GetXaxis()->SetTitleSize (0.060);
    hf->GetYaxis()->SetTitleSize (0.060);
    hf->GetXaxis()->SetLabelSize (0.055);
    hf->GetYaxis()->SetLabelSize (0.055);
    hf->GetXaxis()->SetTitleOffset(1.00);
    hf->GetYaxis()->SetTitleOffset(1.20);
    hf->GetXaxis()->SetNdivisions(510);
    hf->Draw("AXIS");

    // Run 3 (0.5-3 GeV/c): blue filled circle
    gRun3_05to3->SetMarkerStyle(20);
    gRun3_05to3->SetMarkerSize (1.4);
    gRun3_05to3->SetMarkerColor(kBlue+1);
    gRun3_05to3->SetLineColor  (kBlue+1);
    gRun3_05to3->SetLineWidth  (1);
    gRun3_05to3->Draw("P SAME");

    // Run 3 (1-3 GeV/c): green filled triangle
    gRun3_1to3->SetMarkerStyle(22);
    gRun3_1to3->SetMarkerSize (1.4);
    gRun3_1to3->SetMarkerColor(kGreen+2);
    gRun3_1to3->SetLineColor  (kGreen+2);
    gRun3_1to3->SetLineWidth  (1);
    gRun3_1to3->Draw("P SAME");

    // Run 2 (v2 only): red open square
    if (gRun2) {
        gRun2->SetMarkerStyle(25);
        gRun2->SetMarkerSize (1.4);
        gRun2->SetMarkerColor(kRed+1);
        gRun2->SetLineColor  (kRed+1);
        gRun2->SetLineWidth  (1);
        gRun2->Draw("P SAME");
    }

    // Legend — bottom right (default) or top right for 0-10%
    double leg_x1, leg_y1, leg_x2, leg_y2;
    if (legendTopRight) {
        leg_x1 = 0.30; leg_y1 = 0.57; leg_x2 = 0.96; leg_y2 = 0.80;
    } else {
        leg_x1 = 0.32; leg_y1 = 0.20; leg_x2 = 0.96; leg_y2 = 0.43;
    }
    TLegend* leg = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle (0);
    leg->SetTextFont  (42);
    leg->SetTextSize  (0.050);
    leg->AddEntry(gRun3_05to3, "Run 3  0.5 < p_{T} < 3 GeV/c", "p");
    leg->AddEntry(gRun3_1to3,  "Run 3  1.0 < p_{T} < 3 GeV/c", "p");
    if (gRun2)
        leg->AddEntry(gRun2,   "Run 2  1.0 < p_{T} < 3 GeV/c", "p");
    leg->Draw();

    // Labels
    TLatex lat;
    lat.SetNDC();

    // CMS bold, Preliminary normal-weight italic — above the top axis line.
    // Base font must be regular (42, not bold 62/60): #bf{} then bolds only
    // "CMS", and #it{} italicizes "Preliminary" at normal weight; if the base
    // font were already bold, #it{Preliminary} would render bold-italic instead.
    lat.SetTextFont(42);
    lat.SetTextSize(0.065);
    lat.SetTextAlign(11);
    lat.DrawLatex(0.16, 0.915, "#bf{CMS} #it{Preliminary}");

    // PbPb 5.36 TeV — top right, aligned with CMS Preliminary
    lat.SetTextFont(42);
    lat.SetTextSize(0.052);
    lat.SetTextAlign(31);
    lat.DrawLatex(0.96, 0.915, "PbPb 5.36 TeV");

    // Centrality and eta — inside frame, top left
    lat.SetTextFont(42);
    lat.SetTextSize(0.052);
    lat.SetTextAlign(11);
    lat.DrawLatex(0.19, 0.83, Form("Cent: %s", centLabel));
    lat.DrawLatex(0.19, 0.76, "|#eta| < 1.0");
}

// ── main ─────────────────────────────────────────────────────────────────────
void plot_charge_vnVsqbin_comparison_Run2vsRun3()
{
    //TFile* fRun3 = TFile::Open("ROOT/flow_Analysis_chg_out_combined_May11.root");
    TFile* fRun3 = TFile::Open("Charge_flow_Analysis_12NUQbin_diffq2q3_out_combined_Aug24.root");
    bool includeRun2 = true;
    
    setStyle_comp();
    gROOT->ForceStyle();


    TString dateTag = "Sept17";
    TString run2Tag = includeRun2 ? "Run2vsRun3" : "Run3only";
    const double Y_MIN_V2 = 0.0,  Y_MAX_V2 = 0.27;
    const double Y_MIN_V3 = 0.0,  Y_MAX_V3 = 0.12;

    // ─── Open files ──────────────────────────────────────────────────────────
    
    if (!fRun3 || fRun3->IsZombie()) {
        Printf("ERROR: cannot open ROOT file"); return;
    }
    TDirectoryFile* dRun3 = (TDirectoryFile*)fRun3->Get("vsQbin_TProfile");
    if (!dRun3) {
        Printf("ERROR: vsQbin_TProfile directory not found"); return;
    }

    // Run2 input files are only needed (and only opened) when includeRun2 is true.
    TFile* fRun2  = nullptr;
    TFile* fRun2_ = nullptr;
    if (includeRun2) {
        fRun2 = TFile::Open("ch_v2_vs_q2.root");
        if (!fRun2 || fRun2->IsZombie()) {
            Printf("ERROR: cannot open ch_v2_vs_q2.root"); return;
        }
        fRun2_ = TFile::Open("ch_v3_vs_q3_hf_total_new_hist.root");
        if (!fRun2_ || fRun2_->IsZombie()) {
            Printf("ERROR: cannot open ch_v3_vs_q3_hf_total_new_hist.root"); return;
        }
    }

    // ─── Centrality definitions ───────────────────────────────────────────────
    const int NCENT = 5;
    const char* centLabel[NCENT] = {
        "0-10%", "10-20%", "20-30%", "30-40%", "40-50%"
    };
    const char* centTagR3[NCENT] = {
        "cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"
    };
    const char* centTagR2[NCENT] = {
        "_new_cen_0_10", "_new_cen_10_20",
        "_new_cen_20_30", "_new_cen_30_40", "_new_cen_40_50"
    };

    // =========================================================================
    // v2 comparison: Run3 (0.5-3) + Run3 (1-3) + Run2 (1-3)
    // =========================================================================
    {
        TCanvas* cv = new TCanvas("cv_v2_comp", "v2 Run2 vs Run3", 1800, 900);
        cv->Divide(3, 2, 0.002, 0.002);

        for (int ic = 0; ic < NCENT; ic++) {
            TProfile* pR3a = (TProfile*)dRun3->Get(
                Form("hp_v2_vsq2_pt0p5to3_%s", centTagR3[ic]));
            if (!pR3a) { Printf("WARNING: v2 pt0p5to3 %s not found", centTagR3[ic]); continue; }
            TGraphErrors* gR3a = profileToGraph_comp(pR3a, Form("gR3a_v2_%d", ic));

            TProfile* pR3b = (TProfile*)dRun3->Get(
                Form("hp_v2_vsq2_pt1to3_%s", centTagR3[ic]));
            if (!pR3b) { Printf("WARNING: v2 pt1to3 %s not found", centTagR3[ic]); continue; }
            TGraphErrors* gR3b = profileToGraph_comp(pR3b, Form("gR3b_v2_%d", ic));

            TGraphErrors* gR2 = nullptr;
            if (includeRun2) {
                TH1D* hR2 = (TH1D*)fRun2->Get(centTagR2[ic]);
                if (hR2) gR2 = histToGraph_comp(hR2, Form("gR2_v2_%d", ic));
                else     Printf("WARNING: Run2 %s not found", centTagR2[ic]);
            }

            drawCompPad(cv->GetPad(ic + 1), gR3a, gR3b, gR2, 2,
                        centLabel[ic], Y_MIN_V2, Y_MAX_V2,
                        /*legendTopRight=*/ (ic == 0 || ic == 1));
        }

        TString outName = Form("plot_v2_vsq2_%s_%s.pdf", run2Tag.Data(), dateTag.Data());
        cv->SaveAs(outName);
        Printf("Saved: %s", outName.Data());
        delete cv;
    }

    // =========================================================================
    // v3 comparison: Run3 (0.5-3) + Run3 (1-3) + Run2 (1-3)
    // =========================================================================
    {
        TCanvas* cv = new TCanvas("cv_v3_comp", "v3 Run3", 1800, 900);
        cv->Divide(3, 2, 0.002, 0.002);

        for (int ic = 0; ic < NCENT; ic++) {
            TProfile* pR3a = (TProfile*)dRun3->Get(
                Form("hp_v3_vsq3_pt0p5to3_%s", centTagR3[ic]));
            if (!pR3a) { Printf("WARNING: v3 pt0p5to3 %s not found", centTagR3[ic]); continue; }
            TGraphErrors* gR3a = profileToGraph_comp(pR3a, Form("gR3a_v3_%d", ic));

            TProfile* pR3b = (TProfile*)dRun3->Get(
                Form("hp_v3_vsq3_pt1to3_%s", centTagR3[ic]));
            if (!pR3b) { Printf("WARNING: v3 pt1to3 %s not found", centTagR3[ic]); continue; }
            TGraphErrors* gR3b = profileToGraph_comp(pR3b, Form("gR3b_v3_%d", ic));

            TGraphErrors* gR2 = nullptr;
            if (includeRun2) {
                TH1D* hR2_ = (TH1D*)fRun2_->Get(centTagR2[ic]);
                if (hR2_) gR2 = histToGraph_comp(hR2_, Form("gR2_v3_%d", ic), /*shiftHalf=*/true);
                else      Printf("WARNING: Run2 %s not found", centTagR2[ic]);
            }

            drawCompPad(cv->GetPad(ic + 1), gR3a, gR3b, gR2, 3,
                        centLabel[ic], Y_MIN_V2, Y_MAX_V2,
                        /*legendTopRight=*/ (ic == 0 || ic == 1 || ic == 2 || ic == 3 || ic == 4 )); // Show legend on all pads since no Run2 comparison
        }

        TString outName = Form("plot_v3_vsq3_%s_%s.pdf", run2Tag.Data(), dateTag.Data());
        cv->SaveAs(outName);
        Printf("Saved: %s", outName.Data());
        delete cv;
    }

    // =========================================================================
    // v2 correlation: v2(0.5-3 GeV/c) vs v2(1-3 GeV/c), one point per q2 bin
    // =========================================================================
    {
        TCanvas* cv = new TCanvas("cv_v2_corr", "v2 pT correlation", 1800, 900);
        cv->Divide(3, 2, 0.002, 0.002);

        // colour palette cycling through q2 bins (10 bins, indices 0-9)
        const int Q2COLS[10] = {
            kBlack, kRed+1, kOrange+1, kYellow+2, kGreen+2,
            kCyan+2, kBlue+1, kViolet+1, kMagenta+1, kGray+2
        };

        for (int ic = 0; ic < NCENT; ic++) {
            TProfile* pA = (TProfile*)dRun3->Get(
                Form("hp_v2_vsq2_pt0p5to3_%s", centTagR3[ic]));
            if (!pA) { Printf("WARNING: corr v2 pt0p5to3 %s not found", centTagR3[ic]); continue; }
            TProfile* pB = (TProfile*)dRun3->Get(
                Form("hp_v2_vsq2_pt1to3_%s", centTagR3[ic]));
            if (!pB) { Printf("WARNING: corr v2 pt1to3 %s not found", centTagR3[ic]); continue; }

            // build one TGraphErrors per q2 bin so we can colour each separately
            // both profiles have the same binning (bin i = q2-bin i-1)
            int nb = pA->GetNbinsX();

            TVirtualPad* pad = cv->GetPad(ic + 1);
            pad->SetLeftMargin  (0.16);
            pad->SetRightMargin (0.04);
            pad->SetBottomMargin(0.14);
            pad->SetTopMargin   (0.10);
            pad->cd();

            // axis frame
            TH2F* hf = new TH2F(Form("hfcorr_%d", ic), "",
                                  100, 0.0, Y_MAX_V2, 100, 0.0, Y_MAX_V2);
            hf->GetXaxis()->SetTitle("v_{2}(h^{#pm})  1.0 < p_{T} < 3 GeV/c");
            hf->GetYaxis()->SetTitle("v_{2}(h^{#pm})  0.5 < p_{T} < 3 GeV/c");
            hf->GetXaxis()->SetTitleSize (0.055);
            hf->GetYaxis()->SetTitleSize (0.055);
            hf->GetXaxis()->SetLabelSize (0.050);
            hf->GetYaxis()->SetLabelSize (0.050);
            hf->GetXaxis()->SetTitleOffset(1.10);
            hf->GetYaxis()->SetTitleOffset(1.40);
            hf->GetXaxis()->SetNdivisions(505);
            hf->Draw("AXIS");

            // y = x reference line
            TLine* diag = new TLine(0.0, 0.0, Y_MAX_V2, Y_MAX_V2);
            diag->SetLineColor(kGray+1);
            diag->SetLineWidth(1);
            diag->SetLineStyle(2);
            diag->Draw("SAME");

            TLegend* leg = new TLegend(0.64, 0.13, 0.97, 0.43);
            leg->SetBorderSize(0);
            leg->SetFillStyle (0);
            leg->SetTextFont  (42);
            leg->SetTextSize  (0.038);
            leg->SetNColumns  (2);

            for (int iq = 0; iq < nb; iq++) {
                if (pA->GetBinEntries(iq+1) <= 0 || pB->GetBinEntries(iq+1) <= 0) continue;
                double vA  = pA->GetBinContent(iq+1);
                double eA  = pA->GetBinError  (iq+1);
                double vB  = pB->GetBinContent(iq+1);
                double eB  = pB->GetBinError  (iq+1);
                TGraphErrors* gpt = new TGraphErrors(1);
                gpt->SetName(Form("gcorr_%d_%d", ic, iq));
                gpt->SetPoint     (0, vB, vA);
                gpt->SetPointError(0, eB, eA);
                int col = Q2COLS[iq % 10];
                gpt->SetMarkerStyle(20);
                gpt->SetMarkerSize (1.3);
                gpt->SetMarkerColor(col);
                gpt->SetLineColor  (col);
                gpt->SetLineWidth  (1);
                gpt->Draw("P SAME");
                leg->AddEntry(gpt, Form("q_{2} bin %d", iq), "p");
            }
            leg->Draw();

            TLatex lat;
            lat.SetNDC();
            lat.SetTextFont(62);
            lat.SetTextSize(0.065);
            lat.SetTextAlign(11);
            lat.DrawLatex(0.16, 0.915, "#bf{CMS} #it{Preliminary}");
            lat.SetTextFont(42);
            lat.SetTextSize(0.052);
            lat.DrawLatex(0.19, 0.83, Form("Cent: %s", centLabel[ic]));
            lat.DrawLatex(0.19, 0.76, "|#eta| < 1.0");
        }

        TString outName = Form("plot_v2_pTcorr_Run3_%s.pdf", dateTag.Data());
        cv->SaveAs(outName);
        Printf("Saved: %s", outName.Data());
        delete cv;
    }

    // =========================================================================
    // v3 correlation: v3(0.5-3 GeV/c) vs v3(1-3 GeV/c), one point per q3 bin
    // =========================================================================
    {
        TCanvas* cv = new TCanvas("cv_v3_corr", "v3 pT correlation", 1800, 900);
        cv->Divide(3, 2, 0.002, 0.002);

        const int Q3COLS[10] = {
            kBlack, kRed+1, kOrange+1, kYellow+2, kGreen+2,
            kCyan+2, kBlue+1, kViolet+1, kMagenta+1, kGray+2
        };

        for (int ic = 0; ic < NCENT; ic++) {
            TProfile* pA = (TProfile*)dRun3->Get(
                Form("hp_v3_vsq3_pt0p5to3_%s", centTagR3[ic]));
            if (!pA) { Printf("WARNING: corr v3 pt0p5to3 %s not found", centTagR3[ic]); continue; }
            TProfile* pB = (TProfile*)dRun3->Get(
                Form("hp_v3_vsq3_pt1to3_%s", centTagR3[ic]));
            if (!pB) { Printf("WARNING: corr v3 pt1to3 %s not found", centTagR3[ic]); continue; }

            int nb = pA->GetNbinsX();

            TVirtualPad* pad = cv->GetPad(ic + 1);
            pad->SetLeftMargin  (0.16);
            pad->SetRightMargin (0.04);
            pad->SetBottomMargin(0.14);
            pad->SetTopMargin   (0.10);
            pad->cd();

            TH2F* hf = new TH2F(Form("hfcorrv3_%d", ic), "",
                                  100, 0.0, Y_MAX_V3, 100, 0.0, Y_MAX_V3);
            hf->GetXaxis()->SetTitle("v_{3}(h^{#pm})  1.0 < p_{T} < 3 GeV/c");
            hf->GetYaxis()->SetTitle("v_{3}(h^{#pm})  0.5 < p_{T} < 3 GeV/c");
            hf->GetXaxis()->SetTitleSize (0.055);
            hf->GetYaxis()->SetTitleSize (0.055);
            hf->GetXaxis()->SetLabelSize (0.050);
            hf->GetYaxis()->SetLabelSize (0.050);
            hf->GetXaxis()->SetTitleOffset(1.10);
            hf->GetYaxis()->SetTitleOffset(1.40);
            hf->GetXaxis()->SetNdivisions(505);
            hf->Draw("AXIS");

            TLine* diag = new TLine(0.0, 0.0, Y_MAX_V3, Y_MAX_V3);
            diag->SetLineColor(kGray+1);
            diag->SetLineWidth(1);
            diag->SetLineStyle(2);
            diag->Draw("SAME");

            TLegend* leg = new TLegend(0.64, 0.13, 0.97, 0.43);
            leg->SetBorderSize(0);
            leg->SetFillStyle (0);
            leg->SetTextFont  (42);
            leg->SetTextSize  (0.038);
            leg->SetNColumns  (2);

            for (int iq = 0; iq < nb; iq++) {
                if (pA->GetBinEntries(iq+1) <= 0 || pB->GetBinEntries(iq+1) <= 0) continue;
                double vA  = pA->GetBinContent(iq+1);
                double eA  = pA->GetBinError  (iq+1);
                double vB  = pB->GetBinContent(iq+1);
                double eB  = pB->GetBinError  (iq+1);
                TGraphErrors* gpt = new TGraphErrors(1);
                gpt->SetName(Form("gcorrv3_%d_%d", ic, iq));
                gpt->SetPoint     (0, vB, vA);
                gpt->SetPointError(0, eB, eA);
                int col = Q3COLS[iq % 10];
                gpt->SetMarkerStyle(20);
                gpt->SetMarkerSize (1.3);
                gpt->SetMarkerColor(col);
                gpt->SetLineColor  (col);
                gpt->SetLineWidth  (1);
                gpt->Draw("P SAME");
                leg->AddEntry(gpt, Form("q_{3} bin %d", iq), "p");
            }
            leg->Draw();

            TLatex lat;
            lat.SetNDC();
            lat.SetTextFont(62);
            lat.SetTextSize(0.065);
            lat.SetTextAlign(11);
            lat.DrawLatex(0.16, 0.915, "#bf{CMS} #it{Preliminary}");
            lat.SetTextFont(42);
            lat.SetTextSize(0.052);
            lat.DrawLatex(0.19, 0.83, Form("Cent: %s", centLabel[ic]));
            lat.DrawLatex(0.19, 0.76, "|#eta| < 1.0");
        }

        TString outName = Form("plot_v3_pTcorr_Run3_%s.pdf", dateTag.Data());
        cv->SaveAs(outName);
        Printf("Saved: %s", outName.Data());
        delete cv;
    }

    fRun3->Close();
    if (fRun2)  fRun2->Close();
    if (fRun2_) fRun2_->Close();
    Printf("Done.");
}
