// plot_vn_comparison_Run2vsRun3.C
//
// Comparison of charged-particle vn vs qn-bin:
//   Run 3  (PbPb 5.36 TeV, 0.5 < pT < 3 GeV/c) — flow_Analysis_chg_out_combined_May11.root
//   Run 3  (PbPb 5.36 TeV, 1.0 < pT < 3 GeV/c) — same file
//   Run 2  (PbPb 5.02 TeV, 1.0 < pT < 3 GeV/c) — ch_v2_vs_q2.root  (v2 only)
//
// Centrality bins: 0-10, 10-20, 20-30, 30-40, 40-50 %
//
// Output: plot_v2_vsq2_Run2vsRun3_May11.pdf
//         plot_v3_vsq3_Run3_May11.pdf
//
// Usage:
//   root -l -q plot_vn_comparison_Run2vsRun3.C

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

// ── Convert TH1D → TGraphErrors, bin centres shifted by +0.5 ─────────────────────
TGraphErrors* histToGraph_comp(TH1D* h, const char* name)
{
    int n = 0;
    for (int i = 1; i <= h->GetNbinsX(); i++)
        if (h->GetBinContent(i) != 0.0 || h->GetBinError(i) != 0.0) n++;

    TGraphErrors* g = new TGraphErrors(n);
    g->SetName(name);
    int pt = 0;
    for (int i = 1; i <= h->GetNbinsX(); i++) {
        if (h->GetBinContent(i) == 0.0 && h->GetBinError(i) == 0.0) continue;
        // Run2 centres already at 0,1,...,9 — no shift needed
        g->SetPoint    (pt, h->GetBinCenter(i), h->GetBinContent(i));
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

    // Axis frame: all datasets at 0,1,...,9
    TH2F* hf = new TH2F(Form("hf_%s", gRun3_05to3->GetName()), "",
                         10, -0.5, 9.5, 100, yMin, yMax);
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

    // CMS bold — above the top axis line (top margin area)
    lat.SetTextFont(62);
    lat.SetTextSize(0.065);
    lat.SetTextAlign(11);
    lat.DrawLatex(0.16, 0.915, "CMS");

    // Centrality and eta — inside frame, top left
    lat.SetTextFont(42);
    lat.SetTextSize(0.052);
    lat.SetTextAlign(11);
    lat.DrawLatex(0.19, 0.83, Form("Cent: %s", centLabel));
    lat.DrawLatex(0.19, 0.76, "|#eta| < 1.0");
}

// ── main ─────────────────────────────────────────────────────────────────────
void plot_vn_comparison_Run2vsRun3()
{
    setStyle_comp();
    gROOT->ForceStyle();

    // =========================================================================
    // Y-axis ranges
    // =========================================================================
    const double Y_MIN_V2 = 0.0,  Y_MAX_V2 = 0.27;
    const double Y_MIN_V3 = 0.0,  Y_MAX_V3 = 0.12;

    // ─── Open files ──────────────────────────────────────────────────────────
    TFile* fRun3 = TFile::Open("ROOT/flow_Analysis_chg_out_combined_May11.root");
    if (!fRun3 || fRun3->IsZombie()) {
        Printf("ERROR: cannot open ROOT/flow_Analysis_chg_out_combined_May11.root"); return;
    }
    TDirectoryFile* dRun3 = (TDirectoryFile*)fRun3->Get("vsQbin_TProfile");
    if (!dRun3) {
        Printf("ERROR: vsQbin_TProfile directory not found"); return;
    }

    TFile* fRun2 = TFile::Open("ROOT/ch_v2_vs_q2.root");
    if (!fRun2 || fRun2->IsZombie()) {
        Printf("ERROR: cannot open ROOT/ch_v2_vs_q2.root"); return;
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

            TH1D* hR2 = (TH1D*)fRun2->Get(centTagR2[ic]);
            if (!hR2) { Printf("WARNING: Run2 %s not found", centTagR2[ic]); continue; }
            TGraphErrors* gR2 = histToGraph_comp(hR2, Form("gR2_v2_%d", ic));

            drawCompPad(cv->GetPad(ic + 1), gR3a, gR3b, gR2, 2,
                        centLabel[ic], Y_MIN_V2, Y_MAX_V2,
                        /*legendTopRight=*/ (ic == 0 || ic == 1));
        }

        cv->SaveAs("plot_v2_vsq2_Run2vsRun3_May11.pdf");
        Printf("Saved: plot_v2_vsq2_Run2vsRun3_May11.pdf");
        delete cv;
    }

    // =========================================================================
    // v3 comparison: Run3 (0.5-3) + Run3 (1-3)  — no Run2 v3
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

            drawCompPad(cv->GetPad(ic + 1), gR3a, gR3b, nullptr, 3,
                        centLabel[ic], Y_MIN_V3, Y_MAX_V3,
                        /*legendTopRight=*/ (ic == 0 || ic == 1));
        }

        cv->SaveAs("plot_v3_vsq3_Run3_May11.pdf");
        Printf("Saved: plot_v3_vsq3_Run3_May11.pdf");
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
            lat.DrawLatex(0.16, 0.915, "CMS");
            lat.SetTextFont(42);
            lat.SetTextSize(0.052);
            lat.DrawLatex(0.19, 0.83, Form("Cent: %s", centLabel[ic]));
            lat.DrawLatex(0.19, 0.76, "|#eta| < 1.0");
        }

        cv->SaveAs("plot_v2_pTcorr_Run3_May11.pdf");
        Printf("Saved: plot_v2_pTcorr_Run3_May11.pdf");
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
            lat.DrawLatex(0.16, 0.915, "CMS");
            lat.SetTextFont(42);
            lat.SetTextSize(0.052);
            lat.DrawLatex(0.19, 0.83, Form("Cent: %s", centLabel[ic]));
            lat.DrawLatex(0.19, 0.76, "|#eta| < 1.0");
        }

        cv->SaveAs("plot_v3_pTcorr_Run3_May11.pdf");
        Printf("Saved: plot_v3_pTcorr_Run3_May11.pdf");
        delete cv;
    }

    fRun3->Close();
    fRun2->Close();
    Printf("Done.");
}
