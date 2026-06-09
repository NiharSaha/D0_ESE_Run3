// plot_D0diag_and_yield_combined.C
//
// For each (centrality, pT, q-bin) one combined panel:
//   TOP:    chi2/ndf (blue) + Significance (red)  vs SP-bin index
//   BOTTOM: D0 yield                              vs SP-bin index
//   Shared x-axis: SP bin index
//
// One page per (cent, pT) → 10 q-bin cells per page (5 cols × 2 rows)
//
// Input:  ROOT/Flow_MB0to31_May18_out_combined.root
//   directory "diagnostics_chi2_sigma"
//   directory "yield_vs_spbin"
//
// Output: 2 PDFs — one for v2 (45 pages), one for v3 (35 pages)
//
// Usage:
//   root -l -q plot_D0diag_and_yield_combined.C
//   root -l -q 'plot_D0diag_and_yield_combined.C("ROOT/myfile.root","v2.pdf","v3.pdf")'

#include "TFile.h"
#include "TDirectoryFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <utility>

// ============================================================================
// Shared: centrality definitions
// ============================================================================
static const int N_CENT = 5;
static const char* gCentTag  [N_CENT] = {
    "cent0to10","cent10to20","cent20to30","cent30to40","cent40to50"
};
static const char* gCentLabel[N_CENT] = {
    "0-10%","10-20%","20-30%","30-40%","40-50%"
};
static const int gMinCent[N_CENT] = { 0, 10, 20, 30, 40};
static const int gMaxCent[N_CENT] = {10, 20, 30, 40, 50};

// ============================================================================
// Shared: pT bin definitions
// ============================================================================
// --- v2: 9 bins ---
static const int NPT_V2 = 9;
// string tags used by diagnostics
static const char* gPtTagV2  [NPT_V2] = {
    "pT2to3","pT3to4","pT4to5","pT5to6",
    "pT6to8","pT8to10","pT10to15","pT15to30","pT30to100"
};
static const char* gPtLabelV2[NPT_V2] = {
    "2-3","3-4","4-5","5-6",
    "6-8","8-10","10-15","15-30","30-100"
};
// numeric edges used by yield plots
static const double gPtEdgesV2[NPT_V2 + 1] = {
    2, 3, 4, 5, 6, 8, 10, 15, 30, 100
};

// --- v3: 7 bins ---
static const int NPT_V3 = 7;
static const char* gPtTagV3  [NPT_V3] = {
    "pT2to4","pT4to6","pT6to8",
    "pT8to10","pT10to20","pT20to50","pT50to100"
};
static const char* gPtLabelV3[NPT_V3] = {
    "2-4","4-6","6-8",
    "8-10","10-20","20-50","50-100"
};
static const double gPtEdgesV3[NPT_V3 + 1] = {
    2, 4, 6, 8, 10, 20, 50, 100
};

// ============================================================================
// Shared: q-bin count
// ============================================================================
static const int N_QBINS = 10;   // indices 0–9

// ============================================================================
// Global counter for unique object names
// ============================================================================
static int gUID = 0;

void setStyle()
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
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(0);
    gStyle->SetGridWidth(1);
    gStyle->SetGridColor(kGray);
    gStyle->SetGridStyle(3);
}

// ============================================================================
// Draw one combined cell:
//   TOP sub-pad  (upper ~58%): chi2/ndf (blue) + Significance (red) dual axis
//   BOTTOM sub-pad (lower ~42%): D0 yield
//   Shared x-axis range; x-axis labels/title only on bottom panel.
// ============================================================================
void drawCombinedCell(TVirtualPad* cell,
                      TGraph* gChi2, TGraph* gSigma, TH1* hYield,
                      int vn, int iq)
{
    if (!cell) return;
    cell->cd();
    int uid = gUID++;

    // ── Shared x-range ──────────────────────────────────────────────────────
    double xlo = -0.5, xhi = 41.5;
    if (gChi2 && gChi2->GetN() > 0) {
        double xmn = 1e30, xmx = -1e30;
        for (int i = 0; i < gChi2->GetN(); i++) {
            double x, y; gChi2->GetPoint(i, x, y);
            if (y > -5.0) { xmn = std::min(xmn, x); xmx = std::max(xmx, x); }
        }
        if (xmn < xmx) { xlo = xmn - 0.5; xhi = xmx + 0.5; }
    } else if (hYield) {
        xlo = hYield->GetXaxis()->GetXmin();
        xhi = hYield->GetXaxis()->GetXmax();
    }

    const double splitY = 0.43;

    // ── Sub-pads ─────────────────────────────────────────────────────────────
    TPad* pT = new TPad(Form("ptop_%d", uid), "", 0, splitY, 1, 1.0);
    pT->SetFillColor(0); pT->SetBorderMode(0);
    pT->SetLeftMargin(0.17); pT->SetRightMargin(0.17);
    pT->SetTopMargin(0.20);  pT->SetBottomMargin(0.02);
    pT->SetGridx(1); pT->SetGridy(1);
    pT->Draw();

    TPad* pB = new TPad(Form("pbot_%d", uid), "", 0, 0.0, 1, splitY);
    pB->SetFillColor(0); pB->SetBorderMode(0);
    pB->SetLeftMargin(0.17); pB->SetRightMargin(0.17);
    pB->SetTopMargin(0.02);  pB->SetBottomMargin(0.30);
    pB->SetGridx(1);
    pB->Draw();

    // ── TOP PANEL: chi2/ndf + Significance ──────────────────────────────────
    pT->cd();

    bool hasC2  = (gChi2  && gChi2 ->GetN() > 0);
    bool hasSig = (gSigma && gSigma->GetN() > 0);

    if (!hasC2 && !hasSig) {
        TLatex lat; lat.SetNDC(); lat.SetTextFont(42);
        lat.SetTextSize(0.09); lat.SetTextAlign(22);
        lat.DrawLatex(0.5, 0.50, "no data");
    } else {
        int N  = hasC2  ? gChi2 ->GetN() : 0;
        int Ns = hasSig ? gSigma->GetN() : 0;
        std::vector<std::pair<double,double>> chi2v(N), sigv(Ns);
        for (int i = 0; i < N;  i++) { double x,y; gChi2 ->GetPoint(i,x,y); chi2v[i]={x,y}; }
        for (int i = 0; i < Ns; i++) { double x,y; gSigma->GetPoint(i,x,y); sigv [i]={x,y}; }
        std::sort(chi2v.begin(), chi2v.end());
        std::sort(sigv .begin(), sigv .end());

        // chi2/ndf y-range
        double c2lo = -1.0, c2hi = -1e30; bool anyC2 = false;
        for (auto& p : chi2v) if (p.second > -5.0) { c2hi = std::max(c2hi, p.second); anyC2 = true; }
        if (!anyC2) c2hi = c2lo + 1.0;
        else        c2hi += std::max((c2hi - c2lo) * 0.20, 0.30);
        if (c2hi - c2lo < 1e-9) c2hi = c2lo + 1.0;

        // significance y-range
        double slo = 0.0, shi = -1e30; bool anySig2 = false;
        for (auto& p : sigv) if (p.second > -5.0) { shi = std::max(shi, p.second); anySig2 = true; }
        if (!anySig2) shi = slo + 1.0;
        else          shi += std::max((shi - slo) * 0.20, 0.50);
        if (shi - slo < 1e-9) shi = slo + 1.0;

        // axis frame — x-axis labels/ticks suppressed (shared with bottom)
        TH2F* hf = new TH2F(Form("hfcomb_%d", uid), "",
                             100, xlo, xhi, 100, c2lo, c2hi);
        hf->GetXaxis()->SetLabelSize(0); hf->GetXaxis()->SetTickLength(0);
        hf->GetXaxis()->SetNdivisions(510);
        hf->GetYaxis()->SetTitle("#chi^{2}/ndf");
        hf->GetYaxis()->SetLabelSize(0.075); hf->GetYaxis()->SetTitleSize(0.075);
        hf->GetYaxis()->SetTitleOffset(0.90); hf->GetYaxis()->SetNdivisions(505);
        hf->GetYaxis()->SetAxisColor (kBlue+1);
        hf->GetYaxis()->SetTitleColor(kBlue+1);
        hf->GetYaxis()->SetLabelColor(kBlue+1);
        hf->Draw("AXIS");

        // chi2 segments (gaps at sentinel y < -5)
        { std::vector<double> sx, sy; int seg = 0;
          for (int i = 0; i <= N; i++) {
              bool ok = (i < N) && (chi2v[i].second > -5.0);
              if (ok) { sx.push_back(chi2v[i].first); sy.push_back(chi2v[i].second); }
              if (!ok && !sx.empty()) {
                  TGraph* g = new TGraph((int)sx.size());
                  g->SetName(Form("gc2_%d_%d", uid, seg++));
                  for (int j = 0; j < (int)sx.size(); j++) g->SetPoint(j, sx[j], sy[j]);
                  g->SetMarkerStyle(20); g->SetMarkerSize(0.90);
                  g->SetMarkerColor(kBlue+1); g->SetLineColor(kBlue-9); g->SetLineWidth(1);
                  g->Draw((int)sx.size()>1 ? "PL SAME" : "P SAME");
                  sx.clear(); sy.clear();
              }
          }
        }

        pT->Update();
        double uxmax = pT->GetUxmax();
        double uymin = pT->GetUymin();
        double uymax = pT->GetUymax();

        // significance segments (scaled to chi2 y-axis)
        if (anySig2) {
            double rng = shi - slo;
            std::vector<double> sx, sy; int seg = 0;
            for (int i = 0; i <= Ns; i++) {
                bool ok = (i < Ns) && (sigv[i].second > -5.0);
                if (ok) {
                    sx.push_back(sigv[i].first);
                    sy.push_back(uymin + (sigv[i].second - slo) / rng * (uymax - uymin));
                }
                if (!ok && !sx.empty()) {
                    TGraph* g = new TGraph((int)sx.size());
                    g->SetName(Form("gsig_%d_%d", uid, seg++));
                    for (int j = 0; j < (int)sx.size(); j++) g->SetPoint(j, sx[j], sy[j]);
                    g->SetMarkerStyle(22); g->SetMarkerSize(0.90);
                    g->SetMarkerColor(kRed+1); g->SetLineColor(kRed-9); g->SetLineWidth(1);
                    g->Draw((int)sx.size()>1 ? "PL SAME" : "P SAME");
                    sx.clear(); sy.clear();
                }
            }
            TGaxis* rax = new TGaxis(uxmax, uymin, uxmax, uymax, slo, shi, 505, "+L");
            rax->SetLabelFont(42); rax->SetLabelSize(0.075); rax->SetLabelColor(kRed+1);
            rax->SetLineColor(kRed+1); rax->SetTitleFont(42); rax->SetTitleSize(0.075);
            rax->SetTitleColor(kRed+1); rax->SetTitleOffset(0.80);
            rax->SetTitle("Significance"); rax->SetNdivisions(505);
            rax->Draw();
        }
    }

    // q-bin label at top-middle of upper x-axis (drawn last so it sits on top)
    pT->cd();
    TLatex qlat; qlat.SetNDC(); qlat.SetTextFont(42);
    qlat.SetTextSize(0.085); qlat.SetTextAlign(21);
    qlat.DrawLatex(0.50, 0.81, Form("q_{%d} bin = %d", vn, iq));

    // ── BOTTOM PANEL: D0 yield ───────────────────────────────────────────────
    pB->cd();

    if (hYield) {
        TH1* hd = (TH1*)hYield->Clone(Form("hy_%d", uid));
        hd->SetDirectory(nullptr);
        hd->GetXaxis()->SetRangeUser(xlo, xhi);
        hd->SetMarkerStyle(20); hd->SetMarkerSize(0.70);
        hd->SetMarkerColor(kBlack); hd->SetLineColor(kBlack); hd->SetLineWidth(1);
        hd->GetXaxis()->SetTitle(Form("SP_{%d} bin index", vn));
        hd->GetYaxis()->SetTitle("Yield");
        hd->GetXaxis()->SetTitleSize(0.085); hd->GetYaxis()->SetTitleSize(0.085);
        hd->GetXaxis()->SetLabelSize(0.080); hd->GetYaxis()->SetLabelSize(0.080);
        hd->GetXaxis()->SetTitleOffset(0.85); hd->GetYaxis()->SetTitleOffset(0.95);
        hd->GetXaxis()->SetNdivisions(505);  hd->GetYaxis()->SetNdivisions(505);
        hd->GetXaxis()->SetTickLength(0.06);
        hd->Draw("E");
    } else {
        TLatex warn; warn.SetNDC();
        warn.SetTextSize(0.09); warn.SetTextColor(kRed); warn.SetTextAlign(22);
        warn.DrawLatex(0.5, 0.50, "yield not found");
    }
}

// ============================================================================
// Draw one page: (cent ic, pT ip) — 10 q-bin combined cells, 5 cols × 2 rows
// ============================================================================
void drawCombinedPage(TCanvas* cv,
                      TDirectoryFile* dDiag, TDirectoryFile* dYield,
                      int vn, int ic, int ip, const double* pt_edges)
{
    cv->Clear();

    const char* ct     = gCentTag  [ic];
    const char* cl     = gCentLabel[ic];
    const char* pt     = (vn == 2) ? gPtTagV2  [ip] : gPtTagV3  [ip];
    const char* pl     = (vn == 2) ? gPtLabelV2[ip] : gPtLabelV3[ip];
    const char* qTag   = (vn == 2) ? "q2bin"         : "q3bin";
    const char* sigPfx = (vn == 2) ? "sigma_v2_graph" : "sigma_v3_graph";

    // Page header
    TPad* hdrPad = new TPad("hdrpad", "", 0.0, 0.965, 1.0, 1.0);
    hdrPad->SetFillColor(0); hdrPad->SetBorderMode(0); hdrPad->Draw(); hdrPad->cd();
    TLatex ttl; ttl.SetNDC(); ttl.SetTextFont(62);
    ttl.SetTextSize(0.50); ttl.SetTextAlign(22);
    ttl.DrawLatex(0.50, 0.50,
        Form("v_{%d}   Centrality: %s   p_{T}: %s GeV/c", vn, cl, pl));
    cv->cd();

    // Grid pad: 5 cols × 2 rows = 10 cells
    TPad* gridPad = new TPad("gridpad", "", 0.0, 0.0, 1.0, 0.965);
    gridPad->SetFillColor(0); gridPad->SetBorderMode(0); gridPad->Draw();
    gridPad->Divide(5, 2, 0.001, 0.001);

    // pT tag for yield directory lookup (numeric: e.g. "pT10to15")
    TString ptYield = Form("pT%dto%d", (int)pt_edges[ip], (int)pt_edges[ip+1]);

    for (int iq = 0; iq < N_QBINS; iq++) {
        TGraph* gChi2  = nullptr;
        TGraph* gSigma = nullptr;
        if (dDiag) {
            gChi2  = (TGraph*)dDiag->Get(
                Form("chi2ndf_graph_%s_%s_%s_%d", pt, ct, qTag, iq));
            gSigma = (TGraph*)dDiag->Get(
                Form("%s_%s_%s_%s_%d", sigPfx, pt, ct, qTag, iq));
        }
        TH1* hYield = nullptr;
        if (dYield) {
            TString hname = Form("hIdx_v%d_%s_%s_q%dbin%d",
                                 vn, ct, ptYield.Data(), vn, iq);
            hYield = dynamic_cast<TH1*>(dYield->Get(hname));
        }
        drawCombinedCell(gridPad->GetPad(iq + 1), gChi2, gSigma, hYield, vn, iq);
    }

    cv->Modified(); cv->Update();
}

// ============================================================================
// Main
// ============================================================================
void plot_D0diag_and_yield_combined(
    const char* in_file = "ROOT/Flow_MB0to31_May19_out_combined_v3.root",
    const char* out_v2  = "D0_diag_yield_v2_May20.pdf",
    const char* out_v3  = "D0_diag_yield_v3_May20.pdf")
{
    setStyle();
    gROOT->SetBatch(kTRUE);
    gROOT->ForceStyle();

    TFile* f = TFile::Open(in_file);
    if (!f || f->IsZombie()) { Printf("ERROR: cannot open %s", in_file); return; }

    TDirectoryFile* dDiag =
        dynamic_cast<TDirectoryFile*>(f->Get("diagnostics_chi2_sigma"));
    if (!dDiag)
        Printf("WARNING: diagnostics_chi2_sigma not found — diag panels will be blank");

    TDirectoryFile* dYield =
        dynamic_cast<TDirectoryFile*>(f->Get("yield_vs_spbin"));
    if (!dYield)
        Printf("WARNING: yield_vs_spbin not found — yield panels will be blank");

    // Canvas: 5 cols × 2 rows of combined cells
    const int CW = 5 * 340;
    const int CH = 2 * 500 + 40;
    TCanvas* cv = new TCanvas("cv", "D0 diag + yield", CW, CH);

    // v2: 5 cent × 9 pT = 45 pages
    cv->Print(Form("%s[", out_v2));
    for (int ic = 0; ic < N_CENT; ic++)
        for (int ip = 0; ip < NPT_V2; ip++) {
            Printf("[v2]  cent %s  pT %s", gCentLabel[ic], gPtLabelV2[ip]);
            drawCombinedPage(cv, dDiag, dYield, 2, ic, ip, gPtEdgesV2);
            cv->Print(out_v2);
        }
    cv->Print(Form("%s]", out_v2));
    Printf("Saved %s  (%d pages)", out_v2, N_CENT * NPT_V2);

    // v3: 5 cent × 7 pT = 35 pages
    cv->Print(Form("%s[", out_v3));
    for (int ic = 0; ic < N_CENT; ic++)
        for (int ip = 0; ip < NPT_V3; ip++) {
            Printf("[v3]  cent %s  pT %s", gCentLabel[ic], gPtLabelV3[ip]);
            drawCombinedPage(cv, dDiag, dYield, 3, ic, ip, gPtEdgesV3);
            cv->Print(out_v3);
        }
    cv->Print(Form("%s]", out_v3));
    Printf("Saved %s  (%d pages)", out_v3, N_CENT * NPT_V3);

    delete cv;
    f->Close();

    Printf("\n========================================");
    Printf("  Output files");
    Printf("========================================");
    Printf("  [v2]  %s  (%d pages)", out_v2, N_CENT * NPT_V2);
    Printf("  [v3]  %s  (%d pages)", out_v3, N_CENT * NPT_V3);
    Printf("========================================\n");
    Printf("Done.");
}

// ---------------------------------------------------------------------------
// Stand-alone compiled entry point
// ---------------------------------------------------------------------------
#ifndef __CINT__
#ifndef __CLING__
int main() { plot_D0diag_and_yield_combined(); return 0; }
#endif
#endif
