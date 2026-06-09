// plot_scatter_Run2vsRun3_compare.C
//
// Compares D0 v2 vs charged-particle v2 scatter plots (Run2 vs Run3)
// on the same canvas, for each (centrality × D0 pT bin).
//
// Run2 inputs:
//   ROOT/v2_vs_pT_for_q2_bins_5_pT_bins.root
//     v2_graph_cen{min}_{max}_q2bin_{iq}  — TGraphErrors, 5 pT pts per q-bin
//     D0 pT bins: 2-4, 4-6, 6-10, 10-15, 15-30  (midpoints at indices 0-4)
//   ROOT/ch_v2_vs_q2.root
//     _new_cen_{min}_{max}  — TH1D, 10 bins (q2-bin 0-9), chg pT 1-3 GeV/c
//
// Run3 inputs:
//   ROOT/Flow_MB0to31_May19_out_combined_v3.root
//     vn_vs_qbin/v2_vs_q2bin_cent{X}to{Y}_pT{lo}to{hi}  — TGraphErrors
//     D0 pT bins: 2-3, 3-4, 4-5, 5-6, 6-8, 8-10, 10-15, 15-30
//   ROOT/flow_Analysis_chg_out_combined_May11.root
//     vsQbin_TProfile/hp_v2_vsq2_{pt_tag}_{cent}  — TProfile (10 q2-bins)
//     pt_tag: pt0p5to3  or  pt1to3
//
// Mapping Run3 → Run2 pT bins (for overlay):
//   R3 2-3, 3-4 → R2 2-4 (index 0)
//   R3 4-5, 5-6 → R2 4-6 (index 1)
//   R3 6-8, 8-10 → R2 6-10 (index 2)
//   R3 10-15    → R2 10-15 (index 3)
//   R3 15-30    → R2 15-30 (index 4)
//
// Produces two merged PDFs (5 pages each, 4 cols × 2 rows per page):
//   scatter_Run2vsRun3_pt0p5to3_{tag}.pdf  — Run3 chg pT 0.5-3
//   scatter_Run2vsRun3_pt1to3_{tag}.pdf    — Run3 chg pT 1-3
//
// Usage:
//   root -l -b -q 'plot_scatter_Run2vsRun3_compare.C+'

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVirtualFitter.h"
#include "TROOT.h"

void plot_scatter_Run2vsRun3_compare(
    const char* run2_d0_file  = "ROOT/v2_vs_pT_for_q2_bins_5_pT_bins.root",
    const char* run2_chg_file = "ROOT/ch_v2_vs_q2.root",
    const char* run3_d0_file  = "ROOT/Flow_MB0to31_May19_out_combined_v3.root",
    const char* run3_chg_file = "ROOT/flow_Analysis_chg_out_combined_May11.root",
    const char* tag           = "May25_v2",
    
    // If true, use the Run3 file with exactly Run2-like 5 pT bins for comparison instead of the original Run3 file
    bool use_run3_run2ptbin   = false, 
    const char* run3_run2ptbin_file = "ROOT/Flow_MB0to31_ForRun2Comparison_out_combined.root"
    )
{
    // =========================================================================
    // Constants
    // =========================================================================
    const int N_QBINS = 10;
    const int N_CENT  = 5;
    const int minCent[N_CENT] = {  0, 10, 20, 30, 40 };
    const int maxCent[N_CENT] = { 10, 20, 30, 40, 50 };
    const char* cent_r3[N_CENT] = {
        "cent0to10","cent10to20","cent20to30","cent30to40","cent40to50"
    };


    // --- pT binning logic ---
    // Run2 binning is always 5 bins
    const int NPT_R2 = 5;
    const char* r2_ptTag[5] = {"pT2to4","pT4to6","pT6to10","pT10to15","pT15to30"};
    double r2_ptLo[5] = {2,4,6,10,15};
    double r2_ptHi[5] = {4,6,10,15,30};
    TFile* fR2d0 = TFile::Open(run2_d0_file);

    // Run3 binning: 8 bins by default, 5 bins if use_run3_run2ptbin
    int NPT_R3 = 8;
    const char* r3_ptTag_default[8] = {"pT2to3","pT3to4","pT4to5","pT5to6","pT6to8","pT8to10","pT10to15","pT15to30"};
    double r3_ptLo_default[8] = {2,3,4,5,6,8,10,15};
    double r3_ptHi_default[8] = {3,4,5,6,8,10,15,30};
    const char* r3_ptTag_run2bins[5] = {"pT2to4","pT4to6","pT6to10","pT10to15","pT15to30"};
    double r3_ptLo_run2bins[5] = {2,4,6,10,15};
    double r3_ptHi_run2bins[5] = {4,6,10,15,30};
    const char** r3_ptTag = r3_ptTag_default;
    double* r3_ptLo = r3_ptLo_default;
    double* r3_ptHi = r3_ptHi_default;
    TFile* fR3d0 = nullptr;
    if (use_run3_run2ptbin) {
        NPT_R3 = 5;
        r3_ptTag = r3_ptTag_run2bins;
        r3_ptLo = r3_ptLo_run2bins;
        r3_ptHi = r3_ptHi_run2bins;
        fR3d0 = TFile::Open(run3_run2ptbin_file);
    } else {
        fR3d0 = TFile::Open(run3_d0_file);
    }

    // Charged pT configurations
    struct ChgCfg { const char* tag; const char* label; };
    const ChgCfg chg_cfgs[2] = {
        { "pt0p5to3", "0.5-3" },
        { "pt1to3",   "1-3"   }
    };

    // =========================================================================
    // Style
    // =========================================================================
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // =========================================================================
    // Open files
    // =========================================================================
    // fR2d0 and fR3d0 are now opened above
    TFile* fR2chg = TFile::Open(run2_chg_file);
    TFile* fR3chg = TFile::Open(run3_chg_file);

    const char* used_run3_d0_file = use_run3_run2ptbin ? run3_run2ptbin_file : run3_d0_file;
    if (!fR2d0  || fR2d0->IsZombie())  { std::cerr << "ERROR: " << run2_d0_file  << "\n"; return; }
    if (!fR2chg || fR2chg->IsZombie()) { std::cerr << "ERROR: " << run2_chg_file << "\n"; return; }
    if (!fR3d0  || fR3d0->IsZombie())  { std::cerr << "ERROR: " << used_run3_d0_file  << "\n"; return; }
    if (!fR3chg || fR3chg->IsZombie()) { std::cerr << "ERROR: " << run3_chg_file << "\n"; return; }

    auto* dir_r3_qbin = (TDirectory*)fR3d0->Get("vn_vs_qbin");
    auto* dir_r3_chg  = (TDirectory*)fR3chg->Get("vsQbin_TProfile");
    if (!dir_r3_qbin) { std::cerr << "ERROR: vn_vs_qbin dir not found in Run3 D0 file\n";  return; }
    if (!dir_r3_chg)  { std::cerr << "ERROR: vsQbin_TProfile dir not found in Run3 chg file\n"; return; }

    // Strict verification: when using the Run3 file prepared with Run2-like pT bins,
    // ensure it actually contains the expected 5 pT-tagged objects. Abort if missing.
    if (use_run3_run2ptbin) {
        bool missing = false;
        for (int ip = 0; ip < NPT_R3; ++ip) {
            // check for the first centrality as a quick existence probe
            TString chk = Form("v2_vs_q2bin_%s_%s", cent_r3[0], r3_ptTag[ip]);
            if (!dir_r3_qbin->Get(chk)) {
                std::cerr << "ERROR: Expected object not found in " << used_run3_d0_file << ": " << chk << "\n";
                missing = true;
            }
        }
        if (missing) {
            std::cerr << "ERROR: Run3 (Run2-ptbins) file missing expected objects — aborting.\n";
            return;
        }
    }

    // =========================================================================
    // Pre-load Run2 data
    // =========================================================================
    // gr_r2q[ic][iq]: TGraphErrors with NPT_R2 points (D0 v2 vs pT) for q-bin iq
    TGraphErrors* gr_r2q[N_CENT][N_QBINS] = {};
    for (int ic = 0; ic < N_CENT; ++ic) {
        TString ctag = Form("cen%d_%d", minCent[ic], maxCent[ic]);
        for (int iq = 0; iq < N_QBINS; ++iq) {
            // Run2 graphs always have the same naming: v2_graph_<cen>_q2bin_<iq>
            TString nm = Form("v2_graph_%s_q2bin_%d", ctag.Data(), iq);
            gr_r2q[ic][iq] = (TGraphErrors*)fR2d0->Get(nm);
            if (!gr_r2q[ic][iq])
                std::cerr << "  WARN: " << nm << " not found\n";
        }
    }

    // h_r2chg[ic]: TH1D with 10 bins, bin i+1 = charged v2 at q2-bin i, chg pT 1-3
    TH1D* h_r2chg[N_CENT] = {};
    for (int ic = 0; ic < N_CENT; ++ic) {
        TString nm = Form("_new_cen_%d_%d", minCent[ic], maxCent[ic]);
        h_r2chg[ic] = (TH1D*)fR2chg->Get(nm);
        if (!h_r2chg[ic])
            std::cerr << "  WARN: " << nm << " not found\n";
    }

    // =========================================================================
    // Pre-load Run3 data
    // =========================================================================
    // gr_r3d0[ic][ip]: D0 v2 vs q2-bin (N_QBINS points)
    std::vector<std::vector<TGraphErrors*>> gr_r3d0(N_CENT);
    for (int ic = 0; ic < N_CENT; ++ic) {
        gr_r3d0[ic].resize(NPT_R3, nullptr);
        for (int ip = 0; ip < NPT_R3; ++ip) {
            TString nm = Form("v2_vs_q2bin_%s_%s", cent_r3[ic], r3_ptTag[ip]);
            gr_r3d0[ic][ip] = (TGraphErrors*)dir_r3_qbin->Get(nm);
            if (!gr_r3d0[ic][ip])
                std::cerr << "  WARN: " << nm << " not found\n";
        }
    }

    // hp_r3chg[ichg][ic]: charged v2 TProfile (N_QBINS bins), ichg=0→pt0p5to3, 1→pt1to3
    TProfile* hp_r3chg[2][N_CENT] = {};
    for (int ichg = 0; ichg < 2; ++ichg) {
        for (int ic = 0; ic < N_CENT; ++ic) {
            TString nm = Form("hp_v2_vsq2_%s_%s", chg_cfgs[ichg].tag, cent_r3[ic]);
            hp_r3chg[ichg][ic] = (TProfile*)dir_r3_chg->Get(nm);
            if (!hp_r3chg[ichg][ic])
                std::cerr << "  WARN: " << nm << " not found\n";
        }
    }

    // =========================================================================
    // Output directory
    // =========================================================================
    TString plotDir = "scatter_compare";
    gSystem->mkdir(plotDir, kTRUE);

    TLatex ltx;
    ltx.SetNDC();
    ltx.SetTextFont(42);

    // =========================================================================
    // PRE-SCAN: per-centrality y-axis range (scan all R3 pT bins + Run2 matches)
    // D0 v2 is independent of the charged pT range, so compute once.
    // =========================================================================
    // Fixed y-axis range for all plots
    const double Y_LO = -0.15, Y_HI = 0.5;

    // =========================================================================
    // Main loop: one canvas per (chg config × centrality × R3 pT bin)
    // =========================================================================
    // Fixed x-axis range for all plots
    const double X_LO = 0.0, X_HI = 0.30;
    for (int ichg = 0; ichg < 2; ++ichg)
    {
        const char* chg_tag   = chg_cfgs[ichg].tag;
        const char* chg_label = chg_cfgs[ichg].label;
        std::cout << "\n=== Scatter comparison: Run3 chg pT " << chg_label << " GeV/c ===\n";

        for (int ic = 0; ic < N_CENT; ++ic)
        {
            // Extract Run2 charged v2 per q-bin (pT 1-3, fixed)
            double cv2_r2[N_QBINS] = {}, ce2_r2[N_QBINS] = {};
            if (h_r2chg[ic]) {
                for (int iq = 0; iq < N_QBINS; ++iq) {
                    cv2_r2[iq] = h_r2chg[ic]->GetBinContent(iq + 1);
                    ce2_r2[iq] = h_r2chg[ic]->GetBinError(iq + 1);
                }
            }

            // Extract Run3 charged v2 per q-bin for this chg config
            double cv2_r3[N_QBINS] = {}, ce2_r3[N_QBINS] = {};
            if (hp_r3chg[ichg][ic]) {
                for (int iq = 0; iq < N_QBINS; ++iq) {
                    if (hp_r3chg[ichg][ic]->GetBinEntries(iq + 1) > 0) {
                        cv2_r3[iq] = hp_r3chg[ichg][ic]->GetBinContent(iq + 1);
                        ce2_r3[iq] = hp_r3chg[ichg][ic]->GetBinError(iq + 1);
                    }
                }
            }



            // iterate pT bins according to the Run3 pT definition
            // (when use_run3_run2ptbin==true, NPT_R3 == NPT_R2 == 5)
            int nPtBins = NPT_R3;
            for (int ip = 0; ip < nPtBins; ++ip)
            {
                int ir2 = 0;
                if (use_run3_run2ptbin) {
                    // Run3 uses same 5 Run2-like pT bins — direct index mapping
                    ir2 = ip;
                } else {
                    // Explicit 8->5 mapping table for clarity and robustness
                    static const int map_r3_to_r2[8] = {0,0,1,1,2,2,3,4};
                    if (ip >= 0 && ip < 8) ir2 = map_r3_to_r2[ip];
                    else ir2 = 0;
                }

                // Build Run2 scatter: (x=chg_v2, y=D0_v2) per q-bin
                std::vector<double> x2, ex2, y2, ey2;
                for (int iq = 0; iq < N_QBINS; ++iq) {
                    if (!gr_r2q[ic][iq] || ir2 >= gr_r2q[ic][iq]->GetN()) continue;
                    double dv = gr_r2q[ic][iq]->GetPointY(ir2);
                    double de = gr_r2q[ic][iq]->GetErrorY(ir2);
                    double cv = cv2_r2[iq], ce = ce2_r2[iq];
                    if (de <= 0 || ce <= 0) continue;
                    x2.push_back(cv); ex2.push_back(ce);
                    y2.push_back(dv); ey2.push_back(de);
                }

                // Build Run3 scatter: (x=chg_v2, y=D0_v2) per q-bin
                std::vector<double> x3, ex3, y3, ey3;
                if (gr_r3d0[ic][ip]) {
                    for (int iq = 0; iq < N_QBINS && iq < gr_r3d0[ic][ip]->GetN(); ++iq) {
                        double dv = gr_r3d0[ic][ip]->GetPointY(iq);
                        double de = gr_r3d0[ic][ip]->GetErrorY(iq);
                        double cv = cv2_r3[iq], ce = ce2_r3[iq];
                        if (de <= 0 || ce <= 0) continue;
                        x3.push_back(cv); ex3.push_back(ce);
                        y3.push_back(dv); ey3.push_back(de);
                    }
                }

                if ((int)x2.size() < 2 && (int)x3.size() < 2) continue;

                // Axis ranges (x: 0-0.14 for 0-10% centrality, else 0-0.3; y: fixed)
                double xlo = X_LO, xhi = X_HI;
                if (minCent[ic] == 0 && maxCent[ic] == 10) {
                    xhi = 0.14;
                }
                const double ylo = Y_LO, yhi = Y_HI;

                // Canvas
                TString baseName = Form("scatter_comp_cent%dto%d_pT%dto%d_%s",
                    minCent[ic], maxCent[ic],
                    (int)r3_ptLo[ip], (int)r3_ptHi[ip], chg_tag);
                TCanvas* cv = new TCanvas(baseName, baseName, 600, 600);
                cv->SetLeftMargin(0.14);
                cv->SetRightMargin(0.05);
                cv->SetBottomMargin(0.14);
                cv->SetTopMargin(0.08);
                cv->SetGridx();
                cv->SetGridy();

                // Frame
                TH2F* hf = new TH2F(Form("hf_%s", baseName.Data()), "",
                    10, xlo, xhi, 10, ylo, yhi);
                hf->GetXaxis()->SetTitle("charged particle v_{2}");
                hf->GetYaxis()->SetTitle("D^{0} v_{2}");
                hf->GetXaxis()->SetTitleSize(0.050);
                hf->GetYaxis()->SetTitleSize(0.050);
                hf->GetXaxis()->SetTitleOffset(1.10);
                hf->GetYaxis()->SetTitleOffset(1.25);
                hf->GetXaxis()->SetLabelSize(0.040);
                hf->GetYaxis()->SetLabelSize(0.040);
                hf->Draw();

                // y = 0 reference line
                TLine* zero = new TLine(xlo, 0.0, xhi, 0.0);
                zero->SetLineColor(kGray + 1);
                zero->SetLineStyle(2);
                zero->Draw();

                // Run2 scatter + linear fit
                TGraphErrors* gr2  = nullptr;
                TF1*          f2   = nullptr;
                TH1D*         hci2 = nullptr;
                if ((int)x2.size() >= 2) {
                    gr2 = new TGraphErrors((int)x2.size(),
                        x2.data(), y2.data(), ex2.data(), ey2.data());
                    gr2->SetName(Form("gr2_%s", baseName.Data()));

                    f2 = new TF1(Form("f2_%s", baseName.Data()), "[0]+[1]*x", xlo, xhi);
                    f2->SetParameters(0.0, 1.0);
                    f2->SetLineColor(kBlue + 1);
                    f2->SetLineWidth(2);
                    gr2->Fit(f2, "QR");

                    hci2 = new TH1D(Form("hci2_%s", baseName.Data()), "", 400, xlo, xhi);
                    hci2->SetDirectory(nullptr);
                    TVirtualFitter* vf = TVirtualFitter::GetFitter();
                    if (vf) vf->GetConfidenceIntervals(hci2, 0.68);
                    hci2->SetFillColorAlpha(kBlue + 1, 0.20);
                    hci2->SetFillStyle(1001);
                    hci2->SetMarkerSize(0);
                    hci2->Draw("E3 SAME");
                    f2->Draw("SAME");

                    gr2->SetMarkerStyle(20);
                    gr2->SetMarkerSize(1.0);
                    gr2->SetMarkerColor(kBlue + 1);
                    gr2->SetLineColor(kBlue + 1);
                    gr2->Draw("P SAME");
                }

                // Run3 scatter + linear fit
                TGraphErrors* gr3  = nullptr;
                TF1*          f3   = nullptr;
                TH1D*         hci3 = nullptr;
                if ((int)x3.size() >= 2) {
                    gr3 = new TGraphErrors((int)x3.size(),
                        x3.data(), y3.data(), ex3.data(), ey3.data());
                    gr3->SetName(Form("gr3_%s", baseName.Data()));

                    f3 = new TF1(Form("f3_%s", baseName.Data()), "[0]+[1]*x", xlo, xhi);
                    f3->SetParameters(0.0, 1.0);
                    f3->SetLineColor(kRed + 1);
                    f3->SetLineWidth(2);
                    gr3->Fit(f3, "QR");

                    hci3 = new TH1D(Form("hci3_%s", baseName.Data()), "", 400, xlo, xhi);
                    hci3->SetDirectory(nullptr);
                    TVirtualFitter* vf = TVirtualFitter::GetFitter();
                    if (vf) vf->GetConfidenceIntervals(hci3, 0.68);
                    hci3->SetFillColorAlpha(kRed + 1, 0.20);
                    hci3->SetFillStyle(1001);
                    hci3->SetMarkerSize(0);
                    hci3->Draw("E3 SAME");
                    f3->Draw("SAME");

                    gr3->SetMarkerStyle(21);
                    gr3->SetMarkerSize(1.0);
                    gr3->SetMarkerColor(kRed + 1);
                    gr3->SetLineColor(kRed + 1);
                    gr3->Draw("P SAME");
                }

                // CMS label
                ltx.SetTextAlign(11);
                ltx.SetTextFont(62);
                ltx.SetTextSize(0.055);
                ltx.DrawLatex(0.14, 0.935, "CMS");
                ltx.SetTextFont(42);
                ltx.SetTextSize(0.040);
                ltx.SetTextAlign(31);
                ltx.DrawLatex(0.95, 0.935, "PbPb");

                // Cent label (top-left)
                ltx.SetTextAlign(11);
                ltx.SetTextSize(0.040);
                ltx.DrawLatex(0.16, 0.87,
                    Form("%d < Cent. < %d%%", minCent[ic], maxCent[ic]));

                // Legend: Run2/Run3 entries below cent label (top-left)
                TLegend* leg = new TLegend(0.14, 0.73, 0.68, 0.85);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.036);
                // Add Run3 first, then Run2 (match slope/intercept order)
                if (gr3) leg->AddEntry(gr3,
                    Form("Run3 5.36 TeV, D^{0} p_{T} %.0f-%.0f GeV/c",
                         r3_ptLo[ip], r3_ptHi[ip]), "pe");
                if (gr2) leg->AddEntry(gr2,
                    Form("Run2 5.02 TeV, D^{0} p_{T} %.0f-%.0f GeV/c",
                         r2_ptLo[ir2], r2_ptHi[ir2]), "pe");
                leg->Draw();

                // Fit parameters (right-aligned, bottom-right, colour-coded)
                ltx.SetTextAlign(31);
                ltx.SetTextSize(0.034);
                if (f3) {
                    ltx.SetTextColor(kRed + 1);
                    ltx.DrawLatex(0.93, 0.24,
                        Form("Run3: slope=%.3f#pm%.3f, int=%.3f#pm%.3f",
                             f3->GetParameter(1), f3->GetParError(1),
                             f3->GetParameter(0), f3->GetParError(0)));
                    ltx.SetTextColor(kBlack);
                }
                if (f2) {
                    ltx.SetTextColor(kBlue + 1);
                    ltx.DrawLatex(0.93, 0.17,
                        Form("Run2: slope=%.3f#pm%.3f, int=%.3f#pm%.3f",
                             f2->GetParameter(1), f2->GetParError(1),
                             f2->GetParameter(0), f2->GetParError(0)));
                    ltx.SetTextColor(kBlack);
                }

                TString outPath = Form("%s/%s.pdf", plotDir.Data(), baseName.Data());
                cv->SaveAs(outPath);

                delete leg; delete hf; delete zero;
                if (gr2)  delete gr2;
                if (f2)   delete f2;
                if (hci2) delete hci2;
                if (gr3)  delete gr3;
                if (f3)   delete f3;
                if (hci3) delete hci3;
                delete cv;

            } // ip (pT bin)

            std::cout << "  Centrality " << minCent[ic] << "-" << maxCent[ic]
                      << "%: " << NPT_R3 << " plots saved\n";
        } // ic (centrality)

        // =====================================================================
        // Tile individual PDFs into 5-page summary (4 cols × 2 rows)
        // =====================================================================
        TString out_pdf = Form("scatter_Run2vsRun3_%s_%s.pdf", chg_tag, tag);
        std::cout << "\nTiling PDFs -> " << out_pdf << " ...\n";

        TString py_tmp = Form("_tile_scatter_%s.py", chg_tag);
        FILE* fpy = fopen(py_tmp.Data(), "w");
        if (fpy) {
            fprintf(fpy,
                "import sys, os\n"
                "try:\n"
                "    import fitz\n"
                "except ModuleNotFoundError:\n"
                "    sys.exit('ERROR: PyMuPDF not found. pip install pymupdf')\n"
                "\n"
                "plot_dir  = '%s'\n"
                "chg_tag   = '%s'\n"
                "out_pdf   = '%s'\n"
                "page_title = 'chg p_{T}: 1-3 GeV/c (Run2)  vs  %s GeV/c (Run3)'\n"
                "CENT_ORDER = ['0to10','10to20','20to30','30to40','40to50']\n"
                "NCOLS = 4\n"
                "TITLE_H = 36  # points for title bar\n"
                "\n"
                "def pt_key(f):\n"
                "    return int(f.split('_pT')[1].split('to')[0])\n"
                "\n"
                "files = sorted(\n"
                "    [f for f in os.listdir(plot_dir) if f.endswith('.pdf') and chg_tag in f],\n"
                "    key=pt_key)\n"
                "\n"
                "groups = {}\n"
                "for f in files:\n"
                "    try: c = f.split('_cent')[1].split('_')[0]\n"
                "    except: continue\n"
                "    groups.setdefault(c, []).append(f)\n"
                "\n"
                "if not files:\n"
                "    print('No PDFs found'); sys.exit(0)\n"
                "\n"
                "with fitz.open(os.path.join(plot_dir, files[0])) as s:\n"
                "    w, h = s[0].rect.width, s[0].rect.height\n"
                "\n"
                "doc = fitz.open()\n"
                "for cent in [c for c in CENT_ORDER if c in groups]:\n"
                "    fl = sorted(groups[cent], key=pt_key)\n"
                "    rows = -(-len(fl) // NCOLS)\n"
                "    page = doc.new_page(width=w*NCOLS, height=h*rows + TITLE_H)\n"
                "    # Title bar\n"
                "    page.draw_rect(fitz.Rect(0, 0, w*NCOLS, TITLE_H),\n"
                "                  color=(1,1,1), fill=(1,1,1))\n"
                "    label = page_title + '    Cent. ' + cent.replace('to','-') + '%%'\n"
                "    page.insert_text((8, TITLE_H - 8), label,\n"
                "                    fontsize=16, color=(0,0,0))\n"
                "    # Plot grid (offset below title bar)\n"
                "    for i, fn in enumerate(fl[:NCOLS*rows]):\n"
                "        ci, ri = i%%NCOLS, i//NCOLS\n"
                "        with fitz.open(os.path.join(plot_dir, fn)) as src:\n"
                "            pix = src[0].get_pixmap(dpi=150)\n"
                "        page.insert_image(\n"
                "            fitz.Rect(ci*w, TITLE_H + ri*h, (ci+1)*w, TITLE_H + (ri+1)*h),\n"
                "            pixmap=pix)\n"
                "    print(f'  cent{cent}: {len(fl)} plots')\n"
                "doc.save(out_pdf, garbage=4, deflate=True)\n"
                "doc.close()\n"
                "print('Saved ->', out_pdf)\n",
                plotDir.Data(), chg_tag, out_pdf.Data(), chg_label);
            fclose(fpy);

            // Find Python with PyMuPDF
            TString py_exe;
            const char* cands[] = {
                "/usr/local/bin/python3.11",
                "/opt/homebrew/bin/python3",
                "/usr/bin/python3",
                "python3", "python", nullptr
            };
            for (int k = 0; cands[k] && py_exe.IsNull(); ++k)
                if (gSystem->Exec(Form("%s -c 'import fitz' 2>/dev/null", cands[k])) == 0)
                    py_exe = cands[k];

            if (!py_exe.IsNull())
                gSystem->Exec(Form("%s '%s'", py_exe.Data(), py_tmp.Data()));
            else
                std::cerr << "WARNING: PyMuPDF not found; merged PDF skipped.\n";

            gSystem->Unlink(py_tmp.Data());
        }

    } // ichg

    fR2d0->Close();
    fR2chg->Close();
    fR3d0->Close();
    fR3chg->Close();

    std::cout << "\nDone!\n"
              << "  Individual plots: " << plotDir << "/\n"
              << Form("  scatter_Run2vsRun3_pt0p5to3_%s.pdf  (5 pages)\n", tag)
              << Form("  scatter_Run2vsRun3_pt1to3_%s.pdf    (5 pages)\n", tag);
}

#ifndef __CINT__
#ifndef __CLING__
int main() { plot_scatter_Run2vsRun3_compare(); return 0; }
#endif
#endif
