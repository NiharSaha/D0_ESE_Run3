//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// plot_ESE_effW_comparison.C
//
// Standalone macro comparing efficiency-corrected ("with eff") vs
// uncorrected ("without eff") results, overlaid with different colors on
// the same axes. Two comparison categories:
//
//   1) D0 vn vs q-bin-index, with-eff vs without-eff
//      -- D0 has NO shared "_effW" directory: the two variants live in two
//         separate input files (INPUT_D0_WOEFF_FILE / INPUT_D0_WEFF_FILE).
//         D0 depends on its own pT bin, so this is produced per
//         (centrality x D0 pT).
//   2) Charged-particle vn vs q-bin-index, with-eff vs without-eff
//      -- read from the two directories ("vsQbin_TProfile" vs
//         "vsQbin_TProfile_effW") of a single shared charged-particle file
//         (INPUT_CHG_FILE). The charged-particle vn is pT-integrated
//         (0.5-3 or 1-3 GeV/c) and does not depend on D0's pT bin, so this
//         is produced per centrality only (no pT dependence to tile over).
//   3) D0-vs-charged-particle scatter, with-eff vs without-eff
//      -- two full scatter point sets per (centrality x D0 pT): one built
//         entirely from the "without eff" file pair, one entirely from the
//         "with eff" file pair (D0 y from D0_WEFF/WOEFF, charged x from
//         vsQbin_TProfile[_effW] of the same shared charged file).
//
// All 6 centralities are covered (5 differential + the inclusive 0-50% bin,
// "cent0to50"), matching the naming already used in
// plot_ESE_vnVsqn_scatter_and_pearson_combined.C.
//
// Individual per-(cent[,pT]) PDFs are tiled into multi-page/1-page summary
// PDFs via tile_effW_comparison.py (PyMuPDF), following the same tiling
// approach as tile_vn_vs_qbin.py.
//
// Outputs, all under a single tagged output directory:
//   {outdir}/1_vn_vs_qbin_D0_comparison/*.pdf
//   {outdir}/1_vn_vs_qbin_D0_comparison/summary_D0comp_v{2,3}.pdf   (tiled)
//   {outdir}/2_vn_vs_qbin_Charged_comparison/*.pdf
//   {outdir}/2_vn_vs_qbin_Charged_comparison/summary_ChgComp_v{2,3}.pdf (tiled, 1 page)
//   {outdir}/3_scatter_comparison/*.pdf
//   {outdir}/3_scatter_comparison/summary_scatterComp_v{2,3}.pdf   (tiled)
//
// Usage:
//   root -l -b -q plot_ESE_effW_comparison.C
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

// #########################################################################
// ####  INPUT FILES / ARGUMENTS  —  edit these, no need to scroll down  ####
// #########################################################################

static const char* INPUT_D0_WOEFF_FILE = "D0_Flow_woEff_out_combined_Sept25.root";
static const char* INPUT_D0_WEFF_FILE  = "D0_Flow_wEff_combined_Sept29.root";
// Single file containing both "vsQbin_TProfile" (w/o eff) and
// "vsQbin_TProfile_effW" (with eff) charged-particle directories.
static const char* INPUT_CHG_FILE      = "Charge_Flow_wANDwo_Eff_combined_Sept25.root";

// true -> charged particle 0.5 < pT < 3 GeV/c (Run3); false -> 1 < pT < 3 GeV/c (Run2)
static const bool  INPUT_IS_RUN3 = true;

static const std::string INPUT_TAG = "Sept29_v1";
static const std::string INPUT_OUTDIR_TAG_STR = "ESE_effW_comparison_" + INPUT_TAG;
static const char* INPUT_OUTDIR_TAG = INPUT_OUTDIR_TAG_STR.c_str();

// Inclusive centrality (0-50%) tag, matching the naming already used in
// plot_ESE_vnVsqn_scatter_and_pearson_combined.C.
static const char* INCL_CENT_TAG = "cent0to50";

// =========================================================================
// Bin definitions (shared by all 3 comparison categories)
// =========================================================================
static const int   N_CENTBINS = 6;   // 5 differential + 1 inclusive
static const char* CENT_TAG  [N_CENTBINS] = {
    "cent0to10","cent10to20","cent20to30","cent30to40","cent40to50", INCL_CENT_TAG
};
static const int   MIN_CENT[N_CENTBINS] = { 0, 10, 20, 30, 40,  0};
static const int   MAX_CENT[N_CENTBINS] = {10, 20, 30, 40, 50, 50};

static const int    N_PTBINS_V2 = 9;
static const double PT_EDGES_V2[N_PTBINS_V2 + 1] = { 2, 3, 4, 5, 6, 8, 10, 15, 30, 100 };
static const int    N_PTBINS_V3 = 6;
static const double PT_EDGES_V3[N_PTBINS_V3 + 1] = { 2, 4, 6, 8, 10, 20, 50 };

static const int N_QBINS = 12;

// Marker/color style: black circle = without eff, red square = with eff
static const int NOEFF_COLOR  = kBlack;
static const int NOEFF_MARKER = 20;
static const int EFF_COLOR    = kRed + 1;
static const int EFF_MARKER   = 21;

// =========================================================================
// Helper: locate a Python interpreter that has PyMuPDF (fitz) installed.
// =========================================================================
static TString FindPythonWithFitz()
{
    const char* cands[] = {
        "/usr/local/bin/python3.11",
        "/opt/homebrew/bin/python3",
        "/usr/bin/python3",
        "python3",
        "python",
        nullptr
    };
    for (int k = 0; cands[k]; ++k) {
        TString probe = Form("%s -c 'import fitz' 2>/dev/null", cands[k]);
        if (gSystem->Exec(probe.Data()) == 0) return TString(cands[k]);
    }
    return TString();
}

static void setCompStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    TGaxis::SetMaxDigits(2);
}

static void DrawCMSHeader(double leftMargin, double rightMargin, double y)
{
    TLatex cmsLat;
    cmsLat.SetNDC();
    cmsLat.SetTextFont(42);
    cmsLat.SetTextSize(0.048);
    cmsLat.SetTextAlign(11);
    cmsLat.DrawLatex(leftMargin + 0.08, y, "#font[61]{CMS} #font[52]{Preliminary}");

    TLatex enLat;
    enLat.SetNDC();
    enLat.SetTextFont(42);
    enLat.SetTextSize(0.048);
    enLat.SetTextAlign(31);
    enLat.DrawLatex(1.0 - rightMargin, y, "PbPb #sqrt{s_{NN}} = 5.36 TeV");
}

// =========================================================================
// SECTION 1: D0 vn vs q-bin-index, with-eff vs without-eff
// Per (centrality x D0 pT). Reads "vn_vs_qbin/v{ivn}_vs_q{ivn}bin_{cent}_{pT}"
// from the two separate D0 files (no shared _effW directory for D0).
// =========================================================================
static void RunD0Comparison(const TString& outbase, TFile* fNo, TFile* fEf)
{
    TDirectory* dNo = fNo ? (TDirectory*)fNo->Get("vn_vs_qbin") : nullptr;
    TDirectory* dEf = fEf ? (TDirectory*)fEf->Get("vn_vs_qbin") : nullptr;
    if (!dNo) std::cerr << "WARNING: D0 (w/o eff) vn_vs_qbin dir not found\n";
    if (!dEf) std::cerr << "WARNING: D0 (with eff) vn_vs_qbin dir not found\n";

    TString baseDir = outbase + "/1_vn_vs_qbin_D0_comparison";
    gSystem->mkdir(baseDir, kTRUE);

    TCanvas* cv = new TCanvas("cv_d0comp", "", 700, 550);
    TLatex ltx; ltx.SetNDC(); ltx.SetTextFont(42);

    for (int ivn = 2; ivn <= 3; ++ivn) {
        int           n_pt     = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
        const double* pt_edges = (ivn == 2) ? PT_EDGES_V2 : PT_EDGES_V3;

        for (int ic = 0; ic < N_CENTBINS; ++ic) {
            for (int ip = 0; ip < n_pt; ++ip) {
                const double pt_lo = pt_edges[ip];
                const double pt_hi = pt_edges[ip + 1];
                TString ptTag = Form("pT%dto%d", (int)pt_lo, (int)pt_hi);
                TString grName = Form("v%d_vs_q%dbin_%s_%s", ivn, ivn, CENT_TAG[ic], ptTag.Data());

                TGraphErrors* grNo = dNo ? (TGraphErrors*)dNo->Get(grName) : nullptr;
                TGraphErrors* grEf = dEf ? (TGraphErrors*)dEf->Get(grName) : nullptr;
                if ((!grNo || grNo->GetN() == 0) && (!grEf || grEf->GetN() == 0)) continue;

                double xNo[N_QBINS] = {}, yNo[N_QBINS] = {}, eyNo[N_QBINS] = {};
                double xEf[N_QBINS] = {}, yEf[N_QBINS] = {}, eyEf[N_QBINS] = {};
                int nNo = 0, nEf = 0;
                double ylo =  1e9, yhi = -1e9;

                if (grNo) {
                    for (int iq = 0; iq < N_QBINS && iq < grNo->GetN(); ++iq) {
                        double v = grNo->GetPointY(iq), e = grNo->GetErrorY(iq);
                        if (e <= 0) continue;
                        xNo[nNo] = iq; yNo[nNo] = v; eyNo[nNo] = e; ++nNo;
                        ylo = std::min(ylo, v - e); yhi = std::max(yhi, v + e);
                    }
                }
                if (grEf) {
                    for (int iq = 0; iq < N_QBINS && iq < grEf->GetN(); ++iq) {
                        double v = grEf->GetPointY(iq), e = grEf->GetErrorY(iq);
                        if (e <= 0) continue;
                        xEf[nEf] = iq; yEf[nEf] = v; eyEf[nEf] = e; ++nEf;
                        ylo = std::min(ylo, v - e); yhi = std::max(yhi, v + e);
                    }
                }
                if (nNo == 0 && nEf == 0) continue;
                if (ylo > yhi) { ylo = 0.0; yhi = 0.3; }
                double marg = 0.25 * std::max(yhi - ylo, 1e-6);
                ylo -= marg; yhi += marg;

                cv->Clear();
                cv->SetLeftMargin(0.14); cv->SetRightMargin(0.05);
                cv->SetBottomMargin(0.14); cv->SetTopMargin(0.08);
                cv->SetGridx(); cv->SetGridy();

                TH2F* hf = new TH2F("hf_d0comp", "", 10, -0.5, N_QBINS - 0.5, 10, ylo, yhi);
                hf->GetXaxis()->SetTitle(Form("q_{%d} bin index", ivn));
                hf->GetYaxis()->SetTitle(ivn == 2 ? "v_{2}" : "v_{3}");
                hf->GetXaxis()->SetTitleSize(0.060); hf->GetYaxis()->SetTitleSize(0.060);
                hf->GetXaxis()->SetTitleOffset(0.95); hf->GetYaxis()->SetTitleOffset(0.90);
                hf->GetXaxis()->SetLabelSize(0.038); hf->GetYaxis()->SetLabelSize(0.038);
                hf->GetXaxis()->SetNdivisions(N_QBINS);
                hf->Draw();

                TGraphErrors* gNo = nullptr, * gEf = nullptr;
                if (nNo) {
                    gNo = new TGraphErrors(nNo, xNo, yNo, nullptr, eyNo);
                    gNo->SetMarkerStyle(NOEFF_MARKER); gNo->SetMarkerSize(1.2);
                    gNo->SetMarkerColor(NOEFF_COLOR); gNo->SetLineColor(NOEFF_COLOR); gNo->SetLineWidth(2);
                    gNo->Draw("P SAME");
                }
                if (nEf) {
                    gEf = new TGraphErrors(nEf, xEf, yEf, nullptr, eyEf);
                    gEf->SetMarkerStyle(EFF_MARKER); gEf->SetMarkerSize(1.2);
                    gEf->SetMarkerColor(EFF_COLOR); gEf->SetLineColor(EFF_COLOR); gEf->SetLineWidth(2);
                    gEf->Draw("P SAME");
                }

                TLegend* leg = new TLegend(0.50, 0.78, 0.86, 0.90);
                leg->SetBorderSize(0); leg->SetFillStyle(0);
                leg->SetTextFont(42); leg->SetTextSize(0.038);
                if (gNo) leg->AddEntry(gNo, "D^{0}, without eff", "pe");
                if (gEf) leg->AddEntry(gEf, "D^{0}, with eff", "pe");
                leg->Draw();

                ltx.SetTextAlign(11); ltx.SetTextSize(0.038);
                ltx.DrawLatex(0.17, 0.86, Form("Centrality: %d-%d%%", MIN_CENT[ic], MAX_CENT[ic]));
                ltx.DrawLatex(0.17, 0.79, Form("%.0f < p_{T}^{D^{0}} < %.0f GeV/c", pt_lo, pt_hi));

                DrawCMSHeader(0.14, 0.05, 0.945);

                cv->Modified(); cv->Update();
                TString base = Form("vn_vs_qbin_D0comp_v%d_%s_%s", ivn, CENT_TAG[ic], ptTag.Data());
                cv->SaveAs(Form("%s/%s.pdf", baseDir.Data(), base.Data()));

                delete leg; delete hf; delete gNo; delete gEf;
            } // ip
        } // ic
        std::cout << "D0 comparison plots complete for v" << ivn << "\n";
    } // ivn

    delete cv;
}

// =========================================================================
// SECTION 2: Charged-particle vn vs q-bin-index, with-eff vs without-eff
// Per centrality only (charged vn is pT-integrated, no D0-pT dependence).
// Reads "hp_v{ivn}_vsq{ivn}_{chgPtTag}_{cent}" from vsQbin_TProfile[_effW]
// in the single shared charged-particle file.
// =========================================================================
static void RunChargedComparison(const TString& outbase, TFile* fChg, bool isRun3)
{
    TDirectory* dNo = fChg ? (TDirectory*)fChg->Get("vsQbin_TProfile")      : nullptr;
    TDirectory* dEf = fChg ? (TDirectory*)fChg->Get("vsQbin_TProfile_effW") : nullptr;
    if (!dNo) std::cerr << "WARNING: CHG vsQbin_TProfile dir not found\n";
    if (!dEf) std::cerr << "WARNING: CHG vsQbin_TProfile_effW dir not found\n";

    const char* chgPtTag   = isRun3 ? "pt0p5to3" : "pt1to3";
    const char* chgPtRange = isRun3 ? "0.5-3"    : "1-3";

    TString baseDir = outbase + "/2_vn_vs_qbin_Charged_comparison";
    gSystem->mkdir(baseDir, kTRUE);

    TCanvas* cv = new TCanvas("cv_chgcomp", "", 700, 550);
    TLatex ltx; ltx.SetNDC(); ltx.SetTextFont(42);

    for (int ivn = 2; ivn <= 3; ++ivn) {
        for (int ic = 0; ic < N_CENTBINS; ++ic) {
            TString hname = Form("hp_v%d_vsq%d_%s_%s", ivn, ivn, chgPtTag, CENT_TAG[ic]);
            TProfile* pNo = dNo ? (TProfile*)dNo->Get(hname) : nullptr;
            TProfile* pEf = dEf ? (TProfile*)dEf->Get(hname) : nullptr;
            if (!pNo && !pEf) continue;

            double xNo[N_QBINS] = {}, yNo[N_QBINS] = {}, eyNo[N_QBINS] = {};
            double xEf[N_QBINS] = {}, yEf[N_QBINS] = {}, eyEf[N_QBINS] = {};
            int nNo = 0, nEf = 0;
            double ylo =  1e9, yhi = -1e9;

            if (pNo) {
                for (int iq = 0; iq < N_QBINS; ++iq) {
                    if (pNo->GetBinEntries(iq + 1) <= 0) continue;
                    double v = pNo->GetBinContent(iq + 1), e = pNo->GetBinError(iq + 1);
                    if (e <= 0) continue;
                    xNo[nNo] = iq; yNo[nNo] = v; eyNo[nNo] = e; ++nNo;
                    ylo = std::min(ylo, v - e); yhi = std::max(yhi, v + e);
                }
            }
            if (pEf) {
                for (int iq = 0; iq < N_QBINS; ++iq) {
                    if (pEf->GetBinEntries(iq + 1) <= 0) continue;
                    double v = pEf->GetBinContent(iq + 1), e = pEf->GetBinError(iq + 1);
                    if (e <= 0) continue;
                    xEf[nEf] = iq; yEf[nEf] = v; eyEf[nEf] = e; ++nEf;
                    ylo = std::min(ylo, v - e); yhi = std::max(yhi, v + e);
                }
            }
            if (nNo == 0 && nEf == 0) continue;
            if (ylo > yhi) { ylo = 0.0; yhi = 0.3; }
            double marg = 0.25 * std::max(yhi - ylo, 1e-6);
            ylo = std::max(0.0, ylo - marg); yhi += marg;

            cv->Clear();
            cv->SetLeftMargin(0.14); cv->SetRightMargin(0.05);
            cv->SetBottomMargin(0.14); cv->SetTopMargin(0.08);
            cv->SetGridx(); cv->SetGridy();

            TH2F* hf = new TH2F("hf_chgcomp", "", 10, -0.5, N_QBINS - 0.5, 10, ylo, yhi);
            hf->GetXaxis()->SetTitle(Form("q_{%d} bin index", ivn));
            hf->GetYaxis()->SetTitle(ivn == 2 ? "v_{2}" : "v_{3}");
            hf->GetXaxis()->SetTitleSize(0.060); hf->GetYaxis()->SetTitleSize(0.060);
            hf->GetXaxis()->SetTitleOffset(0.95); hf->GetYaxis()->SetTitleOffset(0.90);
            hf->GetXaxis()->SetLabelSize(0.038); hf->GetYaxis()->SetLabelSize(0.038);
            hf->GetXaxis()->SetNdivisions(N_QBINS);
            hf->Draw();

            TGraphErrors* gNo = nullptr, * gEf = nullptr;
            if (nNo) {
                gNo = new TGraphErrors(nNo, xNo, yNo, nullptr, eyNo);
                gNo->SetMarkerStyle(NOEFF_MARKER); gNo->SetMarkerSize(1.2);
                gNo->SetMarkerColor(NOEFF_COLOR); gNo->SetLineColor(NOEFF_COLOR); gNo->SetLineWidth(2);
                gNo->Draw("P SAME");
            }
            if (nEf) {
                gEf = new TGraphErrors(nEf, xEf, yEf, nullptr, eyEf);
                gEf->SetMarkerStyle(EFF_MARKER); gEf->SetMarkerSize(1.2);
                gEf->SetMarkerColor(EFF_COLOR); gEf->SetLineColor(EFF_COLOR); gEf->SetLineWidth(2);
                gEf->Draw("P SAME");
            }

            TLegend* leg = new TLegend(0.50, 0.78, 0.86, 0.90);
            leg->SetBorderSize(0); leg->SetFillStyle(0);
            leg->SetTextFont(42); leg->SetTextSize(0.038);
            if (gNo) leg->AddEntry(gNo, "charged, without eff", "pe");
            if (gEf) leg->AddEntry(gEf, "charged, with eff", "pe");
            leg->Draw();

            ltx.SetTextAlign(11); ltx.SetTextSize(0.038);
            ltx.DrawLatex(0.17, 0.86, Form("Centrality: %d-%d%%", MIN_CENT[ic], MAX_CENT[ic]));
            ltx.DrawLatex(0.17, 0.79, Form("p_{T} : %s GeV/c", chgPtRange));

            DrawCMSHeader(0.14, 0.05, 0.945);

            cv->Modified(); cv->Update();
            TString base = Form("vn_vs_qbin_ChgComp_v%d_%s", ivn, CENT_TAG[ic]);
            cv->SaveAs(Form("%s/%s.pdf", baseDir.Data(), base.Data()));

            delete leg; delete hf; delete gNo; delete gEf;
        } // ic
        std::cout << "Charged-particle comparison plots complete for v" << ivn << "\n";
    } // ivn

    delete cv;
}

// =========================================================================
// SECTION 3: D0-vs-charged-particle scatter, with-eff vs without-eff
// Per (centrality x D0 pT). Two full point sets: "without eff" built
// entirely from D0_WOEFF + vsQbin_TProfile, "with eff" entirely from
// D0_WEFF + vsQbin_TProfile_effW.
// =========================================================================
static void RunScatterComparison(const TString& outbase, TFile* fD0No, TFile* fD0Ef,
                                  TFile* fChg, bool isRun3)
{
    TDirectory* dD0No  = fD0No ? (TDirectory*)fD0No->Get("vn_vs_qbin")        : nullptr;
    TDirectory* dD0Ef  = fD0Ef ? (TDirectory*)fD0Ef->Get("vn_vs_qbin")        : nullptr;
    TDirectory* dChgNo = fChg  ? (TDirectory*)fChg->Get("vsQbin_TProfile")      : nullptr;
    TDirectory* dChgEf = fChg  ? (TDirectory*)fChg->Get("vsQbin_TProfile_effW") : nullptr;

    const char* chgPtTag = isRun3 ? "pt0p5to3" : "pt1to3";

    TString baseDir = outbase + "/3_scatter_comparison";
    gSystem->mkdir(baseDir, kTRUE);

    TCanvas* cv = new TCanvas("cv_sccomp", "", 700, 650);
    TLatex ltx; ltx.SetNDC(); ltx.SetTextFont(42);

    for (int ivn = 2; ivn <= 3; ++ivn) {
        int           n_pt     = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
        const double* pt_edges = (ivn == 2) ? PT_EDGES_V2 : PT_EDGES_V3;

        for (int ic = 0; ic < N_CENTBINS; ++ic) {
            TString hname = Form("hp_v%d_vsq%d_%s_%s", ivn, ivn, chgPtTag, CENT_TAG[ic]);
            TProfile* pNo = dChgNo ? (TProfile*)dChgNo->Get(hname) : nullptr;
            TProfile* pEf = dChgEf ? (TProfile*)dChgEf->Get(hname) : nullptr;

            double chgNo_v[N_QBINS] = {}, chgNo_e[N_QBINS] = {};
            double chgEf_v[N_QBINS] = {}, chgEf_e[N_QBINS] = {};
            for (int iq = 0; iq < N_QBINS; ++iq) {
                if (pNo && pNo->GetBinEntries(iq + 1) > 0) {
                    chgNo_v[iq] = pNo->GetBinContent(iq + 1); chgNo_e[iq] = pNo->GetBinError(iq + 1);
                }
                if (pEf && pEf->GetBinEntries(iq + 1) > 0) {
                    chgEf_v[iq] = pEf->GetBinContent(iq + 1); chgEf_e[iq] = pEf->GetBinError(iq + 1);
                }
            }

            for (int ip = 0; ip < n_pt; ++ip) {
                const double pt_lo = pt_edges[ip];
                const double pt_hi = pt_edges[ip + 1];
                TString ptTag = Form("pT%dto%d", (int)pt_lo, (int)pt_hi);
                TString grName = Form("v%d_vs_q%dbin_%s_%s", ivn, ivn, CENT_TAG[ic], ptTag.Data());

                TGraphErrors* grD0No = dD0No ? (TGraphErrors*)dD0No->Get(grName) : nullptr;
                TGraphErrors* grD0Ef = dD0Ef ? (TGraphErrors*)dD0Ef->Get(grName) : nullptr;
                if ((!grD0No || grD0No->GetN() == 0) && (!grD0Ef || grD0Ef->GetN() == 0)) continue;

                std::vector<double> xNo, yNo, exNo, eyNo;
                std::vector<double> xEf, yEf, exEf, eyEf;
                double xlo = 1e9, xhi = -1e9, ylo = 1e9, yhi = -1e9;

                if (grD0No) {
                    for (int iq = 0; iq < N_QBINS && iq < grD0No->GetN(); ++iq) {
                        double cv_ = chgNo_v[iq], ce = chgNo_e[iq];
                        double dv = grD0No->GetPointY(iq), de = grD0No->GetErrorY(iq);
                        if (ce <= 0 || de <= 0) continue;
                        xNo.push_back(cv_); exNo.push_back(ce);
                        yNo.push_back(dv);  eyNo.push_back(de);
                        xlo = std::min(xlo, cv_ - ce); xhi = std::max(xhi, cv_ + ce);
                        ylo = std::min(ylo, dv - de);  yhi = std::max(yhi, dv + de);
                    }
                }
                if (grD0Ef) {
                    for (int iq = 0; iq < N_QBINS && iq < grD0Ef->GetN(); ++iq) {
                        double cv_ = chgEf_v[iq], ce = chgEf_e[iq];
                        double dv = grD0Ef->GetPointY(iq), de = grD0Ef->GetErrorY(iq);
                        if (ce <= 0 || de <= 0) continue;
                        xEf.push_back(cv_); exEf.push_back(ce);
                        yEf.push_back(dv);  eyEf.push_back(de);
                        xlo = std::min(xlo, cv_ - ce); xhi = std::max(xhi, cv_ + ce);
                        ylo = std::min(ylo, dv - de);  yhi = std::max(yhi, dv + de);
                    }
                }
                if (xNo.size() < 1 && xEf.size() < 1) continue;
                if (xlo > xhi) { xlo = 0.0; xhi = 0.3; }
                if (ylo > yhi) { ylo = 0.0; yhi = 0.3; }
                double xmarg = 0.15 * std::max(xhi - xlo, 1e-6);
                double ymarg = 0.15 * std::max(yhi - ylo, 1e-6);
                xlo -= xmarg; xhi += xmarg; ylo -= ymarg; yhi += ymarg;

                cv->Clear();
                cv->SetLeftMargin(0.14); cv->SetRightMargin(0.05);
                cv->SetBottomMargin(0.14); cv->SetTopMargin(0.07);
                cv->SetGridx(); cv->SetGridy();

                TH2F* hf = new TH2F("hf_sccomp", "", 10, xlo, xhi, 10, ylo, yhi);
                hf->GetXaxis()->SetTitle(Form("charged particle v_{%d}", ivn));
                hf->GetYaxis()->SetTitle(Form("D^{0} v_{%d}", ivn));
                hf->GetXaxis()->SetTitleSize(0.050); hf->GetYaxis()->SetTitleSize(0.050);
                hf->GetXaxis()->SetTitleOffset(1.10); hf->GetYaxis()->SetTitleOffset(1.25);
                hf->GetXaxis()->SetLabelSize(0.040); hf->GetYaxis()->SetLabelSize(0.040);
                hf->Draw();

                TGraphErrors* gNo = nullptr, * gEf = nullptr;
                if (!xNo.empty()) {
                    gNo = new TGraphErrors((int)xNo.size(), xNo.data(), yNo.data(), exNo.data(), eyNo.data());
                    gNo->SetMarkerStyle(NOEFF_MARKER); gNo->SetMarkerSize(1.0);
                    gNo->SetMarkerColor(NOEFF_COLOR); gNo->SetLineColor(NOEFF_COLOR);
                    gNo->Draw("P SAME");
                }
                if (!xEf.empty()) {
                    gEf = new TGraphErrors((int)xEf.size(), xEf.data(), yEf.data(), exEf.data(), eyEf.data());
                    gEf->SetMarkerStyle(EFF_MARKER); gEf->SetMarkerSize(1.0);
                    gEf->SetMarkerColor(EFF_COLOR); gEf->SetLineColor(EFF_COLOR);
                    gEf->Draw("P SAME");
                }

                TLegend* leg = new TLegend(0.16, 0.71, 0.58, 0.92);
                leg->SetBorderSize(0); leg->SetFillStyle(0);
                leg->SetTextFont(42); leg->SetTextSize(0.036);
                if (gNo) leg->AddEntry(gNo, "without eff", "pe");
                if (gEf) leg->AddEntry(gEf, "with eff", "pe");
                leg->AddEntry((TObject*)nullptr, Form("Centrality: %d-%d%%", MIN_CENT[ic], MAX_CENT[ic]), "");
                leg->AddEntry((TObject*)nullptr, Form("%.0f < p_{T}^{D^{0}} < %.0f GeV/c", pt_lo, pt_hi), "");
                leg->Draw();

                DrawCMSHeader(0.14, 0.05, 0.945);

                cv->Modified(); cv->Update();
                TString base = Form("scatter_comp_v%d_%s_%s", ivn, CENT_TAG[ic], ptTag.Data());
                cv->SaveAs(Form("%s/%s.pdf", baseDir.Data(), base.Data()));

                delete leg; delete hf; delete gNo; delete gEf;
            } // ip
        } // ic
        std::cout << "Scatter comparison plots complete for v" << ivn << "\n";
    } // ivn

    delete cv;
}

// =========================================================================
// Main function
// =========================================================================
void plot_ESE_effW_comparison(
    const char* d0_woeff_file = INPUT_D0_WOEFF_FILE,
    const char* d0_weff_file  = INPUT_D0_WEFF_FILE,
    const char* chg_file      = INPUT_CHG_FILE,
    const char* outdir_tag    = INPUT_OUTDIR_TAG,
    bool        isRun3        = INPUT_IS_RUN3)
{
    setCompStyle();
    gROOT->ForceStyle();
    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);

    TString outbase = outdir_tag;
    std::cout << "\n=== Output directory: " << outbase << " ===\n";
    gSystem->mkdir(outbase, kTRUE);

    TFile* fD0No = TFile::Open(d0_woeff_file);
    TFile* fD0Ef = TFile::Open(d0_weff_file);
    TFile* fChg  = TFile::Open(chg_file);
    if (!fD0No || fD0No->IsZombie()) { std::cerr << "ERROR: cannot open D0 (w/o eff) file: " << d0_woeff_file << "\n"; return; }
    if (!fD0Ef || fD0Ef->IsZombie()) { std::cerr << "ERROR: cannot open D0 (with eff) file: " << d0_weff_file << "\n"; return; }
    if (!fChg  || fChg->IsZombie())  { std::cerr << "ERROR: cannot open charged-particle file: " << chg_file << "\n"; return; }

    std::cout << "\n=== SECTION 1: D0 vn-vs-qbin comparison (with vs without eff) ===\n";
    RunD0Comparison(outbase, fD0No, fD0Ef);

    std::cout << "\n=== SECTION 2: Charged-particle vn-vs-qbin comparison (with vs without eff) ===\n";
    RunChargedComparison(outbase, fChg, isRun3);

    std::cout << "\n=== SECTION 3: D0-vs-charged scatter comparison (with vs without eff) ===\n";
    RunScatterComparison(outbase, fD0No, fD0Ef, fChg, isRun3);

    fD0No->Close();
    fD0Ef->Close();
    fChg->Close();

    // =========================================================================
    // Tile individual PDFs into summary PDFs via Python/PyMuPDF
    // =========================================================================
    TString py_exe = FindPythonWithFitz();
    if (py_exe.IsNull()) {
        std::cerr << "WARNING: No Python interpreter with PyMuPDF found; "
                     "summary PDFs will be skipped. Install: pip install pymupdf\n";
    } else {
        std::cout << "\n=== Tiling summary PDFs (Python/PyMuPDF) ===\n";
        gSystem->Exec(Form("%s tile_effW_comparison.py d0 %s/1_vn_vs_qbin_D0_comparison %s/1_vn_vs_qbin_D0_comparison/summary_D0comp",
                            py_exe.Data(), outbase.Data(), outbase.Data()));
        gSystem->Exec(Form("%s tile_effW_comparison.py chg %s/2_vn_vs_qbin_Charged_comparison %s/2_vn_vs_qbin_Charged_comparison/summary_ChgComp",
                            py_exe.Data(), outbase.Data(), outbase.Data()));
        gSystem->Exec(Form("%s tile_effW_comparison.py scatter %s/3_scatter_comparison %s/3_scatter_comparison/summary_scatterComp",
                            py_exe.Data(), outbase.Data(), outbase.Data()));
    }

    std::cout << "\n=== Done ===\n"
              << "  Output directory   -> " << outbase << "/\n"
              << "  D0 comparison      -> " << outbase << "/1_vn_vs_qbin_D0_comparison/  (summary_D0comp_v{2,3}.pdf)\n"
              << "  Charged comparison -> " << outbase << "/2_vn_vs_qbin_Charged_comparison/  (summary_ChgComp_v{2,3}.pdf)\n"
              << "  Scatter comparison -> " << outbase << "/3_scatter_comparison/  (summary_scatterComp_v{2,3}.pdf)\n";
}
