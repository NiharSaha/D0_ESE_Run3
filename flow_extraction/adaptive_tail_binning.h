#ifndef ADAPTIVE_TAIL_BINNING_H
#define ADAPTIVE_TAIL_BINNING_H

// ============================================================
// adaptive_tail_binning.h
//
// Significance-driven adaptive tail binning for D0 v2/v3 ESE.
// Algorithm:
//   Pass 1 - fit inclusive mass histogram (no vn cut) per
//            (i_cen, i_pt) with a double-Gaussian + linear bkg
//            to extract stable reference signal shape.
//   Pass 2 - for each initial tail bin, fix the shape from
//            Pass 1 and float only the overall normalization p0.
//            Significance = p0 / dp0  (same definition as
//            fit_mass_and_flow.C: yield/yield_error).
//            Merge adjacent bins until sig >= minSig or
//            maxMerge bins have been consumed.
//
// Produces per-(i_cen,i_pt):
//   - Significance-vs-bin-center PDF
//   - Mass-fit multi-pad PDF (one pad per final bin)
//   - Text log (edges, per-bin significance, merge flags)
// ============================================================

#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TString.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TMath.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

// ============================================================
// Reference signal shape extracted from inclusive Pass-1 fit
// ============================================================
struct MassShapeParams {
    double mean   = 1.8648; // D0 mass (GeV)
    double sigma1 = 0.030;  // wide Gaussian sigma
    double sigma2 = 0.010;  // narrow Gaussian sigma
    double frac   = 0.5;    // fraction of narrow Gaussian ([4] in fit)
    bool   valid  = false;
};

// ============================================================
// Per-bin diagnostic record filled during adaptive merging
// ============================================================
struct TailBinDiagnostic {
    double lo          = 0;   // vn left edge of (merged) bin
    double hi          = 0;   // vn right edge of (merged) bin
    double entries     = 0;   // raw entries in mass histogram
    double norm        = 0;   // fitted p0
    double norm_err    = 0;   // fitted dp0
    double significance = 0;  // p0 / dp0
    bool   merged      = false; // was merging required?
    int    n_merged    = 1;   // how many initial bins were merged
};

// ============================================================
// Internal counter helpers to avoid ROOT name collisions
// ============================================================
namespace ATB_internal {
    static int ref_counter  = 0;
    static int fit_counter  = 0;
}

// ============================================================
// PASS 1 - inclusive double-Gaussian + linear-bkg mass fit
//
// baseCuts : all cuts EXCEPT mass window (e.g. dca, y, cent, pT)
//            NO vn cut — uses full statistics for shape
// Returns  : MassShapeParams with valid=true if fit is good
// ============================================================
static MassShapeParams get_reference_shape_DG(TNtuple *nt,
                                               TString  baseCuts,
                                               bool     verbose = false)
{
    MassShapeParams p;
    const double fit_lo = 1.74, fit_hi = 2.0;

    TString hname = TString::Format("href_dg_%d",
                                    ATB_internal::ref_counter++);
    TH1D *href = new TH1D(hname, "", 130, fit_lo, fit_hi);
    href->SetDirectory(nullptr);
    nt->Draw(TString::Format("mass>>%s", hname.Data()), baseCuts, "goff");

    if (href->GetEntries() < 50) {
        std::cerr << "[ATB WARN] get_reference_shape_DG: only "
                  << href->GetEntries()
                  << " entries — using defaults\n";
        delete href;
        return p; // valid = false
    }

    // [0] = overall norm
    // [1] = mean  (D0 mass)
    // [2] = sigma1  (wide Gaussian)
    // [3] = sigma2  (narrow Gaussian)
    // [4] = fraction of narrow Gaussian
    // [5] = bkg p0 (const)
    // [6] = bkg p1 (linear)
    // [7] = bkg p2 (quadratic)
    // [8] = bkg p3 (cubic)
    TF1 *ftot = new TF1(TString::Format("ftot_ref_%d",
                            ATB_internal::ref_counter),
        "[0]*( [4]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])"
        "    +(1-[4])*TMath::Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]) )"
        "+ [5] + [6]*x + [7]*x*x + [8]*x*x*x",
        fit_lo, fit_hi);

    ftot->SetParameters(href->GetMaximum() * 0.5,
                        1.8648,
                        0.030, 0.010, 0.5,
                        href->GetMinimum(), 0.0, 0.0, 0.0);
    ftot->SetParLimits(0, 0,     1e9);
    ftot->SetParLimits(2, 0.010, 0.100);  // sigma1 wide
    ftot->SetParLimits(3, 0.003, 0.025);  // sigma2 narrow
    ftot->SetParLimits(4, 0,     1);      // Gaussian fraction in [0,1]

    // two-pass fit: chi2 seed then likelihood polish
    href->Fit(ftot, "Q",   "", fit_lo, fit_hi);
    href->Fit(ftot, "LQM", "", fit_lo, fit_hi);

    p.mean   = ftot->GetParameter(1);
    p.sigma1 = ftot->GetParameter(2);
    p.sigma2 = ftot->GetParameter(3);
    p.frac   = ftot->GetParameter(4);

    // basic sanity checks
    p.valid = (p.mean   > 1.840 && p.mean   < 1.890)
           && (p.sigma1 > 0.005 && p.sigma1 < 0.100)
           && (p.sigma2 > 0.003 && p.sigma2 < 0.025)
           && (p.frac   > 0.0   && p.frac   < 1.0);

    if (verbose || !p.valid)
        std::cout << "[ATB REF FIT]"
                  << " mean="  << p.mean
                  << " s1="    << p.sigma1
                  << " s2="    << p.sigma2
                  << " frac="  << p.frac
                  << (p.valid ? "  [OK]" : "  [FAILED - using defaults]")
                  << "\n";

    delete href;
    delete ftot;
    return p;
}

// ============================================================
// PASS 2 - single-bin significance fit
//
// Fixes signal shape from Pass 1; floats only p0 (norm) and
// the two linear background parameters.
// Significance = p0 / dp0  — consistent with fit_mass_and_flow.C
//
// diag   : optional pointer filled with fit diagnostics
// hout   : if non-null, caller takes ownership of the histogram
//          with the fit function attached (for plotting)
// ============================================================
static double compute_tail_significance(TNtuple              *nt,
                                         TString               baseCuts,
                                         const char           *vnBranch,
                                         double                lo,
                                         double                hi,
                                         const MassShapeParams &shape,
                                         TailBinDiagnostic    *diag    = nullptr,
                                         TH1D                **hout    = nullptr,
                                         bool                  verbose = false)
{
    const double fit_lo = 1.74, fit_hi = 2.0;

    // build full selection string
    TString fullCuts = TString::Format(
        "%s && %s>=%g && %s<%g",
        baseCuts.Data(), vnBranch, lo, vnBranch, hi);

    int idx = ATB_internal::fit_counter++;
    TString hname = TString::Format("htmp_sig_%d", idx);
    TH1D *htmp = new TH1D(hname, "", 130, fit_lo, fit_hi);
    htmp->SetDirectory(nullptr);
    nt->Draw(TString::Format("mass>>%s", hname.Data()), fullCuts, "goff");

    double entries = htmp->GetEntries();

    if (entries < 5) {
        if (verbose)
            std::cout << "[ATB ADAPT] " << vnBranch
                      << " [" << lo << "," << hi << "]"
                      << " entries=" << (int)entries << " -> sig=0\n";
        if (diag) {
            diag->lo = lo; diag->hi = hi;
            diag->entries = entries;
            diag->norm = diag->norm_err = diag->significance = 0;
        }
        if (hout) *hout = htmp;
        else      delete htmp;
        return 0.0;
    }

    // [0] = norm  (floating)
    // [1] = mean  (FIXED from Pass 1)
    // [2] = sigma1 wide  (FIXED)
    // [3] = sigma2 narrow (FIXED)
    // [4] = narrow fraction (FIXED)
    // [5] = bkg p0 const   (floating)
    // [6] = bkg p1 linear  (floating)
    // [7] = bkg p2 quadratic (floating)
    // [8] = bkg p3 cubic   (floating)
    TF1 *ffit = new TF1(TString::Format("ffit_sig_%d", idx),
        "[0]*( [4]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])"
        "    +(1-[4])*TMath::Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]) )"
        "+ [5] + [6]*x + [7]*x*x + [8]*x*x*x",
        fit_lo, fit_hi);

    ffit->SetParameter (0, htmp->GetMaximum() * 0.5);
    ffit->FixParameter (1, shape.mean);
    ffit->FixParameter (2, shape.sigma1);
    ffit->FixParameter (3, shape.sigma2);
    ffit->FixParameter (4, shape.frac);
    ffit->SetParameter (5, 0.0);
    ffit->SetParameter (6, 0.0);
    ffit->SetParameter (7, 0.0);
    ffit->SetParameter (8, 0.0);
    ffit->SetParLimits (0, 0, 1e9);

    // two-pass fit: chi2 seed then likelihood polish
    htmp->Fit(ffit, "Q",   "", fit_lo, fit_hi);
    htmp->Fit(ffit, "LQM", "", fit_lo, fit_hi);

    double norm  = ffit->GetParameter(0);
    double dnorm = ffit->GetParError(0);

    // significance: identical to fit_mass_and_flow.C definition
    // (p5 and bin-width cancel since shape is fixed -> reduces to p0/dp0)
    double sig = (dnorm > 0 && norm > 0) ? norm / dnorm : 0.0;

    if (verbose)
        std::cout << "[ATB ADAPT] " << vnBranch
                  << " [" << lo << ", " << hi << "]"
                  << "  N=" << (int)entries
                  << "  p0=" << norm << " +/- " << dnorm
                  << "  sig=" << sig << "\n";

    if (diag) {
        diag->lo           = lo;
        diag->hi           = hi;
        diag->entries      = entries;
        diag->norm         = norm;
        diag->norm_err     = dnorm;
        diag->significance = sig;
    }

    if (hout) {
        // attach fit to histogram so caller can draw it
        htmp->GetListOfFunctions()->Add(ffit);
        *hout = htmp;
    } else {
        delete htmp;
        delete ffit;
    }

    return sig;
}

// ============================================================
// Adaptive tail edge builder
//
// init_edges : initial guess edges (e.g. vnbinning_side_v2[0..11])
// diags      : optional vector filled with per-final-bin diagnostics
// minSig     : minimum significance threshold (default 3.0)
// maxMerge   : max consecutive initial bins to merge before giving up
//
// Returns a (possibly shorter) edge list where every bin achieves
// significance >= minSig.
// ============================================================
static std::vector<double> adaptive_tail_edges(
    TNtuple                        *nt,
    TString                         baseCuts,
    const char                     *vnBranch,
    const std::vector<double>      &init_edges,
    const MassShapeParams          &shape,
    std::vector<TailBinDiagnostic> *diags    = nullptr,
    std::vector<TH1D*>             *hcache   = nullptr,  // cached histos for plotting
    double                          minSig   = 3.0,
    int                             maxMerge = 4,
    bool                            verbose  = false)
{
    int nEdges = (int)init_edges.size();
    if (nEdges < 2) return init_edges;

    if (diags) diags->clear();

    std::vector<double> result;
    result.push_back(init_edges.front());

    int i = 0;
    while (i < nEdges - 1) {
        int    j          = i + 1; // current right-edge candidate
        int    mergeCount = 0;
        double sig        = 0.0;

        while (j < nEdges) {
            TailBinDiagnostic d;
            TH1D *htmp = nullptr;
            sig = compute_tail_significance(nt, baseCuts, vnBranch,
                                            init_edges[i], init_edges[j],
                                            shape, &d,
                                            hcache ? &htmp : nullptr,
                                            verbose);
            d.merged   = (mergeCount > 0);
            d.n_merged = mergeCount + 1;

            if (sig >= minSig || mergeCount >= maxMerge - 1) {
                if (diags)  diags->push_back(d);
                if (hcache) hcache->push_back(htmp); // transfer ownership
                break;
            }
            // merging further: discard intermediate histogram
            delete htmp;
            ++j;
            ++mergeCount;
        }

        int next = (j < nEdges) ? j : nEdges - 1;
        result.push_back(init_edges[next]);
        i = next;
    }

    return result;
}

// ============================================================
// Diagnostic plot 1: significance vs vn bin center
//
// Blue points  = bins that passed without merging
// Red  points  = bins that required merging
// Green dashed = minSig threshold line
// ============================================================
static void plot_tail_significance(
    const std::vector<TailBinDiagnostic> &diags,
    const char *vnBranch,
    const char *cen_label,
    const char *pt_label,
    const char *outpath,
    double      minSig = 3.0)
{
    if (diags.empty()) return;

    int n = (int)diags.size();
    std::vector<double> x(n), y(n), ex(n), ey(n);
    for (int i = 0; i < n; ++i) {
        x[i]  = 0.5 * (diags[i].lo + diags[i].hi);
        ex[i] = 0.5 * (diags[i].hi - diags[i].lo);
        y[i]  = diags[i].significance;
        ey[i] = 0;
    }

    // split into normal (blue) and merged (red) point sets
    std::vector<double> xn, yn, xm, ym;
    for (int i = 0; i < n; ++i) {
        if (diags[i].merged) { xm.push_back(x[i]); ym.push_back(y[i]); }
        else                 { xn.push_back(x[i]); yn.push_back(y[i]); }
    }

    TCanvas *c = new TCanvas("c_sig_tmp", "", 900, 500);
    c->SetLeftMargin(0.12); c->SetBottomMargin(0.14);

    double ymax = *std::max_element(y.begin(), y.end());
    TGraphErrors *gr = new TGraphErrors(n, x.data(), y.data(),
                                        ex.data(), ey.data());
    gr->SetTitle(TString::Format(
        "Significance %s  %s  p_{T}=%s GeV;"
        "v_{n} bin center;Significance (#sigma)",
        vnBranch, cen_label, pt_label));
    gr->SetMarkerStyle(20); gr->SetMarkerSize(0.9);
    gr->SetMarkerColor(kBlue+1); gr->SetLineColor(kBlue+1);
    gr->Draw("AP");
    gr->GetYaxis()->SetRangeUser(0, std::max(ymax * 1.3, minSig * 2.0));
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleSize(0.05);

    // highlight merged bins in red
    if (!xm.empty()) {
        TGraph *grm = new TGraph((int)xm.size(), xm.data(), ym.data());
        grm->SetMarkerStyle(20); grm->SetMarkerSize(0.9);
        grm->SetMarkerColor(kRed+1);
        grm->Draw("P SAME");
    }

    // threshold reference line using actual data range
    double xlo = diags.front().lo;
    double xhi = diags.back().hi;
    TLine *lref = new TLine(xlo, minSig, xhi, minSig);
    lref->SetLineColor(kGreen+2); lref->SetLineWidth(2);
    lref->SetLineStyle(2);
    lref->Draw("SAME");

    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
    lat.DrawLatex(0.14, 0.92, TString::Format(
        "%s  %s  p_{T}=%s GeV  (red = merged bin)",
        vnBranch, cen_label, pt_label));
    lat.SetTextColor(kGreen+2);
    lat.DrawLatex(0.14, 0.86, TString::Format(
        "-- %.0f#sigma threshold", minSig));

    c->SaveAs(outpath);
    std::cout << "[ATB PLOT] Significance plot: " << outpath << "\n";

    delete lref;
    delete gr;
    delete c;
}

// ============================================================
// Diagnostic plot 2: mass fit per final tail bin  (CACHED version)
//
// Uses histograms cached during adaptive_tail_edges — zero extra
// tree scans.  Caller owns hcache entries; delete after this call.
// ============================================================
static void plot_tail_mass_fits(
    const std::vector<TH1D*>             &hcache,
    const std::vector<TailBinDiagnostic> &diags,
    const char                           *vnBranch,
    const char                           *cen_label,
    const char                           *pt_label,
    const char                           *outpath,
    double                                minSig = 3.0)
{
    int n = (int)diags.size();
    if (n == 0 || (int)hcache.size() < n) return;

    int ncols = std::min(n, 4);
    int nrows = (n + ncols - 1) / ncols;
    TCanvas *c = new TCanvas("c_mfit_tmp", "",
                             ncols * 300, nrows * 280);
    c->Divide(ncols, nrows, 0.001, 0.001);

    for (int i = 0; i < n; ++i) {
        c->cd(i + 1);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.16);
        gPad->SetTopMargin(0.05);

        TH1D *htmp = hcache[i];
        if (!htmp) continue;

        htmp->SetMarkerStyle(20); htmp->SetMarkerSize(0.5);
        htmp->SetLineWidth(1);
        htmp->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
        htmp->GetXaxis()->SetTitleSize(0.07);
        htmp->GetYaxis()->SetTitleSize(0.07);
        htmp->GetXaxis()->SetLabelSize(0.06);
        htmp->GetYaxis()->SetLabelSize(0.06);
        htmp->GetXaxis()->SetRangeUser(1.74, 2.0);
        htmp->SetStats(kFALSE);
        htmp->Draw("EP");

        TF1 *fdrawn = (TF1*)htmp->GetListOfFunctions()->Last();
        if (fdrawn) {
            bool passing = (diags[i].significance >= minSig);
            fdrawn->SetLineColor(passing ? kGreen+2 : kRed+1);
            fdrawn->SetLineWidth(2);
            fdrawn->Draw("SAME");
        }

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.065);
        lat.SetTextColor(kBlack);
        lat.DrawLatex(0.20, 0.88, TString::Format(
            "[%.1f, %.1f]", diags[i].lo, diags[i].hi));
        lat.DrawLatex(0.20, 0.80, TString::Format(
            "sig = %.1f", diags[i].significance));
        lat.DrawLatex(0.20, 0.72, TString::Format(
            "N = %.0f", diags[i].entries));
        if (diags[i].merged) {
            lat.SetTextColor(kRed+1);
            lat.DrawLatex(0.20, 0.64, TString::Format(
                "merged %d bins", diags[i].n_merged));
            lat.SetTextColor(kBlack);
        }
        // htmp ownership stays with caller — do NOT delete here
    }

    c->cd(1);
    TLatex tlat; tlat.SetNDC(); tlat.SetTextSize(0.07);
    tlat.DrawLatex(0.20, 0.96, TString::Format(
        "%s  %s  p_{T}=%s GeV", vnBranch, cen_label, pt_label));

    c->SaveAs(outpath);
    std::cout << "[ATB PLOT] Mass fit plots: " << outpath << "\n";
    delete c;
}

// ============================================================
// Diagnostic plot 2 (re-scan version — kept for reference)
// ============================================================
static void plot_tail_mass_fits(
    TNtuple                              *nt,
    TString                               baseCuts,
    const char                           *vnBranch,
    const std::vector<TailBinDiagnostic> &diags,
    const MassShapeParams                &shape,
    const char                           *cen_label,
    const char                           *pt_label,
    const char                           *outpath,
    double                                minSig = 3.0)
{
    int n = (int)diags.size();
    if (n == 0) return;

    int ncols = std::min(n, 4);
    int nrows = (n + ncols - 1) / ncols;
    TCanvas *c = new TCanvas("c_mfit_tmp", "",
                             ncols * 300, nrows * 280);
    c->Divide(ncols, nrows, 0.001, 0.001);

    for (int i = 0; i < n; ++i) {
        c->cd(i + 1);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.16);
        gPad->SetTopMargin(0.05);

        TH1D *htmp = nullptr;
        compute_tail_significance(nt, baseCuts, vnBranch,
                                  diags[i].lo, diags[i].hi,
                                  shape, nullptr, &htmp, false);
        if (!htmp) continue;

        htmp->SetMarkerStyle(20); htmp->SetMarkerSize(0.5);
        htmp->SetLineWidth(1);
        htmp->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
        htmp->GetXaxis()->SetTitleSize(0.07);
        htmp->GetYaxis()->SetTitleSize(0.07);
        htmp->GetXaxis()->SetLabelSize(0.06);
        htmp->GetYaxis()->SetLabelSize(0.06);
        htmp->GetXaxis()->SetRangeUser(1.74, 2.0);
        htmp->SetStats(kFALSE);
        htmp->Draw("EP");

        TF1 *fdrawn = (TF1*)htmp->GetListOfFunctions()->Last();
        if (fdrawn) {
            bool passing = (diags[i].significance >= minSig);
            fdrawn->SetLineColor(passing ? kGreen+2 : kRed+1);
            fdrawn->SetLineWidth(2);
            fdrawn->Draw("SAME");
        }

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.065);
        lat.SetTextColor(kBlack);
        lat.DrawLatex(0.20, 0.88, TString::Format(
            "[%.1f, %.1f]", diags[i].lo, diags[i].hi));
        lat.DrawLatex(0.20, 0.80, TString::Format(
            "sig = %.1f", diags[i].significance));
        lat.DrawLatex(0.20, 0.72, TString::Format(
            "N = %.0f", diags[i].entries));
        if (diags[i].merged) {
            lat.SetTextColor(kRed+1);
            lat.DrawLatex(0.20, 0.64, TString::Format(
                "merged %d bins", diags[i].n_merged));
            lat.SetTextColor(kBlack);
        }

        delete htmp; // ffit is owned by htmp's function list, deleted with it
    }

    // overall title on first pad
    c->cd(1);
    TLatex tlat; tlat.SetNDC(); tlat.SetTextSize(0.07);
    tlat.DrawLatex(0.20, 0.96, TString::Format(
        "%s  %s  p_{T}=%s GeV", vnBranch, cen_label, pt_label));

    c->SaveAs(outpath);
    std::cout << "[ATB PLOT] Mass fit plots: " << outpath << "\n";
    delete c;
}

// ============================================================
// Write adaptive edge summary to a text log file
// ============================================================
static void write_adaptive_log(
    std::ofstream                        &log,
    const char                           *vnBranch,
    const char                           *cen_label,
    const char                           *pt_label,
    const std::vector<double>            &init_edges,
    const std::vector<double>            &final_edges,
    const std::vector<TailBinDiagnostic> &diags)
{
    if (!log.is_open()) return;

    log << "=== " << vnBranch << "  " << cen_label
        << "  pT=" << pt_label << " ===\n";
    log << "  Initial bins : " << (init_edges.size()  - 1) << "\n";
    log << "  Final bins   : " << (final_edges.size() - 1) << "\n";
    log << "  Final edges  :";
    for (double e : final_edges) log << "  " << e;
    log << "\n";
    log << "  Per-bin significance:\n";
    for (auto &d : diags)
        log << TString::Format(
            "    [%7.2f, %7.2f]  N=%7.0f  sig=%6.2f%s\n",
            d.lo, d.hi, d.entries, d.significance,
            d.merged ? TString::Format("  [MERGED x%d]",
                                       d.n_merged).Data() : "").Data();
    log << "\n";
}

#endif // ADAPTIVE_TAIL_BINNING_H
