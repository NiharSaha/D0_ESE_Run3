// PlotMassDistributions.C
//
// Plots all mass distributions stored in the "mass_fitted" directory of
// D0_Flow_12NUQbin_diffq2q3_out_combined_Jun26.root into multi-page PDFs.
//
// DATA:
//   Histogram naming:  h_mass_cent<CENT>_q2bin<Q>_pT<PT>
//   5 centralities x 9 pT bins x 12 q2 bins = 540 histograms
//   -> 45 PDF pages, 12 histos/page (one page per cent-pT combination,
//      4x3 grid of the 12 q2 bins)
//   Each histogram already carries, in its fFunctions list:
//     - several TF1 fit components (total fit, signal, swap, bkg, KK refl.)
//     - a TLegend (unused/zero-size)
//     - three TLatex annotations: "Yield = ...", "Y/DeltaY = ..." (significance),
//       "chi2/ndf = ..." positioned at NDC (0.14, 0.80/0.75/0.70)
//   We add a new "cent / pT / q2bin" label above these, and push the
//   Yield/significance/chi2 annotations down so nothing overlaps.
//
// MC:
//   Histogram naming:  h_mc_sig_v2_<centIdx>_<ptIdx>_<ptIdx>_cent<CENT>
//   5 centralities x 9 pT bins = 45 histograms
//   -> 5 PDF pages, 9 histos/page (3x3 grid per centrality)
//   These are saved in plain "hist" style with a red total-fit TF1;
//   we redraw as black data points with a blue fit line, no stats/fit box,
//   and add a "cent / pT" label in the top-left corner.
//
// Usage:
//   root -l -b -q 'PlotMassDistributions.C("D0_Flow_12NUQbin_diffq2q3_out_combined_Jun26.root")'

#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TF1.h>
#include <TLatex.h>
#include <TList.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

// Split a "<lo>to<hi>" string (e.g. "10to20") into its lo/hi pieces.
static void ParseRange(const TString& s, TString& lo, TString& hi)
{
    Int_t idx = s.Index("to");
    if (idx < 0) { lo = s; hi = ""; return; }
    lo = s(0, idx);
    hi = s(idx + 2, s.Length() - idx - 2);
}

// Draw a small "cent / pT [/ q2bin]" label in the top-left of the current pad.
// Returns the number of lines drawn (so callers can decide how much room to reserve).
static void DrawBinLabel(const TString& centRange, const TString& ptRange, Int_t q2bin = -1)
{
    TString centLo, centHi, ptLo, ptHi;
    ParseRange(centRange, centLo, centHi);
    ParseRange(ptRange, ptLo, ptHi);

    TLatex* line1 = new TLatex(0.15, 0.88, Form("cent: %s-%s%%", centLo.Data(), centHi.Data()));
    line1->SetNDC();
    line1->SetTextFont(42);
    line1->SetTextSize(0.055);
    line1->Draw();

    TString line2text = (q2bin >= 0)
        ? Form("%s<p_{T}<%s GeV, q2bin%d", ptLo.Data(), ptHi.Data(), q2bin)
        : Form("%s<p_{T}<%s GeV", ptLo.Data(), ptHi.Data());

    TLatex* line2 = new TLatex(0.15, 0.80, line2text);
    line2->SetNDC();
    line2->SetTextFont(42);
    line2->SetTextSize(0.055);
    line2->Draw();
}

void PlotMassDistributions(const char* infile = "D0_Flow_12NUQbin_diffq2q3_out_combined_Jun26.root",
                            const char* outdir = ".")
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);   // neither data nor MC needs the automatic stats/fit box
    gStyle->SetOptTitle(1);

    TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: could not open file " << infile << std::endl;
        return;
    }

    TDirectory* dir = (TDirectory*)f->Get("mass_fitted");
    if (!dir) {
        std::cerr << "ERROR: could not find directory 'mass_fitted' in " << infile << std::endl;
        return;
    }

    // ---- Bin definitions (must match the ROOT file's naming convention) ----
    const std::vector<TString> centBins = {"0to10", "10to20", "20to30", "30to40", "40to50"};
    const std::vector<TString> ptBins   = {"2to3", "3to4", "4to5", "5to6", "6to8",
                                            "8to10", "10to15", "15to30", "30to100"};
    const int nQ2Bins = 12;

    // =====================================================================
    // 1) DATA: 45 pages (one per cent x pT), 12 q2-bin mass plots per page
    // =====================================================================
    {
        TString outPdfData = TString::Format("%s/DataMass_AllQ2Bins.pdf", outdir);
        TCanvas* cData = new TCanvas("cData", "Data Mass Distributions", 1600, 1200);

        bool firstPage = true;
        int pageCount = 0;

        for (const auto& cent : centBins) {
            for (const auto& pt : ptBins) {
                cData->Clear();
                cData->Divide(4, 3); // 12 pads: 4 columns x 3 rows

                for (int q = 0; q < nQ2Bins; ++q) {
                    TString hname = TString::Format("h_mass_cent%s_q2bin%d_pT%s",
                                                      cent.Data(), q, pt.Data());
                    TH1* h = (TH1*)dir->Get(hname);

                    cData->cd(q + 1);
                    if (!h) {
                        std::cerr << "WARNING: missing histogram " << hname << std::endl;
                        TPaveText* pt_missing = new TPaveText(0.1, 0.4, 0.9, 0.6, "NDC");
                        pt_missing->AddText(Form("MISSING: %s", hname.Data()));
                        pt_missing->Draw();
                        continue;
                    }

                    h->SetStats(0);
                    h->SetTitle("");
                    h->GetXaxis()->SetTitle("Mass (GeV)");
                    h->GetYaxis()->SetTitle("Entries/5 MeV");
                    h->GetXaxis()->SetTitleSize(0.055);
                    h->GetYaxis()->SetTitleSize(0.055);
                    h->GetXaxis()->SetLabelSize(0.05);
                    h->GetYaxis()->SetLabelSize(0.05);
                    h->GetYaxis()->SetTitleOffset(1.3);

                    // Push the existing Yield / significance / chi2 TLatex
                    // annotations (saved inside the histogram's fFunctions)
                    // down, to make room for the new cent/pT/q2bin label.
                    const double shiftY = 0.20;
                    TIter nextFunc(h->GetListOfFunctions());
                    TObject* fobj;
                    while ((fobj = nextFunc())) {
                        if (fobj->InheritsFrom("TLatex")) {
                            TLatex* lat = (TLatex*)fobj;
                            lat->SetY(lat->GetY() - shiftY);
                        }
                    }

                    h->Draw("PE"); // draws histogram + all attached fit curves/annotations

                    DrawBinLabel(cent, pt, q);
                }

                TString pageOpt = firstPage ? (outPdfData + "(") : outPdfData;
                cData->Print(pageOpt, "pdf");
                firstPage = false;
                pageCount++;
            }
        }
        // Close the multi-page PDF
        cData->Print(outPdfData + ")", "pdf");
        std::cout << "Wrote " << pageCount << " pages -> " << outPdfData << std::endl;
        delete cData;
    }

    // =====================================================================
    // 2) MC: 5 pages (one per centrality), 9 pT-bin h_mc_sig plots per page
    // =====================================================================
    {
        TString outPdfMC = TString::Format("%s/MCMass_Sig_AllPtBins.pdf", outdir);
        TCanvas* cMC = new TCanvas("cMC", "MC Signal Mass Distributions", 1600, 1200);

        bool firstPage = true;
        int pageCount = 0;

        for (size_t ic = 0; ic < centBins.size(); ++ic) {
            cMC->Clear();
            cMC->Divide(3, 3); // 9 pads: 3 columns x 3 rows

            for (size_t ip = 0; ip < ptBins.size(); ++ip) {
                TString hname = TString::Format("h_mc_sig_v2_%d_%d_%d_cent%s",
                                                  (int)ic, (int)ip, (int)ip, centBins[ic].Data());
                TH1* h = (TH1*)dir->Get(hname);

                cMC->cd(ip + 1);
                if (!h) {
                    std::cerr << "WARNING: missing histogram " << hname << std::endl;
                    TPaveText* pt_missing = new TPaveText(0.1, 0.4, 0.9, 0.6, "NDC");
                    pt_missing->AddText(Form("MISSING: %s", hname.Data()));
                    pt_missing->Draw();
                    continue;
                }

                h->SetStats(0);
                h->SetTitle("");

                // Black data points instead of the saved "hist" style.
                h->SetMarkerStyle(20);
                h->SetMarkerColor(kBlack);
                h->SetMarkerSize(0.7);
                h->SetLineColor(kBlack);

                h->GetXaxis()->SetTitle("Mass (GeV)");
                h->GetYaxis()->SetTitle("Entries");
                h->GetXaxis()->SetTitleSize(0.055);
                h->GetYaxis()->SetTitleSize(0.055);
                h->GetXaxis()->SetLabelSize(0.05);
                h->GetYaxis()->SetLabelSize(0.05);
                h->GetYaxis()->SetTitleOffset(1.3);

                // Fit curve in blue.
                TIter nextFunc(h->GetListOfFunctions());
                TObject* fobj;
                while ((fobj = nextFunc())) {
                    if (fobj->InheritsFrom("TF1")) {
                        TF1* fit = (TF1*)fobj;
                        fit->SetLineColor(kBlue);
                        fit->SetLineWidth(2);
                    }
                }

                h->Draw("PE"); // black points + blue fit curve, no stats/fit box

                DrawBinLabel(centBins[ic], ptBins[ip]);
            }

            TString pageOpt = firstPage ? (outPdfMC + "(") : outPdfMC;
            cMC->Print(pageOpt, "pdf");
            firstPage = false;
            pageCount++;
        }
        cMC->Print(outPdfMC + ")", "pdf");
        std::cout << "Wrote " << pageCount << " pages -> " << outPdfMC << std::endl;
        delete cMC;
    }

    f->Close();
}