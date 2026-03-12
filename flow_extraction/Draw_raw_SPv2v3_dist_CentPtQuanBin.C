#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include <vector>
#include <iostream>

#include "Analysis_bin.h"

void Draw_raw_SPv2v3_dist_CentPtQuanBin() {

    const char* inputFile = "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Flow_output_Mar9_v0/output/Flow_out_combined.root";
    TFile *inf = TFile::Open(inputFile);
    if (!inf || inf->IsZombie()) {
        std::cerr << "Cannot open input file." << std::endl;
        return;
    }

    TDirectory *dir = (TDirectory*)inf->Get("histograms_v2v3dist");
    if (!dir) {
        std::cerr << "Cannot find directory histograms_v2v3dist." << std::endl;
        return;
    }

     //TGaxis::SetMaxDigits(1); 

    const int col_v2 = kRed+1;
    const int col_v3 = kBlue+1;

    bool firstPage = true;

    for (int i_cen = 0; i_cen < N_CENTBINS; i_cen++) {
        for (int i_pt = 0; i_pt < N_PTBINS; i_pt++) {

            // recreate canvas fresh every page
            TCanvas *c = new TCanvas("c", "c", 1600, 750);
            c->SetFillStyle(0);

            // ---- title pad: top 6% ----
            TPad *pad_title = new TPad("pad_title", "pad_title", 0.0, 0.94, 1.0, 1.0);
            pad_title->SetFillStyle(0);
            pad_title->SetBorderSize(0);
            pad_title->Draw();
            pad_title->cd();

            TLatex tex_title;
            tex_title.SetNDC(kTRUE);
            tex_title.SetTextFont(62);
            tex_title.SetTextSize(0.60);
            tex_title.SetTextAlign(22);
            tex_title.DrawLatex(0.50, 0.45, Form("Cent: %d-%d%%,  %.1f < p_{T} < %.1f GeV/c",
                                cen_edges[i_cen], cen_edges[i_cen+1],
                                pt_edges[i_pt],   pt_edges[i_pt+1]));

            // ---- main pad: bottom 94%, divided 5x2 ----
            c->cd();
            TPad *pad_main = new TPad("pad_main", "pad_main", 0.0, 0.0, 1.0, 0.94);
            pad_main->SetFillStyle(0);
            pad_main->SetBorderSize(0);
            pad_main->Draw();
            pad_main->Divide(5, 2, 0.001, 0.001);

            for (int i_q = 0; i_q < N_QBINS; i_q++) {

                TPad *pad = (TPad*)pad_main->GetPad(i_q + 1);
                if (!pad) {
                    std::cerr << "Null pad for i_q=" << i_q << std::endl;
                    continue;
                }
                pad->cd();
                pad->SetLeftMargin(0.16);
                pad->SetRightMargin(0.04);
                pad->SetTopMargin(0.10);
                pad->SetBottomMargin(0.16);
                pad->SetLogy();

                // ---- fetch histograms ----
                TString hname_v2 = Form("h_v2dist_%s_q2bin%d_%s",
                                        cen_name[i_cen], i_q, pt_name[i_pt]);
                TString hname_v3 = Form("h_v3dist_%s_q3bin%d_%s",
                                        cen_name[i_cen], i_q, pt_name[i_pt]);

                TH1D *hv2 = (TH1D*)dir->Get(hname_v2);
                TH1D *hv3 = (TH1D*)dir->Get(hname_v3);

                if (!hv2) std::cerr << "Missing: " << hname_v2 << std::endl;
                if (!hv3) std::cerr << "Missing: " << hname_v3 << std::endl;

                // ---- style v2 ----
                if (hv2) {
                    hv2->SetLineColor(col_v2);
                    hv2->SetLineWidth(1);
                    hv2->SetFillStyle(0);
                    hv2->SetStats(kFALSE);
                    hv2->SetTitle("");
                    hv2->GetXaxis()->SetTitle("v_{n}^{SP}");
                    hv2->GetXaxis()->SetTitleSize(0.07);
                    hv2->GetXaxis()->SetLabelSize(0.06);
                    hv2->GetYaxis()->SetTitle("Entries");
                    hv2->GetYaxis()->SetTitleSize(0.07);
                    hv2->GetYaxis()->SetLabelSize(0.06);
                    hv2->GetYaxis()->SetTitleOffset(0.95);
                    hv2->GetYaxis()->SetMaxDigits(1);
                }

                // ---- style v3 ----
                if (hv3) {
                    hv3->SetLineColor(col_v3);
                    hv3->SetLineWidth(1);
                    hv3->SetFillStyle(0);
                    hv3->SetStats(kFALSE);
                    hv3->SetTitle("");
                    hv3->GetXaxis()->SetTitle("v_{n}^{SP}");
                    hv3->GetXaxis()->SetTitleSize(0.07);
                    hv3->GetXaxis()->SetLabelSize(0.06);
                    hv3->GetYaxis()->SetTitle("Entries");
                    hv3->GetYaxis()->SetTitleSize(0.07);
                    hv3->GetYaxis()->SetLabelSize(0.06);
                    hv3->GetYaxis()->SetTitleOffset(0.95);
                    hv3->GetYaxis()->SetMaxDigits(1);
                }

                // ---- global y max ----
                double maxY = 0;
                if (hv2) maxY = std::max(maxY, hv2->GetMaximum());
                if (hv3) maxY = std::max(maxY, hv3->GetMaximum());

                // ---- draw ----
                bool drawn = false;
                if (hv2) {
                    hv2->GetYaxis()->SetRangeUser(0.1, 10.0 * maxY);
                    hv2->Draw("hist");
                    drawn = true;
                }
                if (hv3) {
                    hv3->GetYaxis()->SetRangeUser(0.1, 10.0 * maxY);
                    if (drawn) hv3->Draw("hist SAME");
                    else       hv3->Draw("hist");
                }

                // ---- q bin label: top-left inside each plot ----
                TLatex tex_pad;
                tex_pad.SetNDC(kTRUE);
                tex_pad.SetTextFont(62);
                tex_pad.SetTextSize(0.075);
                tex_pad.SetTextAlign(11);
                tex_pad.DrawLatex(0.20, 0.84, Form("q bin %d", i_q));

                // ---- legend: top-right inside each pad ----
                TLegend *leg = new TLegend(0.68, 0.65, 0.94, 0.82);
                leg->SetBorderSize(0);
                leg->SetTextSize(0.075);
                leg->SetTextFont(42);
                leg->SetFillStyle(0);
                if (hv2) leg->AddEntry(hv2, "v_{2}^{SP}", "l");
                if (hv3) leg->AddEntry(hv3, "v_{3}^{SP}", "l");
                leg->Draw();

            } // i_q

            c->Update();

            // ---- print to PDF ----
            bool isLast = (i_cen == N_CENTBINS-1 && i_pt == N_PTBINS-1);
            if      (firstPage && isLast) c->Print("vndist_all.pdf");
            else if (firstPage)           { c->Print("vndist_all.pdf("); firstPage = false; }
            else if (isLast)              c->Print("vndist_all.pdf)");
            else                          c->Print("vndist_all.pdf");

            delete c; // delete every page to free memory

        } // i_pt
    } // i_cen

    inf->Close();

    std::cout << "Done! Output: vndist_all.pdf  ("
              << N_CENTBINS * N_PTBINS << " pages)" << std::endl;
}

