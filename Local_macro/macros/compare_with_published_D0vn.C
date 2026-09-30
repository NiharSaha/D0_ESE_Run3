// compare_with_published.C
// Compares D0 vn results from Flow_out_combined.root (summary_pt directory)
// with published CMS results from HEPData-ins1819820-v1-Figure_2[a-f].root
//
// Layout: 2 rows x 3 columns
//   Row 1: v2 for 0-10%, 10-30%, 30-50%  (HEPData Fig 2a, 2b, 2c)
//   Row 2: v3 for 0-10%, 10-30%, 30-50%  (HEPData Fig 2d, 2e, 2f)

#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"

void compare_with_published_D0vn() {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // -----------------------------------------------------------------------
    // Open the user's file and retrieve graphs from summary_pt
    // -----------------------------------------------------------------------
    TFile *fUser = TFile::Open("Flow_out_combined.root", "READ");
    if (!fUser || fUser->IsZombie()) {
        printf("ERROR: Cannot open Flow_out_combined.root\n");
        return;
    }

    const char *centLabels[3]  = {"Cent: 0-10%", " Cent: 10-30%", " Cent: 30-50%"};
    const char *centNames[3]   = {"cent0to10", "cent10to30", "cent30to50"};
    const char *hepLetters_v2[3] = {"a", "b", "c"};
    const char *hepLetters_v3[3] = {"d", "e", "f"};

    TGraphErrors      *gUser[2][3];   // [harmonic(0=v2,1=v3)][centrality]
    TGraphAsymmErrors *gHEP[2][3];

    // Load user graphs
    fUser->cd("summary_pt");
    for (int ic = 0; ic < 3; ic++) {
        TString gname_v2 = Form("v2_graph_%s", centNames[ic]);
        TString gname_v3 = Form("v3_graph_%s", centNames[ic]);
        gUser[0][ic] = (TGraphErrors*)gDirectory->Get(gname_v2);
        gUser[1][ic] = (TGraphErrors*)gDirectory->Get(gname_v3);
        if (!gUser[0][ic]) printf("WARNING: %s not found\n", gname_v2.Data());
        if (!gUser[1][ic]) printf("WARNING: %s not found\n", gname_v3.Data());
    }

    // -----------------------------------------------------------------------
    // Load HEPData graphs
    // -----------------------------------------------------------------------
    for (int ic = 0; ic < 3; ic++) {
        // v2 files: Figure_2a, 2b, 2c
        TString fname_v2 = Form("HEPData-ins1819820-v1-Figure_2%s.root", hepLetters_v2[ic]);
        TFile *fv2 = TFile::Open(fname_v2, "READ");
        if (!fv2 || fv2->IsZombie()) {
            printf("ERROR: Cannot open %s\n", fname_v2.Data());
            return;
        }
        TString dirnam_v2 = Form("Figure 2%s", hepLetters_v2[ic]);
        fv2->cd(dirnam_v2);
        TGraphAsymmErrors *tmp_v2 = (TGraphAsymmErrors*)gDirectory->Get("Graph1D_y1");
        gHEP[0][ic] = tmp_v2 ? (TGraphAsymmErrors*)tmp_v2->Clone() : nullptr;
        fv2->Close();

        // v3 files: Figure_2d, 2e, 2f
        TString fname_v3 = Form("HEPData-ins1819820-v1-Figure_2%s.root", hepLetters_v3[ic]);
        TFile *fv3 = TFile::Open(fname_v3, "READ");
        if (!fv3 || fv3->IsZombie()) {
            printf("ERROR: Cannot open %s\n", fname_v3.Data());
            return;
        }
        TString dirnam_v3 = Form("Figure 2%s", hepLetters_v3[ic]);
        fv3->cd(dirnam_v3);
        TGraphAsymmErrors *tmp_v3 = (TGraphAsymmErrors*)gDirectory->Get("Graph1D_y1");
        gHEP[1][ic] = tmp_v3 ? (TGraphAsymmErrors*)tmp_v3->Clone() : nullptr;
        fv3->Close();
    }

    // -----------------------------------------------------------------------
    // Style settings
    // -----------------------------------------------------------------------
    // HEP published: open red circles
    // User result  : filled blue squares

    for (int ih = 0; ih < 2; ih++) {
        for (int ic = 0; ic < 3; ic++) {
            if (gHEP[ih][ic]) {
                gHEP[ih][ic]->SetMarkerStyle(24);   // open circle
                gHEP[ih][ic]->SetMarkerColor(kRed+1);
                gHEP[ih][ic]->SetLineColor(kRed+1);
                gHEP[ih][ic]->SetMarkerSize(1.0);
                gHEP[ih][ic]->SetFillColorAlpha(kRed+1, 0.25);
            }
            if (gUser[ih][ic]) {
                gUser[ih][ic]->SetMarkerStyle(21);  // filled square
                gUser[ih][ic]->SetMarkerColor(kBlue+1);
                gUser[ih][ic]->SetLineColor(kBlue+1);
                gUser[ih][ic]->SetMarkerSize(0.9);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Canvas: 2 rows x 3 columns
    // -----------------------------------------------------------------------
    TCanvas *c = new TCanvas("c_compare", "D^{0} v_{n} Comparison with Published", 1400, 850);
    c->Divide(3, 2, 0.003, 0.003);

    const char *yTitles[2] = {"v_{2}", "v_{3}"};
    // y-axis ranges tuned to data
    double yMin[2] = {-0.02, -0.05};
    double yMax[2] = { 0.25,  0.12};

    for (int ih = 0; ih < 2; ih++) {          // harmonic: 0=v2, 1=v3
        for (int ic = 0; ic < 3; ic++) {      // centrality
            int padIdx = ih * 3 + ic + 1;     // pad numbering: 1..6
            c->cd(padIdx);

            gPad->SetLeftMargin(0.14);
            gPad->SetRightMargin(0.04);
            gPad->SetBottomMargin(0.14);
            gPad->SetTopMargin(0.07);

            // Draw a frame
            TH1F *frame = gPad->DrawFrame(0.5, yMin[ih], 60.0, yMax[ih]);
            frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            frame->GetYaxis()->SetTitle(yTitles[ih]);
            frame->GetXaxis()->SetTitleSize(0.055);
            frame->GetYaxis()->SetTitleSize(0.060);
            frame->GetXaxis()->SetLabelSize(0.048);
            frame->GetYaxis()->SetLabelSize(0.048);
            frame->GetXaxis()->SetTitleOffset(1.05);
            frame->GetYaxis()->SetTitleOffset(1.10);

            // Zero line
            TLine *zline = new TLine(0.5, 0.0, 60.0, 0.0);
            zline->SetLineStyle(2);
            zline->SetLineColor(kGray+1);
            zline->Draw();

            // Draw HEPData (published) first, then user on top
            if (gHEP[ih][ic])  gHEP[ih][ic]->Draw("P2 SAME");   // P2: points + error band
            if (gHEP[ih][ic])  gHEP[ih][ic]->Draw("P SAME");
            if (gUser[ih][ic]) gUser[ih][ic]->Draw("P SAME");

            // Centrality label
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.052);
            latex.SetTextFont(42);
            latex.DrawLatex(0.18, 0.80, centLabels[ic]);

            // Harmonic label on left column only
            /*if (ic == 0) {
                latex.SetTextSize(0.056);
                latex.SetTextFont(62);
                latex.DrawLatex(0.18, 0.78, yTitles[ih]);
            }*/

            // Legend only in first pad (top-left)
            //if (ih == 0 && ic == 0) {
                TLegend *leg = new TLegend(0.33, 0.63, 0.86, 0.84);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.044);
                if (gHEP[ih][ic])  leg->AddEntry(gHEP[ih][ic],  "HIN-19-008", "pe");
                if (gUser[ih][ic]) leg->AddEntry(gUser[ih][ic], "PbPb 2023",           "pe");
                leg->Draw();
            //}
        }
    }

    c->cd();
    // Overall title
    TLatex title;
    title.SetNDC();
    title.SetTextSize(0.022);
    title.SetTextFont(42);
    title.SetTextAlign(22);
    title.DrawLatex(0.5, 0.985, "D^{0} Meson v_{n} vs p_{T} - Pb+Pb #sqrt{s_{NN}} = 5.02 TeV");

    c->Update();
    c->SaveAs("compare_with_published.pdf");
    c->SaveAs("compare_with_published.png");

    printf("\nDone! Saved: compare_with_published.pdf and compare_with_published.png\n");
}
