// compare_chg_v2.C
//
// Compares charged-particle v2{SP} vs pT between:
//   - Published CMS data (HEPData-ins1511868-v1-Table_*.root)
//   - Our new results  (flow_Analysis_chg_out_combined.root)
//
// Centrality bins: 0-5%, 5-10%, 10-20%, 20-30%, 30-40%, 40-50%
//
// Usage:
//   root -l -q compare_charge_vn.C

#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TPad.h"
#include "TLine.h"
#include "TROOT.h"
#include "TH1F.h"

void setStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetTitleFont(42,"xyz");
}

// Format a TProfile for our new results
void styleNewProfile(TProfile *p, int color, int marker)
{
    p->SetMarkerStyle(marker);
    p->SetMarkerSize(1.1);
    p->SetMarkerColor(color);
    p->SetLineColor(color);
    p->SetLineWidth(2);
}

// Format a TGraphAsymmErrors for published (HEPData) results
void styleRefGraph(TGraphAsymmErrors *g, int color, int marker)
{
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.1);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
    // no fill — markers only for published data
}

// ─────────────────────────────────────────────────────────────────────────────
void compare_with_published_chargeVn()
{
    setStyle();
    gROOT->ForceStyle();

    // ── centrality bin definitions ───────────────────────────────────────────
    const int NCENT = 6;
    const char* centLabel[NCENT] = {
        "0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%"
    };
    // HEPData table index (1-based) for each centrality bin
    // Table 1 → 0-5% (smallest v2), ..., Table 6 → 40-50%
    const int hepTable[NCENT] = {1, 2, 3, 4, 5, 6};
    // TProfile names in our output file
    const char* ourProfName[NCENT] = {
        "hp_v2_incl_cent0to5",
        "hp_v2_incl_cent5to10",
        "hp_v2_incl_cent10to20",
        "hp_v2_incl_cent20to30",
        "hp_v2_incl_cent30to40",
        "hp_v2_incl_cent40to50"
    };

    // Colors and markers
    const int colNew = kRed+1;
    const int colRef = kBlue+1;

    // ── open files ───────────────────────────────────────────────────────────
    TFile *fNew = TFile::Open("flow_Analysis_chg_out_combined.root");
    if (!fNew || fNew->IsZombie()) {
        Printf("ERROR: cannot open flow_Analysis_chg_out_combined.root");
        return;
    }
    TDirectoryFile *dNew = (TDirectoryFile*)fNew->Get("Inclusive_TProfile");
    if (!dNew) {
        Printf("ERROR: directory Inclusive_TProfile not found");
        return;
    }

    // ── create output PDF ────────────────────────────────────────────────────
    TCanvas *cAll = new TCanvas("cAll","chg v2 comparison",1800,900);
    cAll->Divide(3,2,0.002,0.002);

    for (int ic = 0; ic < NCENT; ic++) {

        // --- get our new result ---
        TProfile *pNew = (TProfile*)dNew->Get(ourProfName[ic]);
        if (!pNew) {
            Printf("WARNING: %s not found in new file", ourProfName[ic]);
            continue;
        }
        pNew = (TProfile*)pNew->Clone(Form("pNew_cent%d",ic));
        styleNewProfile(pNew, colNew, 20);

        // --- get HEPData reference ---
        TString refName = Form("HEPData-ins1511868-v1-Table_%d.root",
                               hepTable[ic]);
        TFile *fRef = TFile::Open(refName);
        if (!fRef || fRef->IsZombie()) {
            Printf("WARNING: cannot open %s", refName.Data());
            continue;
        }
        TString dirName = Form("Table %d", hepTable[ic]);
        TDirectoryFile *dRef = (TDirectoryFile*)fRef->Get(dirName);
        TGraphAsymmErrors *gRef = (TGraphAsymmErrors*)dRef->Get("Graph1D_y1");
        if (!gRef) {
            Printf("WARNING: Graph1D_y1 not found in %s", refName.Data());
            fRef->Close();
            continue;
        }
        // Clone so we can close the file safely
        gRef = (TGraphAsymmErrors*)gRef->Clone(Form("gRef_cent%d",ic));
        fRef->Close();
        styleRefGraph(gRef, colRef, 24);

        // --- pad ---
        TVirtualPad *pad = cAll->cd(ic + 1);
        pad->SetLeftMargin(0.15);
        pad->SetBottomMargin(0.14);
        pad->SetTopMargin(0.07);
        pad->SetRightMargin(0.04);

        // dummy frame to control axes
        double xMax = 20.0;
        double yMin =  0.0;
        double yMax =  0.35;
        TH1F *hFrame = new TH1F(Form("hf%d",ic),"",100,0.,xMax);
        hFrame->SetMinimum(yMin);
        hFrame->SetMaximum(yMax);
        hFrame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hFrame->GetYaxis()->SetTitle("v_{2}{SP}");
        hFrame->GetXaxis()->SetTitleSize(0.055);
        hFrame->GetYaxis()->SetTitleSize(0.055);
        hFrame->GetXaxis()->SetLabelSize(0.048);
        hFrame->GetYaxis()->SetLabelSize(0.048);
        hFrame->GetXaxis()->SetTitleOffset(1.05);
        hFrame->GetYaxis()->SetTitleOffset(1.20);
        hFrame->GetXaxis()->SetNdivisions(505);
        hFrame->Draw("AXIS");

        // draw reference markers only (no shaded band)
        gRef->Draw("P same");

        // dotted line at y=0
        TLine *zline = new TLine(0., 0., xMax, 0.);
        zline->SetLineStyle(3);
        zline->SetLineWidth(1);
        zline->SetLineColor(kBlack);
        zline->Draw("same");

        // draw our new result
        // Note: statistical errors on TProfile are very small (~1e-4 to 1e-3)
        // due to large track statistics; E1 draws them correctly but they are
        // sub-pixel at this scale.
        pNew->Draw("E1 X0 same");

        // --- centrality label (upper-left) ---
        TLatex lat;
        lat.SetNDC();
        lat.SetTextFont(42);
        lat.SetTextSize(0.055);
        lat.DrawLatex(0.18, 0.85, Form("Cent: %s", centLabel[ic]));

        // --- data legend (right side) ---
        TLegend *leg = new TLegend(0.45, 0.72, 0.94, 0.92);
        leg->SetTextSize(0.040);
        leg->AddEntry(gRef, "HIN-15-014 (5.02 TeV)", "P");
        leg->AddEntry(pNew, "CMS PbPb 5.36 TeV", "P");
        leg->Draw();
    }

    cAll->SaveAs("compare_chg_v2.pdf");
    Printf("Saved compare_chg_v2.pdf");

    // =========================================================================
    // v3 comparison
    // =========================================================================
    // HEPData table indices for v3 (Tables 8-13: 0-5% ... 40-50%)
    const int hepTableV3[NCENT] = {8, 9, 10, 11, 12, 13};
    // TProfile names for v3 in our output file
    const char* ourProfNameV3[NCENT] = {
        "hp_v3_incl_cent0to5",
        "hp_v3_incl_cent5to10",
        "hp_v3_incl_cent10to20",
        "hp_v3_incl_cent20to30",
        "hp_v3_incl_cent30to40",
        "hp_v3_incl_cent40to50"
    };

    TCanvas *cAllV3 = new TCanvas("cAllV3","chg v3 comparison",1800,900);
    cAllV3->Divide(3,2,0.002,0.002);

    for (int ic = 0; ic < NCENT; ic++) {

        // --- get our new v3 result ---
        TProfile *pNewV3 = (TProfile*)dNew->Get(ourProfNameV3[ic]);
        if (!pNewV3) {
            Printf("WARNING: %s not found in new file", ourProfNameV3[ic]);
            continue;
        }
        pNewV3 = (TProfile*)pNewV3->Clone(Form("pNewV3_cent%d",ic));
        styleNewProfile(pNewV3, colNew, 20);

        // --- get HEPData v3 reference ---
        TString refNameV3 = Form("HEPData-ins1511868-v1-Table_%d.root",
                                 hepTableV3[ic]);
        TFile *fRefV3 = TFile::Open(refNameV3);
        if (!fRefV3 || fRefV3->IsZombie()) {
            Printf("WARNING: cannot open %s", refNameV3.Data());
            continue;
        }
        TString dirNameV3 = Form("Table %d", hepTableV3[ic]);
        TDirectoryFile *dRefV3 = (TDirectoryFile*)fRefV3->Get(dirNameV3);
        TGraphAsymmErrors *gRefV3 = (TGraphAsymmErrors*)dRefV3->Get("Graph1D_y1");
        if (!gRefV3) {
            Printf("WARNING: Graph1D_y1 not found in %s", refNameV3.Data());
            fRefV3->Close();
            continue;
        }
        gRefV3 = (TGraphAsymmErrors*)gRefV3->Clone(Form("gRefV3_cent%d",ic));
        fRefV3->Close();
        styleRefGraph(gRefV3, colRef, 24);

        // --- pad ---
        TVirtualPad *padV3 = cAllV3->cd(ic + 1);
        padV3->SetLeftMargin(0.15);
        padV3->SetBottomMargin(0.14);
        padV3->SetTopMargin(0.07);
        padV3->SetRightMargin(0.04);

        // dummy frame
        double xMaxV3 = 20.0;
        double yMinV3 =  -0.02;
        double yMaxV3 =  0.15;
        TH1F *hFrameV3 = new TH1F(Form("hfv3_%d",ic),"",100,0.,xMaxV3);
        hFrameV3->SetMinimum(yMinV3);
        hFrameV3->SetMaximum(yMaxV3);
        hFrameV3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hFrameV3->GetYaxis()->SetTitle("v_{3}{SP}");
        hFrameV3->GetXaxis()->SetTitleSize(0.055);
        hFrameV3->GetYaxis()->SetTitleSize(0.055);
        hFrameV3->GetXaxis()->SetLabelSize(0.048);
        hFrameV3->GetYaxis()->SetLabelSize(0.048);
        hFrameV3->GetXaxis()->SetTitleOffset(1.05);
        hFrameV3->GetYaxis()->SetTitleOffset(1.20);
        hFrameV3->GetXaxis()->SetNdivisions(505);
        hFrameV3->Draw("AXIS");

        // draw reference markers only
        gRefV3->Draw("P same");

        // dotted line at y=0
        TLine *zlineV3 = new TLine(0., 0., xMaxV3, 0.);
        zlineV3->SetLineStyle(3);
        zlineV3->SetLineWidth(1);
        zlineV3->SetLineColor(kBlack);
        zlineV3->Draw("same");

        // draw our new result
        pNewV3->Draw("E1 X0 same");

        // --- centrality label (upper-left) ---
        TLatex latV3;
        latV3.SetNDC();
        latV3.SetTextFont(42);
        latV3.SetTextSize(0.055);
        latV3.DrawLatex(0.18, 0.85, Form("Cent: %s", centLabel[ic]));

        // --- data legend (right side) ---
        TLegend *legV3 = new TLegend(0.45, 0.72, 0.94, 0.92);
        legV3->SetTextSize(0.040);
        legV3->AddEntry(gRefV3, "HIN-15-014 (5.02 TeV)", "P");
        legV3->AddEntry(pNewV3, "CMS PbPb 5.36 TeV", "P");
        legV3->Draw();
    }

    cAllV3->SaveAs("compare_chg_v3.pdf");
    Printf("Saved compare_chg_v3.pdf");

    fNew->Close();
    Printf("Done.");
}
