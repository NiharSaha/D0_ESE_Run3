// plot_vn_vs_qbin_D0Chg_combined.C
//
// Generates both Q-Bin Index and Mean Q plots sequentially.
// Includes boolean toggles to select the x-axis error type:
//   - showMeanErr: plots the statistical error on the mean (RMS/sqrt(N))
//   - showStdDev: plots the physical standard deviation (width) of the q-bin

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <numeric>

#include "TFile.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TPad.h"

static double PearsonR(const std::vector<double>& x, const std::vector<double>& y)
{
    int n = (int)x.size();
    if (n < 2) return 0.0;
    double mx = 0, my = 0;
    for (int i = 0; i < n; ++i) { mx += x[i]; my += y[i]; }
    mx /= n; my /= n;
    double num = 0, dx2 = 0, dy2 = 0;
    for (int i = 0; i < n; ++i) {
        double dx = x[i] - mx, dy = y[i] - my;
        num += dx * dy; dx2 += dx * dx; dy2 += dy * dy;
    }
    double denom = std::sqrt(dx2 * dy2);
    return (denom > 0) ? num / denom : 0.0;
}

void plot_vn_vs_qbin_D0Chg_combined(
				    const char* d0_file     = "D0_Flow_12NUQbin_diffq2q3_out_combined_Aug24.root",
				    const char* chg_file    = "Charge_flow_Analysis_12NUQbin_diffq2q3_out_combined_Aug24.root",
				    const char* tag         = "Sept14_v0",
				    bool        isRun3      = true,
				    bool        showMeanErr = true,
				    bool        showStdDev  = false) // NEW: StdDev boolean added
{
    const int N_QBINS = 12;
    const int N_PTBINS_V2 = 9;
    const double pt_edges_v2[N_PTBINS_V2 + 1] = { 2, 3, 4, 5, 6, 8, 10, 15, 30, 100 };
    const int N_PTBINS_V3 = 6;
    const double pt_edges_v3[N_PTBINS_V3 + 1] = { 2, 4, 6, 8, 10, 20, 50};

    const int N_CENTBINS_CHG = 5;
    const int    min_cent_chg[N_CENTBINS_CHG] = {  0, 10, 20, 30, 40 };
    const int    max_cent_chg[N_CENTBINS_CHG] = { 10, 20, 30, 40, 50 };
    const char* cent_lbl_chg[N_CENTBINS_CHG] = {
        "cent0to10","cent10to20","cent20to30","cent30to40","cent40to50"
    };

    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    // =========================================================================
    // Two-Pass System: Pass 0 = Q-Bin Index | Pass 1 = Mean Q
    // =========================================================================
    for (int pass = 0; pass < 2; ++pass)
    {
        bool doMeanQ = (pass == 1);
        TString typeTag = doMeanQ ? "meanq" : "qbin";

        std::cout << "\n====================================================================\n";
        std::cout << ">>> STARTING PASS: " << (doMeanQ ? "Mean Q Value Plots" : "Q-Bin Index Plots") << "\n";
        std::cout << "====================================================================\n";

        TFile* fD0  = TFile::Open(d0_file);
        TFile* fChg = TFile::Open(chg_file);
        if (!fD0  || fD0->IsZombie())  { std::cerr << "ERROR: Cannot open D0 file.\n"; return; }
        if (!fChg || fChg->IsZombie()) { std::cerr << "ERROR: Cannot open chg file.\n"; return; }

        TString d0_dir_name  = "vn_vs_qbin";
        TString chg_dir_name = "vsQbin_TProfile"; 
        
        TString d0_qdist_dir  = "q_distributions"; 
        TString chg_qdist_dir = "Qdist_perBin"; 

        auto* dir_d0_qbin   = (TDirectory*)fD0->Get(d0_dir_name);
        auto* dir_chg_vsq   = (TDirectory*)fChg->Get(chg_dir_name);
        auto* dir_d0_qdist  = (TDirectory*)fD0->Get(d0_qdist_dir);
        auto* dir_chg_qdist = (TDirectory*)fChg->Get(chg_qdist_dir);

        if (!dir_d0_qbin) std::cerr << "WARNING: D0 " << d0_dir_name << " dir not found\n";
        if (!dir_chg_vsq) std::cerr << "WARNING: CHG " << chg_dir_name << " dir not found\n";
        if (doMeanQ && !dir_d0_qdist)  std::cerr << "WARNING: D0 " << d0_qdist_dir << " dir not found\n";
        if (doMeanQ && !dir_chg_qdist) std::cerr << "WARNING: CHG " << chg_qdist_dir << " dir not found\n";

        TString plot_base_dir = Form("vn_vs_%s_plots", typeTag.Data());
        gSystem->mkdir(plot_base_dir, kTRUE);
        gSystem->mkdir(plot_base_dir + "/v2", kTRUE);
        gSystem->mkdir(plot_base_dir + "/v3", kTRUE);

        TString out_file = Form("ESE_vn_vs_%s_%s.root", typeTag.Data(), tag);
        TFile* fOut   = new TFile(out_file.Data(), "RECREATE");
        auto* dir_v2 = fOut->mkdir("vn_vs_q_v2");
        auto* dir_v3 = fOut->mkdir("vn_vs_q_v3");

        TLatex ltx;
        ltx.SetNDC();
        ltx.SetTextFont(42);

        const char* chg_pt_tag = isRun3 ? "pt0p5to3" : "pt1to3";
        TString chg_prof_base2 = "hp_v2_vsq2";
        TString chg_prof_base3 = "hp_v3_vsq3";

        TProfile* hp_chg_v2[N_CENTBINS_CHG] = {};
        TProfile* hp_chg_v3[N_CENTBINS_CHG] = {};
        
        TH1D* hq2_dist_d0[N_CENTBINS_CHG][N_QBINS]  = {};
        TH1D* hq3_dist_d0[N_CENTBINS_CHG][N_QBINS]  = {};
        TH1D* hq2_dist_chg[N_CENTBINS_CHG][N_QBINS] = {};
        TH1D* hq3_dist_chg[N_CENTBINS_CHG][N_QBINS] = {};

        for (int ic = 0; ic < N_CENTBINS_CHG; ++ic) {
            hp_chg_v2[ic] = dir_chg_vsq ? (TProfile*)dir_chg_vsq->Get(Form("%s_%s_%s", chg_prof_base2.Data(), chg_pt_tag, cent_lbl_chg[ic])) : nullptr;
            hp_chg_v3[ic] = dir_chg_vsq ? (TProfile*)dir_chg_vsq->Get(Form("%s_%s_%s", chg_prof_base3.Data(), chg_pt_tag, cent_lbl_chg[ic])) : nullptr;

            if (doMeanQ) {
                for (int iq = 0; iq < N_QBINS; ++iq) {
                    if (dir_d0_qdist) {
                        hq2_dist_d0[ic][iq] = (TH1D*)dir_d0_qdist->Get(Form("hq2_dist_%s_q2bin%d", cent_lbl_chg[ic], iq));
                        hq3_dist_d0[ic][iq] = (TH1D*)dir_d0_qdist->Get(Form("hq3_dist_%s_q3bin%d", cent_lbl_chg[ic], iq));
                    }
                    if (dir_chg_qdist) {
                        hq2_dist_chg[ic][iq] = (TH1D*)dir_chg_qdist->Get(Form("hq2_dist_cen%d_q2bin%d", ic, iq));
                        hq3_dist_chg[ic][iq] = (TH1D*)dir_chg_qdist->Get(Form("hq3_dist_cen%d_q3bin%d", ic, iq));
                    }
                }
            }
        }

        double cent_d0_ylo[2][N_CENTBINS_CHG],  cent_d0_yhi[2][N_CENTBINS_CHG];
        double cent_chg_ylo[2][N_CENTBINS_CHG], cent_chg_yhi[2][N_CENTBINS_CHG];

        for (int ivn = 2; ivn <= 3; ++ivn) {
            int vi             = ivn - 2;
            int n_pt_scan      = (ivn == 2) ? N_PTBINS_V2 : N_PTBINS_V3;
            const double* pt_e = (ivn == 2) ? pt_edges_v2 : pt_edges_v3;
            TProfile** hp_scan = (ivn == 2) ? hp_chg_v2   : hp_chg_v3;
            TString d0_gr_base = Form("v%d_vs_q%dbin", ivn, ivn);

            for (int ic = 0; ic < N_CENTBINS_CHG; ++ic) {
                double gd0_min =  1e9, gd0_max = -1e9;
                double gch_min =  1e9, gch_max = -1e9;
                
                for (int iq = 0; iq < N_QBINS; ++iq) {
                    if (hp_scan[ic] && hp_scan[ic]->GetBinEntries(iq+1) > 0) {
                        double cv = hp_scan[ic]->GetBinContent(iq+1);
                        double ce = hp_scan[ic]->GetBinError(iq+1);
                        if (ce <= 0) continue;
                        gch_min = std::min(gch_min, cv - ce);
                        gch_max = std::max(gch_max, cv + ce);
                    }
                }
                
                for (int ip = 0; ip < n_pt_scan; ++ip) {
                    if (ivn==2) {if (pt_e[ip] >= 30) continue;}
                    if (ivn==3) {if ( pt_e[ip] >= 20) continue;}
                    TString ptTag  = Form("pT%dto%d", (int)pt_e[ip], (int)pt_e[ip+1]);
                    TString grName = Form("%s_%s_%s", d0_gr_base.Data(), cent_lbl_chg[ic], ptTag.Data());
                    TGraphErrors* gr = dir_d0_qbin ? (TGraphErrors*)dir_d0_qbin->Get(grName) : nullptr;
                    
                    if (!gr || gr->GetN() == 0) continue;
                    for (int iq = 0; iq < N_QBINS && iq < gr->GetN(); ++iq) {
                        double dv = gr->GetPointY(iq), de = gr->GetErrorY(iq);
                        double ce = (hp_scan[ic] && hp_scan[ic]->GetBinEntries(iq+1) > 0) ? hp_scan[ic]->GetBinError(iq+1) : 0.0;
                        if (de <= 0 || ce <= 0) continue;
                        gd0_min = std::min(gd0_min, dv - de);
                        gd0_max = std::max(gd0_max, dv + de);
                    }
                }
                if (gd0_min > gd0_max) { gd0_min = 0.0; gd0_max = 0.3; }
                if (gch_min > gch_max) { gch_min = 0.0; gch_max = 0.3; }
                double d0_marg  = 0.25 * std::max(gd0_max - gd0_min, 1e-6);
                double chg_marg = 0.25 * std::max(gch_max - gch_min, 1e-6);
                cent_d0_ylo[vi][ic]  = gd0_min  - d0_marg;
                cent_d0_yhi[vi][ic]  = gd0_max  + d0_marg;
                cent_chg_ylo[vi][ic] = std::max(0.0, gch_min - chg_marg);
                cent_chg_yhi[vi][ic] = gch_max  + chg_marg;
                
                if (ivn == 2) {
                    double unified_lo = std::min(cent_d0_ylo[vi][ic], cent_chg_ylo[vi][ic]);
                    double unified_hi = std::max(cent_d0_yhi[vi][ic], cent_chg_yhi[vi][ic]);
                    cent_d0_ylo[vi][ic]  = unified_lo;
                    cent_d0_yhi[vi][ic]  = unified_hi;
                    cent_chg_ylo[vi][ic] = unified_lo;
                    cent_chg_yhi[vi][ic] = unified_hi;
                }
            }
        }

        // --- Plot Generation ---
        for (int ivn = 2; ivn <= 3; ++ivn)
        {
            int           n_pt       = (ivn == 2) ? N_PTBINS_V2  : N_PTBINS_V3;
            const double* pt_edges   = (ivn == 2) ? pt_edges_v2  : pt_edges_v3;
            TProfile** hp_chg        = (ivn == 2) ? hp_chg_v2    : hp_chg_v3;
            TDirectory* dir_out      = (ivn == 2) ? dir_v2       : dir_v3;
            TString d0_gr_base       = Form("v%d_vs_q%dbin", ivn, ivn);

            for (int ic = 0; ic < N_CENTBINS_CHG; ++ic)
            {
                double chg_vn_q[N_QBINS] = {}, chg_vn_q_e[N_QBINS] = {};
                for (int iq = 0; iq < N_QBINS; ++iq) {
                    if (hp_chg[ic] && hp_chg[ic]->GetBinEntries(iq + 1) > 0) {
                        chg_vn_q[iq]   = hp_chg[ic]->GetBinContent(iq + 1);
                        chg_vn_q_e[iq] = hp_chg[ic]->GetBinError(iq + 1);
                    }
                }

                for (int ip = 0; ip < n_pt; ++ip)
                {
                    const double pt_lo = pt_edges[ip];
                    const double pt_hi = pt_edges[ip + 1];
                    TString ptTag  = Form("pT%dto%d", (int)pt_lo, (int)pt_hi);
                    TString grName = Form("%s_%s_%s", d0_gr_base.Data(), cent_lbl_chg[ic], ptTag.Data());

                    TGraphErrors* gr_d0 = dir_d0_qbin ? (TGraphErrors*)dir_d0_qbin->Get(grName) : nullptr;
                    if (!gr_d0 || gr_d0->GetN() == 0) continue;

                    std::vector<double> d0_v, chg_v;
                    double d0_q[N_QBINS]  = {}, d0_q_e[N_QBINS]  = {};
                    double chg_q[N_QBINS] = {}, chg_q_e[N_QBINS] = {};
                    
                    double qx_d0[N_QBINS] = {}, qex_d0[N_QBINS] = {};
                    double qx_chg[N_QBINS] = {}, qex_chg[N_QBINS] = {};
                    int npts = 0;

                    for (int iq = 0; iq < N_QBINS && iq < gr_d0->GetN(); ++iq) {
                        double dv  = gr_d0->GetPointY(iq);
                        double de  = gr_d0->GetErrorY(iq);
                        double cv  = chg_vn_q[iq];
                        double ce  = chg_vn_q_e[iq];
                        if (de <= 0 || ce <= 0) continue;

                        double current_qx_d0, current_qex_d0;
                        double current_qx_chg, current_qex_chg;

                        if (doMeanQ) {
                            // Extract D0 mean and error
                            TH1D* hdist_d0 = (ivn == 2) ? hq2_dist_d0[ic][iq] : hq3_dist_d0[ic][iq];
                            if (hdist_d0 && hdist_d0->GetEntries() > 0) {
                                current_qx_d0  = hdist_d0->GetMean();
                                // NEW: Priority logic for x-axis error selection
                                if (showStdDev)       current_qex_d0 = hdist_d0->GetStdDev();
                                else if (showMeanErr) current_qex_d0 = hdist_d0->GetMeanError();
                                else                  current_qex_d0 = 0.0;
                            } else {
                                current_qx_d0  = gr_d0->GetPointX(iq);
                                current_qex_d0 = 0.0;
                            }
                            
                            // Extract Charged Particle mean and error
                            TH1D* hdist_chg = (ivn == 2) ? hq2_dist_chg[ic][iq] : hq3_dist_chg[ic][iq];
                            if (hdist_chg && hdist_chg->GetEntries() > 0) {
                                current_qx_chg  = hdist_chg->GetMean();
                                // NEW: Priority logic for x-axis error selection
                                if (showStdDev)       current_qex_chg = hdist_chg->GetStdDev();
                                else if (showMeanErr) current_qex_chg = hdist_chg->GetMeanError();
                                else                  current_qex_chg = 0.0;
                            } else {
                                current_qx_chg  = current_qx_d0;
                                current_qex_chg = 0.0;
                            }
                        } else {
                            current_qx_d0  = (double)iq;  current_qex_d0  = 0.0;
                            current_qx_chg = (double)iq;  current_qex_chg = 0.0;
                        }

                        d0_q[npts]    = dv;  d0_q_e[npts]   = de;
                        chg_q[npts]   = cv;  chg_q_e[npts]  = ce;
                        
                        qx_d0[npts]   = current_qx_d0;   qex_d0[npts]  = current_qex_d0;
                        qx_chg[npts]  = current_qx_chg;  qex_chg[npts] = current_qex_chg;
                        
                        d0_v.push_back(dv);
                        chg_v.push_back(cv);
                        ++npts;
                    }
                    if (npts < 2) continue;

                    double rho = PearsonR(d0_v, chg_v);
                    (void)rho;

                    double d0_ylo  = cent_d0_ylo[ivn - 2][ic];
                    double d0_yhi  = cent_d0_yhi[ivn - 2][ic];

                    double actual_xmin = std::min(*std::min_element(qx_d0, qx_d0 + npts), *std::min_element(qx_chg, qx_chg + npts));
                    double actual_xmax = std::max(*std::max_element(qx_d0, qx_d0 + npts), *std::max_element(qx_chg, qx_chg + npts));
                    double x_range = actual_xmax - actual_xmin;
                    
                    double frame_xmin = doMeanQ ? std::max(0.0, actual_xmin - 0.15 * x_range) : -0.5;
                    double frame_xmax = doMeanQ ? (actual_xmax + 0.15 * x_range) : ((double)N_QBINS - 0.5);

                    TString base = Form("vn_vs_%s_v%d_%s_%s", typeTag.Data(), ivn, cent_lbl_chg[ic], ptTag.Data());

                    TGraphErrors* gr_d0_plot = new TGraphErrors(npts, qx_d0, d0_q, qex_d0, d0_q_e);
                    gr_d0_plot->SetName(Form("gr_d0_%s",  base.Data()));
                    gr_d0_plot->SetMarkerStyle(20);
                    gr_d0_plot->SetMarkerSize(1.4);
                    gr_d0_plot->SetMarkerColor(kBlue + 1);
                    gr_d0_plot->SetLineColor(kBlue + 1);
                    gr_d0_plot->SetLineWidth(2);

                    TGraphErrors* gr_chg_plot = new TGraphErrors(npts, qx_chg, chg_q, qex_chg, chg_q_e);
                    gr_chg_plot->SetName(Form("gr_chg_%s", base.Data()));
                    gr_chg_plot->SetMarkerStyle(21);
                    gr_chg_plot->SetMarkerSize(1.4);
                    gr_chg_plot->SetMarkerColor(kRed + 1);
                    gr_chg_plot->SetLineColor(kRed + 1);
                    gr_chg_plot->SetLineWidth(2);

                    TString cvName = Form("cv_%s", base.Data());
                    TCanvas* cv = new TCanvas(cvName, base, 700, 550);
                    cv->SetLeftMargin(0.14);
                    cv->SetRightMargin(0.14);
                    cv->SetBottomMargin(0.14);
                    cv->SetTopMargin(0.08);
                    cv->SetGridx();
                    cv->SetGridy();

                    TH2F* hf = new TH2F(Form("hf_%s", base.Data()), "",
                                         10, frame_xmin, frame_xmax,
                                         10, d0_ylo, d0_yhi);
                    hf->GetXaxis()->SetTitle(doMeanQ ? Form("<q_{%d}>", ivn) : Form("q_{%d} bin index", ivn));
                    hf->GetYaxis()->SetTitle(ivn == 2 ? "v_{2}" : "v_{3}");
                    hf->GetXaxis()->SetTitleSize(0.060);
                    hf->GetYaxis()->SetTitleSize(0.060);
                    hf->GetXaxis()->SetTitleOffset(0.95);
                    hf->GetYaxis()->SetTitleOffset(0.85);
                    hf->GetXaxis()->SetLabelSize(0.038);
                    hf->GetYaxis()->SetLabelSize(0.038);
                    if (!doMeanQ) hf->GetXaxis()->SetNdivisions(10);
                    hf->GetYaxis()->SetTitleColor(kBlack);
                    hf->GetYaxis()->SetLabelColor(kBlack);
                    hf->Draw();

                    gr_d0_plot->Draw("P SAME");
                    gr_chg_plot->Draw("P SAME");

                    TLegend* leg = new TLegend(0.52, 0.78, 0.86, 0.90);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->SetTextFont(42);
                    leg->SetTextSize(0.038);
                    leg->AddEntry(gr_d0_plot,  Form("D^{0} v_{%d}", ivn), "pe");
                    leg->AddEntry(gr_chg_plot, Form("charged particle v_{%d}", ivn), "pe");
                    leg->Draw();

                    ltx.SetTextAlign(11);
                    ltx.SetTextFont(42);
                    ltx.SetTextSize(0.038);
                    ltx.DrawLatex(0.17, 0.86, Form("Centrality: %d-%d%%", min_cent_chg[ic], max_cent_chg[ic]));
                    ltx.DrawLatex(0.17, 0.79, Form("%.0f < p_{T}^{D^{0}} < %.0f GeV/c", pt_lo, pt_hi));

                    ltx.SetTextAlign(11);
                    ltx.SetTextFont(42);
                    ltx.SetTextSize(0.060);
                    ltx.DrawLatex(0.14, 0.945, "#bf{CMS} #it{Preliminary}");
                    ltx.SetTextAlign(31);
                    ltx.SetTextFont(42);
                    ltx.SetTextSize(0.052);
                    ltx.DrawLatex(0.86, 0.945, "PbPb 5.36 TeV");
                    ltx.SetTextAlign(11);

                    cv->SaveAs(Form("%s/v%d/%s.pdf", plot_base_dir.Data(), ivn, base.Data()));

                    dir_out->cd();
                    gr_d0_plot->Write();
                    gr_chg_plot->Write();
                    cv->Write();

                    delete leg;
                    delete hf;
                    delete gr_d0_plot;
                    delete gr_chg_plot;
                    delete cv;
                }
            }
        }

        fOut->Write();
        fOut->Close();
        fD0->Close();
        fChg->Close();

        TString sum_v2 = Form("summary_v2_%s_%s", typeTag.Data(), tag);
        TString sum_v3 = Form("summary_v3_%s_%s", typeTag.Data(), tag);

        std::cout << "\nTiling " << typeTag << " PDFs...\n";
        gSystem->Exec(Form("/usr/local/bin/python3.11 tile_vn_vs_qbin.py 2 %s %s", plot_base_dir.Data(), sum_v2.Data()));
        gSystem->Exec(Form("/usr/local/bin/python3.11 tile_vn_vs_qbin.py 3 %s %s", plot_base_dir.Data(), sum_v3.Data()));
    }
    
    std::cout << "\n====================================================================\n";
    std::cout << "DONE! Generated 4 merged PDF files total.\n";
}
