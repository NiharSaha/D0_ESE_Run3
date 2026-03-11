#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFileCollection.h"
#include "TLegend.h"
#include "TMath.h"
#include "TComplex.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include <TNtuple.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TParameter.h>

#include "Analysis_bin.h"

using namespace std;




int fit_mass_and_flow(int target_cent = -1)
{
    auto nt_Data = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/SP_MB11to21_Mar5_PubResults_v0/ROOT/flow_Analysis_out_combined.root");
    if (!nt_Data || nt_Data->IsZombie()) { cerr << "Cannot open ntuple file\n"; return 1; }
    auto input_nt = (TNtuple*)nt_Data->Get("nt");
    if (!input_nt) { cerr << "Cannot find TNtuple 'nt'\n"; return 1; }

    auto inf_MC = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/MC_template_MC2023_Feb16_v2/ROOT/D0_MCtemplate_out_combined.root");
    if (!inf_MC || inf_MC->IsZombie()) { cerr << "Cannot open MC file\n"; return 1; }

    // needed_cen_group: -1 = all coarse bins, 0/1/2 = specific coarse bin
    const int needed_cen_group = target_cent; // target_cent IS the coarse bin index directly

    TH1::AddDirectory(kFALSE);
    gROOT->cd();

    Float_t cen_val, pT_val, mass_val, v2_val, v3_val, dca_val, y_val;
    input_nt->SetBranchAddress("cent", &cen_val);
    input_nt->SetBranchAddress("pT",   &pT_val);
    input_nt->SetBranchAddress("mass", &mass_val);
    input_nt->SetBranchAddress("dca",  &dca_val);
    input_nt->SetBranchAddress("y",    &y_val);
    input_nt->SetBranchAddress("v2",   &v2_val);
    input_nt->SetBranchAddress("v3",   &v3_val);


    TH1D *h_mass_default[N_CENTBINS][N_PTBINS]          = {};
    TH1D *h_mass_v2_fit[N_CENTBINS][N_PTBINS][N_VBINS_V2]       = {};  // was N_VBINS
    TH1D *h_mass_v3_fit[N_CENTBINS][N_PTBINS][N_VBINS_V3]       = {};


    // histogram names produced by save_mass_distributions_forPubResults.C:
    // h_mass_cent<label>_pT<label>_def
    // h_mass_v2_cent<label>_pT<label>_v2bin<iv>
    // h_mass_v3_cent<label>_pT<label>_v3bin<iv>
    /*for (int i_cen = 0; i_cen < N_CENTBINS; ++i_cen) {
        if (needed_cen_group >= 0 && i_cen != needed_cen_group) continue;
        for (int i_pt = 0; i_pt < N_PTBINS; ++i_pt) {

            TString name_def = TString::Format("h_mass_cent%s_%s_def",
                                               cen_name[i_cen], pt_name[i_pt]);
            h_mass_default[i_cen][i_pt] = (TH1D*)massFile->Get(name_def);
            if (!h_mass_default[i_cen][i_pt])
                cerr << "[WARN] Missing: " << name_def << "\n";
            else
                h_mass_default[i_cen][i_pt]->SetDirectory(nullptr);

            for (int iv = 0; iv < N_VBINS; ++iv) {
                TString name_v2 = TString::Format("h_mass_v2_cent%s_%s_v2bin%d",
                                                  cen_name[i_cen], pt_name[i_pt], iv);
                h_mass_v2_fit[i_cen][i_pt][iv] = (TH1D*)massFile->Get(name_v2);
                if (h_mass_v2_fit[i_cen][i_pt][iv])
                    h_mass_v2_fit[i_cen][i_pt][iv]->SetDirectory(nullptr);

                TString name_v3 = TString::Format("h_mass_v3_cent%s_%s_v3bin%d",
                                                  cen_name[i_cen], pt_name[i_pt], iv);
                h_mass_v3_fit[i_cen][i_pt][iv] = (TH1D*)massFile->Get(name_v3);
                if (h_mass_v3_fit[i_cen][i_pt][iv])
                    h_mass_v3_fit[i_cen][i_pt][iv]->SetDirectory(nullptr);
            }
        }
    }*/


    static Double_t vnbinning_v2[N_CENTBINS][N_PTBINS][N_VBINS_V2+1];  // was N_VBINS+1
    static Double_t vnbinning_v3[N_CENTBINS][N_PTBINS][N_VBINS_V3+1];  // was N_VBINS+1
    for (int ic = 0; ic < N_CENTBINS; ++ic) {
        for (int ip = 0; ip < N_PTBINS; ++ip) {
            for (int iv = 0; iv <= N_VBINS_V2; ++iv)
                vnbinning_v2[ic][ip][iv] = vnb_v2_table[ic][ip][iv];
            for (int iv = 0; iv <= N_VBINS_V3; ++iv)
                vnbinning_v3[ic][ip][iv] = vnb_v3_table[ic][ip][iv];
        }
    }

     for (int i_cen = 0; i_cen < N_CENTBINS; ++i_cen) {
        if (needed_cen_group >= 0 && i_cen != needed_cen_group) continue;
        for (int i_pt = 0; i_pt < N_PTBINS; ++i_pt) {

            TString name_def = TString::Format("h_mass_cent%s_%s_def", cen_name[i_cen], pt_name[i_pt]);
            h_mass_default[i_cen][i_pt] = new TH1D(name_def, name_def, N_MASSBINS, fit_range_low, fit_range_high);
            h_mass_default[i_cen][i_pt]->SetDirectory(nullptr);

            // ── 3) separate loop limits for v2 and v3 ──
            for (int iv = 0; iv < N_VBINS_V2; ++iv) {       // was N_VBINS
                TString name_v2 = TString::Format("h_mass_v2_cent%s_%s_v2bin%d", cen_name[i_cen], pt_name[i_pt], iv);
                h_mass_v2_fit[i_cen][i_pt][iv] = new TH1D(name_v2, name_v2, N_MASSBINS, fit_range_low, fit_range_high);
                h_mass_v2_fit[i_cen][i_pt][iv]->SetDirectory(nullptr);
            }
            for (int iv = 0; iv < N_VBINS_V3; ++iv) {       // was N_VBINS
                TString name_v3 = TString::Format("h_mass_v3_cent%s_%s_v3bin%d", cen_name[i_cen], pt_name[i_pt], iv);
                h_mass_v3_fit[i_cen][i_pt][iv] = new TH1D(name_v3, name_v3, N_MASSBINS, fit_range_low, fit_range_high);
                h_mass_v3_fit[i_cen][i_pt][iv]->SetDirectory(nullptr);
            }
        }
    }

    cout << "Filling histograms from ntuple..." << endl;
    Long64_t N_ENTRIES = input_nt->GetEntriesFast();
    for (Long64_t i_entry = 0; i_entry < N_ENTRIES; ++i_entry) {
        input_nt->GetEntry(i_entry);
        if (i_entry % 1000000 == 0)
            cout << i_entry << " / " << N_ENTRIES << "  " << 100*i_entry/N_ENTRIES << "%" << endl;

        if (dca_val >= 0.0085) continue;
        if (TMath::Abs(y_val) >= 1.0) continue;

        for (int i_cen = 0; i_cen < N_CENTBINS; ++i_cen) {
            if (needed_cen_group >= 0 && i_cen != needed_cen_group) continue;
            if (cen_val < cen_coarse_edges[i_cen] || cen_val >= cen_coarse_edges[i_cen+1]) continue;

            for (int i_pt = 0; i_pt < N_PTBINS; ++i_pt) {
                if (pT_val < pt_edges[i_pt] || pT_val >= pt_edges[i_pt+1]) continue;

                if (h_mass_default[i_cen][i_pt])
                    h_mass_default[i_cen][i_pt]->Fill(mass_val);

                for (int iv = 0; iv < N_VBINS_V2; ++iv) {   // was N_VBINS
                    if (v2_val >= vnbinning_v2[i_cen][i_pt][iv] && v2_val < vnbinning_v2[i_cen][i_pt][iv+1])
                        if (h_mass_v2_fit[i_cen][i_pt][iv]) h_mass_v2_fit[i_cen][i_pt][iv]->Fill(mass_val);
                }
                for (int iv = 0; iv < N_VBINS_V3; ++iv) {   // was N_VBINS
                    if (v3_val >= vnbinning_v3[i_cen][i_pt][iv] && v3_val < vnbinning_v3[i_cen][i_pt][iv+1])
                        if (h_mass_v3_fit[i_cen][i_pt][iv]) h_mass_v3_fit[i_cen][i_pt][iv]->Fill(mass_val);
                }
            }
        }
    }
    cout << "Ntuple loop done." << endl;
    nt_Data->Close();



    Double_t yield_v2[N_VBINS_V2]={}, yield_error_v2[N_VBINS_V2]={};
    Double_t v2_x[N_VBINS_V2]={}, v2_x_err[N_VBINS_V2]={};
    Double_t yield_v3[N_VBINS_V3]={}, yield_error_v3[N_VBINS_V3]={};
    Double_t v3_x[N_VBINS_V3]={}, v3_x_err[N_VBINS_V3]={};
    Double_t pt_x[N_PTBINS]={}, pt_x_error[N_PTBINS]={};
    Double_t mean_val_v2[N_PTBINS]={}, mean_error_v2[N_PTBINS]={};
    Double_t mean_val_v3[N_PTBINS]={}, mean_error_v3[N_PTBINS]={};

    Double_t v2_x_store[N_CENTBINS][N_PTBINS][N_VBINS_V2] = {};
    Double_t v2_x_err_store[N_CENTBINS][N_PTBINS][N_VBINS_V2] = {};
    Double_t v3_x_store[N_CENTBINS][N_PTBINS][N_VBINS_V3] = {};
    Double_t v3_x_err_store[N_CENTBINS][N_PTBINS][N_VBINS_V3] = {};

    double chi2_ndf_v2[N_CENTBINS][N_PTBINS][N_VBINS_V2] = {};
    double chi2_ndf_v3[N_CENTBINS][N_PTBINS][N_VBINS_V3] = {};
    double sigma_v2[N_CENTBINS][N_PTBINS][N_VBINS_V2] = {};
    double sigma_v3[N_CENTBINS][N_PTBINS][N_VBINS_V3] = {};
    double good_v2[N_CENTBINS][N_PTBINS][N_VBINS_V2] = {};
    double good_v3[N_CENTBINS][N_PTBINS][N_VBINS_V3] = {};

    TGraphErrors *yield_v2_graph[N_CENTBINS][N_PTBINS] = {};
    TGraphErrors *yield_v3_graph[N_CENTBINS][N_PTBINS] = {};
    TGraph       *sigma_v2_fit[N_CENTBINS][N_PTBINS]   = {};
    TGraph       *sigma_v3_fit[N_CENTBINS][N_PTBINS]   = {};
    TGraph       *chi2_v2_graph[N_CENTBINS][N_PTBINS]  = {};
    TGraph       *chi2_v3_graph[N_CENTBINS][N_PTBINS]  = {};
    TH1D         *h_v2_hist[N_CENTBINS][N_PTBINS]      = {};
    TH1D         *h_v3_hist[N_CENTBINS][N_PTBINS]      = {};
    TGraphErrors *v2_cen[N_CENTBINS]                   = {};
    TGraphErrors *v3_cen[N_CENTBINS]                   = {};

    // book yield histograms with variable bin edges from vnbinning_manual
    for (int i_cen = 0; i_cen < N_CENTBINS; ++i_cen) {
        if (needed_cen_group >= 0 && i_cen != needed_cen_group) continue;
        for (int i_pt = 0; i_pt < N_PTBINS; ++i_pt) {
            TString nv2 = TString::Format("h_v2_hist_%s_%s", cen_name[i_cen], pt_name[i_pt]);
            h_v2_hist[i_cen][i_pt] = new TH1D(nv2, nv2, N_VBINS_V2, vnbinning_v2[i_cen][i_pt]);
            h_v2_hist[i_cen][i_pt]->Sumw2();
            h_v2_hist[i_cen][i_pt]->SetDirectory(nullptr);

            TString nv3 = TString::Format("h_v3_hist_%s_%s", cen_name[i_cen], pt_name[i_pt]);
            h_v3_hist[i_cen][i_pt] = new TH1D(nv3, nv3, N_VBINS_V3, vnbinning_v3[i_cen][i_pt]);
            h_v3_hist[i_cen][i_pt]->Sumw2();
            h_v3_hist[i_cen][i_pt]->SetDirectory(nullptr);
        }
    }


    std::string outname = "Flow_output_pubcomp";
    if (target_cent >= 0) outname += std::string("_cen") + cen_name[target_cent];
    const char *env_job  = std::getenv("SLURM_ARRAY_JOB_ID");
    if (!env_job) env_job = std::getenv("SLURM_JOB_ID");
    const char *env_task = std::getenv("SLURM_ARRAY_TASK_ID");
    if (env_job)  outname += std::string("_job")  + env_job;
    if (env_task) outname += std::string("_task") + env_task;
    outname += ".root";

    // ---- open log file with same base name ----
    std::string logname = outname;
    logname.replace(logname.size()-5, 5, "_yield_check.txt");  // replace .root -> _yield_check.txt
    std::ofstream logf(logname);
    if (!logf.is_open()) cerr << "[WARN] Cannot open log file: " << logname << endl;
    else cout << "Log file: " << logname << endl;

    TFile *outf = new TFile(outname.c_str(), "recreate");
    if (!outf || outf->IsZombie()) { cerr << "Cannot create output file\n"; return 1; }

    TDirectory *dir_yield_vn       = outf->mkdir("yields_vN");
    TDirectory *dir_hist_sp        = outf->mkdir("histograms_SP");
    TDirectory *dir_summary        = outf->mkdir("summary_pt");
    TDirectory *dir_diag           = outf->mkdir("diagnostics_chi2_sigma");
    TDirectory *dir_mass_filled = outf->mkdir("mass_filled");
    TDirectory *dir_mass_fitted    = outf->mkdir("mass_fitted");

    dir_mass_filled->cd();
    for (int i_cen = 0; i_cen < N_CENTBINS; ++i_cen) {
        if (needed_cen_group >= 0 && i_cen != needed_cen_group) continue;
        for (int i_pt = 0; i_pt < N_PTBINS; ++i_pt) {
            if (h_mass_default[i_cen][i_pt]) h_mass_default[i_cen][i_pt]->Write();
            for (int iv = 0; iv < N_VBINS_V2; ++iv) {
                if (h_mass_v2_fit[i_cen][i_pt][iv]) h_mass_v2_fit[i_cen][i_pt][iv]->Write();
            }
            for (int iv = 0; iv < N_VBINS_V3; ++iv) {
                if (h_mass_v3_fit[i_cen][i_pt][iv]) h_mass_v3_fit[i_cen][i_pt][iv]->Write();
            }
        }
    }

    TLatex *tex = new TLatex;
    tex->SetNDC(); tex->SetTextFont(42); tex->SetTextSize(0.045); tex->SetLineWidth(2);
    auto c1 = new TCanvas("c1","c1",600,600);
    gSystem->mkdir("prompt_mass_plot_withchi2_sigma", true);


    for (int i_cen = 0; i_cen < N_CENTBINS; ++i_cen) {
        if (needed_cen_group >= 0 && i_cen != needed_cen_group) continue;

        for (int i_pt = 0; i_pt < N_PTBINS; ++i_pt) {

            TH1D *h_mc_match_signal;
            TH1D *h_mc_match_all;


            h_mc_match_signal = (TH1D*)inf_MC->Get(Form("hMass_Signal_%s", pt_name[i_pt]));
            h_mc_match_all    = (TH1D*)inf_MC->Get(Form("hMass_Signal_Plus_Swap_%s", pt_name[i_pt]));

            if (h_mc_match_signal) h_mc_match_signal->SetTitle(Form("hMass_Signal_%s", pt_name[i_pt]));
            if (h_mc_match_all)    h_mc_match_all->SetTitle(Form("hMass_Signal_Plus_Swap_%s", pt_name[i_pt]));

        
            // for high pT bins fall back to pT bin 9 (20-40 GeV)
            if (i_pt == 10 || i_pt == 11) {
                h_mc_match_signal = (TH1D*)inf_MC->Get(Form("hMass_Signal_%s", pt_name[9]));
                h_mc_match_all    = (TH1D*)inf_MC->Get(Form("hMass_Signal_Plus_Swap_%s", pt_name[9]));

                if (h_mc_match_signal) h_mc_match_signal->SetTitle(Form("hMass_Signal_%s", pt_name[9]));
                if (h_mc_match_all)    h_mc_match_all->SetTitle(Form("hMass_Signal_Plus_Swap_%s", pt_name[9]));
            }

            //TH1D *h_def  = h_mass_default[i_cen][i_pt];
            TH1D *hist_mass = h_mass_default[i_cen][i_pt];  
            
            
        

            int ymax = hist_mass->GetBinContent(hist_mass->GetMaximumBin());


            TF1 *f = new TF1(Form("f_ptbin_%d_%d", i_cen, i_pt),
                "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0+[6]))/(sqrt(2*3.14159)*[2]*(1.0+[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0+[6]))/(sqrt(2*3.14159)*[3]*(1.0+[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0+[6]))/(sqrt(2*3.14159)*[7]*(1.0+[6])))+ [8]+[9]*x+[10]*x*x+ [11]*ROOT::Math::crystalball_function(x,2.2,17,0.0267*(1+[6]),1.96*(1+[12]))+ 4*[11]*(ROOT::Math::crystalball_function(x,0.34,5,0.0146*(1+[6]),1.7734*(1+[13])))",
                fit_range_low, fit_range_high);

            f->SetLineColor(2); f->SetLineWidth(1);
            f->SetParameter(0, 100.); 
            f->SetParameter(1, 1.8648);
            f->SetParameter(2, 0.03); 
            f->SetParameter(3, 0.005);
            f->SetParameter(4, 0.1);
            f->FixParameter(5, 1);   
            f->FixParameter(6, 0);
            f->FixParameter(7, 0.1); 
            f->FixParameter(8, 0);
            f->FixParameter(9, 0);   
            f->FixParameter(10, 0);
            f->FixParameter(11, 0);  
            f->FixParameter(12, 0);
            f->FixParameter(13, 0);

            if (i_pt < 5) { f->SetParLimits(2,0.01,0.5);  f->SetParLimits(3,0.001,0.25); }
            else { f->SetParLimits(2,0.005,0.15); f->SetParLimits(3,0.001,0.08); }
            f->SetParLimits(4,0,1); f->SetParLimits(5,0,1);


            if (h_mc_match_signal) {
                f->FixParameter(1, 1.8648);
                h_mc_match_signal->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "q",  "", fit_range_low, fit_range_high);
                h_mc_match_signal->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "q",  "", fit_range_low, fit_range_high);
                f->ReleaseParameter(1);
                h_mc_match_signal->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "L q", "", fit_range_low, fit_range_high);
                h_mc_match_signal->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "L q", "", fit_range_low, fit_range_high);
                h_mc_match_signal->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "L q m","", fit_range_low, fit_range_high);
                f->FixParameter(1, f->GetParameter(1));
                f->FixParameter(2, f->GetParameter(2));
                f->FixParameter(3, f->GetParameter(3));
                f->FixParameter(4, f->GetParameter(4));
            }

            if (h_mc_match_all) {
                f->ReleaseParameter(5); 
                f->ReleaseParameter(7);
                f->SetParameter(7, 0.1);
                h_mc_match_all->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt),"L q","",fit_range_low,fit_range_high);
                h_mc_match_all->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt),"L q","",fit_range_low,fit_range_high);
                h_mc_match_all->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt),"L q","",fit_range_low,fit_range_high);
                h_mc_match_all->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt),"L q","",fit_range_low,fit_range_high);
                h_mc_match_all->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt),"L q","",fit_range_low,fit_range_high);  
            }
            f->FixParameter(5, f->GetParameter(5));
            f->FixParameter(7, f->GetParameter(7));
            f->FixParameter(1, f->GetParameter(1));
            f->FixParameter(2, f->GetParameter(2));
            f->FixParameter(3, f->GetParameter(3));
            f->FixParameter(4, f->GetParameter(4));
            //f->FixParameter(6, 0);
            f->ReleaseParameter(8); f->ReleaseParameter(9); f->ReleaseParameter(10);


            /*if (h_def) {
                f->FixParameter(1, 1.8648); f->FixParameter(6, 0);
                h_def->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "Lq",  "", fit_range_low, fit_range_high);
                h_def->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "Lq",  "", fit_range_low, fit_range_high);
                f->ReleaseParameter(1); f->ReleaseParameter(6);
                
                h_def->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "Lq", "", fit_range_low, fit_range_high);
                h_def->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "Lq", "", fit_range_low, fit_range_high);
                h_def->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "Lq", "", fit_range_low, fit_range_high);

                f->FixParameter(6, f->GetParameter(6));
                f->FixParameter(1, f->GetParameter(1));
            }*/

            f->FixParameter(1, 1.8648); f->FixParameter(6, 0);
            hist_mass->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "q",   "", fit_range_low, fit_range_high);
            hist_mass->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "q",   "", fit_range_low, fit_range_high);
            hist_mass->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "q",   "", fit_range_low, fit_range_high);
            f->ReleaseParameter(1);
            f->ReleaseParameter(6);
            hist_mass->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "L q", "", fit_range_low, fit_range_high);
            hist_mass->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "L q", "", fit_range_low, fit_range_high);
            f->FixParameter(1, f->GetParameter(1));
            hist_mass->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "L q", "", fit_range_low, fit_range_high);
            hist_mass->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "L q m", "", fit_range_low, fit_range_high);

            if (i_pt != 0) { f->ReleaseParameter(11); 
                f->SetParLimits(11,0,f->GetParameter(0)*(1-f->GetParameter(5))); }
            else f->FixParameter(11, 0);


            hist_mass->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "L q m", "", fit_range_low, fit_range_high);
            hist_mass->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "L q m", "", fit_range_low, fit_range_high);
            f->SetParLimits(12,-0.03,0.03); 
            f->SetParLimits(13,-0.03,0.03);
            f->ReleaseParameter(12); 
            f->ReleaseParameter(13);
            
            
            hist_mass->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "L q m", "", fit_range_low, fit_range_high);
            hist_mass->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "L q m", "", fit_range_low, fit_range_high);
            hist_mass->Fit(Form("f_ptbin_%d_%d",i_cen,i_pt), "L q m", "", fit_range_low, fit_range_high);


            hist_mass->SetMinimum(0);
            hist_mass->SetMarkerSize(0.5); hist_mass->SetMarkerStyle(20);
            hist_mass->SetLineWidth(1); hist_mass->SetOption("e");
            hist_mass->SetStats(kFALSE);
            hist_mass->GetXaxis()->SetRangeUser(fit_range_low, fit_range_high);
            hist_mass->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
            hist_mass->GetYaxis()->SetTitle("Entries / 5 MeV");
            hist_mass->GetYaxis()->SetRangeUser(0, 1.35*ymax);
            hist_mass->SetTitle(Form("Mass_%s_%s", cen_name[i_cen], pt_name[i_pt]));
            hist_mass->Draw("ep");

            // component functions
            TF1 *f1 = new TF1(Form("f_sig_%d_%d",i_cen,i_pt),
                "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0+[6]))/(sqrt(2*3.14159)*[2]*(1.0+[6]))"
                "+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0+[6]))/(sqrt(2*3.14159)*[3]*(1.0+[6]))))", fit_range_low, fit_range_high);
            f1->SetLineColor(kOrange-3); f1->SetLineWidth(1); f1->SetLineStyle(2);
            f1->SetFillColorAlpha(kOrange-3,0.3); f1->SetFillStyle(1001);
            for (int p=0; p<=6; ++p) f1->FixParameter(p, f->GetParameter(p));
            f1->Draw("LSAME");

            TF1 *f2 = new TF1(Form("f_swap_%d_%d",i_cen,i_pt), "[0]*((1-[1])*TMath::Gaus(x,[4],[3]*(1.0+[2]))/(sqrt(2*3.14159)*[3]*(1.0+[2])))", fit_range_low, fit_range_high);
            f2->SetLineColor(kGreen+4); f2->SetLineWidth(1); f2->SetLineStyle(1);
            f2->SetFillColorAlpha(kGreen+4,0.3); f2->SetFillStyle(1001);
            f2->FixParameter(0,f->GetParameter(0)); f2->FixParameter(1,f->GetParameter(5));
            f2->FixParameter(2,f->GetParameter(6)); f2->FixParameter(3,f->GetParameter(7));
            f2->FixParameter(4,f->GetParameter(1));
            f2->Draw("LSAME");

            TF1 *f3 = new TF1(Form("f_bkg_%d_%d",i_cen,i_pt),"[0]+[1]*x+[2]*x*x", fit_range_low, fit_range_high);
            f3->SetLineColor(4); 
            f3->SetLineWidth(1); 
            f3->SetLineStyle(2);
            f3->FixParameter(0,f->GetParameter(8)); 
            f3->FixParameter(1,f->GetParameter(9));
            f3->FixParameter(2,f->GetParameter(10));
            f3->Draw("LSAME");

            TF1 *f4 = new TF1(Form("f_kk_%d_%d",i_cen,i_pt),
                "[0]*ROOT::Math::crystalball_function(x,2.2,17,0.0267*(1+[1]),1.96*(1+[2]))"
                "+ 4*[0]*(ROOT::Math::crystalball_function(x,0.34,5,0.0146*(1+[1]),1.7734*(1+[3])))",
                fit_range_low, fit_range_high);
            f4->SetLineColor(kViolet-4); f4->SetLineWidth(1); f4->SetLineStyle(1);
            f4->SetFillColorAlpha(kViolet-4,0.3); f4->SetFillStyle(1001);
            f4->FixParameter(0,f->GetParameter(11)); f4->FixParameter(1,f->GetParameter(6));
            f4->FixParameter(2,f->GetParameter(12)); f4->FixParameter(3,f->GetParameter(13));
            f4->Draw("LSAME");

            double yield_val = f->GetParameter(0)*f->GetParameter(5)/width;
            double yield_err = (f->GetParError(0)*f->GetParameter(5)+f->GetParError(5)*f->GetParameter(0))/width;
            double significance = (yield_err>0) ? fabs(yield_val/yield_err) : 0.;
            double chi2ndf_mass = (f->GetNDF()>0) ? f->GetChisquare()/f->GetNDF() : -1.;

            TLegend *leg = new TLegend(0.65,0.6,0.81,0.9,NULL,"brNDC");
            leg->SetBorderSize(0); leg->SetTextSize(0.035); leg->SetTextFont(42); leg->SetFillStyle(0);
            leg->AddEntry(hist_mass," Data","lep");
            leg->AddEntry(f," Fit","L");
            leg->AddEntry(f1," D^{0}+#bar{D^{#lower[0.2]{0}}} Signal","l");
            leg->AddEntry(f2," K-#pi swap","l");
            leg->AddEntry(f4," K-#bar{K}, #pi-#bar{#pi}","l");
            leg->AddEntry(f3," Combinatorial","l");
            leg->Draw("SAME");

            TLatex *texYield = new TLatex(0.14,0.80,Form("Yield = %.0f #pm %.0f",yield_val,yield_err));
            texYield->SetNDC(); texYield->SetTextFont(42); texYield->SetTextSize(0.035); texYield->Draw();
            TLatex *texSign  = new TLatex(0.14,0.75,Form("Y/#DeltaY = %.2f",significance));
            texSign->SetNDC();  texSign->SetTextFont(42); texSign->SetTextSize(0.035);  texSign->Draw();
            TLatex *texChi   = new TLatex(0.14,0.70,Form("#chi^{2}/ndf = %.2f",chi2ndf_mass));
            texChi->SetNDC();   texChi->SetTextFont(42); texChi->SetTextSize(0.035);   texChi->Draw();
            TLatex *texCENT  = new TLatex(0.14,0.86,Form("Cent %d-%d%%",
                                          cen_coarse_edges[i_cen], cen_coarse_edges[i_cen+1]));
            texCENT->SetNDC();  texCENT->SetTextFont(42); texCENT->SetTextSize(0.035); texCENT->Draw();

            //c1->SaveAs(Form("prompt_mass_plot_withchi2_sigma/hmassfit_def_%s_%s.pdf", cen_name[i_cen], pt_name[i_pt]));
            //c1->Clear();
            hist_mass->GetListOfFunctions()->Add(f);
            hist_mass->GetListOfFunctions()->Add(f1);
            hist_mass->GetListOfFunctions()->Add(f2);
            hist_mass->GetListOfFunctions()->Add(f3);
            hist_mass->GetListOfFunctions()->Add(f4);
            hist_mass->GetListOfFunctions()->Add(leg);
            
            // write mass histograms
            dir_mass_fitted->cd();
            // write MC clones once per (cen, pt)
         if (h_mc_match_signal) {
                TH1D *cl = (TH1D*)h_mc_match_signal->Clone(Form("hMass_Signal_cen%s_%s",cen_name[i_cen],pt_name[i_pt]));
                cl->Write(); delete cl;
            }
            if (h_mc_match_all) {
                TH1D *cl = (TH1D*)h_mc_match_all->Clone(Form("hMass_Signal_Plus_Swap_cen%s_%s",cen_name[i_cen],pt_name[i_pt]));
                cl->Write(); delete cl;
            }
            
            hist_mass->Write();


            for (int i_v = 0; i_v < N_VBINS_V2; ++i_v) {

                v2_x[i_v]     = 0.5*(vnbinning_v2[i_cen][i_pt][i_v] + vnbinning_v2[i_cen][i_pt][i_v+1]);
                v2_x_err[i_v] = 0.5*fabs(vnbinning_v2[i_cen][i_pt][i_v+1] - vnbinning_v2[i_cen][i_pt][i_v]);
               v2_x_store[i_cen][i_pt][i_v]     = v2_x[i_v];
                v2_x_err_store[i_cen][i_pt][i_v] = v2_x_err[i_v];


                TH1D *h_v2 = h_mass_v2_fit[i_cen][i_pt][i_v];
                if (!h_v2 || h_v2->GetEntries() < 10) {
                    yield_v2[i_v]=0; yield_error_v2[i_v]=0;
                    chi2_ndf_v2[i_cen][i_pt][i_v]=-10; sigma_v2[i_cen][i_pt][i_v]=-10;
                    continue;
                }

                 TF1 *fitFcn_v2 = new TF1(Form("fit_v2_%d_%d_%d",i_cen,i_pt,i_v),
                    "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0+[6]))/(sqrt(2*3.14159)*[2]*(1.0+[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0+[6]))/(sqrt(2*3.14159)*[3]*(1.0+[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0+[6]))/(sqrt(2*3.14159)*[7]*(1.0+[6])))+ [8]+[9]*x+[10]*x*x+ [11]*ROOT::Math::crystalball_function(x,2.2,17,0.0267*(1+[6]),1.96*(1+[12]))+ 4*[11]*(ROOT::Math::crystalball_function(x,0.34,5,0.0146*(1+[6]),1.7734*(1+[13])))",
                    fit_range_low, fit_range_high);

                fitFcn_v2->SetParameter(0, 100);
                for (int p=1; p<=7;  ++p) fitFcn_v2->FixParameter(p, f->GetParameter(p));
                fitFcn_v2->FixParameter(11, f->GetParameter(11));
                fitFcn_v2->FixParameter(12, f->GetParameter(12));
                fitFcn_v2->FixParameter(13, f->GetParameter(13));

                // Float: normalization [0] and background [8],[9],[10]
                fitFcn_v2->SetParameter(0, f->GetParameter(0));  // start from inclusive yield
                fitFcn_v2->SetParameter(8, f->GetParameter(8));
                fitFcn_v2->SetParameter(9, f->GetParameter(9));
                fitFcn_v2->SetParameter(10, f->GetParameter(10));

                h_v2->Fit(fitFcn_v2,"M",  "", fit_range_low, fit_range_high);
                h_v2->Fit(fitFcn_v2,"LQ", "", fit_range_low, fit_range_high);
                h_v2->Fit(fitFcn_v2,"LQ", "", fit_range_low, fit_range_high);
                h_v2->Fit(fitFcn_v2,"LM", "", fit_range_low, fit_range_high);

                chi2_ndf_v2[i_cen][i_pt][i_v] = (fitFcn_v2->GetNDF()>0) ?
                    fitFcn_v2->GetChisquare()/fitFcn_v2->GetNDF() : -1.;

                yield_v2[i_v]       = fitFcn_v2->GetParameter(0)*fitFcn_v2->GetParameter(5)/width;
                yield_error_v2[i_v] = fitFcn_v2->GetParError(0)*fitFcn_v2->GetParameter(5)/width;
                sigma_v2[i_cen][i_pt][i_v] = (yield_error_v2[i_v]>0) ? yield_v2[i_v]/yield_error_v2[i_v] : -10;

                /*double chi2_v2   = chi2_ndf_v2[i_cen][i_pt][i_v];
                double sig_v2    = sigma_v2[i_cen][i_pt][i_v];        // yield/error
                double norm_v2   = fitFcn_v2->GetParameter(0);
                double norm_err_v2 = fitFcn_v2->GetParError(0);

                
                bool chi2_ok_v2  = (chi2_v2 > 0.0 && chi2_v2 < 6.0);
                bool sig_ok_v2   = (sig_v2 > 1.0);
                bool norm_ok_v2  = (norm_v2 > 0 && norm_err_v2 > 0 && norm_v2 > norm_err_v2);
                bool ndf_ok_v2   = (fitFcn_v2->GetNDF() > 0);
                double inclusive_yield = f->GetParameter(0)*f->GetParameter(5)/width;
                bool range_ok_v2 = (yield_v2[i_v] < 3.0 * inclusive_yield);

                good_v2[i_cen][i_pt][i_v] = chi2_ok_v2 && sig_ok_v2 && norm_ok_v2 && ndf_ok_v2 && range_ok_v2;

                if (!good_v2[i_cen][i_pt][i_v]) {
                    yield_v2[i_v]                  = 0;
                    yield_error_v2[i_v]            = 0;
                    chi2_ndf_v2[i_cen][i_pt][i_v]  = -10;
                    sigma_v2[i_cen][i_pt][i_v]     = -10;
                }*/


                /*if (yield_v2[i_v] <= 0) {
                        yield_v2[i_v]=0; yield_error_v2[i_v]=0;
                        chi2_ndf_v2[i_cen][i_pt][i_v]=-10;
                        sigma_v2[i_cen][i_pt][i_v]=-10;
                    }*/

                    if (yield_v2[i_v] <= 0) {
                    bool in_sp_range = (v2_x[i_v] >= -5.0 && v2_x[i_v] <= 5.0);
                    if (in_sp_range && i_v > 0 && yield_v2[i_v-1] > 0) {
                        logf << "[WARN] yield_v2 <= 0 at cen=" << cen_name[i_cen]
                             << " pt=" << pt_name[i_pt]
                             << " i_v=" << i_v
                             << " (v2_x=" << v2_x[i_v] << ", in SP range)"
                             << " => using prev bin values (yield=" << yield_v2[i_v-1]
                             << " +/- " << yield_error_v2[i_v-1] << ")\n";
                        yield_v2[i_v]       = yield_v2[i_v-1];
                        yield_error_v2[i_v] = yield_error_v2[i_v-1];
                        chi2_ndf_v2[i_cen][i_pt][i_v] = -10;
                        sigma_v2[i_cen][i_pt][i_v]    = -10;
                    } else {
                        logf << "[WARN] yield_v2 <= 0 at cen=" << cen_name[i_cen]
                             << " pt=" << pt_name[i_pt]
                             << " i_v=" << i_v
                             << " (v2_x=" << v2_x[i_v] << ")"
                             << (in_sp_range ? " [SP range but no valid prev bin]" : " [outside SP range]")
                             << " => zeroing out.\n";
                        yield_v2[i_v]                      = 0;
                        yield_error_v2[i_v]                = 0;
                        chi2_ndf_v2[i_cen][i_pt][i_v]      = -10;
                        sigma_v2[i_cen][i_pt][i_v]         = -10;
                    }
                }
                

                h_v2_hist[i_cen][i_pt]->SetBinContent(i_v+1, yield_v2[i_v]);
                h_v2_hist[i_cen][i_pt]->SetBinError  (i_v+1, yield_error_v2[i_v]);

                // draw v2 bin mass fit
                c1->cd();
                h_v2->SetStats(kFALSE);
                h_v2->GetXaxis()->SetRangeUser(fit_range_low,fit_range_high);
                h_v2->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
                h_v2->GetYaxis()->SetTitle("Entries / 5 MeV");
                h_v2->SetMarkerStyle(20); h_v2->SetMarkerSize(0.5); h_v2->SetLineWidth(1);
                h_v2->Draw("ep");
                fitFcn_v2->SetLineColor(2); 
                fitFcn_v2->SetLineWidth(1); 
                fitFcn_v2->Draw("LSAME");

                TLatex *texCENT_v2 = new TLatex(0.14,0.85,Form("%d < Cent < %d",
                                      cen_coarse_edges[i_cen],cen_coarse_edges[i_cen+1]));
                texCENT_v2->SetNDC(); texCENT_v2->SetTextFont(42); texCENT_v2->SetTextSize(0.035); texCENT_v2->Draw();
                TLatex *texPT_v2   = new TLatex(0.14,0.80,Form("%.1f < p_{T} < %.1f GeV/c",
                                      pt_edges[i_pt],pt_edges[i_pt+1]));
                texPT_v2->SetNDC();   texPT_v2->SetTextFont(42);   texPT_v2->SetTextSize(0.035); texPT_v2->Draw();
                TLatex *texVbin_v2 = new TLatex(0.14, 0.75, Form("%.3f < v_{2}^{i} < %.3f",
                          vnbinning_v2[i_cen][i_pt][i_v], vnbinning_v2[i_cen][i_pt][i_v+1]));
                texVbin_v2->SetNDC(); texVbin_v2->SetTextFont(42); texVbin_v2->SetTextSize(0.035); texVbin_v2->Draw();
                TLatex *texY_v2    = new TLatex(0.14,0.70,Form("Yield = %.0f #pm %.0f",yield_v2[i_v],yield_error_v2[i_v]));
                texY_v2->SetNDC();    texY_v2->SetTextFont(42);    texY_v2->SetTextSize(0.035); texY_v2->Draw();
                TLatex *texC_v2    = new TLatex(0.14,0.65,Form("#chi^{2}/ndf = %.2f",chi2_ndf_v2[i_cen][i_pt][i_v]));
                texC_v2->SetNDC();    texC_v2->SetTextFont(42);    texC_v2->SetTextSize(0.035); texC_v2->Draw();

                //c1->SaveAs(Form("prompt_mass_plot_withchi2_sigma/hmassfit_%s_%s_v2bin_%d.pdf", cen_name[i_cen], pt_name[i_pt], i_v));
                c1->Clear();
                }// ---i_v loop for v2

                for (int i_v = 0; i_v < N_VBINS_V3; ++i_v) {

                  v3_x[i_v]     = 0.5*(vnbinning_v3[i_cen][i_pt][i_v] + vnbinning_v3[i_cen][i_pt][i_v+1]);
                  v3_x_err[i_v] = 0.5*fabs(vnbinning_v3[i_cen][i_pt][i_v+1] - vnbinning_v3[i_cen][i_pt][i_v]);
                  v3_x_store[i_cen][i_pt][i_v]     = v3_x[i_v];
                  v3_x_err_store[i_cen][i_pt][i_v] = v3_x_err[i_v];


                TH1D *h_v3 = h_mass_v3_fit[i_cen][i_pt][i_v];
                if (!h_v3 || h_v3->GetEntries() < 10) {
                    yield_v3[i_v]=0; yield_error_v3[i_v]=0;
                    chi2_ndf_v3[i_cen][i_pt][i_v]=-10; sigma_v3[i_cen][i_pt][i_v]=-10;
                    continue;
                }

                TF1 *fitFcn_v3 = new TF1(Form("fit_v3_%d_%d_%d",i_cen,i_pt,i_v),
                    "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0+[6]))/(sqrt(2*3.14159)*[2]*(1.0+[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0+[6]))/(sqrt(2*3.14159)*[3]*(1.0+[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0+[6]))/(sqrt(2*3.14159)*[7]*(1.0+[6])))+ [8]+[9]*x+[10]*x*x+ [11]*ROOT::Math::crystalball_function(x,2.2,17,0.0267*(1+[6]),1.96*(1+[12]))+ 4*[11]*(ROOT::Math::crystalball_function(x,0.34,5,0.0146*(1+[6]),1.7734*(1+[13])))",
                    fit_range_low, fit_range_high);

                fitFcn_v3->SetParameter(0,100);
                for (int p=1; p<=7; ++p) fitFcn_v3->FixParameter(p, f->GetParameter(p));
                //fitFcn_v3->SetParameter(8,1); fitFcn_v3->SetParameter(9,1); fitFcn_v3->SetParameter(10,1);
                //fitFcn_v3->FixParameter(11,0);
                fitFcn_v3->FixParameter(11, f->GetParameter(11));
                fitFcn_v3->FixParameter(12, f->GetParameter(12));
                fitFcn_v3->FixParameter(13, f->GetParameter(13));

                // Float: normalization [0] and background [8],[9],[10]
                fitFcn_v3->SetParameter(0,  f->GetParameter(0));
                fitFcn_v3->SetParameter(8,  f->GetParameter(8));
                fitFcn_v3->SetParameter(9,  f->GetParameter(9));
                fitFcn_v3->SetParameter(10, f->GetParameter(10));


                h_v3->Fit(fitFcn_v3,"M",  "", fit_range_low, fit_range_high);
                h_v3->Fit(fitFcn_v3,"LQ", "", fit_range_low, fit_range_high);
                h_v3->Fit(fitFcn_v3,"LQ", "", fit_range_low, fit_range_high);
                h_v3->Fit(fitFcn_v3,"LM", "", fit_range_low, fit_range_high);

                chi2_ndf_v3[i_cen][i_pt][i_v] = (fitFcn_v3->GetNDF()>0) ?
                    fitFcn_v3->GetChisquare()/fitFcn_v3->GetNDF() : -1.;

                yield_v3[i_v]       = fitFcn_v3->GetParameter(0)*fitFcn_v3->GetParameter(5)/width;
                yield_error_v3[i_v] = fitFcn_v3->GetParError(0)*fitFcn_v3->GetParameter(5)/width;
                sigma_v3[i_cen][i_pt][i_v] = (yield_error_v3[i_v]>0) ? yield_v3[i_v]/yield_error_v3[i_v] : -10;

                /*double chi2_v3_val = chi2_ndf_v3[i_cen][i_pt][i_v];
                double sig_v3      = sigma_v3[i_cen][i_pt][i_v];
                double norm_v3     = fitFcn_v3->GetParameter(0);
                double norm_err_v3 = fitFcn_v3->GetParError(0);

                bool chi2_ok_v3  = (chi2_v3_val > 0.0 && chi2_v3_val < 6.0);
                bool sig_ok_v3   = (sig_v3 > 1.0);
                bool norm_ok_v3  = (norm_v3 > 0 && norm_err_v3 > 0 && norm_v3 > norm_err_v3);
                bool ndf_ok_v3   = (fitFcn_v3->GetNDF() > 0);
                double inclusive_yield = f->GetParameter(0)*f->GetParameter(5)/width;
                bool range_ok_v3 = (yield_v3[i_v] < 3.0 * inclusive_yield);

                good_v3[i_cen][i_pt][i_v] = chi2_ok_v3 && sig_ok_v3 && norm_ok_v3 && ndf_ok_v3 && range_ok_v3;

                if (!good_v3[i_cen][i_pt][i_v]) {
                    yield_v3[i_v]                  = 0;
                    yield_error_v3[i_v]            = 0;
                    chi2_ndf_v3[i_cen][i_pt][i_v]  = -10;
                    sigma_v3[i_cen][i_pt][i_v]     = -10;
                }*/

                /*if (yield_v3[i_v] <= 0) {
                        yield_v3[i_v]=0; yield_error_v3[i_v]=0;
                        chi2_ndf_v3[i_cen][i_pt][i_v]=-10;
                        sigma_v3[i_cen][i_pt][i_v]=-10;
                    }*/

                    if (yield_v3[i_v] <= 0) {
                    bool in_sp_range = (v3_x[i_v] >= -5.0 && v3_x[i_v] <= 5.0);
                    if (in_sp_range && i_v > 0 && yield_v3[i_v-1] > 0) {
                        logf << "[WARN] yield_v3 <= 0 at cen=" << cen_name[i_cen]
                             << " pt=" << pt_name[i_pt]
                             << " i_v=" << i_v
                             << " (v3_x=" << v3_x[i_v] << ", in SP range)"
                             << " => using prev bin values (yield=" << yield_v3[i_v-1]
                             << " +/- " << yield_error_v3[i_v-1] << ")\n";
                        yield_v3[i_v]       = yield_v3[i_v-1];
                        yield_error_v3[i_v] = yield_error_v3[i_v-1];
                        chi2_ndf_v3[i_cen][i_pt][i_v] = -10;
                        sigma_v3[i_cen][i_pt][i_v]    = -10;
                    } else {
                        logf << "[WARN] yield_v3 <= 0 at cen=" << cen_name[i_cen]
                             << " pt=" << pt_name[i_pt]
                             << " i_v=" << i_v
                             << " (v3_x=" << v3_x[i_v] << ")"
                             << (in_sp_range ? " [SP range but no valid prev bin]" : " [outside SP range]")
                             << " => zeroing out.\n";
                        yield_v3[i_v]                      = 0;
                        yield_error_v3[i_v]                = 0;
                        chi2_ndf_v3[i_cen][i_pt][i_v]      = -10;
                        sigma_v3[i_cen][i_pt][i_v]         = -10;
                    }
                }
                

                h_v3_hist[i_cen][i_pt]->SetBinContent(i_v+1, yield_v3[i_v]);
                h_v3_hist[i_cen][i_pt]->SetBinError  (i_v+1, yield_error_v3[i_v]);

                // draw v3 bin mass fit
                c1->cd();
                h_v3->SetStats(kFALSE);
                h_v3->GetXaxis()->SetRangeUser(fit_range_low,fit_range_high);
                h_v3->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
                h_v3->GetYaxis()->SetTitle("Entries / 5 MeV");
                h_v3->SetMarkerStyle(20); h_v3->SetMarkerSize(0.5); h_v3->SetLineWidth(1);
                h_v3->Draw("ep");
                fitFcn_v3->SetLineColor(2); fitFcn_v3->SetLineWidth(1); fitFcn_v3->Draw("LSAME");

                TLatex *texCENT_v3 = new TLatex(0.14,0.85,Form("%d < Cent < %d",
                                      cen_coarse_edges[i_cen],cen_coarse_edges[i_cen+1]));
                texCENT_v3->SetNDC(); texCENT_v3->SetTextFont(42); texCENT_v3->SetTextSize(0.035); texCENT_v3->Draw();
                TLatex *texPT_v3   = new TLatex(0.14,0.80,Form("%.1f < p_{T} < %.1f GeV/c",
                                      pt_edges[i_pt],pt_edges[i_pt+1]));
                texPT_v3->SetNDC();   texPT_v3->SetTextFont(42);   texPT_v3->SetTextSize(0.035); texPT_v3->Draw();
                TLatex *texVbin_v3 = new TLatex(0.14,0.75,Form("%.3f < v_{3}^{i} < %.3f",
                                      vnbinning_v3[i_cen][i_pt][i_v],vnbinning_v3[i_cen][i_pt][i_v+1]));
                texVbin_v3->SetNDC(); texVbin_v3->SetTextFont(42); texVbin_v3->SetTextSize(0.035); texVbin_v3->Draw();
                TLatex *texY_v3    = new TLatex(0.14,0.70,Form("Yield = %.0f #pm %.0f",yield_v3[i_v],yield_error_v3[i_v]));
                texY_v3->SetNDC();    texY_v3->SetTextFont(42);    texY_v3->SetTextSize(0.035); texY_v3->Draw();
                TLatex *texC_v3    = new TLatex(0.14,0.65,Form("#chi^{2}/ndf = %.2f",chi2_ndf_v3[i_cen][i_pt][i_v]));
                texC_v3->SetNDC();    texC_v3->SetTextFont(42);    texC_v3->SetTextSize(0.035); texC_v3->Draw();

                //c1->SaveAs(Form("prompt_mass_plot_withchi2_sigma/hmassfit_%s_%s_v3bin_%d.pdf", cen_name[i_cen], pt_name[i_pt], i_v));
                c1->Clear();

            }

            // build yield graphs per (cen, pt)
            yield_v2_graph[i_cen][i_pt] = new TGraphErrors(N_VBINS_V2,v2_x,yield_v2,v2_x_err,yield_error_v2);
            yield_v2_graph[i_cen][i_pt]->SetNameTitle(
                Form("v2_yield_graph_%s_%s",cen_name[i_cen],pt_name[i_pt]),
                Form("v2_yield_graph_%s_%s",cen_name[i_cen],pt_name[i_pt]));

            yield_v3_graph[i_cen][i_pt] = new TGraphErrors(N_VBINS_V3,v3_x,yield_v3,v3_x_err,yield_error_v3);
            yield_v3_graph[i_cen][i_pt]->SetNameTitle(
                Form("v3_yield_graph_%s_%s",cen_name[i_cen],pt_name[i_pt]),
                Form("v3_yield_graph_%s_%s",cen_name[i_cen],pt_name[i_pt]));

        } 

    } 

    
    for (int i_cen = 0; i_cen < N_CENTBINS; ++i_cen) {
        if (needed_cen_group >= 0 && i_cen != needed_cen_group) continue;

        for (int i_pt = 0; i_pt < N_PTBINS; ++i_pt) {
            pt_x[i_pt]       = 0.5*(pt_edges[i_pt]+pt_edges[i_pt+1]);
            pt_x_error[i_pt] = 0.5*(pt_edges[i_pt+1]-pt_edges[i_pt]);

            mean_val_v2[i_pt]   = h_v2_hist[i_cen][i_pt] ? h_v2_hist[i_cen][i_pt]->GetMean()      : 0;
            mean_error_v2[i_pt] = h_v2_hist[i_cen][i_pt] ? h_v2_hist[i_cen][i_pt]->GetMeanError() : 0;
            mean_val_v3[i_pt]   = h_v3_hist[i_cen][i_pt] ? h_v3_hist[i_cen][i_pt]->GetMean()      : 0;
            mean_error_v3[i_pt] = h_v3_hist[i_cen][i_pt] ? h_v3_hist[i_cen][i_pt]->GetMeanError() : 0;

            chi2_v2_graph[i_cen][i_pt] = new TGraph(N_VBINS_V2,v2_x_store[i_cen][i_pt],chi2_ndf_v2[i_cen][i_pt]);
            chi2_v2_graph[i_cen][i_pt]->SetNameTitle(
                Form("chi2ndf_v2_%s_%s",cen_name[i_cen],pt_name[i_pt]),
                Form("chi2ndf_v2_%s_%s",cen_name[i_cen],pt_name[i_pt]));
            chi2_v2_graph[i_cen][i_pt]->GetXaxis()->SetTitle("v_{2}^{i}");
            chi2_v2_graph[i_cen][i_pt]->GetYaxis()->SetTitle("#chi^{2}/ndf");
            chi2_v2_graph[i_cen][i_pt]->SetMarkerStyle(20); chi2_v2_graph[i_cen][i_pt]->SetMarkerSize(0.8);

            chi2_v3_graph[i_cen][i_pt] = new TGraph(N_VBINS_V3,v3_x_store[i_cen][i_pt],chi2_ndf_v3[i_cen][i_pt]);
            chi2_v3_graph[i_cen][i_pt]->SetNameTitle(
                Form("chi2ndf_v3_%s_%s",cen_name[i_cen],pt_name[i_pt]),
                Form("chi2ndf_v3_%s_%s",cen_name[i_cen],pt_name[i_pt]));
            chi2_v3_graph[i_cen][i_pt]->GetXaxis()->SetTitle("v_{3}^{i}");
            chi2_v3_graph[i_cen][i_pt]->GetYaxis()->SetTitle("#chi^{2}/ndf");
            chi2_v3_graph[i_cen][i_pt]->SetMarkerStyle(20); chi2_v3_graph[i_cen][i_pt]->SetMarkerSize(0.8);

            sigma_v2_fit[i_cen][i_pt] = new TGraph(N_VBINS_V2,v2_x_store[i_cen][i_pt],sigma_v2[i_cen][i_pt]);
            sigma_v2_fit[i_cen][i_pt]->SetNameTitle(
                Form("sigma_v2_%s_%s",cen_name[i_cen],pt_name[i_pt]),
                Form("sigma_v2_%s_%s",cen_name[i_cen],pt_name[i_pt]));
            sigma_v2_fit[i_cen][i_pt]->SetMarkerStyle(20); sigma_v2_fit[i_cen][i_pt]->SetMarkerSize(0.8);

            sigma_v3_fit[i_cen][i_pt] = new TGraph(N_VBINS_V3,v3_x_store[i_cen][i_pt],sigma_v3[i_cen][i_pt]);
            sigma_v3_fit[i_cen][i_pt]->SetNameTitle(
                Form("sigma_v3_%s_%s",cen_name[i_cen],pt_name[i_pt]),
                Form("sigma_v3_%s_%s",cen_name[i_cen],pt_name[i_pt]));
            sigma_v3_fit[i_cen][i_pt]->SetMarkerStyle(20); sigma_v3_fit[i_cen][i_pt]->SetMarkerSize(0.8);
        }

        v2_cen[i_cen] = new TGraphErrors(N_PTBINS,pt_x,mean_val_v2,pt_x_error,mean_error_v2);
        v2_cen[i_cen]->SetNameTitle(Form("v2_graph_%s",cen_name[i_cen]),
                                    Form("v2_graph_%s",cen_name[i_cen]));
        v3_cen[i_cen] = new TGraphErrors(N_PTBINS,pt_x,mean_val_v3,pt_x_error,mean_error_v3);
        v3_cen[i_cen]->SetNameTitle(Form("v3_graph_%s",cen_name[i_cen]),
                                    Form("v3_graph_%s",cen_name[i_cen]));
    }

    
    for (int i_cen = 0; i_cen < N_CENTBINS; ++i_cen) {
        if (needed_cen_group >= 0 && i_cen != needed_cen_group) continue;
        for (int i_pt = 0; i_pt < N_PTBINS; ++i_pt) {
            dir_yield_vn->cd();
            if (yield_v2_graph[i_cen][i_pt]) yield_v2_graph[i_cen][i_pt]->Write();
            if (yield_v3_graph[i_cen][i_pt]) yield_v3_graph[i_cen][i_pt]->Write();

            dir_hist_sp->cd();
            if (h_v2_hist[i_cen][i_pt]) h_v2_hist[i_cen][i_pt]->Write();
            if (h_v3_hist[i_cen][i_pt]) h_v3_hist[i_cen][i_pt]->Write();

            dir_diag->cd();
            if (chi2_v2_graph[i_cen][i_pt])  chi2_v2_graph[i_cen][i_pt]->Write();
            if (chi2_v3_graph[i_cen][i_pt])  chi2_v3_graph[i_cen][i_pt]->Write();
            if (sigma_v2_fit[i_cen][i_pt])   sigma_v2_fit[i_cen][i_pt]->Write();
            if (sigma_v3_fit[i_cen][i_pt])   sigma_v3_fit[i_cen][i_pt]->Write();
            
            // ---- index-based versions (x = 0,1,2,...) ----
            {
                Double_t idx_v2[N_VBINS_V2], idx_v3[N_VBINS_V3];
                for (int iv = 0; iv < N_VBINS_V2; ++iv) idx_v2[iv] = iv;
                for (int iv = 0; iv < N_VBINS_V3; ++iv) idx_v3[iv] = iv;

                // chi2/ndf v2 vs bin index
                TGraph *chi2_v2_idx = new TGraph(N_VBINS_V2, idx_v2, chi2_ndf_v2[i_cen][i_pt]);
                chi2_v2_idx->SetNameTitle(
                    Form("chi2ndf_v2_idx_%s_%s", cen_name[i_cen], pt_name[i_pt]),
                    Form("chi2ndf_v2_idx_%s_%s", cen_name[i_cen], pt_name[i_pt]));
                chi2_v2_idx->GetXaxis()->SetTitle("v_{2} bin index");
                chi2_v2_idx->GetYaxis()->SetTitle("#chi^{2}/ndf");
                chi2_v2_idx->SetMarkerStyle(20); chi2_v2_idx->SetMarkerSize(0.8);
                chi2_v2_idx->Write();
                delete chi2_v2_idx;

                // chi2/ndf v3 vs bin index
                TGraph *chi2_v3_idx = new TGraph(N_VBINS_V3, idx_v3, chi2_ndf_v3[i_cen][i_pt]);
                chi2_v3_idx->SetNameTitle(
                    Form("chi2ndf_v3_idx_%s_%s", cen_name[i_cen], pt_name[i_pt]),
                    Form("chi2ndf_v3_idx_%s_%s", cen_name[i_cen], pt_name[i_pt]));
                chi2_v3_idx->GetXaxis()->SetTitle("v_{3} bin index");
                chi2_v3_idx->GetYaxis()->SetTitle("#chi^{2}/ndf");
                chi2_v3_idx->SetMarkerStyle(20); chi2_v3_idx->SetMarkerSize(0.8);
                chi2_v3_idx->Write();
                delete chi2_v3_idx;

                // sigma (yield/error) v2 vs bin index
                TGraph *sigma_v2_idx = new TGraph(N_VBINS_V2, idx_v2, sigma_v2[i_cen][i_pt]);
                sigma_v2_idx->SetNameTitle(
                    Form("sigma_v2_idx_%s_%s", cen_name[i_cen], pt_name[i_pt]),
                    Form("sigma_v2_idx_%s_%s", cen_name[i_cen], pt_name[i_pt]));
                sigma_v2_idx->GetXaxis()->SetTitle("v_{2} bin index");
                sigma_v2_idx->GetYaxis()->SetTitle("Yield / #DeltaYield");
                sigma_v2_idx->SetMarkerStyle(20); sigma_v2_idx->SetMarkerSize(0.8);
                sigma_v2_idx->Write();
                delete sigma_v2_idx;

                // sigma (yield/error) v3 vs bin index
                TGraph *sigma_v3_idx = new TGraph(N_VBINS_V3, idx_v3, sigma_v3[i_cen][i_pt]);
                sigma_v3_idx->SetNameTitle(
                    Form("sigma_v3_idx_%s_%s", cen_name[i_cen], pt_name[i_pt]),
                    Form("sigma_v3_idx_%s_%s", cen_name[i_cen], pt_name[i_pt]));
                sigma_v3_idx->GetXaxis()->SetTitle("v_{3} bin index");
                sigma_v3_idx->GetYaxis()->SetTitle("Yield / #DeltaYield");
                sigma_v3_idx->SetMarkerStyle(20); sigma_v3_idx->SetMarkerSize(0.8);
                sigma_v3_idx->Write();
                delete sigma_v3_idx;
            }





        }
        dir_summary->cd();
        if (v2_cen[i_cen]) v2_cen[i_cen]->Write();
        if (v3_cen[i_cen]) v3_cen[i_cen]->Write();
    }

    // Note: outf->Write() intentionally removed — all objects written individually
    outf->Close();
    inf_MC->Close();

    delete tex; tex = nullptr;
    delete c1;  c1  = nullptr;

    cout << "Saved " << outname << endl;
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc == 1) {
        fit_mass_and_flow();
    } else if (argc == 2) {
        int idx = atoi(argv[1]);
        if (idx < 0 || idx >= N_CENTBINS) {
            cerr << "ERROR: target_cent must be 0-" << N_CENTBINS-1
                 << " (0=0-10%, 1=10-30%, 2=30-50%)\n"; return 1;
        }
        fit_mass_and_flow(idx);
    } else {
        cout << "Usage: ./fit_mass_and_flow [target_cent (0=0-10,1=10-30,2=30-50)]\n";
        return 1;
    }
    return 0;
}
