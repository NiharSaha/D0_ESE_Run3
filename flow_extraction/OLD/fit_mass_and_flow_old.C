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

//#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/ESE_Cuts_MB11to21_Jan22.h"

using namespace std;


int fit_mass_and_flow(int target1pct = -1, int target_pt = -1)
{
  const int N_QBINS = 10;
  const int N_CENTBINS = 4;
  const int N_PTBINS = 10;
  const int N_VBINS = 44;
  const int N_MASSBINS = 52;
  Float_t fit_range_low = 1.74;
  Float_t fit_range_high = 2.00;
  Float_t width = (fit_range_high-fit_range_low)/N_MASSBINS;
  
  TString cen_label[N_CENTBINS]={"cent0to10","cent10to30","cent30to50","cent50to90"};
  std::string cen_name[N_CENTBINS] = {"cent0to10","cent10to30","cent30to50","cent50to90"};
  std::string pt_name[N_PTBINS] = {"pT1to2","pT2to3","pT3to4","pT4to5","pT5to6","pT6to8","pT8to10","pT10to15","pT15to20","pT20to40"};
  double pt_edges[N_PTBINS+1] = {1,2,3,4,5,6,8,10,15,20,40};



  
  auto inf_MC = TFile::Open("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/D0mass_MC_5p36TeV_OfficialMC_Jan26/ROOT/D0_Mass_MC_out_combined.root");

  TFile *massFile = TFile::Open("mass_distributions.root","READ");
  if (!massFile || massFile->IsZombie()) {
    std::cerr<<"Cannot open mass_distributions.root. Run save_mass_distributions first."<<std::endl;
    return 1;
  }





  TH1D* h_mass_default[N_CENTBINS][N_PTBINS] = {};
  TH1D* h_mass[N_CENTBINS][N_QBINS][N_PTBINS] = {};
  TH1D* h_mass_v2_fit[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS] = {};
  TH1D* h_mass_v3_fit[N_CENTBINS][N_QBINS][N_PTBINS][N_VBINS] = {};
  
  // initialize empty grouped histograms (same binning as originals)
  for (int ic=0; ic<N_CENTBINS; ++ic) {
    for (int iq=0; iq<N_QBINS; ++iq) {
      for (int ip=0; ip<N_PTBINS; ++ip) {
        TString hm_name = Form("h_mass_%s_q2bin%d_%s", cen_name[ic].c_str(), iq, pt_name[ip].c_str());
        h_mass[ic][iq][ip] = new TH1D(hm_name, hm_name, N_MASSBINS, fit_range_low, fit_range_high);
        h_mass[ic][iq][ip]->Sumw2();
        for (int iv=0; iv<N_VBINS; ++iv) {
          TString hv2 = Form("hist_mass_v2_%s_q2bin%d_pt_%s_in_v2bin_idx_%d", cen_name[ic].c_str(), iq, pt_name[ip].c_str(), iv);
          h_mass_v2_fit[ic][iq][ip][iv] = new TH1D(hv2, hv2, N_MASSBINS, fit_range_low, fit_range_high);
          h_mass_v2_fit[ic][iq][ip][iv]->Sumw2();
          TString hv3 = Form("hist_mass_v3_%s_q2bin%d_pt_%s_in_v3bin_idx_%d", cen_name[ic].c_str(), iq, pt_name[ip].c_str(), iv);
          h_mass_v3_fit[ic][iq][ip][iv] = new TH1D(hv3, hv3, N_MASSBINS, fit_range_low, fit_range_high);
          h_mass_v3_fit[ic][iq][ip][iv]->Sumw2();
        }
        TString hdef = Form("h_mass_%s_%s_def", cen_name[ic].c_str(), pt_name[ip].c_str());
        h_mass_default[ic][ip] = new TH1D(hdef, hdef, N_MASSBINS, fit_range_low, fit_range_high);
        h_mass_default[ic][ip]->Sumw2();
      }
    }
  }
  
  const int N_CENTBINS_1 = N_CENTBINS*10; // 40 if N_CENTBINS==4
  for (int i_pt=0; i_pt<N_PTBINS; ++i_pt) {
    //if (target_pt >= 0 && i_pt != target_pt) continue;
    for (int i_cen=0; i_cen<N_CENTBINS_1; ++i_cen) {
      //if (target1pct >= 0 && i_cen != target1pct) continue;
      int cen_group = i_cen/10; // map 1% sub-bin -> coarse centrality group
      if (cen_group < 0 || cen_group >= N_CENTBINS) continue;
      
      // default 1% histogram
      std::string name_def = std::string("h_mass_cen") + std::to_string(i_cen+1) + "_" + pt_name[i_pt] + "_def";
      TH1D* h_def = (TH1D*)massFile->Get(name_def.c_str());
      if (h_def && h_mass_default[cen_group][i_pt]) h_mass_default[cen_group][i_pt]->Add(h_def);

      for (int i_q2=0; i_q2<N_QBINS; ++i_q2) {
        // q2 histogram
        std::string name_q2 = std::string("h_mass_cen") + std::to_string(i_cen+1) + "_q2bin" + std::to_string(i_q2) + "_" + pt_name[i_pt];
        TH1D* h_q2 = (TH1D*)massFile->Get(name_q2.c_str());
        if (h_q2 && h_mass[cen_group][i_q2][i_pt]) h_mass[cen_group][i_q2][i_pt]->Add(h_q2);

        // per-v2-bin histograms
        for (int i_v=0; i_v<N_VBINS; ++i_v) {
          std::string name_v2 = std::string("hist_mass_v2_cen") + std::to_string(i_cen+1) + "_q2bin" + std::to_string(i_q2) + "_" + pt_name[i_pt] + "_in_v2bin_idx_" + std::to_string(i_v);
          TH1D* h_v2 = (TH1D*)massFile->Get(name_v2.c_str());
          if (h_v2 && h_mass_v2_fit[cen_group][i_q2][i_pt][i_v]) h_mass_v2_fit[cen_group][i_q2][i_pt][i_v]->Add(h_v2);

	  std::string name_v3 = std::string("hist_mass_v3_cen") + std::to_string(i_cen+1) + "_q2bin" + std::to_string(i_q2) + "_" + pt_name[i_pt] + "_in_v3bin_idx_" + std::to_string(i_v);
	  TH1D* h_v3 = (TH1D*)massFile->Get(name_v3.c_str());
          if (h_v3 && h_mass_v3_fit[cen_group][i_q2][i_pt][i_v]) h_mass_v3_fit[cen_group][i_q2][i_pt][i_v]->Add(h_v3);
	} // i_v
	
      } // i_q2
    } // i_cen
  } // i_pt
  
  
  
  
  // Prepare output structures
  Double_t yield_v2[N_VBINS], yield_error_v2[N_VBINS], v2_x[N_VBINS], v2_x_err[N_VBINS], mean_val_v2[N_PTBINS], mean_error_v2[N_PTBINS];
  Double_t yield_v3[N_VBINS], yield_error_v3[N_VBINS], v3_x[N_VBINS], v3_x_err[N_VBINS], mean_val_v3[N_PTBINS], mean_error_v3[N_PTBINS];
  Double_t pt_x[N_PTBINS], pt_x_error[N_PTBINS];
 
  TGraphErrors* yield_v2_graph[N_CENTBINS][N_QBINS][N_PTBINS];
  TGraphErrors* yield_v3_graph[N_CENTBINS][N_QBINS][N_PTBINS];
  TGraph* sigma_v2_fit[N_CENTBINS][N_QBINS][N_PTBINS];
  TGraph* sigma_v3_fit[N_CENTBINS][N_QBINS][N_PTBINS];
  TH1D* h_v2_hist[N_CENTBINS][N_QBINS][N_PTBINS] = {};
  TH1D* h_v3_hist[N_CENTBINS][N_QBINS][N_PTBINS] = {};

  TGraphErrors* v2_cen[N_CENTBINS][N_Q2BINS];
  TGraph* chi2_ndf_for_q2_v2_graph[N_CENTBINS][N_Q2BINS][N_PTBINS];

  TGraphErrors* v3_cen[N_CENTBINS][N_Q2BINS];
  TGraph* chi2_ndf_for_q3_v3_graph[N_CENTBINS][N_Q2BINS][N_PTBINS];

  

  TH1D *h_mass = nullptr;

  TLatex* tex = new TLatex;
  tex->SetNDC(); tex->SetTextFont(42); tex->SetTextSize(0.045); tex->SetLineWidth(2);
  auto c1= new TCanvas("c1","c1",600,600);

  // Loop and fit (retrieve histograms from saved file)
  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {

    for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {

      TH1D* h_mc_match_signal = (TH1D*)inf_MC->Get(Form("hMass_Signal_%s", pt_name[i_pt].c_str()));
      TH1D* h_mc_match_all = (TH1D*)inf_MC->Get(Form("hMass_Swap_%s", pt_name[i_pt].c_str()));



      for (int i_q=0; i_q<N_QBINS; i_q++) {


	h_mass = h_mass[i_cen][i_q][i_pt];
        if (!h_mass) {
          std::cerr << "[WARN] Missing grouped histogram for cen="<<cen_name[i_cen]<<", q2="<<i_q<<", pt="<<pt_name[i_pt]<<")\n";
          continue;
        }
	
        int ymax = h_mass->GetBinContent(h_mass->GetMaximumBin());

        h_mass->SetMinimum(0);
        h_mass->SetMarkerSize(0.5);
        h_mass->SetTitle("");
        h_mass->SetMarkerStyle(20);
        h_mass->SetLineWidth(1);
        h_mass->SetOption("e");
        h_mass->GetXaxis()->SetRangeUser(fit_range_low,fit_range_high);
        h_mass->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
        h_mass->GetYaxis()->SetTitle("Entries / 2.5 MeV");
        h_mass->SetStats(kFALSE);

        TF1* f = new TF1(Form("f_ptbin_%d",i_pt),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])) ) + [8] + [9]*x + [10]*x*x + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12])) + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])) )",fit_range_low,fit_range_high);
        f->SetLineColor(2);
	f->SetLineWidth(1);
        f->SetParameter(0,100.);
	f->SetParameter(1,1.8648);
	f->SetParameter(2,0.03);
	f->SetParameter(3,0.005);
	f->SetParameter(4,0.1);
        f->FixParameter(5,1);
	f->FixParameter(6,0);
	f->FixParameter(7,0.1);
	f->FixParameter(8,0);
	f->FixParameter(9,0);
	f->FixParameter(10,0);
	f->FixParameter(11,0);
	f->FixParameter(12,0);
	f->FixParameter(13,0);
        f->SetParLimits(2,0.01,0.1);
	f->SetParLimits(3,0.001,0.05);
	f->SetParLimits(4,0,1);
	f->SetParLimits(5,0,1);

        // fit MC templates to get fixed params

        if (h_mc_match_signal) {
          f->FixParameter(1,1.8648);
          h_mc_match_signal->Fit(Form("f_ptbin_%d",i_pt),"q","",fit_range_low,fit_range_high);
          h_mc_match_signal->Fit(Form("f_ptbin_%d",i_pt),"q","",fit_range_low,fit_range_high);
          f->ReleaseParameter(1);
          h_mc_match_signal->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
          h_mc_match_signal->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
          h_mc_match_signal->Fit(Form("f_ptbin_%d",i_pt),"L q m","",fit_range_low,fit_range_high);
          f->FixParameter(1,f->GetParameter(1)); f->FixParameter(2,f->GetParameter(2)); f->FixParameter(3,f->GetParameter(3)); f->FixParameter(4,f->GetParameter(4));
        }
        if (h_mc_match_all) {
          f->ReleaseParameter(5);
	  f->ReleaseParameter(7);
          f->SetParameter(7,0.1);
          h_mc_match_all->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
          h_mc_match_all->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
          h_mc_match_all->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
          f->FixParameter(5,f->GetParameter(5)); f->FixParameter(7,f->GetParameter(7));
        }

        // Fit data
        h_mass->Fit(Form("f_ptbin_%d",i_pt),"q","",fit_range_low,fit_range_high);
        h_mass->Fit(Form("f_ptbin_%d",i_pt),"q","",fit_range_low,fit_range_high);
        f->ReleaseParameter(1); f->ReleaseParameter(6); f->SetParameter(6,0); f->SetParLimits(6,-1,1);
        h_mass->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
        h_mass->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
        h_mass->Fit(Form("f_ptbin_%d",i_pt),"L q","",fit_range_low,fit_range_high);
        h_mass->Fit(Form("f_ptbin_%d",i_pt),"L q m","",fit_range_low,fit_range_high);
        if (i_pt!=0) { f->ReleaseParameter(11); f->SetParLimits(11,0,f->GetParameter(0)*(1-f->GetParameter(5))); } else { f->FixParameter(11,0); }
        h_mass->Fit(Form("f_ptbin_%d",i_pt),"L q m","",fit_range_low,fit_range_high);
        h_mass->Fit(Form("f_ptbin_%d",i_pt),"L q m","",fit_range_low,fit_range_high);
        f->SetParLimits(12,-0.03,0.03); f->SetParLimits(13,-0.03,0.03);
        f->ReleaseParameter(12); f->ReleaseParameter(13);
        h_mass->Fit(Form("f_ptbin_%d",i_pt),"L q m","",fit_range_low,fit_range_high);
        h_mass->Fit(Form("f_ptbin_%d",i_pt),"L q m","",fit_range_low,fit_range_high);

        h_mass->GetYaxis()->SetRangeUser(0,1.35*ymax);
        h_mass->Draw("ep");


	//draw D0 signal separately
        TF1* f1 = new TF1(Form("f_sig_%d",i_pt),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
        f1->SetLineColor(kOrange-3);
        f1->SetLineWidth(1);
        f1->SetLineStyle(2);
        f1->SetFillColorAlpha(kOrange-3,0.3);
        f1->SetFillStyle(1001);
        f1->FixParameter(0,f->GetParameter(0));
        f1->FixParameter(1,f->GetParameter(1));
        f1->FixParameter(2,f->GetParameter(2));
        f1->FixParameter(3,f->GetParameter(3));
        f1->FixParameter(4,f->GetParameter(4));
        f1->FixParameter(5,f->GetParameter(5));
        f1->FixParameter(6,f->GetParameter(6));

        f1->Draw("LSAME");

        //write Yield uncertainty         
        TF1* f1e = new TF1(Form("f_e_%d",i_pt),"[0]", fit_range_low, fit_range_high);
        f1e->SetLineColor(kOrange-3);
        f1e->SetLineWidth(1);
        f1e->SetLineStyle(2);
        f1e->SetFillColorAlpha(kOrange-3,0.3);
        f1e->SetFillStyle(1001);
        f1e->FixParameter(0,(f->GetParError(0)*f->GetParameter(5)+f->GetParError(5)*f->GetParameter(0))/(f->GetParameter(0)*f->GetParameter(5)));

        //draw swap separately
        TF1* f2 = new TF1(Form("f_swap_%d",i_pt),"[0]*((1-[1])*TMath::Gaus(x,[4],[3]*(1.0 +[2]))/(sqrt(2*3.14159)*[3]*(1.0 +[2])))", fit_range_low, fit_range_high);
        f2->SetLineColor(kGreen+4);
        f2->SetLineWidth(1);
        f2->SetLineStyle(1);
        f2->SetFillColorAlpha(kGreen+4,0.3);
        f2->SetFillStyle(1001);
        f2->FixParameter(0,f->GetParameter(0));
        f2->FixParameter(1,f->GetParameter(5));
        f2->FixParameter(2,f->GetParameter(6));
        f2->FixParameter(3,f->GetParameter(7));
        f2->FixParameter(4,f->GetParameter(1));

        f2->Draw("LSAME");

        //draw KK branch separately 
        TF1* f4 = new TF1(Form("f_swap_%d",i_pt),"[0]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[1]), 1.96*(1+[2])) + 4*[0]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[1]), 1.7734*(1+[3])) )", fit_range_low, fit_range_high);
        f4->SetLineColor(kViolet-4);
        f4->SetLineWidth(1);
        f4->SetLineStyle(1);
        f4->SetFillColorAlpha(kViolet-4,0.3);
        f4->SetFillStyle(1001);
        f4->FixParameter(0,f->GetParameter(11));
        f4->FixParameter(1,f->GetParameter(6));
        f4->FixParameter(2,f->GetParameter(12));
        f4->FixParameter(3,f->GetParameter(13));

        f4->Draw("LSAME");

        //draw poly bkg separately
        TF1* f3 = new TF1(Form("f_bkg_%d",i_pt),"[0] + [1]*x + [2]*x*x ", fit_range_low, fit_range_high);
        f3->SetLineColor(4);
        f3->SetLineWidth(1);
        f3->SetLineStyle(2);
        f3->FixParameter(0,f->GetParameter(8));
        f3->FixParameter(1,f->GetParameter(9));
        f3->FixParameter(2,f->GetParameter(10));
        //f3->FixParameter(3,0);
        f3->Draw("LSAME");

        tex->DrawLatex(0.14,0.86,Form("cent %s pT %s GeV/c q2 %d",cen_label[i_cen].Data(),pt_label[i_pt].Data(),i_q));
        //tex->DrawLatex(0.37,0.86,Form("%s GeV/c",pt_label[i_pt].Data()));
        //texCMS->DrawLatex(.18,.97,"#font[61]{CMS} #it{Preliminary}");
        //texCMS->DrawLatex(0.62,0.97, "#scale[0.8]{PbPb #sqrt{s_{NN}} = 5.02 TeV}");

        TLegend* leg = new TLegend(0.65,0.6,0.81,0.9,NULL,"brNDC");
        leg->SetBorderSize(0);
        leg->SetTextSize(0.035);
        leg->SetTextFont(42);
        leg->SetFillStyle(0);
        //leg->AddEntry(h_mass[i_cen][i_q][i_pt]," Data","lep"); //for full q2
        leg->AddEntry(h_mass[i_cen][i_q][i_pt]," Data","lep");
        leg->AddEntry(f," Fit","L");
        leg->AddEntry(f1," D^{0}+#bar{D^{#lower[0.2]{0}}} Signal","f");
        leg->AddEntry(f2," K-#pi swap","f");
        leg->AddEntry(f4," K-#bar{K}, #pi-#bar{#pi}","f");
        leg->AddEntry(f3," Combinatorial","l");
        leg->Draw("SAME");

        int Yield = f->GetParameter(0)*f->GetParameter(5)/width;
        int YieldEr = (f->GetParError(0)*f->GetParameter(5)+f->GetParError(5)*f->GetParameter(0))/width;
        tex->DrawLatex(0.14,0.75,Form("Yield = %d #pm %d",Yield,YieldEr));
        tex->DrawLatex(0.14,0.80,Form("Significance: %.2f",fabs(1.0*Yield/YieldEr)));

        //c1->SaveAs(Form("prompt_default_fit_plots/massplot_cen_%s_q2bin_%d_pt_%s.pdf",cen_name[i_cen].Data(),i_q,pt_name[i_pt].Data()));



	for (int i_v=0; i_v<N_VBINS; i_v++) {


	  //=====================
	  // For v3 calculations
	  //=====================
	  
	  if (h_mass_v2_fit[i_cen][i_q][i_pt][i_v]->GetEntries()<10)
	    {
	      yield_v2[i_v] = 0; 
	      yield_error_v2[i_v] = 0;
	      chi2_ndf_for_q2_v2[i_cen][i_q][i_pt][i_v] = -10;
	      sigma_v2[i_cen][i_q][i_pt][i_v] = -10;
	      continue;
	    }

          TF1* fitFcn_v2 = new TF1(Form("fit_v2_%d_%d_%d", i_pt, i_q, iv),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])) ) + [8] + [9]*x + [10]*x*x + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12])) + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])) )", fit_range_low, fit_range_high);

	  fitFcn_v2->SetParameter(0,100);
          fitFcn_v2->FixParameter(1,f->GetParameter(1));
          fitFcn_v2->FixParameter(2,f->GetParameter(2));
          fitFcn_v2->FixParameter(3,f->GetParameter(3));
          fitFcn_v2->FixParameter(4,f->GetParameter(4));
          fitFcn_v2->FixParameter(5,f->GetParameter(5));
          fitFcn_v2->FixParameter(6,f->GetParameter(6));
          fitFcn_v2->FixParameter(7,f->GetParameter(7));
          fitFcn_v2->SetParameter(8,1);
          fitFcn_v2->SetParameter(9,1);
          fitFcn_v2->SetParameter(10,1);
          fitFcn_v2->FixParameter(11,0);

          h_mass_v2_fit[i_cen][i_q][i_pt][i_v]->Draw("AEP");
          h_mass_v2_fit[i_cen][i_q][i_pt][i_v]->Fit(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v),"M","",fit_range_low,fit_range_high);
          h_mass_v2_fit[i_cen][i_q][i_pt][i_v]->Fit(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v),"L q","",fit_range_low,fit_range_high);
          h_mass_v2_fit[i_cen][i_q][i_pt][i_v]->Fit(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v),"L q","",fit_range_low,fit_range_high);
          h_mass_v2_fit[i_cen][i_q][i_pt][i_v]->Fit(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v),"L m","",fit_range_low,fit_range_high);

          chi2_ndf_for_q2_v2[i_cen][i_q][i_pt][i_v] = 1.0*fitFcn_v2->GetChisquare()/fitFcn_v2->GetNDF();

          c1->SaveAs(Form("prompt_mass_plot_withchi2_sigma/hmassfit_pt_%s_cen%s_q2bin_%d_v2bin_%d.pdf",pt_name[i_pt].Data(),cen_name[i_cen].Data(),i_q,i_v));

          yield_v2[i_v] = fitFcn_v2->GetParameter(0)*fitFcn_v2->GetParameter(5)/width;
          yield_error_v2[i_v] = fitFcn_v2->GetParError(0)*fitFcn_v2->GetParameter(5)/width;
          sigma_v2[i_cen][i_q][i_pt][i_v] = yield_v2[i_v]/yield_error_v2[i_v];

          if (yield_v2[i_v] <= 0) {
            yield_v2[i_v] = 0; 
            yield_error_v2[i_v] = 0;
            chi2_ndf_for_q2_v2[i_cen][i_q][i_pt][i_v] = -10;
            sigma_v2[i_cen][i_q][i_pt][i_v] = -10;
          }

	  //NEED TO CHANGE
	  v2_x[i_v] = (vnbinning[i_cen][i_pt][i_v]+vnbinning[i_cen][i_pt][i_v+1])/2;
          v2_x_err[i_v] = fabs(vnbinning[i_cen][i_pt][i_v+1]-vnbinning[i_cen][i_pt][i_v])/2;
	  
	  h_v2_hist[i_cen][i_q][i_pt]->SetBinContent(i_v+1,yield_v2[i_v]);
          h_v2_hist[i_cen][i_q][i_pt]->SetBinError(i_v+1,yield_error_v2[i_v]);



	  //=====================
	  // For v3 calculations
	  //=====================

	  if (h_mass_v3_fit[i_cen][i_q][i_pt][i_v]->GetEntries()<10)
	    {
	      yield_v3[i_v] = 0; 
	      yield_error_v3[i_v] = 0;
	      chi2_ndf_for_q3_v3[i_cen][i_q][i_pt][i_v] = -10;
	      sigma_v3[i_cen][i_q][i_pt][i_v] = -10;
	      continue;
	    }


	  TF1* fitFcn_v3 = new TF1(Form("fit_v3_%d_%d_%d", i_pt, i_q, iv),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])) ) + [8] + [9]*x + [10]*x*x + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12])) + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])) )", fit_range_low, fit_range_high);
          fitFcn_v3->SetParameter(0,100);
          fitFcn_v3->FixParameter(1,f->GetParameter(1));
	  fitFcn_v3->FixParameter(2,f->GetParameter(2));
	  fitFcn_v3->FixParameter(3,f->GetParameter(3));
	  fitFcn_v3->FixParameter(4,f->GetParameter(4));
	  fitFcn_v3->FixParameter(5,f->GetParameter(5));
	  fitFcn_v3->FixParameter(6,f->GetParameter(6));
	  fitFcn_v3->FixParameter(7,f->GetParameter(7));
          fitFcn_v3->SetParameter(8,1);
	  fitFcn_v3->SetParameter(9,1);
	  fitFcn_v3->SetParameter(10,1);
	  fitFcn_v3->FixParameter(11,0);
	  

	  h_mass_v3_fit[i_cen][i_q][i_pt][i_v]->Draw("AEP");
          h_mass_v3_fit[i_cen][i_q][i_pt][i_v]->Fit(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v),"M","",fit_range_low,fit_range_high);
          h_mass_v3_fit[i_cen][i_q][i_pt][i_v]->Fit(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v),"L q","",fit_range_low,fit_range_high);
          h_mass_v3_fit[i_cen][i_q][i_pt][i_v]->Fit(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v),"L q","",fit_range_low,fit_range_high);
          h_mass_v3_fit[i_cen][i_q][i_pt][i_v]->Fit(Form("f_ptbin_%d_v2bin_%d",i_pt,i_v),"L m","",fit_range_low,fit_range_high);

          chi2_ndf_for_q3_v3[i_cen][i_q][i_pt][i_v] = 1.0*fitFcn_v3->GetChisquare()/fitFcn_v3->GetNDF();

          c1->SaveAs(Form("prompt_mass_plot_withchi2_sigma/hmassfit_pt_%s_cen%s_q2bin_%d_v3bin_%d.pdf",pt_name[i_pt].Data(),cen_name[i_cen].Data(),i_q,i_v));

          yield_v3[i_v] = fitFcn_v3->GetParameter(0)*fitFcn_v3->GetParameter(5)/width;
          yield_error_v3[i_v] = fitFcn_v3->GetParError(0)*fitFcn_v3->GetParameter(5)/width;
          sigma_v3[i_cen][i_q][i_pt][i_v] = yield_v3[i_v]/yield_error_v3[i_v];

          if (yield_v3[i_v] <= 0) {
            yield_v3[i_v] = 0; 
            yield_error_v3[i_v] = 0;
            chi2_ndf_for_q3_v3[i_cen][i_q][i_pt][i_v] = -10;
            sigma_v3[i_cen][i_q][i_pt][i_v] = -10;
          }

	  //NEED TO CHANGE
	  v3_x[i_v] = (vnbinning[i_cen][i_pt][i_v]+vnbinning[i_cen][i_pt][i_v+1])/2;
          v3_x_err[i_v] = fabs(vnbinning[i_cen][i_pt][i_v+1]-vnbinning[i_cen][i_pt][i_v])/2;
	  
	  h_v3_hist[i_cen][i_q][i_pt]->SetBinContent(i_v+1,yield_v3[i_v]);
          h_v3_hist[i_cen][i_q][i_pt]->SetBinError(i_v+1,yield_error_v3[i_v]);

	}// ---- i_v loop

	yield_v2_graph[i_cen][i_q][i_pt] = new TGraphErrors(N_VBINS,v2_x,yield_v2,v2_x_err,yield_error_v2);
        yield_v2_graph[i_cen][i_q][i_pt]->SetNameTitle("v2_graph_"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q),"v2_graph_"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q));

	yield_v3_graph[i_cen][i_q][i_pt] = new TGraphErrors(N_VBINS,v3_x,yield_v3,v3_x_err,yield_error_v3);
        yield_v3_graph[i_cen][i_q][i_pt]->SetNameTitle("v3_graph_"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q3bin_"+to_string(i_q),"v3_graph_"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q3bin_"+to_string(i_q));
	
      }//-- NQBINS --
    }// -- PTBIN ---
  }//-- CENTBIN --
	  



  /*
	// Now per-v2-bin fits: find all histograms matching this cen/q2/pt with suffix _in_v2bin_
        std::vector<std::tuple<double,double,TH1D*>> v2_hists; // (center,width,hist)
        for (int iv=0; iv<N_VBINS; ++iv) {
          TH1D* hh = h_mass_v2_fit[i_cen][i_q2][i_pt][iv];
          if (!hh) continue;
          double low, high;
          if (parse_bin_edges_from_name(std::string(hh->GetName()), low, high)) {
            double center = (low+high)/2.0;
            double halfwidth = fabs(high-low)/2.0;
            v2_hists.emplace_back(center, halfwidth, hh);
          }
        }
        sort(v2_hists.begin(), v2_hists.end(), [](auto &a, auto &b){ return get<0>(a) < get<0>(b); });
        int iv=0;
        for (auto &t : v2_hists) {
          double center = get<0>(t); double halfw = get<1>(t); TH1D* hh = get<2>(t);
          if (!hh) continue;
          if (hh->GetEntries()<10) { yield_v2[iv]=0; yield_error_v2[iv]=0; iv++; continue; }
          TF1* fitFcn_v2 = new TF1(Form("fit_v2_%d_%d_%d", i_pt, i_q2, iv),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])) ) + [8] + [9]*x + [10]*x*x + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12])) + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])) )", fit_range_low, fit_range_high);
          fitFcn_v2->SetParameter(0,100);
          fitFcn_v2->FixParameter(1,f->GetParameter(1)); fitFcn_v2->FixParameter(2,f->GetParameter(2)); fitFcn_v2->FixParameter(3,f->GetParameter(3)); fitFcn_v2->FixParameter(4,f->GetParameter(4)); fitFcn_v2->FixParameter(5,f->GetParameter(5)); fitFcn_v2->FixParameter(6,f->GetParameter(6)); fitFcn_v2->FixParameter(7,f->GetParameter(7));
          fitFcn_v2->SetParameter(8,1); fitFcn_v2->SetParameter(9,1); fitFcn_v2->SetParameter(10,1); fitFcn_v2->FixParameter(11,0);

	  hh->Fit(fitFcn_v2, "L Q M", "", fit_range_low, fit_range_high);
          double chi2ndf = fitFcn_v2->GetNDF()>0 ? fitFcn_v2->GetChisquare()/fitFcn_v2->GetNDF() : -1;
          yield_v2[iv] = fitFcn_v2->GetParameter(0)*fitFcn_v2->GetParameter(5)/width;
          yield_error_v2[iv] = fitFcn_v2->GetParError(0)*fitFcn_v2->GetParameter(5)/width;
          if (yield_v2[iv] <= 0) { yield_v2[iv]=0; yield_error_v2[iv]=0; }
          v2_x[iv] = center;
	  v2_x_err[iv] = halfw;
          iv++;
          delete fitFcn_v2;
        }
        if (iv>0) yield_v2_graph[i_cen][i_q2][i_pt] = new TGraphErrors(iv, v2_x, yield_v2, v2_x_err, yield_error_v2);

	if (iv>0) {
	  std::vector<double> edges(iv+1);
	  for (int b=0;b<iv; ++b) edges[b] = v2_x[b] - v2_x_err[b];
          edges[iv] = v2_x[iv-1] + v2_x_err[iv-1];
	  h_v2_hist[i_cen][i_q2][i_pt] = new TH1D(Form("h_v2_cen%s_q2bin%d_pt_%s", cen_name[i_cen].c_str(), i_q2, pt_name[i_pt].c_str()), Form("yield vs v2;v2;Yield (entries/bin)"), iv, edges.data());
	  for (int b=0;b<iv; ++b) {
            h_v2_hist[i_cen][i_q2][i_pt]->SetBinContent(b+1, yield_v2[b]);
	    h_v2_hist[i_cen][i_q2][i_pt]->SetBinError(b+1, yield_error_v2[b]);
	  }
        }




	// v3: same approach (respect directory layout)
	std::vector<std::tuple<double,double,TH1D*>> v3_hists;
        for (int iv=0; iv<N_VBINS; ++iv) {
          TH1D* hh = h_mass_v3_fit[i_cen][i_q2][i_pt][iv];
          if (!hh) continue;
          double low, high;
          if (parse_bin_edges_from_name(std::string(hh->GetName()), low, high)) {
            double center = (low+high)/2.0; double halfwidth = fabs(high-low)/2.0;
            v3_hists.emplace_back(center, halfwidth, hh);
          }
        }

        sort(v3_hists.begin(), v3_hists.end(), [](auto &a, auto &b){ return get<0>(a) < get<0>(b); });
        iv=0;
        for (auto &t : v3_hists) {
          double center = get<0>(t); double halfw = get<1>(t); TH1D* hh = get<2>(t);
          if (!hh) continue;
          if (hh->GetEntries()<10) { yield_v3[iv]=0; yield_error_v3[iv]=0; iv++; continue; }
          TF1* fitFcn_v3 = new TF1(Form("fit_v3_%d_%d_%d", i_pt, i_q2, iv),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])) ) + [8] + [9]*x + [10]*x*x + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12])) + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])) )", fit_range_low, fit_range_high);
          fitFcn_v3->SetParameter(0,100);
          fitFcn_v3->FixParameter(1,f->GetParameter(1)); fitFcn_v3->FixParameter(2,f->GetParameter(2)); fitFcn_v3->FixParameter(3,f->GetParameter(3)); fitFcn_v3->FixParameter(4,f->GetParameter(4)); fitFcn_v3->FixParameter(5,f->GetParameter(5)); fitFcn_v3->FixParameter(6,f->GetParameter(6)); fitFcn_v3->FixParameter(7,f->GetParameter(7));
          fitFcn_v3->SetParameter(8,1); fitFcn_v3->SetParameter(9,1); fitFcn_v3->SetParameter(10,1); fitFcn_v3->FixParameter(11,0);

	  hh->Fit(fitFcn_v3, "L Q M", "", fit_range_low, fit_range_high);

	  yield_v3[iv] = fitFcn_v3->GetParameter(0)*fitFcn_v3->GetParameter(5)/width;
          yield_error_v3[iv] = fitFcn_v3->GetParError(0)*fitFcn_v3->GetParameter(5)/width;

	  if (yield_v3[iv] <= 0) { yield_v3[iv]=0; yield_error_v3[iv]=0; }
          v3_x[iv] = center; v3_x_err[iv] = halfw;
          iv++;
          delete fitFcn_v3;
        }
        if (iv>0) yield_v3_graph[i_cen][i_q2][i_pt] = new TGraphErrors(iv, v3_x, yield_v3, v3_x_err, yield_error_v3);

	if (iv>0) {
          std::vector<double> edges3(iv+1);
          for (int b=0;b<iv; ++b) edges3[b] = v3_x[b] - v3_x_err[b];
          edges3[iv] = v3_x[iv-1] + v3_x_err[iv-1];
          h_v3_hist[i_cen][i_q2][i_pt] = new TH1D(Form("h_v3_cen%s_q2bin%d_pt_%s", cen_name[i_cen].c_str(), i_q2, pt_name[i_pt].c_str()),
                                                  Form("yield vs v3;v3;Yield (entries/bin)"), iv, edges3.data());
          for (int b=0;b<iv; ++b) {
            h_v3_hist[i_cen][i_q2][i_pt]->SetBinContent(b+1, yield_v3[b]);
            h_v3_hist[i_cen][i_q2][i_pt]->SetBinError(b+1, yield_error_v3[b]);
          }
        }
	
        delete f;

      }//-- QBIN loop --
    }//--pT loop ---
  }// -- CENT loop --

  */

  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    for (int i_q2=0; i_q2<N_QBINS; i_q2++){
      for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
        pt_x[i_pt]=(pt_edges[i_pt]+pt_edges[i_pt+1])/2;
        pt_x_error[i_pt]=(-1.0*pt_edges[i_pt]+pt_edges[i_pt+1])/2;
	
	mean_val_v2[i_pt] = h_v2_hist[i_cen][i_q2][i_pt]->GetMean();
        mean_error_v2[i_pt] = h_v2_hist[i_cen][i_q2][i_pt]->GetMeanError();

	mean_val_v3[i_pt] = h_v3_hist[i_cen][i_q2][i_pt]->GetMean();
        mean_error_v3[i_pt] = h_v3_hist[i_cen][i_q2][i_pt]->GetMeanError();

	chi2_ndf_for_q2_v2_graph[i_cen][i_q2][i_pt] = new TGraph(N_VBINS,v2_x,chi2_ndf_for_q2_v2[i_cen][i_q2][i_pt]);
        chi2_ndf_for_q2_v2_graph[i_cen][i_q2][i_pt]->SetNameTitle("chi2ndf_graph_pt"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2),"chi2ndf_graph_pt"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2));
        chi2_ndf_for_q2_v2_graph[i_cen][i_q2][i_pt]->Write();
        chi2_ndf_for_q2_v2_graph[i_cen][i_q2][i_pt]->GetXaxis()->SetTitle("v2");
        chi2_ndf_for_q2_v2_graph[i_cen][i_q2][i_pt]->GetYaxis()->SetTitle("chi2/ndf");
        chi2_ndf_for_q2_v2_graph[i_cen][i_q2][i_pt]->Draw("ACP");
        c1->SaveAs(Form("prompt_sigma_chi2_plot/chi2_pt_%s_cen%s_q2bin_%d.pdf",pt_name[i_pt].Data(),cen_name[i_cen].Data(),i_q2));

	chi2_ndf_for_q3_v3_graph[i_cen][i_q2][i_pt] = new TGraph(N_VBINS,v3_x,chi2_ndf_for_q3_v3[i_cen][i_q2][i_pt]);
        chi2_ndf_for_q3_v3_graph[i_cen][i_q2][i_pt]->SetNameTitle("chi2ndf_graph_pt"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2),"chi2ndf_graph_pt"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2));
        chi2_ndf_for_q3_v3_graph[i_cen][i_q2][i_pt]->Write();
        chi2_ndf_for_q3_v3_graph[i_cen][i_q2][i_pt]->GetXaxis()->SetTitle("v2");
        chi2_ndf_for_q3_v3_graph[i_cen][i_q2][i_pt]->GetYaxis()->SetTitle("chi2/ndf");
        chi2_ndf_for_q3_v3_graph[i_cen][i_q2][i_pt]->Draw("ACP");
        c1->SaveAs(Form("prompt_sigma_chi2_plot/chi2_pt_%s_cen%s_q3bin_%d.pdf",pt_name[i_pt].Data(),cen_name[i_cen].Data(),i_q2));

	sigma_v2_fit[i_cen][i_q2][i_pt] = new TGraph(N_VBINS,v2_x,sigma_v2[i_cen][i_q2][i_pt]);
        sigma_v2_fit[i_cen][i_q2][i_pt]->SetNameTitle("sigma_v2_graph_pt"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2),"sigma_v2_graph_pt"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2));
        sigma_v2_fit[i_cen][i_q2][i_pt]->Write();
        sigma_v2_fit[i_cen][i_q2][i_pt]->GetXaxis()->SetTitle("v2");
        sigma_v2_fit[i_cen][i_q2][i_pt]->GetYaxis()->SetTitle("sigma");
        sigma_v2_fit[i_cen][i_q2][i_pt]->Draw("ACP");
        c1->SaveAs(Form("prompt_sigma_chi2_plot/sigma_pt_%s_cen%s_q2bin_%d.pdf",pt_name[i_pt].Data(),cen_name[i_cen].Data(),i_q2));

	sigma_v3_fit[i_cen][i_q2][i_pt] = new TGraph(N_VBINS,v3_x,sigma_v3[i_cen][i_q2][i_pt]);
        sigma_v3_fit[i_cen][i_q2][i_pt]->SetNameTitle("sigma_v2_graph_pt"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2),"sigma_v2_graph_pt"+pt_name[i_pt]+"_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2));
        sigma_v3_fit[i_cen][i_q2][i_pt]->Write();
        sigma_v3_fit[i_cen][i_q2][i_pt]->GetXaxis()->SetTitle("v3");
        sigma_v3_fit[i_cen][i_q2][i_pt]->GetYaxis()->SetTitle("sigma");
        sigma_v3_fit[i_cen][i_q2][i_pt]->Draw("ACP");
        c1->SaveAs(Form("prompt_sigma_chi2_plot/sigma_pt_%s_cen%s_q3bin_%d.pdf",pt_name[i_pt].Data(),cen_name[i_cen].Data(),i_q2));     

      }
      v2_cen[i_cen][i_q2] = new TGraphErrors(N_PTBINS,pt_x,mean_val_v2,pt_x_error,mean_error_v2);
      v2_cen[i_cen][i_q2]->SetNameTitle("v2_graph_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2),"v2_graph_cen"+cen_name[i_cen]+"_q2bin_"+to_string(i_q2));

      v3_cen[i_cen][i_q2] = new TGraphErrors(N_PTBINS,pt_x,mean_val_v3,pt_x_error,mean_error_v3);
      v3_cen[i_cen][i_q2]->SetNameTitle("v3_graph_cen"+cen_name[i_cen]+"_q3bin_"+to_string(i_q2),"v3_graph_cen"+cen_name[i_cen]+"_q3bin_"+to_string(i_q2));


    }
  }

  

  
  // write output ROOT file
  TFile *outf = new TFile("prompt_vn_in_q2HF_total_graph_fullrange_withchi2_sigma_vs_v2v3_from_hists.root","recreate");
  outf->cd();
  for (int i_cen=0; i_cen<N_CENTBINS; i_cen++) {
    for (int i_q=0; i_q<N_QBINS; i_q++){
      for (int i_pt=0; i_pt<N_PTBINS; i_pt++) {
        if (yield_v2_graph[i_cen][i_q][i_pt]) yield_v2_graph[i_cen][i_q][i_pt]->Write();
        if (yield_v3_graph[i_cen][i_q][i_pt]) yield_v3_graph[i_cen][i_q][i_pt]->Write();
	if (h_v2_hist[i_cen][i_q][i_pt])h_v2_hist[i_cen][i_q][i_pt]->Write();
	if (h_v3_hist[i_cen][i_q][i_pt])h_v3_hist[i_cen][i_q][i_pt]->Write();
      }
      if (v2_cen[i_cen][i_q])v2_cen[i_cen][i_q]->Write();
      if (v3_cen[i_cen][i_q]) v3_cen[i_cen][i_q]->Write();
    }
  }





  outf->Write();
  outf->Close();
  massFile->Close();
  return 0;
}












int main(int argc, char *argv[])
{
  if (argc==1) {
    fit_mass_and_flow();
  } else if (argc==2) {
    int idx = atoi(argv[1]);
    fit_mass_and_flow(idx);
  } else if (argc==3) {
    int cen_idx = atoi(argv[1]);
    int pt_idx  = atoi(argv[2]);
    fit_mass_and_flow(cen_idx, pt_idx);
  } else {
    std::cout << "Usage: ./fit_mass_and_flow [target1pct_index] [target_pt_index]" << std::endl;
    return 1;
  }
  return 0;
}
