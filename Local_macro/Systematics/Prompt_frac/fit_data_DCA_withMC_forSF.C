#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TCut.h>
#include <TStyle.h>
#include <TF1.h>

using namespace std;



/*TH1F*h_prompt_MC;
TH1F*h_nonprompt_MC;
TH1F *h_data_DCA;*/


const int bins =11;
Double_t DCA[bins+1]={0.0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.016, 0.022, 0.03, 0.05, 0.12};


TH1D* DCA_fit;
Double_t xMin;
Double_t xMax;
Double_t total_yield;
//int fit_lowbin[N_CENTBIN][N_PTBIN];                                                                                                            
//int fit_highbin[N_CENTBIN][N_PTBIN];                                                                                                          
TAxis* xAxis;
TCanvas *c0;
TPad *pad1;
TPad *pad2;
TLatex * tex_1;
TLegend *leg;
TLegendEntry*entry;
TLine* line;
TF1 *ftot;



Double_t ftotal(Double_t *x, Double_t *par){
	Double_t xx = x[0];
	Int_t bin = h_nonprompt_MC->GetXaxis()->FindBin(xx);
	Double_t br = par[1]*(1-par[0])*h_nonprompt_MC->GetBinContent(bin)/h_nonprompt_MC->Integral();
	Double_t sr = par[1]*par[0]*h_prompt_MC->GetBinContent(bin)/h_prompt_MC->Integral();
	return sr+ br;
}

/*
Double_t funNonPrompt(Double_t *x, Double_t *par){
	Double_t xx = x[0];
	Int_t bin = h_nonprompt_MC->GetXaxis()->FindBin(xx);
	Double_t br = par[1]*(1-par[0])*h_nonprompt_MC->GetBinContent(bin)/h_nonprompt_MC->Integral();
	return br;
}

Double_t funPrompt(Double_t *x, Double_t *par){
	Double_t xx = x[0];
	Int_t bin = h_nonprompt_MC->GetXaxis()->FindBin(xx);
	Double_t sr = par[1]*par[0]*h_prompt_MC->GetBinContent(bin)/h_prompt_MC->Integral();
	return sr;
}
*/

void BinWidthNormalization(TH1* hist) {
    Int_t numBins = hist->GetNbinsX();
    for (Int_t i = 0; i < numBins; ++i) {
        hist->SetBinContent(i+1, hist->GetBinContent(i+1)/hist->GetBinWidth(i+1));
        hist->SetBinError(i+1, hist->GetBinError(i+1)/hist->GetBinWidth(i+1));
    }
}

void fit(){
  const char*cent="cent30to50";
  const char*pt="pt810";
  const char*sf[100] = {};
  Int_t numSF = ;

  TFile *promptMC_DCA = new TFile("Hist_MC_Prompt_Official_DCA_Nov28_combined.root");
  TFile *nonpromptMC_DCA = new TFile("Hist_MC_NonPrompt_Official_DCA_Nov28_combined.root");
  TFile *data_DCA = new TFile("Output_DATA_DCA_allCent_Nov21.root");
  Double_t chi2s[100];
  for(int i=0; i<numSF; i++){
    TH1F*h_prompt_MC = (TH1F*)promptMC_DCA->Get(Form("hdca_%s_%s_%s", cent, pt, sf[i]));
    TH1F*h_nonprompt_MC = (TH1F*)nonpromptMC_DCA->Get(Form("hdca_%s_%s_%s", cent, pt, sf[i]));
    TH1F*h_data_DCA = (TH1F*)data_DCA->Get(Form("hdca_%s_%s", cent, pt));
    
    fit_data_DCA_withMC_forSF(h_prompt_MC, h_nonprompt_MC, h_data_DCA, chi2);
    chi2s[i]=chi2;

    std::cout<<"chi2="chi2s[i]<<std::endl;
  }

}



void fit_data_DCA_withMC_forSF(TH1F*h_prompt_MC, TH1F*h_nonprompt_MC, TH1F*h_data_DCA, Double_t chi2){

  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(true);

  
  bool binWidthNorm = true;
  /*
    TFile *promptMC_DCA = new TFile("Hist_MC_Prompt_Official_DCA_Nov28_combined.root");
    TFile *nonpromptMC_DCA = new TFile("Hist_MC_NonPrompt_Official_DCA_Nov28_combined.root");
    TFile *data_DCA = new TFile("Output_DATA_DCA_allCent_Nov21.root");

  
    h_prompt_MC = (TH1F*) promptMC_DCA->Get(Form("hdca_%s_%s_%s", cent, pt, sf));
    h_nonprompt_MC = (TH1F*) nonpromptMC_DCA->Get(Form("hdca_%s_%s_%s", cent, pt, sf));
    h_data_DCA = (TH1F*) data_DCA->Get(Form("hdca_%s_%s", cent, pt));
  */
  //h_data_DCA->Draw();
  
  h_prompt_MC->GetXaxis()->SetRangeUser(0,0.12);
  h_nonprompt_MC->GetXaxis()->SetRangeUser(0.0,0.12);
  h_data_DCA->GetXaxis()->SetRangeUser(0,0.12);
  

if(binWidthNorm){
  BinWidthNormalization(h_prompt_MC);
  BinWidthNormalization(h_nonprompt_MC);
  BinWidthNormalization(h_data_DCA);
 }

 h_data_DCA->Draw();
 
//TH1F *h_nonprompt_d = (TH1F*)h_nonprompt_MC->Clone("h_nonprompt_d");
//TH1F *h_prompt_d = (TH1F*)h_prompt_MC->Clone("h_prompt_d");
 c0 = new TCanvas("c", "c", 0,53,900,700);
 c0->Range(0,0,1,1);
 c0->SetFillColor(0);
 c0->SetBorderMode(0);
 c0->SetBorderSize(2);
 c0->SetFrameBorderMode(0);

 pad1 = new TPad("pad1", "pad1",0,0.03703704,1,1);
 pad1->Range(-0.003478303,-1125.566,0.03161739,3479.021);

 pad1->SetFillColor(0);
 pad1->SetBorderMode(0);
 pad1->SetBorderSize(2);
 pad1->SetBottomMargin(0.2444444);
 pad1->SetFrameFillStyle(0);
 pad1->SetFrameBorderMode(0);
 pad1->SetFrameFillStyle(0);
 pad1->SetFrameBorderMode(0);
 //pad1->SetLogy(1);
  
 pad2 = new TPad("pad2", "pad2",0.001113586,0.004237288,0.9988864,0.2394068);
 pad2->Range(-0.01244011,-0.5486192,0.1142228,1.520914);
 pad2->SetFillColor(0);
 pad2->SetBorderMode(0);
 pad2->SetBorderSize(2);
 pad2->SetTopMargin(0.01056133);
 pad2->SetBottomMargin(0.3025364);
 pad2->SetFrameFillStyle(0);
 pad2->SetFrameBorderMode(0);
 pad2->SetFrameFillStyle(0);
 pad2->SetFrameBorderMode(0);
  
 pad1->Draw();
 pad1->cd();
 
 
 total_yield=h_data_DCA->Integral();
 
 h_data_DCA->SetFillColor(2);
 h_data_DCA->SetFillStyle(3001);
 h_data_DCA->SetLineColor(1);
 h_data_DCA->SetMarkerStyle(20);
 h_data_DCA->SetMarkerSize(1.3);
 h_data_DCA->GetXaxis()->SetTitle("DCA (cm)");
 //h_data_DCA->GetXaxis()->SetRange(1,9);                                                                        
 h_data_DCA->GetXaxis()->CenterTitle(true);
 h_data_DCA->GetXaxis()->SetLabelFont(42);
 h_data_DCA->GetXaxis()->SetTitleOffset(1.17);
 h_data_DCA->GetXaxis()->SetTitleFont(132);
 h_data_DCA->GetYaxis()->SetTitle("dN/dDCA (cm^{-1})");
 h_data_DCA->GetYaxis()->CenterTitle(true);
 h_data_DCA->GetYaxis()->SetLabelFont(42);
 h_data_DCA->GetYaxis()->SetTitleFont(132);
 h_data_DCA->GetZaxis()->SetLabelFont(42);
 h_data_DCA->GetZaxis()->SetTitleOffset(1);
 h_data_DCA->GetZaxis()->SetTitleFont(42);
 
 xAxis = h_data_DCA->GetXaxis();
 xMin = xAxis->GetXmin();
 xMax = xAxis->GetXmax();
 
 ftot = new TF1("ftot",ftotal,xMin,xMax,2);
 ftot->SetLineColor(kRed);
 ftot->SetLineWidth(5);
 
 ftot->SetParameter(0,0.8);
 ftot->SetParLimits(0,0,1);
 ftot->SetParameter(1,total_yield);

 //ftot->SetParLimits(1,0,2*total_yield);
 ftot->SetFillColor(19);
 ftot->SetFillStyle(0);
 ftot->SetLineColor(2);
 ftot->SetLineWidth(2);
 
 //std::cout<<"UPTO THIS CORRECT"<<std::endl;
 
 h_data_DCA->Fit("ftot", "SR","", xMin, xMax);
 h_data_DCA->Fit("ftot", "SR","", xMin, xMax);
 ftot->ReleaseParameter(0);
 ftot->ReleaseParameter(1);
 h_data_DCA->Fit("ftot", "SR","", xMin, xMax);

 //h_data_DCA->Fit("ftot", "SR");
 //h_data_DCA->Fit("ftot","b");                                                                                    
 
 h_data_DCA->Draw("E1");
 
 h_prompt_MC->SetFillColor(42);
 h_prompt_MC->SetLineColor(4);
 h_prompt_MC->SetMarkerColor(42);
 h_prompt_MC->SetMarkerStyle(6);
 h_prompt_MC->SetMarkerSize(1.6);
 //h_prompt_MC->GetXaxis()->SetRange(1,9);                                                                         
 h_prompt_MC->GetXaxis()->SetLabelFont(42);
 h_prompt_MC->GetXaxis()->SetTitleOffset(1);
 h_prompt_MC->GetXaxis()->SetTitleFont(42);
 h_prompt_MC->GetYaxis()->SetLabelFont(42);
 h_prompt_MC->GetYaxis()->SetTitleFont(42);
 h_prompt_MC->GetZaxis()->SetLabelFont(42);
 h_prompt_MC->GetZaxis()->SetTitleOffset(1);
 h_prompt_MC->GetZaxis()->SetTitleFont(42);
 h_prompt_MC->Scale(1.0/h_prompt_MC->Integral()*ftot->GetParameter(0)*ftot->GetParameter(1));
 
 h_prompt_MC->Draw("SAME HIST");
 
 h_nonprompt_MC->SetFillColor(51);
 h_nonprompt_MC->SetFillStyle(3001);
 // h_nonprompt_MC->SetLineColor("#000099");                                                                       
 h_nonprompt_MC->SetMarkerColor(4);
 h_nonprompt_MC->SetMarkerSize(1.4);
 //h_nonprompt_MC->GetXaxis()->SetRange(1,9);                                                                      
 h_nonprompt_MC->GetXaxis()->SetLabelFont(42);
 h_nonprompt_MC->GetXaxis()->SetTitleOffset(1);
 h_nonprompt_MC->GetXaxis()->SetTitleFont(42);
 h_nonprompt_MC->GetYaxis()->SetLabelFont(42);
 h_nonprompt_MC->GetYaxis()->SetTitleFont(42);
 h_nonprompt_MC->GetZaxis()->SetLabelFont(42);
 h_nonprompt_MC->GetZaxis()->SetTitleOffset(1);
 h_nonprompt_MC->GetZaxis()->SetTitleFont(42);
 
 h_nonprompt_MC->Scale(1.0/h_nonprompt_MC->Integral()*(1.0-ftot->GetParameter(0))*ftot->GetParameter(1));
 
 h_nonprompt_MC->Draw("SAME HIST");
 
 cout<<"Prompt fraction is "<<ftot->GetParameter(0)<<"  +/-  "<<ftot->GetParError(0)<<endl;
 cout<<"Non-Prompt fraction is "<<(1-ftot->GetParameter(0))<<endl;
 cout<<"chi2/Ndof = "<<double(ftot->GetChisquare()/ftot->GetNDF())<<endl;

 
 chi2=(ftot->GetChisquare()/ftot->GetNDF());
 
 tex_1 = new TLatex(0.5556793,0.8692308,Form("Prompt fraction =%1.4f #pm %1.4f", ftot->GetParameter(0), ftot->GetParError(0)));
 tex_1->SetNDC();
 tex_1->SetTextAlign(32);
 tex_1->SetTextFont(42);
 tex_1->SetTextSize(0.03);
 tex_1->SetLineWidth(4);
 tex_1->Draw();
 tex_1 = new TLatex(0.4142539,0.8338462, Form("chi2/Ndof =%1.2f ", (ftot->GetChisquare()/ftot->GetNDF())));
 tex_1->SetNDC();
 tex_1->SetTextAlign(32);
 tex_1->SetTextFont(42);
 tex_1->SetTextSize(0.03);
 tex_1->SetLineWidth(4);
 tex_1->Draw();
 
 
 leg = new TLegend(0.6002227,0.6948148,0.889755,0.8918519,NULL,"brNDC");
 leg->SetBorderSize(0);
 leg->SetLineColor(1);
 leg->SetLineStyle(1);
 leg->SetLineWidth(1);
 leg->SetFillColor(0);
 leg->SetFillStyle(1001);
      
 entry=leg->AddEntry(h_data_DCA,"Data","lpf");
 entry->SetFillColor(2);
 entry->SetFillStyle(3001);
 entry->SetLineColor(2);
 entry->SetLineStyle(1);
 entry->SetLineWidth(1);
 entry->SetMarkerColor(1);
 entry->SetMarkerStyle(20);
 entry->SetMarkerSize(1.3);
 entry->SetTextFont(42);
 entry=leg->AddEntry(h_prompt_MC,"MC Prompt Ds^{#pm}","lpf");
 entry->SetFillColor(42);
 entry->SetFillStyle(1001);
 entry->SetLineColor(4);
 entry->SetLineStyle(1);
 entry->SetLineWidth(1);
 entry->SetMarkerColor(42);
 entry->SetMarkerStyle(6);
 entry->SetMarkerSize(1.6);
 entry->SetTextFont(42);
 entry=leg->AddEntry(h_nonprompt_MC,"MC NonPrompt Ds^{#pm}","lpf");
 entry->SetFillColor(51);
 entry->SetFillStyle(3001);
 // entry->SetLineColor("#000099");                                                                                
 entry->SetLineStyle(1);
 entry->SetLineWidth(1);
 entry->SetMarkerColor(4);
 entry->SetMarkerStyle(1);
 entry->SetMarkerSize(1.4);
 entry->SetTextFont(42);
 entry=leg->AddEntry(ftot,"Fit","l");
 leg->Draw();
 
 pad2->Draw();
 pad2->cd();
 
 std::pair<double, double>dca_fit(0.0, 0.12);                                                                                
 
 DCA_fit = new TH1D("DCA_fit","",bins, DCA);
 Double_t fit_lowbin = DCA_fit->GetXaxis()->FindBin(dca_fit.first);
 Double_t fit_highbin = DCA_fit->GetXaxis()->FindBin(dca_fit.second);
 
 for(int k =fit_lowbin; k<=fit_highbin; k++)
   {
     double bin_upedge = h_data_DCA->GetXaxis()->GetBinUpEdge(k);
     double bin_lowedge = h_data_DCA->GetXaxis()->GetBinLowEdge(k);
     
     double bin_center = h_data_DCA->GetXaxis()->GetBinCenter(k);	\
     
     double fits_bincontent = ftot->Eval(bin_center);
     double data_bincontent = h_data_DCA->GetBinContent(k);
     double data_binerror = h_data_DCA->GetBinError(k);
     
     if(data_binerror == 0) continue;
     
     double ratio = data_bincontent/fits_bincontent;
     double ratio_error = data_binerror/fits_bincontent;
     
     DCA_fit->SetBinContent(k+1, ratio);
     DCA_fit->SetBinError(k+1, ratio_error);
     
   }

 DCA_fit->Draw("E1P");
 DCA_fit->SetMarkerStyle(20);
 DCA_fit->SetMarkerSize(1.2);
 DCA_fit->SetLineColor(1);
 DCA_fit->GetXaxis()->SetTitle("DCA (cm)");
 DCA_fit->GetXaxis()->CenterTitle(true);
 DCA_fit->GetXaxis()->SetLabelFont(132);
 DCA_fit->GetXaxis()->SetLabelSize(0.15);
 DCA_fit->GetXaxis()->SetTitleSize(0.14);
 DCA_fit->GetXaxis()->SetTitleOffset(1);
 DCA_fit->GetXaxis()->SetTitleFont(132);
 DCA_fit->GetYaxis()->SetTitle("Data/Fit");
 DCA_fit->GetYaxis()->CenterTitle(true);
 DCA_fit->GetYaxis()->SetNdivisions(408);
 DCA_fit->GetYaxis()->SetLabelFont(132);
 DCA_fit->GetYaxis()->SetLabelSize(0.1);
 DCA_fit->GetYaxis()->SetTitleSize(0.15);
 DCA_fit->GetYaxis()->SetTitleOffset(0.23);
 DCA_fit->GetYaxis()->SetTitleFont(132);
 
 line = new TLine(DCA_fit->GetXaxis()->GetXmin(), 1.0, DCA_fit->GetXaxis()->GetXmax(), 1.0);
 line->SetLineColor(kRed);
 line->SetLineStyle(2);
 line->Draw();
 
 
 
 
 /*
   TF1 *ftot = new TF1("ftot",ftotal,0,0.07,2);
   double total_yield=h_data_DCA->Integral();
   cout<<total_yield<<endl;
   ftot->SetParameter(0,0.7);
   ftot->SetParameter(1,total_yield);
   ftot->SetParLimits(0,0,1);
   ftot->SetParLimits(1,0,2*total_yield);
   
   TCanvas *c0 = new TCanvas("c0","fit",900,700);
   
   h_data_DCA->Fit("ftot","b");
   
   cout<<"Prompt fraction is "<<ftot->GetParameter(0)<<endl;
   cout<<"chi2/Ndof = "<<int(ftot->GetChisquare()/ftot->GetNDF())<<endl;
 */
 
}


