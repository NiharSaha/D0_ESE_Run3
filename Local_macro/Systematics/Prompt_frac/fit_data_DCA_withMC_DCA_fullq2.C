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
#include <TLatex.h>
#include <TGaxis.h>
#include <TLegend.h>
TH1D *prompt_whole_10_20;
TH1D *nonprompt_whole_10_20;

Double_t ftotal(Double_t *x, Double_t *par){
  Double_t xx = x[0];
  Int_t bin = nonprompt_whole_10_20->GetXaxis()->FindBin(xx);
  Double_t br = par[1]*(1-par[0])*nonprompt_whole_10_20->GetBinContent(bin)/(1.0*nonprompt_whole_10_20->Integral(0,-1));
  Double_t sr = par[1]*par[0]*prompt_whole_10_20->GetBinContent(bin)/(1.0*prompt_whole_10_20->Integral(0,-1));
  return sr + br;
}
/*
Double_t funNonPrompt(Double_t *x, Double_t *par){
  Double_t xx = x[0];
  Int_t bin = nonprompt_whole_10_20->GetXaxis()->FindBin(xx);
  Double_t br = par[1]*(1-par[0])*nonprompt_whole_10_20->GetBinContent(bin)/nonprompt_whole_10_20->Integral();
  return br;
}

Double_t funPrompt(Double_t *x, Double_t *par){
  Double_t xx = x[0];
  Int_t bin = prompt_whole_10_20->GetXaxis()->FindBin(xx);
  Double_t sr = par[1]*par[0]*prompt_whole_10_20->GetBinContent(bin)/prompt_whole_10_20->Integral();
  return sr;
}
*/

void fit_data_DCA_withMC_DCA_fullq2(float sc, string cent_name, string pt_name)
{
  string inf_MC_prompt = "MC_DCA_template/MC_prompt_cen"+cent_name+"_dca_hist_with_res_"+to_string(sc)+".root";
  TFile *MC_prompt_DCA = new TFile(inf_MC_prompt.c_str());
  string inf_MC_nonprompt = "MC_DCA_template/MC_nonprompt_cen"+cent_name+"_dca_hist_with_res_"+to_string(sc)+".root" ;
  TFile *MC_nonprompt_DCA = new TFile(inf_MC_nonprompt.c_str());
  //TFile *data_DCA = new TFile("dca_hist_output_pol2_rebin_allq2.root");
  //TFile *data_DCA = new TFile("chi2_minimization_for_each_pT_cen/dca_hist_output_pol3.root");
  TFile *data_DCA = new TFile("data_signal_count_dca.root");

  TCanvas *c1 = new TCanvas("c1","c1",600,800);
  c1->SetLeftMargin(0.155);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gPad->SetLogx();
  //gPad->SetLogy();
  TPad *pad1 = new TPad("pad1", "Top pad",0.0,0.25,1.0,1.0);
  TPad *pad2 = new TPad("pad2", "Bottom pad",0.0,0.,1.0,0.25);
  pad1->SetBottomMargin(0.02);
  pad1->SetLogy();
  //pad1->SetLogx();
  //pad1->Range(2.059422,275625,2.464241,319375);
  pad1->SetLeftMargin(0.17);
  pad1->SetRightMargin(0.02);
  pad2->SetRightMargin(0.02);
  pad1->Draw();
  pad1->cd();

  string MC_prompt_name = "hist_MC_prompt_dca_cen_"+cent_name+"_pt_"+pt_name;
  prompt_whole_10_20 = (TH1D*) MC_prompt_DCA->Get(MC_prompt_name.c_str())->Clone("prompt_whole_10_20");
  string MC_nonprompt_name = "hist_MC_nonprompt_dca_cen_"+cent_name+"_pt_"+pt_name;
  nonprompt_whole_10_20 = (TH1D*) MC_nonprompt_DCA->Get(MC_nonprompt_name.c_str())->Clone("nonprompt_whole_10_20");
  string data_name = "hist_dca_cen_"+cent_name+"_pt_"+pt_name;
  TH1D *h_data_DCA = (TH1D*) data_DCA->Get(data_name.c_str())->Clone("h_data_dca");
  //TF1 *ftot = new TF1("ftot",ftotal,0,0.7,2);

  h_data_DCA->Sumw2();
  prompt_whole_10_20->Sumw2();
  nonprompt_whole_10_20->Sumw2();

  //Add the Edges of the bins here
  //
  const Int_t NBINS = 13;
  double edges[NBINS +1], bin_center[NBINS];

  
  double pull_val_denominator[NBINS];

  for (int k=0; k<NBINS; k++)
  {
    pull_val_denominator[k] = sqrt((h_data_DCA->GetBinError(k+1))*(h_data_DCA->GetBinError(k+1))+(prompt_whole_10_20->GetBinError(k+1))*(prompt_whole_10_20->GetBinError(k+1))+(nonprompt_whole_10_20->GetBinError(k+1))*(nonprompt_whole_10_20->GetBinError(k+1)));
    cout << "Pull error for bin "<<k<<": "<<pull_val_denominator[k]<<endl;
  }


  for(int i_dca=0; i_dca<NBINS+1; i_dca++)
  {
    edges[i_dca]=h_data_DCA->GetBinLowEdge(i_dca+1);
    //cout << "Edge values "<<i_dca<<": "<< edges[i_dca] << endl;
    if (i_dca==NBINS) break;
    bin_center[i_dca] = h_data_DCA->GetBinCenter(i_dca+1);
    //cout << "Bin center "<<i_dca<<": "<<bin_center[i_dca] << endl;
  }
  //edges[0]=0.0003; //Using this for log scale x-axis
  //edges[0]=0.000;

  auto h_data_DCA_empty = new TH1D ("h_data_DCA_empty","h_data_DCA_empty",NBINS,edges);
  //h_data_DCA_empty->SetBinContent(1,h_data_DCA->GetBinContent(1));

  TH1F *h_nonprompt_d = (TH1F*)nonprompt_whole_10_20->Clone("h_nonprompt_d");
  TH1F *h_prompt_d = (TH1F*)prompt_whole_10_20->Clone("h_prompt_d");
  h_nonprompt_d->Sumw2();
  h_prompt_d->Sumw2();

  TF1 *ftot = new TF1("ftot",ftotal,0,0.135,2);

  double total_yield=h_data_DCA->Integral();
  cout<<total_yield<<endl;
  ftot->SetParameter(0,0.7);
  ftot->SetParameter(1,total_yield);
  ftot->SetParLimits(0,0,1);
  ftot->SetParLimits(1,0,2*total_yield);
  //h_data_DCA->Fit("ftot","b","",0,0.135675);
  //h_data_DCA->Fit("ftot","b L","",0,0.135675);
  
  h_data_DCA->Fit("ftot","E R 0","",0,0.08);
  h_data_DCA->Fit("ftot","E R 0","",0,0.08);
  h_data_DCA->Fit("ftot","E R b","",0,0.08);
  

  double errpar0 = ftot->GetParError(0);
  double errpar1 = ftot->GetParError(1);

  double par0 = ftot->GetParameter(0);
  double par1 = ftot->GetParameter(1);

  double ymax = h_data_DCA->GetBinContent(h_data_DCA->GetMaximumBin());

  /*ofstream outfrac;
    outfrac.open ("DCA_fit_merge/chi2_ndf_resolution.txt",fstream::app);
    outfrac << sc*100 << "	" << (ftot->GetChisquare()/ftot->GetNDF()) << "    " << par0 << " +- " << errpar0 << "    " << par1 << " +- " << errpar1 << endl;*/

  cout<<"prompt ratio " <<ftot->GetParameter(0)<<endl;
  cout<<"chi2  =  "<<(ftot->GetChisquare())<< "  ";
  cout<<"ndf = "<<ftot->GetNDF()<< "  ";
  cout<<"chi2/ndf  =  "<<(ftot->GetChisquare()/ftot->GetNDF())<<endl;
   

  TGaxis::SetMaxDigits(3);
  h_data_DCA_empty->SetMarkerColor(1);
  h_data_DCA_empty->SetMarkerStyle(20);
  h_data_DCA_empty->SetMarkerSize(0);
  h_data_DCA->SetMarkerColor(1);
  h_data_DCA->SetMarkerStyle(20);
  h_data_DCA->SetMarkerSize(0.4);
  h_data_DCA_empty->GetXaxis()->CenterTitle();
  h_data_DCA_empty->GetYaxis()->CenterTitle();
  h_data_DCA_empty->GetXaxis()->SetTitleOffset(1.0);
  h_data_DCA_empty->GetYaxis()->SetTitleOffset(1.4);
  h_data_DCA_empty->GetXaxis()->SetLabelOffset(0.007);
  h_data_DCA_empty->GetYaxis()->SetLabelOffset(0.007);
  h_data_DCA_empty->GetXaxis()->SetTitleSize(0.045);
  h_data_DCA_empty->GetYaxis()->SetTitleSize(0.045);
  h_data_DCA_empty->GetXaxis()->SetTitleFont(42);
  h_data_DCA_empty->GetYaxis()->SetTitleFont(42);
  h_data_DCA_empty->GetXaxis()->SetLabelFont(42);
  h_data_DCA_empty->GetYaxis()->SetLabelFont(42);
  h_data_DCA_empty->GetXaxis()->SetLabelSize(0.00);
  h_data_DCA_empty->GetYaxis()->SetLabelSize(0.04);
  h_data_DCA_empty->GetXaxis()->SetTitle("DCA (cm)");
  h_data_DCA_empty->GetYaxis()->SetTitle("dN / d(DCA) (cm^{-1})");
  h_data_DCA_empty->SetAxisRange(10,10*ymax,"Y");
  //h_data_DCA_empty->SetAxisRange(0.0003,edges[NBINS-1],"X");
  //h_data_DCA->SetAxisRange(0.0000001,0.3);

  h_data_DCA_empty->Draw();
  h_data_DCA->Draw("esame");

  //h_data_DCA->Fit("ftot","b");
  //h_data_DCA->Fit("ftot","b L");
  //ftot->Draw("same");
  cout<<ftot->GetParameter(0)<<endl;
  TH1F *h_PNP_ratio =  new TH1F ("h_PNP_ratio","h_PNP_ratio",10,0,1);
  h_PNP_ratio->Sumw2();
  h_PNP_ratio->SetBinContent(1,ftot->GetParameter(0));

  //const Int_t NBINS = 9;
  //double edges[NBINS +1]={0,0.0016,0.0032,0.0048,0.0064,0.008,0.0096,0.0128,0.016,0.16};

  /*const Int_t NBINS = 8;
  double edges[NBINS +1], bin_center[NBINS];

  for(int i_dca=0; i_dca<NBINS+1; i_dca++)
  {
    edges[i_dca]=h_data_DCA->GetBinLowEdge(i_dca+1);
    cout << "Edge values "<<i_dca<<": "<< edges[i_dca] << endl; 
    if (i_dca==NBINS) break;
    bin_center[i_dca] = h_data_DCA->GetBinCenter(i_dca+1);
    cout << "Bin center "<<i_dca<<": "<<bin_center[i_dca] << endl;
  }

  //TH1 *h_data_dca = new TH1D("h_data_dca","h_data_dca",NBINS,edges);

  //TH1 *h_tot = new TH1D("h_tot","h_tot",NBINS,edges);
  //TH1 *h_np = new TH1D("h_np","h_np",NBINS,edges);

    h_tot = ftot->DoCreateHistogram(0,0.13567500,true);
    h_tot->SetLineColor(2);
    h_tot->SetFillStyle(3001);
    h_tot->SetFillColor(2);
    h_tot->Draw("histsame");

    for (int i_bin=0; i_bin<=h_tot->GetNbinsX(); i_bin++)
    {
    cout << "h_tot bin edge: "<<h_tot->GetBinLowEdge(i_bin+1);
    }

    TF1* fNP = new TF1("fNP",&funNonPrompt, 0, 0.13567500, 2);
    fNP->SetParameters(ftot->GetParameter(0),ftot->GetParameter(1));

    h_np = fNP->DoCreateHistogram(0,0.13567500,false);
    h_np->SetLineColor(4);
    h_np->SetFillStyle(3001);
    h_np->SetFillColor(4);
  //h_np->Draw("same");

    int h_np_bins = h_np->GetNbinsX();
    for (int i_bin=0; i_bin<=h_np_bins; i_bin++)
    {
    cout << "h_np bin edge: "<<h_np->GetBinLowEdge(i_bin+1) << endl;
    }

    TF1* fP = new TF1("fP",&funPrompt, 0, 0.13567500, 2);
    fP->SetParameters(ftot->GetParameter(0),ftot->GetParameter(1));
    fP->SetLineColor(5);
    fP->SetFillStyle(1001);
    fP->SetFillColorAlpha(5,0.01);
  //fP->Draw("same");*/

  /*nonprompt_whole_10_20->Scale((1-ftot->GetParameter(0))*ftot->GetParameter(1)/nonprompt_whole_10_20->Integral());
  prompt_whole_10_20->Scale((ftot->GetParameter(0))*ftot->GetParameter(1)/prompt_whole_10_20->Integral());
  prompt_whole_10_20->Add(nonprompt_whole_10_20);

  prompt_whole_10_20->SetLineColor(2);
  prompt_whole_10_20->SetFillStyle(3001);
  prompt_whole_10_20->SetFillColor(2);
  prompt_whole_10_20->Draw("histsame");

  nonprompt_whole_10_20->SetLineColor(4);
  nonprompt_whole_10_20->SetFillStyle(3001);
  nonprompt_whole_10_20->SetFillColor(4);
  nonprompt_whole_10_20->Draw("histsame");*/

  h_nonprompt_d->Scale((1-ftot->GetParameter(0))*ftot->GetParameter(1)/h_nonprompt_d->Integral());
  h_prompt_d->Scale((ftot->GetParameter(0))*ftot->GetParameter(1)/h_prompt_d->Integral());
  h_prompt_d->Add(h_nonprompt_d);
  h_prompt_d->SetMarkerStyle(20);
  h_nonprompt_d->SetMarkerStyle(25);

  h_prompt_d->SetLineColor(2);
  h_prompt_d->SetFillStyle(3004);
  h_prompt_d->SetFillColor(2);
  h_prompt_d->SetMarkerSize(0);
  h_prompt_d->Draw("ehistsame");

  h_nonprompt_d->SetLineColor(4);
  h_nonprompt_d->SetFillStyle(3005);
  h_nonprompt_d->SetFillColor(4);
  h_nonprompt_d->SetMarkerSize(0);
  h_nonprompt_d->Draw("ehistsame");

  TLatex* tex;
  tex = new TLatex(0.20,0.86,Form("p_{T} %s cen %s ",pt_name.c_str(),cent_name.c_str()));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.03);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.715,0.86,Form("resolution scale %0.3f",sc));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.03);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.65,0.81,Form("Prompt fraction %0.2f #pm %0.2f",ftot->GetParameter(0),ftot->GetParError(0)));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.03);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.20,0.81,Form("#chi^{2}/NDF: %f",1.0*ftot->GetChisquare()/ftot->GetNDF()));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.03);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.22,0.75,"|y| < 1");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  //tex->Draw();

  TLatex Tl; 
  Tl.SetNDC();
  Tl.SetTextAlign(12);
  Tl.SetTextSize(0.05);
  Tl.SetTextFont(42);
  //Tl.DrawLatex(0.23,0.92, "#font[61]{CMS}");
  //Tl.DrawLatex(0.56,0.92, "#scale[0.8]{pp 252 nb^{-1} (5.02 TeV)}");

  TLatex Tl2;
  Tl2.SetNDC();
  Tl2.SetTextAlign(12);
  Tl2.SetTextSize(0.05*0.75);
  Tl2.SetTextFont(42);
  //Tl2.DrawLatex(0.35,0.92, "#font[52]{Preliminary}");

  /*auto *leg1 = new TLegend(0.3,0.5,0.59,0.65);
    leg1->SetTextSize(0.037);
    leg1->AddEntry("h_data_DCA","Data","pe");
    leg1->AddEntry("h_np","nonprompt component","f");
    leg1->AddEntry("h_tot","total","f");
    leg1->SetTextFont(42);
    leg1->SetBorderSize(0);
    leg1->Draw();*/

  double pull_val[NBINS];
  TH1D *pull__1 = new TH1D("pull__1","pull",NBINS,edges);
  pull__1 = (TH1D*)h_data_DCA_empty->Clone("pull__1");

  for (int k=0; k<NBINS; k++)
  {
    //pull__1->SetBinContent(0,k+1);
    //pull_val[k]=(ds_fy3001[k]-graph->Eval(ds_fx3001[k]))/ds_fely3001[k];
    //pull_val[k] = ((h_data_DCA->GetBinContent(k+1))-ftot->Eval(h_data_DCA->GetBinCenter(k+1)))/sqrt((h_data_DCA->GetBinError(k+1))*(h_data_DCA->GetBinError(k+1))+(prompt_whole_10_20->GetBinError(k+1))*(prompt_whole_10_20->GetBinError(k+1))+(nonprompt_whole_10_20->GetBinError(k+1))*(nonprompt_whole_10_20->GetBinError(k+1)));
    pull_val[k] = (h_data_DCA->GetBinContent(k+1)-ftot->Eval(h_data_DCA->GetBinCenter(k+1))) / sqrt(h_data_DCA->GetBinError(k+1)*h_data_DCA->GetBinError(k+1) + h_prompt_d->GetBinError(k+1)*h_prompt_d->GetBinError(k+1));
    //pull_val[k] = (h_data_DCA->GetBinContent(k+1)-ftot->Eval(h_data_DCA->GetBinCenter(k+1))) / pull_val_denominator[k];
    cout << "Pull values denominator for bin "<<k<<": "<<sqrt(h_data_DCA->GetBinError(k+1)*h_data_DCA->GetBinError(k+1) + h_prompt_d->GetBinError(k+1)*h_prompt_d->GetBinError(k+1))<<endl;
    cout << "Pull values for bin "<<k<<": "<<pull_val[k]<<endl;

  }

  TGraph *graph_pull = new TGraph(NBINS,bin_center,pull_val);
  Double_t pull_max=3.5;

  pull__1->SetMinimum(-1.0*pull_max);
  pull__1->SetMaximum(1.0*pull_max);
  //pull__1->SetEntries(1);
  pull__1->SetDirectory(0);
  pull__1->SetStats(0);
  pull__1->GetYaxis()->SetLabelSize(0.00);
  pull__1->GetYaxis()->SetTickLength(0.00);
  pull__1->GetYaxis()->SetTitleSize(0.00);

  //pull__1->GetXaxis()->SetRangeUser(0.0000,0.13567500);
  pull__1->GetXaxis()->SetTitle("DCA (cm)");
  pull__1->GetXaxis()->CenterTitle(true);
  //pull__1->GetXaxis()->SetLabelFont(42);
  pull__1->GetXaxis()->SetLabelSize(0.108);
  pull__1->GetXaxis()->SetTitleSize(0.15);
  pull__1->GetXaxis()->SetTitleOffset(1.02);
  pull__1->GetXaxis()->SetTitleFont(42);


  c1->cd(0);
  //pad2->RangeAxis(0.000001,-10.0,0.136,10.0);
  //pad2->SetLogx();
  pad2->SetGridy();
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.45);
  pad2->SetLeftMargin(0.17);
  pad2->Draw();
  pad2->cd();

  //TLine *line = new TLine(0.0003,0,edges[NBINS-1],0);
  //line->SetLineColor(kRed);
  pull__1->Draw();
  //h_data_DCA->Draw("FUNC");
  //pull__1->GetXaxis()->SetRangeUser(0.000001,0.136);
  graph_pull->SetMarkerStyle(20);
  graph_pull->SetMarkerSize(0.4);
  graph_pull->Draw("P");
  //graph_pull->GetXaxis()->SetRangeUser(0.000001,0.136);
  //line->Draw("same");
  //TGaxis *axis = new TGaxis(0.0003,-7,0.0003,7,-7,7,505,""); //Using this for log scale x-axis
  TGaxis *axis = new TGaxis(0.000,-1.0*pull_max,0.000,1.0*pull_max,-1.0*pull_max,1.0*pull_max,505,"");
  axis->SetLabelSize(0.115);

  //axis->SetTitle("pull test");
  axis->SetTitle("Pull");
  axis->CenterTitle(true);
  axis->SetTitleOffset(0.41);
  axis->SetTitleSize(0.15);
  axis->SetTitleFont(42);
  axis->SetLabelFont(42);
  axis->Draw();

  c1->SaveAs(Form("DCA_fit_Sept_3_2025_nolikelihood/DCA_fit_cen_%s_pt_%s_sc_%0.3f.pdf",cent_name.c_str(),pt_name.c_str(),sc));

  fstream f(Form("DCA_fit_Sept_3_2025_nolikelihood/chi2_ndf_cen_%s_pt_%s.txt",cent_name.c_str(),pt_name.c_str()), f.out | f.app);
  f << endl << sc << "  " <<  (ftot->GetChisquare()/ftot->GetNDF()) << "  " << ftot->GetParameter(0) << "  " << ftot->GetParError(0);


  /*TFile * result = new TFile(Form("DCA_fit_merge/%s_%0.2f.root",output.Data(),sc),"RECREATE");
    h_PNP_ratio->Write();
    result->Close();*/

}

