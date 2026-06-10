#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TMath.h"
#include "TNtuple.h"

using namespace std;

void make_q2_slices_PbPb2018(TString input_txt, TString output_path, int istart, int iend) {

    TH1::SetDefaultSumw2();

    ifstream file_stream(input_txt.Data());
    TString outfile_ntup = TString::Format("%s/ROOT/Quantiles_ESE_ntuples_%d_%d.root", output_path.Data(), istart, iend);
    TString outfile_hist = TString::Format("%s/ROOT/Quantiles_ESE_hists_%d_%d.root", output_path.Data(), istart, iend);
    TFile *fout_ntup = new TFile(outfile_ntup, "RECREATE");
    TFile *fout_hist = new TFile(outfile_hist, "RECREATE");
    string filename;
    int ifile = 0;

    const int N_CENTBIN = 4;
    Int_t min_centbin[N_CENTBIN] = {0, 10, 30, 50};
    Int_t max_centbin[N_CENTBIN] = {10, 30, 50, 90};
    const int N_CENTBINS_1PC = 90;
    
    // Define your 1% histograms exactly as before
    TH1D *hist_q2_tot[N_CENTBINS_1PC];
    TH1D *hist_q3_tot[N_CENTBINS_1PC];


    for (int j = 0; j < N_CENTBINS_1PC; ++j) {
      hist_q2_tot[j] = new TH1D(Form("hist_q2_tot_cent%i_%i", j, j+1), Form("hist_q2_tot_cent%i_%i", j, j+1), 7000, 0.0, 0.70);
    }


    //TNtuple *nt_ese_global = new TNtuple("nt_global", "nt_global", "cent:q2_hfp:q2_hfm:q2_hf_total");

    TNtuple* nt[N_CENTBINS_1PC];
    for (int i_cen=0; i_cen<N_CENTBINS_1PC; i_cen++)
      {
	nt[i_cen] = new TNtuple(Form("nt_q2_cent%i_%i",i_cen,i_cen+1),Form("nt_q2_cent%i_%i",i_cen,i_cen+1),"q2_hfp:q2_hfm:q2_hf_total: hfp_w: hfm_w: hf_tot_w");
      }
    
    while (file_stream >> filename) {
      if (ifile < istart) {
	ifile++;
	continue;
      }
      if (ifile >= iend) break;
      
      TFile *fin = TFile::Open(filename.c_str());
      if (!fin || fin->IsZombie()) { ifile++; continue; }

	std::cout << ">>> Processing ifile=" << ifile << " : " << filename << std::endl;
	
        TTree *t_eventinfoana = (TTree*)fin->Get("Dfinder/ntDkpi");
        Int_t centrality;
        Float_t ephfmQ[3], ephfpQ[3], ephfmSumW, ephfpSumW;
        t_eventinfoana->SetBranchAddress("centrality", &centrality);
        t_eventinfoana->SetBranchAddress("ephfmQ", ephfmQ);
        t_eventinfoana->SetBranchAddress("ephfpQ", ephfpQ);
        t_eventinfoana->SetBranchAddress("ephfmSumW", &ephfmSumW);
        t_eventinfoana->SetBranchAddress("ephfpSumW", &ephfpSumW);


	Int_t n_entries = t_eventinfoana->GetEntries();
	std::cout<<"Nevents = "<<n_entries<<std::endl;
	
	for (Long64_t ii = 0; ii < n_entries; ii++) {
            t_eventinfoana->GetEntry(ii);

	    if(ii % 100000 == 0)
	      printf("Current entry of the loop = %lld out of %d : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);

	    Int_t cent = centrality / 2;
            if (cent >= 90) continue;

	    Double_t q2_hfp = (ephfpSumW > 0) ? ephfpQ[1]/ephfpSumW : -999;
	    Double_t q2_hfm = (ephfpSumW > 0) ? ephfmQ[1]/ephfmSumW : 999;
	    Double_t q2_total_norm = (ephfmSumW + ephfpSumW > 0) ? (ephfmQ[1] + ephfpQ[1])/(ephfmSumW + ephfpSumW) : -999;
            //Double_t q3_total_norm = (ephfmSumW[2] + ephfpSumW[2] > 0) ? (ephfmQ[2] + ephfpQ[2])/(ephfmSumW[2] + ephfpSumW[2]) : -999;

	    //nt_ese_global->Fill(cent, q2_hfp, q2_hfm, q2_total_norm);
	    
	    
	    if (cent >= 0 && cent < N_CENTBINS_1PC) {
	      hist_q2_tot[cent]->Fill(q2_total_norm);
	      nt[cent]->Fill(q2_hfp,q2_hfm,q2_total_norm, ephfpSumW, ephfmSumW, (ephfmSumW + ephfpSumW));
	    }
        }
        fin->Close();
	delete fin;
        ifile++;
    }



    fout_ntup->cd();
    fout_ntup->mkdir("Q_ntuple");
    fout_ntup->cd("Q_ntuple");
    //nt_ese_global->Write();
    for(int i_cen=0; i_cen<N_CENTBINS_1PC; i_cen++){
      nt[i_cen]->Write();
    }

    fout_hist->cd();
    fout_hist->mkdir("q2_raw_1pc");
    fout_hist->cd("q2_raw_1pc");

    for (int j = 0; j < N_CENTBINS_1PC; ++j){
      if (hist_q2_tot[j]->GetEntries() > 0) hist_q2_tot[j]->Write();
    }

    /*fout_hist->cd();
    fout_hist->mkdir("q3_raw_1pc");
    fout_hist->cd("q3_raw_1pc");
    for (int i = 0; i < N_CENTBIN; ++i){ 
      for (int j = 0; j < N_CENTBIN_1PC; ++j){
	if (hist_q3_tot[i][j]->GetEntries() > 0) hist_q3_tot[i][j]->Write();
      }
      }*/
    fout_ntup->Close();
    fout_hist->Close();
}



int main(int argc, char *argv[])
{
    if (argc == 5) // Expecting 4 arguments: input_txt, output_path, istart, iend
                                                                                                                                                                                             
    {
        TString input_txt = argv[1];
        TString output_path = argv[2];
        int istart = std::stoi(argv[3]);
        int iend = std::stoi(argv[4]);

        make_q2_slices_PbPb2018(input_txt, output_path, istart, iend);
    }
    else
    {
        std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
        return 1;
    }
    return 0;
}
