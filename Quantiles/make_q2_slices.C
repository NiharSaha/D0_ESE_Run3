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

void make_q2_slices(TString input_txt, TString output_path, int istart, int iend) {

    TH1::SetDefaultSumw2();

    ifstream file_stream(input_txt.Data());
    TString outfile = TString::Format("%s/ROOT/Quantiles_ESE_%d_%d.root", output_path.Data(), istart, iend);
    TFile *fout = new TFile(outfile, "RECREATE");
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
      hist_q3_tot[j] = new TH1D(Form("hist_q3_tot_cent%i_%i", j, j+1), Form("hist_q3_tot_cent%i_%i", j, j+1), 7000, 0.0, 0.70);
    }


    //TNtuple *nt_ese_global = new TNtuple("nt_global", "nt_global", "centrality:q2_hfp:q2_hfm:q3_hfp:q3_hfm:q2_total:q3_total");

    TNtuple* nt[N_CENTBINS_1PC];
    for (int i_cen=0; i_cen<N_CENTBINS_1PC; i_cen++)
      {
        nt[i_cen] = new TNtuple(Form("nt_qn_cent%i_%i",i_cen,i_cen+1),Form("nt_qn_cent%i_%i",i_cen,i_cen+1),"q2_hfp:q2_hfm:q3_hfp:q3_hfm:q2_hf_total:q3_hf_total:q2_hfp_w:q2_hfm_w:q3_hfp_w:q3_hfm_w:q2_hf_total_w:q3_hf_total_w");
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
	
        TTree *t_eventinfoana = (TTree*)fin->Get("eventinfoana/EventInfoNtuple");
        Short_t centrality;
        Float_t ephfmQ[3], ephfpQ[3], ephfmSumW[3], ephfpSumW[3];
        t_eventinfoana->SetBranchAddress("centrality", &centrality);
        t_eventinfoana->SetBranchAddress("ephfmQ", ephfmQ);
        t_eventinfoana->SetBranchAddress("ephfpQ", ephfpQ);
        t_eventinfoana->SetBranchAddress("ephfmSumW", ephfmSumW);
        t_eventinfoana->SetBranchAddress("ephfpSumW", ephfpSumW);


	Int_t n_entries = t_eventinfoana->GetEntries();
	std::cout<<"Nevents = "<<n_entries<<std::endl;
	
	for (Long64_t ii = 0; ii < n_entries; ii++) {
            t_eventinfoana->GetEntry(ii);

	    if(ii % 10000 == 0)
	      printf("Current entry of the loop = %lld out of %d : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);

	    Int_t cent = centrality / 2;
            if (cent >= 90) continue;

	    Double_t q2_hfp = (ephfpSumW[1] > 0) ? ephfpQ[1]/ephfpSumW[1] : -999;
            Double_t q2_hfm = (ephfmSumW[1] > 0) ? ephfmQ[1]/ephfmSumW[1] : -999;
	    Double_t q3_hfp = (ephfpSumW[2] > 0) ? ephfpQ[2]/ephfpSumW[2] : -999;
            Double_t q3_hfm = (ephfmSumW[2] > 0) ? ephfmQ[2]/ephfmSumW[2] : -999;
            Double_t q2_total_norm = (ephfmSumW[1] + ephfpSumW[1] > 0) ? (ephfmQ[1] + ephfpQ[1])/(ephfmSumW[1] + ephfpSumW[1]) : -999;
            Double_t q3_total_norm = (ephfmSumW[2] + ephfpSumW[2] > 0) ? (ephfmQ[2] + ephfpQ[2])/(ephfmSumW[2] + ephfpSumW[2]) : -999;

	    //nt_ese_global->Fill(cent, q2_hfp, q2_hfm, q3_hfp, q3_hfm, q2_total_norm, q3_total_norm);

	    
	    if (cent>=0 && cent<N_CENTBINS_1PC){
	      hist_q2_tot[cent]->Fill(q2_total_norm);
	      hist_q3_tot[cent]->Fill(q3_total_norm);
	      nt[cent]->Fill(q2_hfp,q2_hfm,q3_hfp, q3_hfm, q2_total_norm, q3_total_norm, ephfpSumW[1], ephfmSumW[1], ephfpSumW[2], ephfmSumW[2], (ephfmSumW[1] + ephfpSumW[1]), (ephfmSumW[2] + ephfpSumW[2]));
	    }
            
	    
        }// -- evt loop --
        fin->Close();
	delete fin;
        ifile++;
    }




    fout->cd();
    fout->mkdir("Q_ntuple");
    fout->cd("Q_ntuple");
    //nt_ese_global->Write();
    for(int i_cen=0; i_cen<N_CENTBINS_1PC; i_cen++){
      nt[i_cen]->Write();
    }

    fout->cd();
    fout->mkdir("q2_raw_1pc");
    fout->cd("q2_raw_1pc");

    for (int j = 0; j < N_CENTBINS_1PC; ++j){
      if (hist_q2_tot[j]->GetEntries() > 0) hist_q2_tot[j]->Write();
    }

    fout->cd();
    fout->mkdir("q3_raw_1pc");
    fout->cd("q3_raw_1pc");

    for (int j = 0; j < N_CENTBINS_1PC; ++j){
      if (hist_q3_tot[j]->GetEntries() > 0) hist_q3_tot[j]->Write();
    }

    fout->Close();
}



int main(int argc, char *argv[])
{
    if (argc == 5) // Expecting 4 arguments: input_txt, output_path, istart, iend
                                                                                                                                                                                             
    {
        TString input_txt = argv[1];
        TString output_path = argv[2];
        int istart = std::stoi(argv[3]);
        int iend = std::stoi(argv[4]);

        make_q2_slices(input_txt, output_path, istart, iend);
    }
    else
    {
        std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
        return 1;
    }
    return 0;
}
