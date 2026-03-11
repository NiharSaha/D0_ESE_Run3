#include<iostream>
#include <cstdlib>
#include <fstream>
#include <map>
#include <string>
#include <chrono>
#include <ctime>
#include <cstdlib>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFileCollection.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TComplex.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


  
using namespace std;
using namespace TMVA;
using namespace std::chrono;

const int N_CENTBIN =3;
const char* label_cent[N_CENTBIN]={"cent010","cent1030","cent3050",};
Int_t min_cent[N_CENTBIN]= {0, 20, 60};
Int_t max_cent[N_CENTBIN]= {20, 60, 100};
Double_t centbinning[N_CENTBIN+1]= {0., 20., 60., 100.};


const int N_PTBIN = 10;  
const char* label_pTbin[N_PTBIN]={"pT12", "pT23", "pT34", "pT45", "pT56", "pT68", "pT810", "pT1015", "pT1520", "pT2030"};
Float_t min_pTbin[N_PTBIN]={1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0};
Float_t max_pTbin[N_PTBIN]={2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0};
Double_t ptbinning[N_PTBIN+1]={1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0};


Double_t bdtcut_0_10[N_PTBIN] = {0.29,0.37,0.38,0.32,0.25,0.18,0.15,0.10,0.00,0.00};
Double_t bdtcut_10_30[N_PTBIN] = {0.29,0.36,0.36,0.30,0.20,0.15,0.10,0.06,-0.05,-0.06};
Double_t bdtcut_30_50[N_PTBIN] = {0.29,0.34,0.33,0.24,0.15,0.10,0.05,0.00,-0.05,-0.05};


/*Double_t GetBDTCutPrompt(Double_t pt, Int_t centrality) {
  Double_t bdtcut = -1.;
  Double_t bdtcut_0_10[N_PTBIN] = {0.29,0.37,0.38,0.32,0.25,0.18,0.15,0.10,0.00,0.00};
  Double_t bdtcut_10_30[N_PTBIN] = {0.29,0.36,0.36,0.30,0.20,0.15,0.10,0.06,-0.05,-0.06};
  Double_t bdtcut_30_50[N_PTBIN] = {0.29,0.34,0.33,0.24,0.15,0.10,0.05,0.00,-0.05,-0.05};
  Double_t ptBDTbin[N_PTBIN+1] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0};

  if (centrality < 10) {
    for (int i=0; i<N_PTBIN; i++) if (ptBDTbin[i]<pt && ptBDTbin[i+1]>pt) bdtcut = bdtcut_0_10[i];
  }
  if (centrality >= 10 && centrality < 30) {
    for (int i=0; i<N_PTBIN; i++) if (ptBDTbin[i]<pt && ptBDTbin[i+1]>pt) bdtcut = bdtcut_10_30[i];
  }
  if (centrality >= 30 && centrality < 50) {
    for (int i=0; i<N_PTBIN; i++) if (ptBDTbin[i]<pt && ptBDTbin[i+1]>pt) bdtcut = bdtcut_30_50[i];
  }
  return bdtcut;
  
}
*/

TH1D * h_dmass[N_CENTBIN][N_PTBIN];
TH2D *h_binning = new TH2D("h_binning","h_binning", N_CENTBIN, centbinning, N_PTBIN, ptbinning);





void Dmass(TString input_txt, TString output_path, int istart, int iend){

  TH1::StatOverflows(kTRUE);
  TH2::StatOverflows(kTRUE);
  TH3::StatOverflows(kTRUE);
  TH1::SetDefaultSumw2(); //Added newly                                                                                                     
  TH2::SetDefaultSumw2(); //Added newly
  TH3::SetDefaultSumw2(); //Added newly
  

  //ifstream file_stream("/home/saha115/Purdue/CMSSW_13_2_11/src/PbPb_files_new.txt"); 
  //TString outfile = TString::Format("/scratch/bell/saha115/D0_RAA_new/ROOT/Dmass_file_%d_%d.root", istart, iend);
  //ifstream file_stream = TString::Format("%s", input_txt); 
  ifstream file_stream(input_txt.Data());
  TString outfile = TString::Format("%s/ROOT/Dmass_file_%d_%d.root",output_path.Data(), istart, iend);

  string filename;


  //TFile* fout = new TFile("out.root", "recreate");


  
  //if(istest)
  //{
  //  outfile = "test.root";
  //}
  
  int ifile=0;
  int counts = 0;

  for(int icent=1; icent<2; icent++){  
    for(int ipt=0; ipt<N_PTBIN; ipt++){  
      h_dmass[icent][ipt]=new TH1D(Form("Dmass_%s_%s",label_cent[icent], label_pTbin[ipt]), Form("Dmass_%s_%s",label_cent[icent], label_pTbin[ipt]), 50, 1.75, 2.0);
      
    }
  }


  
  TH1D * hist_cent = new TH1D("Centrality", "Centrality", 200, 0, 200);



  //const float max_entries = (10 * tree->GetEntries()) / 5;

  
  while(true) {
    //if (counts > max_entries) break;

    //if (counts > 2) break;
    file_stream >> filename;
    if(file_stream.eof()) break;
    if (ifile < istart) {
      ifile++;
      continue;
    }
    if (ifile >= iend) break;
    
  
  TFile *fin = TFile::Open(filename.c_str());
  if(fin->IsZombie() || fin->Get("mvaTree")==nullptr) {
     ifile ++;
     continue;
   }	
  cout<<"ifile="<<ifile<<endl;

  TTree* tree = (TTree*)fin->Get("mvaTree");

  
  const Int_t MAX=20000;

  Int_t centrality;
  Int_t candSize;
  
   
  Float_t mass[MAX], pT[MAX], alpha[MAX], y[MAX], BDT_weight[MAX]; 
   
   
  tree->SetBranchAddress("centrality", &centrality);
  tree->SetBranchAddress("candSize", &candSize);
  tree->SetBranchAddress("mass", &mass);
  //tree->SetBranchAddress("Dgen", &Dgen);
  tree->SetBranchAddress("pT", &pT);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("alpha", &alpha);
  tree->SetBranchAddress("BDT_weight", &BDT_weight);
  
  
  tree->SetBranchStatus("*", 0);
  
  for (const auto& p : {"centrality", "candSize", "mass", "pT", "y", "alpha" , "BDT_weight"})
    tree->SetBranchStatus(p, 1);
  
   
   Int_t nevent =tree->GetEntries();


   
   for (int ievt =0; ievt < nevent; ievt++)
     {
       tree->GetEntry(ievt);


       if(ievt % 500000 == 0)
	 printf("current entry of event loop = %d out of %d : %.3f %%\n", ievt, nevent, (Double_t)ievt / nevent * 100);
       
       hist_cent->Fill(centrality);
       
       //if (centrality >= max_cent[N_CENTBIN-1]) continue;
       if ( centrality < 20 || centrality >= 60) continue;
       
       Int_t i_cent = (h_binning->GetXaxis()->FindBin(centrality))-1;
       //cout<<"CandSize="<<candSize<<endl;       

       for (int j =0; j < candSize; j++)
	 {

	   if (pT[j] > max_pTbin[N_PTBIN-1] || pT[j] < min_pTbin[0]) continue;

	   if (y[j]>1.0) continue;
	   	   


	   //std::cout<<"pt="<<pT[j]<<" "<<"cent="<<centrality<<" "<<"BDT="<<GetBDTCutPrompt(pT[j],centrality)<<std::endl;

	   //if(BDT_weight[j] <=GetBDTCutPrompt(pT[j],centrality)) continue;


	   Int_t i_pt=(h_binning->GetYaxis()->FindBin(pT[j]))-1;

	   if(BDT_weight[j] <= bdtcut_10_30[i_pt]) continue;
	   

	   h_dmass[i_cent][i_pt]->Fill(mass[j]); //for hist
	   

	 }//---particle loop--
       
     }//---event loop---


   ifile++;
   fin->Close();
   //delete tree;

  }

  //cout<<"UPTO THIS IS OKAY"<<endl;          
  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();  
  hist_cent->Write();
   for(int icent=1; icent<2; icent++){
     for(int ipt=0; ipt<N_PTBIN; ipt++){
       h_dmass[icent][ipt]->Write();  
     }
   }

   fout->Write();
   fout->Close();
   
    cout << "DONE" << endl;   
    //auto stop = high_resolution_clock::now();
    //auto duration = duration_cast<minutes>(stop - start);
    
    //cout << "Total time taken: "<< duration.count() << "minutes" << endl;
   
  }





int main(int argc, char *argv[])
{
    if (argc == 5) // Expecting 4 arguments: input_txt, output_path, istart, iend
    {
        TString input_txt = argv[1];
        TString output_path = argv[2];
        int istart = std::stoi(argv[3]);
        int iend = std::stoi(argv[4]);

        Dmass(input_txt, output_path, istart, iend);
    }
    else
    {
        std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
        return 1;
    }
    return 0;
}




