#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <iomanip>

void AnalyzeD0MC(TString input_txt, TString output_path, int istart, int iend) {
    
    // 1. Initialize Histograms
    const int N_PTBIN = 10;
    Float_t min_pTbin[N_PTBIN] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0};
    Float_t max_pTbin[N_PTBIN] = {2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 40.0};
    TString label_pTbin[N_PTBIN] = {"pT1to2", "pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8", "pT8to10", "pT10to15", "pT15to20", "pT20to40"};

    TH1F *hMass_Signal[N_PTBIN];
    TH1F *hMass_Signal_Swap[N_PTBIN];

    for (int j = 0; j < N_PTBIN; j++) {
        hMass_Signal[j] = new TH1F(Form("hMass_Signal_%s", label_pTbin[j].Data()), "Signal;mass;Entries", 60, 1.7, 2.0);
        hMass_Signal_Swap[j]   = new TH1F(Form("hMass_Signal_Plus_Swap_%s", label_pTbin[j].Data()), "Signal+Swap;mass;Entries", 60, 1.7, 2.0);
    }

    // 2. Open File List
    std::ifstream file_stream(input_txt.Data());
    if (!file_stream.is_open()) {
        std::cerr << "Error: Could not open input list " << input_txt << std::endl;
        return;
    }

    std::string filename;
    int ifile = 0;

    // 3. Loop over the file list
    while (file_stream >> filename) {
        if (ifile < istart) { ifile++; continue; }
        if (ifile >= iend) break;

        TFile *fin = TFile::Open(filename.c_str());
        if (!fin || fin->IsZombie()) {
            std::cout << "Warning: Skipping bad file: " << filename << std::endl;
            if (fin) { fin->Close(); delete fin; }
            ifile++;
            continue;
        }

        std::cout << ">>> Processing ifile=" << ifile << " : " << filename << std::endl;

        // Get Trees and Add Friend
        TTree *tree = (TTree*)fin->Get("d0ana_newreduced/VertexCompositeNtuple"); // Using your previously mentioned path
        //TTree *t_eventinfoana = (TTree*)fin->Get("eventinfoana/EventInfoNtuple");
        
        if (!tree) {
            std::cout << "Skipping: Tree not found in " << filename << std::endl;
            fin->Close(); delete fin;
            ifile++;
            continue;
        }
        //tree->AddFriend(t_eventinfoana);

        // Define Variables for this tree
        const int MAX_CAND = 10000;
        float pT[MAX_CAND], mass[MAX_CAND];
        int candSize;

	bool matchGEN[MAX_CAND], isSwap[MAX_CAND];

	tree->SetBranchAddress("candSize", &candSize);
        tree->SetBranchAddress("pT", pT);
        tree->SetBranchAddress("mass", mass);
        tree->SetBranchAddress("matchGEN", matchGEN);
        tree->SetBranchAddress("isSwap", isSwap);
        

        for (Long64_t i = 0; i < tree->GetEntries(); i++) {
            tree->GetEntry(i);
            for (int k = 0; k < candSize; k++) {
                for (int p = 0; p < N_PTBIN; p++) {
                    if (pT[k] >= min_pTbin[p] && pT[k] < max_pTbin[p]) {
		      if (matchGEN[k] == 1 && isSwap[k] == 0) { hMass_Signal[p]->Fill(mass[k]);}
		      if (matchGEN[k] == 1 || isSwap[k] == 1) {hMass_Signal_Swap[p]->Fill(mass[k]);}
		      break;
		    }
		    
		}
	    }
	}
	

	fin->Close();
	delete fin;
	ifile++;
    }

    // 5. Save Output
    TString outfile = TString::Format("%s/ROOT/D0_Mass_Analysis_%d_%d.root", output_path.Data(), istart, iend);
    TFile *fOut = new TFile(outfile, "RECREATE");
    for (int j = 0; j < N_PTBIN; j++) {
        hMass_Signal[j]->Write();
        hMass_Signal_Swap[j]->Write();
    }
    fOut->Close();
    std::cout << "Successfully saved to " << outfile << std::endl;
}



int main(int argc, char *argv[])
{
    if (argc == 5) // Expecting 4 arguments: input_txt, output_path, istart, iend                                                                                                            

    {
        TString input_txt = argv[1];
        TString output_path = argv[2];
        int istart = std::stoi(argv[3]);
        int iend = std::stoi(argv[4]);

        AnalyzeD0MC(input_txt, output_path, istart, iend);
    }
    else
    {
        std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
        return 1;
    }
    return 0;
}
