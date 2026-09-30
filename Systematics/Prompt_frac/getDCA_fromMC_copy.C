#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <cmath>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TMath.h"
#include "TVector3.h"

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/Prompt_frac/BDTHandler.h"
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/Prompt_frac/BDTHandler.cc"

using namespace std;
using namespace std::chrono;

void getDCA_fromMC(TString input_txt, TString output_file, int istart, int iend, TString sample_type)
{
    TH1::SetDefaultSumw2();
    auto start = high_resolution_clock::now();

    if (!sample_type.EqualTo("Prompt") && !sample_type.EqualTo("NonPrompt"))
    {
        cout << "Error: sample_type must be 'Prompt' or 'NonPrompt', got '"
             << sample_type << "'" << endl;
        return;
    }

    BDTHandler bdtHandler;

    const double MAX_PT_ANA = 100.0;
    const double MIN_PT_ANA = 2.0;
    const double MAX_Y_ANA = 1.0;

    // const int N_CENTBIN = 5;
    // Double_t centbinning[N_CENTBIN + 1] = {0, 10, 20, 30, 40, 50};
    // const char *label_centbin[N_CENTBIN] = {"cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"};

    const int nCentApprox = 5;
    double centEdges[nCentApprox + 1] = {0, 10, 20, 30, 40, 50};
    double centApprox[nCentApprox] = {10.0, 30.0, 50, 70.0, 90.0}; // 2xcent

    const int N_PTBIN = 9;
    Double_t ptbinning[N_PTBIN + 1] = {2, 3, 4, 5, 6, 8, 10, 15, 30, 100};
    const char *label_pTbin[N_PTBIN] = {"pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8", "pT8to10", "pT10to15", "pT15to30", "pT30to100"};

    const Int_t dca_bins = 18;
    Double_t DCA[dca_bins + 1] = {0,0.0005,0.0011,0.0014,0.0024,0.0029,0.0039,0.0045,0.0059,0.0067,0.0085,0.01179,0.016,0.0214,0.028,0.0366,0.0475,0.079,0.135};

    // Scale factors 0.00 .. 1.50 in steps of 0.01, generated in a loop so there are
    // never missing points. numSF = 151, SF = 1.0 at index 100.
    const Int_t numSF = 151;
    Double_t ScaleFactor[numSF];
    TString sf_label[numSF];
    for (int i = 0; i < numSF; i++)
    {
        ScaleFactor[i] = i / 100.0;           // i is the value in hundredths (0..150)
        int whole = i / 100;
        int frac = i % 100;
        sf_label[i] = (frac % 10 == 0) ? Form("x%dp%d", whole, frac / 10)   // 0.50 -> x0p5
                                       : Form("x%dp%02d", whole, frac);      // 0.51 -> x0p51
    }

    TH3D *h_binning = new TH3D("centbins", "", nCentApprox, centEdges, N_PTBIN, ptbinning, dca_bins, DCA);

    // Create a histogram for every combination of Cent, pT, and SF
    TH1F *h_DCA[nCentApprox][N_PTBIN][numSF];

    for (int icase = 0; icase < nCentApprox; icase++)
    {
        for (int ipt = 0; ipt < N_PTBIN; ipt++)
        {
            for (int isf = 0; isf < numSF; isf++)
            {
                TString hname = Form("hdca_cent%dto%d_%s_%s",
                                     (int)centEdges[icase],
                                     (int)centEdges[icase + 1],
                                     label_pTbin[ipt],
                                     sf_label[isf].Data());
                h_DCA[icase][ipt][isf] = new TH1F(hname, hname, dca_bins, DCA);
                h_DCA[icase][ipt][isf]->SetDirectory(0);
            }
        }
    }

    // Gen 3D decay-length distributions (SF-independent), one per Cent x pT bin.
    // Filled with the same dl3d_gen used to build the DCA templates and written
    // to a dedicated "DecayLength" directory in the output file.
    const Int_t dl_bins = 250;
    const Double_t DL_MIN = 0.0, DL_MAX = 0.5;
    TH1F *h_DL[nCentApprox][N_PTBIN];
    for (int icase = 0; icase < nCentApprox; icase++)
    {
        for (int ipt = 0; ipt < N_PTBIN; ipt++)
        {
            TString hname = Form("hdl_cent%dto%d_%s_%s",
                                 (int)centEdges[icase],
                                 (int)centEdges[icase + 1],
                                 label_pTbin[ipt],
                                 sample_type.Data());
            h_DL[icase][ipt] = new TH1F(hname, hname, dl_bins, DL_MIN, DL_MAX);
            h_DL[icase][ipt]->SetDirectory(0);
        }
    }

    TH1D *hcent = new TH1D(Form("hcent_%s", sample_type.Data()),
                           Form("Centrality %s;Centrality (%%);Events", sample_type.Data()),
                           100, 0.0, 100.0);

    Long64_t nTotalMatched = 0;

    ifstream file_stream(input_txt.Data());
    if (!file_stream.is_open())
    {
        cout << "Error: Could not open input list " << input_txt << endl;
        return;
    }

    string filename;
    int ifile = 0;

    while (file_stream >> filename)
    {
        if (ifile < istart)
        {
            ifile++;
            continue;
        }
        if (ifile >= iend)
            break;

        TFile *fin = TFile::Open(filename.c_str());
        if (!fin || fin->IsZombie())
        {
            cout << "Warning: Skipping bad file: " << filename << endl;
            if (fin)
            {
                fin->Close();
                delete fin;
            }
            ifile++;
            continue;
        }

        cout << ">>> Processing [" << sample_type << "] ifile=" << ifile << " : " << filename << endl;

        TTree *tree = (TTree *)fin->Get("d0Analyzer/VCNtuple_D02kpi");
        TTree *t_eventinfoana = (TTree *)fin->Get("eventinfoana/EventInfoNtuple");

        if (!tree || !t_eventinfoana)
        {
            cout << "Skipping: Tree(s) missing in " << filename << endl;
            fin->Close();
            delete fin;
            ifile++;
            continue;
        }

        tree->AddFriend(t_eventinfoana);
        

        Int_t candSize = 0;
        Int_t centrality = 0;
        

        std::vector<Bool_t> *matchGEN = nullptr;
        std::vector<Bool_t> *isSwap = nullptr;
        std::vector<Float_t> *pT = nullptr;
        std::vector<Float_t> *mass = nullptr;
        std::vector<Float_t> *y = nullptr;
        std::vector<Float_t> *mva = nullptr;
        std::vector<Float_t> *dca = nullptr;
        std::vector<Float_t> *pT_gen = nullptr;
        std::vector<Float_t> *eta_gen = nullptr;
        std::vector<Float_t> *phi_gen = nullptr;
        std::vector<Float_t> *y_gen = nullptr;

        
        std::vector<bool> *Dgen_isPrompt = nullptr; 
        std::vector<Float_t> *gen3DDecayLength = nullptr;
        std::vector<Float_t> *gen3DPointingAngle = nullptr;

        std::vector<Float_t> *DgenprodvtxX = nullptr;
        std::vector<Float_t> *DgenprodvtxY = nullptr;
        std::vector<Float_t> *DgenprodvtxZ = nullptr;
        std::vector<Float_t> *DgendecayvtxX = nullptr;
        std::vector<Float_t> *DgendecayvtxY = nullptr;
        std::vector<Float_t> *DgendecayvtxZ = nullptr;
        std::vector<Float_t> *Dgen_bAncestorVtxX = nullptr;
        std::vector<Float_t> *Dgen_bAncestorVtxY = nullptr;
        std::vector<Float_t> *Dgen_bAncestorVtxZ = nullptr;


        tree->SetBranchStatus("*", 0);
        for (const auto &p : {"candSize", "centrality", "y", "pT", "mass", "mva", "ip3d", 
                              "matchGEN", "isSwap", 
                              "pT_gen", "eta_gen", "phi_gen", "y_gen",
                              "Dgen_isPrompt", "gen3DDecayLength", "gen3DPointingAngle",
                              "DgenprodvtxX", "DgenprodvtxY", "DgenprodvtxZ",
                              "DgendecayvtxX", "DgendecayvtxY", "DgendecayvtxZ",
                              "Dgen_bAncestorVtxX", "Dgen_bAncestorVtxY", "Dgen_bAncestorVtxZ"})
            tree->SetBranchStatus(p, 1);

        tree->SetBranchAddress("candSize", &candSize);
        tree->SetBranchAddress("centrality", &centrality);
        tree->SetBranchAddress("y", &y);    
        tree->SetBranchAddress("pT", &pT);
        tree->SetBranchAddress("mass", &mass);
        tree->SetBranchAddress("mva", &mva);
        tree->SetBranchAddress("ip3d", &dca);
        tree->SetBranchAddress("pT_gen", &pT_gen);
        tree->SetBranchAddress("eta_gen", &eta_gen);
        tree->SetBranchAddress("phi_gen", &phi_gen);
        tree->SetBranchAddress("y_gen", &y_gen);
        tree->SetBranchAddress("DgenprodvtxX", &DgenprodvtxX);
        tree->SetBranchAddress("DgenprodvtxY", &DgenprodvtxY);
        tree->SetBranchAddress("DgenprodvtxZ", &DgenprodvtxZ);
        tree->SetBranchAddress("DgendecayvtxX", &DgendecayvtxX);
        tree->SetBranchAddress("DgendecayvtxY", &DgendecayvtxY);
        tree->SetBranchAddress("DgendecayvtxZ", &DgendecayvtxZ);
        tree->SetBranchAddress("Dgen_bAncestorVtxX", &Dgen_bAncestorVtxX);
        tree->SetBranchAddress("Dgen_bAncestorVtxY", &Dgen_bAncestorVtxY);
        tree->SetBranchAddress("Dgen_bAncestorVtxZ", &Dgen_bAncestorVtxZ);
        tree->SetBranchAddress("matchGEN", &matchGEN);
        tree->SetBranchAddress("isSwap", &isSwap);
        tree->SetBranchAddress("Dgen_isPrompt", &Dgen_isPrompt);
        tree->SetBranchAddress("gen3DDecayLength", &gen3DDecayLength);
        tree->SetBranchAddress("gen3DPointingAngle", &gen3DPointingAngle);

        Int_t nevent = tree->GetEntries();

        for (int ievt = 0; ievt < nevent; ievt++)
        {
            tree->GetEntry(ievt);

            if (centrality < 0 || centrality >= 200)
                continue;

            hcent->Fill(centrality / 2.0);

            if (ievt % 50000 == 0)
            {
                printf("  [%s] entry = %d / %d (%.2f%%)\n",
                       sample_type.Data(), ievt, nevent, (Double_t)ievt / nevent * 100.0);
            }

            // int nCands = std::min(candSize, Max_cand);

            for (int icand = 0; icand < candSize; icand++)
            {

                if (matchGEN->at(icand) != 1 || isSwap->at(icand) != 0) continue;

                float Pt = pT->at(icand);
                float Y = y->at(icand);
                float Mva = mva->at(icand);
                
                
                
                if (Pt >= MAX_PT_ANA || Pt < MIN_PT_ANA || fabs(Y) >= MAX_Y_ANA) continue;
		        

                float dca3d_reco = dca->at(icand);
                bool isPrompt = Dgen_isPrompt->at(icand);
                
                if (sample_type.EqualTo("Prompt")    && !isPrompt) continue;
                if (sample_type.EqualTo("NonPrompt") &&  isPrompt) continue;

                nTotalMatched++;

                // ---- generator DCA definition (flip this flag & recompile;
                //      submit one definition at a time) ----------------------
                //   true  : from ntuple branches,
                //           gen3DDecayLength * sin(gen3DPointingAngle)   (both origin-referenced)
                //   false : self-computed,
                //           |D0 decay vertex - refVtx| * sin(pointing angle), where the
                //           pointing angle is recomputed from the gen D0 momentum w.r.t. the
                //           same refVtx; refVtx = Dgenprodvtx (prompt) / Dgen_bAncestorVtx (non-prompt)
                const bool useBranchDca = false;

                float dl3d_gen = 0.f;
                float pa3d_gen = 0.f;

                if (useBranchDca)
                {
                    dl3d_gen = gen3DDecayLength->at(icand);
                    pa3d_gen = gen3DPointingAngle->at(icand);
                }
                else
                {
                    // flight vector: b-ancestor vertex -> D0 decay vertex
                    TVector3 flight(DgendecayvtxX->at(icand) - Dgen_bAncestorVtxX->at(icand),
                                    DgendecayvtxY->at(icand) - Dgen_bAncestorVtxY->at(icand),
                                    DgendecayvtxZ->at(icand) - Dgen_bAncestorVtxZ->at(icand));

                    // gen D0 momentum from the matched gen kinematics
                    TVector3 p_gen;
                    p_gen.SetPtEtaPhi(pT_gen->at(icand), eta_gen->at(icand), phi_gen->at(icand));

                    dl3d_gen = flight.Mag();

                    // pointing angle: angle between gen D0 momentum and flight vector
                    if (flight.Mag() > 0. && p_gen.Mag() > 0.)
                        pa3d_gen = p_gen.Angle(flight);
                }

                float dca3d_gen = dl3d_gen * sin(pa3d_gen);


                for (int icase = 0; icase < nCentApprox; icase++)
                {
                    int cent_case = (int)centApprox[icase];

                    double bdt_cut = bdtHandler.getBDTCut(Y, cent_case, Pt);
                    if (Mva <= bdt_cut)
                        continue;

                    Int_t i_pTbin = h_binning->GetYaxis()->FindBin(Pt) - 1;
                    if (i_pTbin < 0 || i_pTbin >= N_PTBIN)
                        continue;


                    if (dl3d_gen >= DL_MIN && dl3d_gen < DL_MAX)
                        h_DL[icase][i_pTbin]->Fill(dl3d_gen);

                    for (int isf = 0; isf < numSF; isf++)
                    { 
                        float scaled_Dca = dca3d_gen + ScaleFactor[isf] * (dca3d_reco - dca3d_gen);

                        if (scaled_Dca >= DCA[0] && scaled_Dca < DCA[dca_bins])
                        {
                            h_DCA[icase][i_pTbin][isf]->Fill(scaled_Dca);
                        }
                    }
                }
            } // -- Cand Loop --
        } // -- Evt loop --

        fin->Close();
        delete fin;
        ifile++;
    }

    cout << "[" << sample_type << "] Total matched candidates in acceptance: "
         << nTotalMatched << endl;

    // persist the count alongside histograms for later cross-checks
    TH1D *hSummary = new TH1D(Form("hMatchedSummary_%s", sample_type.Data()), "", 1, 0, 1);
    hSummary->SetBinContent(1, nTotalMatched);
    hSummary->GetXaxis()->SetBinLabel(1, "TotalMatched");

    // Save histograms to output file
    TFile *fout = new TFile(output_file, "RECREATE");
    fout->cd();

    hcent->Write();
    hSummary->Write();
    
    // Gen decay-length distributions in their own directory (raw counts, Sumw2 errors)
    TDirectory *dir_dl = fout->mkdir("DecayLength");
    dir_dl->cd();
    for (int icase = 0; icase < nCentApprox; icase++)
    {
        for (int ipt = 0; ipt < N_PTBIN; ipt++)
        {
            h_DL[icase][ipt]->Write();
        }
    }

    fout->cd();
    for (int icase = 0; icase < nCentApprox; icase++)
    {
        for (int ipt = 0; ipt < N_PTBIN; ipt++)
        {
            for (int isf = 0; isf < numSF; isf++)
            {
                TH1F *hist = h_DCA[icase][ipt][isf];

                // --- APPLY BIN-WIDTH NORMALIZATION ---
                Int_t numBins = hist->GetNbinsX();
                for (Int_t i = 1; i <= numBins; ++i)
                {
                    Double_t binW = hist->GetBinWidth(i);
                    hist->SetBinContent(i, hist->GetBinContent(i) / binW);
                    hist->SetBinError(i, hist->GetBinError(i) / binW);
                }
                // -------------------------------------

                hist->Write();
            }
        }
    }

    
   

    fout->Close();
    delete fout;
    delete h_binning;

    auto stop = high_resolution_clock::now();
    cout << "\nFinished " << sample_type << " in "
         << duration_cast<minutes>(stop - start).count() << " minutes." << endl;
}

int main(int argc, char **argv)
{
    if (argc == 6)
    {
        TString input_txt = argv[1];
        TString output_file = argv[2];
        int istart = std::stoi(argv[3]);
        int iend = std::stoi(argv[4]);
        TString sample_type = argv[5];

        getDCA_fromMC(input_txt, output_file, istart, iend, sample_type);
    }
    else
    {
        std::cout << "Usage: ./getDCA_MC <input_list.txt> <output.root> <istart> <iend> <sample_type (Prompt/NonPrompt)>\n";
        return 1;
    }
    return 0;
}
