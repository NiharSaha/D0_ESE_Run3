#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include <iostream>

void computeEfficiency(TString input_file = "Eff_out_combined_Sept24.root",
		       TString output_file = "Eff_factor_output_Sept24.root")
{

  TH1::StatOverflows(kTRUE);
  TH2::StatOverflows(kTRUE);
  TH3::StatOverflows(kTRUE);
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  

  TFile *fin = TFile::Open(input_file);
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: Cannot open input file: " << input_file << std::endl;
        return;
    }

    // --- Load 2D histograms ---
    TH2D *h_centptreco       = (TH2D *)fin->Get("centptreco");
    TH2D *h_centptreco_ptw   = (TH2D *)fin->Get("centptreco_ptw");
    TH2D *h_centptreco_woBDT = (TH2D *)fin->Get("centptreco_woBDT");
    TH2D *h_centptgen        = (TH2D *)fin->Get("centptgen");


    

    // --- Load 1D histograms per centrality bin ---
    /*TH1D *h_ptgen_cent010  = (TH1D *)fin->Get("ptgen_cent010");
    TH1D *h_ptgen_cent1020 = (TH1D *)fin->Get("ptgen_cent1020");
    TH1D *h_ptgen_cent2030 = (TH1D *)fin->Get("ptgen_cent2030");
    TH1D *h_ptgen_cent3040 = (TH1D *)fin->Get("ptgen_cent3040");
    TH1D *h_ptgen_cent4050 = (TH1D *)fin->Get("ptgen_cent4050");

    TH1D *h_ptreco_cent010  = (TH1D *)fin->Get("ptreco_cent010");
    TH1D *h_ptreco_cent1020 = (TH1D *)fin->Get("ptreco_cent1020");
    TH1D *h_ptreco_cent2030 = (TH1D *)fin->Get("ptreco_cent2030");
    TH1D *h_ptreco_cent3040 = (TH1D *)fin->Get("ptreco_cent3040");
    TH1D *h_ptreco_cent4050 = (TH1D *)fin->Get("ptreco_cent4050");

    TH1D *h_ptreco_cent010_woBDT  = (TH1D *)fin->Get("ptreco_cent010_woBDT");
    TH1D *h_ptreco_cent1020_woBDT = (TH1D *)fin->Get("ptreco_cent1020_woBDT");
    TH1D *h_ptreco_cent2030_woBDT = (TH1D *)fin->Get("ptreco_cent2030_woBDT");
    TH1D *h_ptreco_cent3040_woBDT = (TH1D *)fin->Get("ptreco_cent3040_woBDT");
    TH1D *h_ptreco_cent4050_woBDT = (TH1D *)fin->Get("ptreco_cent4050_woBDT");*/

    // Basic null check
    if (!h_centptreco || !h_centptgen || !h_centptreco_woBDT) {
        std::cerr << "Error: Missing 2D histograms in input file." << std::endl;
        fin->Close();
        return;
    }

    TFile *fout = new TFile(output_file, "RECREATE");
    fout->cd();

    // --------------------------------------------------
    // 2D efficiency: with BDT
    // --------------------------------------------------
    TH2D *h_centpteff = (TH2D *)h_centptreco->Clone("h_centpteff");
    h_centpteff->SetTitle("Reco/Gen efficiency;cent;pT");
    h_centpteff->Divide(h_centptreco, h_centptgen);
    h_centpteff->Write();

    // --------------------------------------------------
    // 2D efficiency: without BDT
    // --------------------------------------------------
    TH2D *h_centpteff_woBDT = (TH2D *)h_centptreco_woBDT->Clone("h_centpteff_woBDT");
    h_centpteff_woBDT->SetTitle("Reco/Gen efficiency (no BDT);cent;pT");
    h_centpteff_woBDT->Divide(h_centptreco_woBDT, h_centptgen);
    h_centpteff_woBDT->Write();

    // --------------------------------------------------
    // 2D efficiency: with BDT, pT-weighted reco / unweighted gen
    // (weighting gen too would cancel the weight bin by bin)
    // --------------------------------------------------
    if (h_centptreco_ptw) {
        TH2D *h_centpteff_ptw = (TH2D *)h_centptreco_ptw->Clone("h_centpteff_ptw");
        h_centpteff_ptw->SetTitle("Reco(pT-weighted)/Gen efficiency;cent;pT");
        h_centpteff_ptw->Divide(h_centptreco_ptw, h_centptgen);
        h_centpteff_ptw->Write();
    } else {
        std::cerr << "Warning: Missing centptreco_ptw, skipping pT-weighted 2D efficiency." << std::endl;
    }

    // --------------------------------------------------
    // 1D efficiency vs pT, integrated over centrality
    // NOTE: project reco and gen separately, then divide.
    // ProjectionY() of an efficiency histogram sums efficiencies over cent bins (not < 1).
    // --------------------------------------------------
    auto makePtEff = [&](TH2D *hreco2D, TH2D *hgen2D, const char *name, const char *title, Option_t *opt = "B") {
        if (!hreco2D || !hgen2D) {
            std::cerr << "Warning: Missing 2D histogram for " << name << ", skipping." << std::endl;
            return;
        }
        TH1D *hreco = hreco2D->ProjectionY(Form("%s_reco", name));
        TH1D *hgen  = hgen2D->ProjectionY(Form("%s_gen", name));
        TH1D *heff  = (TH1D *)hreco->Clone(name);
        heff->SetTitle(title);
        heff->Divide(hreco, hgen, 1.0, 1.0, opt);
        heff->Write();
        delete hreco;
        delete hgen;
    };

    makePtEff(h_centptreco,       h_centptgen,     "h_pteff",       "Eff (0-50%);pT;Eff");
    makePtEff(h_centptreco_woBDT, h_centptgen,     "h_pteff_woBDT", "Eff no BDT (0-50%);pT;Eff");
    // reco is weighted but gen is not -> not a binomial subset, use plain error propagation
    makePtEff(h_centptreco_ptw,   h_centptgen,     "h_pteff_ptw",   "Eff reco pT-weighted (0-50%);pT;Eff", "");

    // --------------------------------------------------
    // Helper lambda for 1D efficiency
    // --------------------------------------------------
    auto make1DEff = [&](TH1D *hreco, TH1D *hgen, const char *name, const char *title) -> TH1D* {
        if (!hreco || !hgen) {
            std::cerr << "Warning: Missing histogram for " << name << ", skipping." << std::endl;
            return nullptr;
        }
        TH1D *heff = (TH1D *)hreco->Clone(name);
        heff->SetTitle(title);
        heff->Divide(hreco, hgen);
        heff->Write();
        return heff;
    };

    // --------------------------------------------------
    // 1D efficiency per centrality bin (with BDT)
    // --------------------------------------------------
    /*make1DEff(h_ptreco_cent010,  h_ptgen_cent010,  "h_pteff_cent010",  "Eff (0-10%);pT;Eff");
    make1DEff(h_ptreco_cent1020, h_ptgen_cent1020, "h_pteff_cent1020", "Eff (10-20%);pT;Eff");
    make1DEff(h_ptreco_cent2030, h_ptgen_cent2030, "h_pteff_cent2030", "Eff (20-30%);pT;Eff");
    make1DEff(h_ptreco_cent3040, h_ptgen_cent3040, "h_pteff_cent3040", "Eff (30-40%);pT;Eff");
    make1DEff(h_ptreco_cent4050, h_ptgen_cent4050, "h_pteff_cent4050", "Eff (40-50%);pT;Eff");

    // --------------------------------------------------
    // 1D efficiency per centrality bin (without BDT)
    // --------------------------------------------------
    make1DEff(h_ptreco_cent010_woBDT,  h_ptgen_cent010,  "h_pteff_cent010_woBDT",  "Eff no BDT (0-10%);pT;Eff");
    make1DEff(h_ptreco_cent1020_woBDT, h_ptgen_cent1020, "h_pteff_cent1020_woBDT", "Eff no BDT (10-20%);pT;Eff");
    make1DEff(h_ptreco_cent2030_woBDT, h_ptgen_cent2030, "h_pteff_cent2030_woBDT", "Eff no BDT (20-30%);pT;Eff");
    make1DEff(h_ptreco_cent3040_woBDT, h_ptgen_cent3040, "h_pteff_cent3040_woBDT", "Eff no BDT (30-40%);pT;Eff");
    make1DEff(h_ptreco_cent4050_woBDT, h_ptgen_cent4050, "h_pteff_cent4050_woBDT", "Eff no BDT (40-50%);pT;Eff");*/

    fout->Close();
    fin->Close();

    std::cout << "Efficiency histograms written to: " << output_file << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc == 3) {
        computeEfficiency(argv[1], argv[2]);
    } else {
        std::cout << "Usage: root -l -b -q 'computeEfficiency.C(\"input.root\", \"output_eff.root\")'" << std::endl;
        std::cout << "  or compile and run: ./computeEfficiency input.root output_eff.root" << std::endl;
        return 1;
    }
    return 0;
}
