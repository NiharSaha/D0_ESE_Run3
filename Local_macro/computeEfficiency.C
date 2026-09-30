#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TEfficiency.h"
#include "TString.h"
#include <iostream>

void computeEfficiency(TString input_file = "ROOT/Eff_out_combined.root",
		       TString output_file = "Eff_output.root")
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
    TH2D *h_centptreco_woBDT = (TH2D *)fin->Get("centptreco_woBDT");
    TH2D *h_centptgen        = (TH2D *)fin->Get("centptgen");

    // --- Load 1D histograms per centrality bin ---
    TH1D *h_ptgen_cent010  = (TH1D *)fin->Get("ptgen_cent010");
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
    TH1D *h_ptreco_cent4050_woBDT = (TH1D *)fin->Get("ptreco_cent4050_woBDT");

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

    TH2D *h_centpteff_Bino = (TH2D *)h_centptreco->Clone("h_centpteff_Bino");
    h_centpteff_Bino->SetTitle("Reco/Gen efficiency (Binomial);cent;pT");
    h_centpteff_Bino->Divide(h_centptreco, h_centptgen, 1.0, 1.0, "B");
    h_centpteff_Bino->Write();

    // --------------------------------------------------
    // 2D efficiency: without BDT
    // --------------------------------------------------
    TH2D *h_centpteff_woBDT = (TH2D *)h_centptreco_woBDT->Clone("h_centpteff_woBDT");
    h_centpteff_woBDT->SetTitle("Reco/Gen efficiency (no BDT);cent;pT");
    h_centpteff_woBDT->Divide(h_centptreco_woBDT, h_centptgen);
    h_centpteff_woBDT->Write();

    TH2D *h_centpteff_woBDT_Bino = (TH2D *)h_centptreco_woBDT->Clone("h_centpteff_woBDT_Bino");
    h_centpteff_woBDT_Bino->SetTitle("Reco/Gen efficiency (no BDT, Binomial);cent;pT");
    h_centpteff_woBDT_Bino->Divide(h_centptreco_woBDT, h_centptgen, 1.0, 1.0, "B");
    h_centpteff_woBDT_Bino->Write();

    // --------------------------------------------------
    // 2D efficiency: TEfficiency with weighted-events statistics.
    // The events here carry the pT-hat stitching weight, which spans many
    // orders of magnitude across the merged samples, so the plain "B"
    // (binomial) Divide() approximation above is unreliable. TEfficiency's
    // exact/Bayesian statistic options assume unweighted integer counts and
    // are not statistically well-defined for weighted entries, so kFNormal
    // (normal/Wald approximation, built from the Sumw2 already tracked on
    // both histograms) is the only option paired with SetUseWeightedEvents().
    // --------------------------------------------------
    if (TEfficiency::CheckConsistency(*h_centptreco, *h_centptgen)) {
        TEfficiency *eff2D = new TEfficiency(*h_centptreco, *h_centptgen);
        eff2D->SetUseWeightedEvents();
        eff2D->SetStatisticOption(TEfficiency::kFNormal);
        eff2D->SetName("h_centpteff_TEff");
        eff2D->SetTitle("Reco/Gen efficiency (weighted TEfficiency);cent;pT");
        eff2D->Write();
    } else {
        std::cerr << "Warning: h_centptreco/h_centptgen inconsistent, skipping TEfficiency." << std::endl;
    }

    if (TEfficiency::CheckConsistency(*h_centptreco_woBDT, *h_centptgen)) {
        TEfficiency *eff2D_woBDT = new TEfficiency(*h_centptreco_woBDT, *h_centptgen);
        eff2D_woBDT->SetUseWeightedEvents();
        eff2D_woBDT->SetStatisticOption(TEfficiency::kFNormal);
        eff2D_woBDT->SetName("h_centpteff_woBDT_TEff");
        eff2D_woBDT->SetTitle("Reco/Gen efficiency (no BDT, weighted TEfficiency);cent;pT");
        eff2D_woBDT->Write();
    } else {
        std::cerr << "Warning: h_centptreco_woBDT/h_centptgen inconsistent, skipping TEfficiency." << std::endl;
    }

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
        heff->Divide(hreco, hgen, 1.0, 1.0, "B");
        heff->Write();
        return heff;
    };

    // --------------------------------------------------
    // Helper lambda for 1D efficiency: TEfficiency with weighted-events
    // statistics (see the 2D block above for why kFNormal is used here).
    // --------------------------------------------------
    auto make1DEffWeighted = [&](TH1D *hreco, TH1D *hgen, const char *name, const char *title) -> TEfficiency* {
        if (!hreco || !hgen) {
            std::cerr << "Warning: Missing histogram for " << name << ", skipping." << std::endl;
            return nullptr;
        }
        if (hgen->GetEntries() == 0) {
            std::cerr << "Warning: " << hgen->GetName() << " has zero entries, skipping " << name << std::endl;
            return nullptr;
        }
        if (!TEfficiency::CheckConsistency(*hreco, *hgen)) {
            std::cerr << "Warning: " << name << " -- passed/total histograms are inconsistent, skipping." << std::endl;
            return nullptr;
        }

        TEfficiency *eff = new TEfficiency(*hreco, *hgen);
        eff->SetUseWeightedEvents();
        eff->SetStatisticOption(TEfficiency::kFNormal);
        eff->SetName(name);
        eff->SetTitle(title);
        eff->Write();
        return eff;
    };

    // --------------------------------------------------
    // 1D efficiency per centrality bin (with BDT)
    // --------------------------------------------------
    make1DEff(h_ptreco_cent010,  h_ptgen_cent010,  "h_pteff_cent010",  "Eff (0-10%);pT;Eff");
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
    make1DEff(h_ptreco_cent4050_woBDT, h_ptgen_cent4050, "h_pteff_cent4050_woBDT", "Eff no BDT (40-50%);pT;Eff");

    // --------------------------------------------------
    // 1D efficiency per centrality bin: weighted TEfficiency (with BDT)
    // --------------------------------------------------
    make1DEffWeighted(h_ptreco_cent010,  h_ptgen_cent010,  "h_pteff_cent010_TEff",  "Eff (0-10%, weighted TEfficiency);pT;Eff");
    make1DEffWeighted(h_ptreco_cent1020, h_ptgen_cent1020, "h_pteff_cent1020_TEff", "Eff (10-20%, weighted TEfficiency);pT;Eff");
    make1DEffWeighted(h_ptreco_cent2030, h_ptgen_cent2030, "h_pteff_cent2030_TEff", "Eff (20-30%, weighted TEfficiency);pT;Eff");
    make1DEffWeighted(h_ptreco_cent3040, h_ptgen_cent3040, "h_pteff_cent3040_TEff", "Eff (30-40%, weighted TEfficiency);pT;Eff");
    make1DEffWeighted(h_ptreco_cent4050, h_ptgen_cent4050, "h_pteff_cent4050_TEff", "Eff (40-50%, weighted TEfficiency);pT;Eff");

    // --------------------------------------------------
    // 1D efficiency per centrality bin: weighted TEfficiency (without BDT)
    // --------------------------------------------------
    make1DEffWeighted(h_ptreco_cent010_woBDT,  h_ptgen_cent010,  "h_pteff_cent010_woBDT_TEff",  "Eff no BDT (0-10%, weighted TEfficiency);pT;Eff");
    make1DEffWeighted(h_ptreco_cent1020_woBDT, h_ptgen_cent1020, "h_pteff_cent1020_woBDT_TEff", "Eff no BDT (10-20%, weighted TEfficiency);pT;Eff");
    make1DEffWeighted(h_ptreco_cent2030_woBDT, h_ptgen_cent2030, "h_pteff_cent2030_woBDT_TEff", "Eff no BDT (20-30%, weighted TEfficiency);pT;Eff");
    make1DEffWeighted(h_ptreco_cent3040_woBDT, h_ptgen_cent3040, "h_pteff_cent3040_woBDT_TEff", "Eff no BDT (30-40%, weighted TEfficiency);pT;Eff");
    make1DEffWeighted(h_ptreco_cent4050_woBDT, h_ptgen_cent4050, "h_pteff_cent4050_woBDT_TEff", "Eff no BDT (40-50%, weighted TEfficiency);pT;Eff");

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
