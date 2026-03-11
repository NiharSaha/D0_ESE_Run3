#include <iostream>
#include <vector>
#include "TFile.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TDirectory.h"

using namespace std;

// --- Global Config ---
const int N_CENT_GROUPS = 4;
int min_cent[] = {0, 10, 30, 50};
int max_cent[] = {10, 30, 50, 90};
const int N_Q_BINS = 10;
const int N_SP_TOTAL = 80;
Double_t sp_edges[81];

void SetupSPBins() {
    int b = 0;
    for(double x = -10.0; x < -4.0;  x += 0.5)  sp_edges[b++] = x;
    for(double x = -4.0;  x < -2.0;  x += 0.25) sp_edges[b++] = x;
    for(double x = -2.0;  x < 0.0;   x += 0.1)  sp_edges[b++] = x;
    for(double x = 0.0;   x < 2.0;   x += 0.1)  sp_edges[b++] = x;
    for(double x = 2.0;   x < 4.0;   x += 0.25) sp_edges[b++] = x;
    for(double x = 4.0;   x <= 10.0; x += 0.5)  sp_edges[b++] = x;
}


void RunPreFit(TH1D* &hInclusive, double* parBuffer) {

  if (!hInclusive) return;

    double fitMin = 1.75;
    double fitMax = 1.98;
    double peak = hInclusive->GetMaximum();
    double bkgInitial = peak * 0.8;;

    //double signalHeight = peak - bkgInitial;
    // 2. Define the complex 14-parameter function
    TString funcName = TString::Format("fGlobal_%s", hInclusive->GetName());
    
    TF1* fGlobal = new TF1(funcName,
        "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0+[6]))/(sqrt(2*3.14159)*[2]*(1.0+[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0+[6]))/(sqrt(2*3.14159)*[3]*(1.0+[6]))) + (1-[5])*TMath::Gaus(x,[1],[7]*(1.0+[6]))/(sqrt(2*3.14159)*[7]*(1.0+[6])))"
        " + [8] + [9]*x + [10]*x*x"
        " + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12]))"
        " + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])))",
        fitMin, fitMax);

    fGlobal->SetParameter(0, (peak - bkgInitial) * 0.05); // Signal Yield Norm
    fGlobal->SetParameter(1, 1.865);                     // Mass Peak
    fGlobal->SetParameter(2, 0.012);                     // Sigma 1
    fGlobal->SetParameter(3, 0.025);                     // Sigma 2
    fGlobal->SetParameter(4, 0.7);                       // Frac 1
    fGlobal->SetParameter(5, 0.8);                       // Frac 2
    fGlobal->SetParameter(6, 0.0);                       // Resolution Scale
    fGlobal->SetParameter(7, 0.045);                     // Sigma 3
    fGlobal->SetParameter(8, bkgInitial);                // Bkg Pol0
    fGlobal->SetParameter(9, -0.1);                      // Bkg Pol1
    fGlobal->SetParameter(10, 0.01);                     // Bkg Pol2
    fGlobal->SetParameter(11, peak * 0.01);              // Reflection Norm
    fGlobal->SetParameter(12, 0.0);                      // Refl Scale 1
    fGlobal->SetParameter(13, 0.0);                      // Refl Scale 2

    // Step-wise fitting strategy for better convergence
    fGlobal->FixParameter(11, 0); // Start with reflections off
    hInclusive->Fit(fGlobal, "RLQN", "", fitMin, fitMax);

    fGlobal->ReleaseParameter(11);
    fGlobal->SetParLimits(1, 1.855, 1.875); // Keep peak near D0 mass
    hInclusive->Fit(fGlobal, "RLQ", "", fitMin, fitMax);
    

    // Fill the buffer to lock shape for the SP bins
    for(int i=0; i<14; i++) parBuffer[i] = fGlobal->GetParameter(i);
    
    hInclusive->GetListOfFunctions()->Add(fGlobal);

}


void ExtractYield(TH1D* hMass, double &yield, double &yieldErr, TString hName, double* globalPars) {
    // 1. Minimum Entry Guard (More strict for 14-par fits)
  if (!hMass || hMass->GetEntries() < 20) {
    yield = 0; yieldErr = 0;
    return;
  }

    double low = 1.75, high = 1.98;
    // 2. Define the Full 14-Parameter Function
    
    TF1 *fTotal = new TF1(Form("fTotal_%s", hName.Data()),
        "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0+[6]))/(sqrt(2*3.14159)*[2]*(1.0+[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0+[6]))/(sqrt(2*3.14159)*[3]*(1.0+[6]))) + (1-[5])*TMath::Gaus(x,[1],[7]*(1.0+[6]))/(sqrt(2*3.14159)*[7]*(1.0+[6])))" // Signal
        " + [8] + [9]*x + [10]*x*x" // Combinatorial Background			  
        " + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12]))" // Reflection 1
        " + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])))", // Reflection 2
        low, high);
    
    // 3. Set Parameters from Pre-Fit (globalPars)
    fTotal->SetParameters(globalPars);
    // 4. FIX SHAPE PARAMETERS (Crucial for stability in Q-bins/SP-bins)
    // We assume the D0 shape doesn't change with Event Shape, only its yield does.
    fTotal->FixParameter(2, globalPars[2]); // Sigma 1
    fTotal->FixParameter(3, globalPars[3]); // Sigma 2
    fTotal->FixParameter(4, globalPars[4]); // Fraction 1
    fTotal->FixParameter(5, globalPars[5]); // Fraction 2
    fTotal->FixParameter(7, globalPars[7]); // Sigma 3
    
    // 5. Constrain Floating Parameters
    fTotal->SetParLimits(1, globalPars[1]-0.003, globalPars[1]+0.003); // Peak position shouldn't move much
    fTotal->SetParLimits(6, -0.05, 0.05);  // Allow resolution to vary by +/- 5%
   

    // 6. Perform the Fit (Likelihood + Range + No Graphics + Quiet)
    
    //hMass->Fit(fTotal, "QNRL");
    hMass->Fit(fTotal, "QNRL", "", low, high);

    // 7. Extract Components for Visualization
    // --- Signal Component ---
    
    TF1 *fSig = new TF1(Form("fSig_%s", hName.Data()), "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0+[6]))/(sqrt(2*3.14159)*[2]*(1.0+[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0+[6]))/(sqrt(2*3.14159)*[3]*(1.0+[6])))+(1-[5])*TMath::Gaus(x,[1],[7]*(1.0+[6]))/(sqrt(2*3.14159)*[7]*(1.0+[6])))", low, high);
    for(int i=0; i<8; i++) fSig->SetParameter(i, fTotal->GetParameter(i));
    fSig->SetLineColor(kRed);
    fSig->SetNpx(500);
    // --- Background Component ---
    TF1 *fBkg = new TF1(Form("fBkg_%s", hName.Data()), "pol2", low, high);
    for(int i=0; i<3; i++) fBkg->SetParameter(i, fTotal->GetParameter(8+i));
    fBkg->SetLineColor(kBlue);
    fBkg->SetLineStyle(2);
    // --- Reflection Component ---
    
    TF1 *fRef = new TF1(Form("fRef_%s", hName.Data()), "[0]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[1]), 1.96*(1+[2])) + 4*[0]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[1]), 1.7734*(1+[3])))", low, high);
    fRef->SetParameters(fTotal->GetParameter(11), fTotal->GetParameter(6), fTotal->GetParameter(12), fTotal->GetParameter(13));
    fRef->SetLineColor(kGreen+2);
    fRef->SetLineStyle(7);


    if (fTotal) hMass->GetListOfFunctions()->Add((TF1*)fTotal->Clone());
    if (fSig)   hMass->GetListOfFunctions()->Add((TF1*)fSig->Clone());
    if (fBkg)   hMass->GetListOfFunctions()->Add((TF1*)fBkg->Clone());
    if (fRef)   hMass->GetListOfFunctions()->Add((TF1*)fRef->Clone());
    // 9. Final Yield Calculation
    
    // Total Yield is the area under the signal function [Parameter 0]
    
    yield = fTotal->GetParameter(0) / hMass->GetBinWidth(1);
    yieldErr = fTotal->GetParError(0) / hMass->GetBinWidth(1);

    delete fTotal;
    delete fSig;
    delete fBkg;
    delete fRef;

}


void SP_new(TString input_file, TString output_dir, int harm, int cIdx, int pIdx) {
    SetupSPBins();
    TFile *fIn = TFile::Open(input_file, "READ");
    THnSparseD *hn = (THnSparseD*)fIn->Get(Form("D0_vn_THn/hn_v%d", harm));

    TString outName = TString::Format("%s/Result_v%d_c%d_p%d.root", output_dir.Data(), harm, cIdx, pIdx);
    TFile *fOut = new TFile(outName, "RECREATE");

    for(int i=3; i<=5; i++) hn->GetAxis(i)->SetRange(0, 0); 
    hn->GetAxis(1)->SetRange(pIdx + 1, pIdx + 1); // pT bin
    hn->GetAxis(2)->SetRange(min_cent[cIdx] + 1, max_cent[cIdx]); // Cent bin
    
    TString incName = Form("hInclusive_v%d_c%d_p%d", harm, cIdx, pIdx);
    TH1D* hInclusive = (TH1D*)hn->Projection(0, incName.Data());
    if (!hInclusive) {
      cout << "Error: Inclusive projection failed!" << endl;
      return;
    }
    hInclusive->SetDirectory(0);
    
    // --- 2. RUN THE PRE-FIT TO GET INITIAL PARAMETERS ---
    double parBuf[14];
    RunPreFit(hInclusive, parBuf); 
    
    fOut->cd();
    hInclusive->Write();



    for (int q = 0; q < N_Q_BINS; q++) {
        TDirectory *qDir = fOut->mkdir(Form("Qbin_%d", q));
        
        vector<double> yields, yieldErrs, sp_means;
        double sum_y = 0, sum_y_vn = 0;

        for (int s = 0; s < N_SP_TOTAL; s++) {
            hn->GetAxis(1)->SetRange(pIdx + 1, pIdx + 1);
            hn->GetAxis(2)->SetRange(min_cent[cIdx] + 1, max_cent[cIdx]);
            hn->GetAxis(3)->SetRange(q + 1, q + 1);
            hn->GetAxis(5)->SetRange(s + 1, s + 1);

            TH1D *hS = (TH1D*)hn->Projection(0);
            double y, yE;
            ExtractYield(hS, y, yE, Form("hM_q%d_s%d", q, s), parBuf);
            
            yields.push_back(y); yieldErrs.push_back(yE);
            sp_means.push_back((sp_edges[s] + sp_edges[s+1]) / 2.0);
            
            if (y > 0) { sum_y += y; sum_y_vn += (y * sp_means.back()); }
            qDir->cd(); hS->Write(); delete hS;
        }

        if (sum_y > 0) {
            double vn_avg = sum_y_vn / sum_y;
            double sum_err_sq = 0;
            for (int i = 0; i < N_SP_TOTAL; i++) {
                if (yields[i] <= 0) continue;
                sum_err_sq += TMath::Power(sp_means[i] - vn_avg, 2) * TMath::Power(yieldErrs[i], 2);
            }
            double vn_err = TMath::Sqrt(sum_err_sq) / sum_y;

            TGraphErrors *gr = new TGraphErrors(1);
            gr->SetName(Form("gr_v%d_c%d_q%d_p%d", harm, cIdx, q, pIdx));
            gr->SetPoint(0, pIdx, vn_avg); // Store index or mid-pT
            gr->SetPointError(0, 0, vn_err);
            fOut->cd(); gr->Write();
        }
    }
    fOut->Close();
}
