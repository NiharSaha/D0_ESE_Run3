#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TDirectory.h"

using namespace std;

//==================================
// Global variables definition
//==================================

const int N_CENT = 4;
const int N_PT = 10;
const int N_Q = 10;
int min_cent[N_CENT] = {0, 10, 30, 50};
int max_cent[N_CENT] = {10, 30, 50, 90};
Float_t min_pt[N_PT] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0};
Float_t max_pt[N_PT] = {2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 40.0};

const int N_SP_TOTAL = 80;
Double_t sp_edges[81];

void SetupSPBins() {
    int b = 0;
    // Negative Side
    for(double x = -10.0; x < -4.0;  x += 0.5)  sp_edges[b++] = x;
    for(double x = -4.0;  x < -2.0;  x += 0.25) sp_edges[b++] = x;
    for(double x = -2.0;  x < 0.0;   x += 0.1)  sp_edges[b++] = x;
    // Positive Side
    for(double x = 0.0;   x < 2.0;   x += 0.1)  sp_edges[b++] = x;
    for(double x = 2.0;   x < 4.0;   x += 0.25) sp_edges[b++] = x;
    for(double x = 4.0;   x <= 10.0; x += 0.5)  sp_edges[b++] = x;
}

//==================================
// Fitting for inclusive q quantiles bin to fix parameters
//==================================

void RunPreFit(THnSparseD* hn, int harm, int centIdx, int ptBin, double* parBuffer, TH1D* &hInclusive) {


  for(int i=1; i<=5; i++) hn->GetAxis(i)->SetRange(0, 0);

  // 1. Project Inclusive Mass (Resetting axes 3, 4, and 5 to cover everything)
    hn->GetAxis(1)->SetRange(ptBin, ptBin);
    hn->GetAxis(2)->SetRange(min_cent[centIdx]+1, max_cent[centIdx]);


    TString name = Form("hInclusive_v%d_c%d_p%d", harm, centIdx, ptBin);
    hInclusive = (TH1D*)hn->Projection(0, name.Data());
    if (!hInclusive) { 
      cout << "Error: Projection failed for " << name << endl; 
      return; 
    }
    hInclusive->SetDirectory(0);
    //hInclusive->SetTitle(Form("Inclusive Mass v%d, CentIdx %d, pT Bin %d", harm, centIdx, ptBin));

    double fitMin = 1.75;
    double fitMax = 1.98;
    double bkgInitial = 10.5e6;
    double peak = hInclusive->GetMaximum();
    double signalHeight = peak - bkgInitial;
    // 2. Define the complex 14-parameter function
    TF1* fGlobal = new TF1(Form("fGlobal_%s", name.Data()), 
        "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0+[6]))/(sqrt(2*3.14159)*[2]*(1.0+[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0+[6]))/(sqrt(2*3.14159)*[3]*(1.0+[6]))) + (1-[5])*TMath::Gaus(x,[1],[7]*(1.0+[6]))/(sqrt(2*3.14159)*[7]*(1.0+[6])))"
        " + [8] + [9]*x + [10]*x*x"
        " + [11]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[6]), 1.96*(1+[12]))"
        " + 4*[11]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[6]), 1.7734*(1+[13])))", 
        fitMin, fitMax);

    // 3. Robust Initial Guesses for Inclusive Peak

    fGlobal->SetParameter(0, signalHeight * 0.02);
    fGlobal->SetParameter(1, 1.865);
    fGlobal->SetParameter(2, 0.012);
    fGlobal->SetParameter(3, 0.025);
    fGlobal->SetParameter(4, 0.7);
    fGlobal->SetParameter(5, 0.8);
    fGlobal->SetParameter(6, 0.0);
    fGlobal->SetParameter(7, 0.045);
    fGlobal->SetParameter(8, bkgInitial);
    fGlobal->SetParameter(9, -peak*0.01);
    fGlobal->SetParameter(10, 0.001);
    fGlobal->SetParameter(11, peak*0.005);
    fGlobal->SetParameter(12, 0.0);
    fGlobal->SetParameter(13, 0.0);


    fGlobal->FixParameter(1, 1.865); 
    fGlobal->FixParameter(11, 0); // Temporary no reflection
    
    // Fit only background and normalization
    hInclusive->Fit(fGlobal, "RLQN", "", fitMin, fitMax);

    fGlobal->ReleaseParameter(1);
    fGlobal->ReleaseParameter(11);
    fGlobal->SetParLimits(1, 1.855, 1.875);
    
    hInclusive->Fit(fGlobal, "RLQ", "", fitMin, fitMax);
    //hInclusive->Fit(fGlobal, "RLQ");

    // 4. Fill parBuffer to pass to ExtractYield
    for(int i=0; i<14; i++) parBuffer[i] = fGlobal->GetParameter(i);
    hInclusive->GetListOfFunctions()->Add(fGlobal);

}


//==================================
// Fitting for q quantiles bin
//==================================
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
    hMass->Fit(fTotal, "QNRL", "", fitMin, fitMax);

    
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

    // 8. Attach functions to the histogram so they save into the .root file
    /*hMass->GetListOfFunctions()->Add(fTotal);
    hMass->GetListOfFunctions()->Add(fSig);
    hMass->GetListOfFunctions()->Add(fBkg);
    hMass->GetListOfFunctions()->Add(fRef);
    */

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



//==================================
// Main function to perform analysis
//==================================


void Extract_Flow_Final_Directories(
				    TString inputFile,
				    int harm,
				    int specificCent,
				    int specificPt,
				    TString outputFilePath
				    ) {

  SetupSPBins();

  double globalPars[14];
  TH1D* hPreFitResult = nullptr;



  TFile *fIn = TFile::Open(inputFile);
    THnSparseD *hn = (THnSparseD*)fIn->Get(Form("D0_vn_THn/hn_v%d", harm));
    TH2D *hRes = (TH2D*)fIn->Get(Form("D0_vn_THn/hRes_v%d", harm));


    // --- Data Storage ---
    TGraphErrors* gVn_vs_Pt_Int[N_CENT];
    TGraphErrors* gVn_vs_Q[N_CENT][N_PT];
    TH1D* hSP_Dist[N_CENT][N_PT][N_Q][2];
    TH1D* hMassFits[N_CENT][N_PT][N_Q][2][N_SP_TOTAL];

    for(int i=0; i<N_CENT; i++) {
        gVn_vs_Pt_Int[i] = new TGraphErrors();
        for(int p=0; p<N_PT; p++) gVn_vs_Q[i][p] = new TGraphErrors();
    }

    for (int i = 0; i < N_CENT; i++) {

      if (specificCent != -1 && i != specificCent) continue;
      cout << "\n>>> Processing v" << harm << " | Centrality: " << min_cent[i] << "-" << max_cent[i] << "%" << endl;

      for (int p = 1; p <= N_PT; p++) {

	if (specificPt != -1 && p != specificPt) continue;
	cout << ">>> Processing v" << harm << " | Cent " << i << " | pT Bin " << p << endl;



	
	if (specificPt == -1 || p == specificPt) {
	  RunPreFit(hn, harm, i, p, globalPars, hPreFitResult);
	  
	  //cout << ">>> Processing RunPreFit <<<<" << endl;

	  //cout << "DEBUG: hRes pointer is valid: " << hRes->GetName() << endl;
	
	double global_sum_YiVi = 0, global_sum_Yi = 0;
	vector<double> glob_vi, glob_err;

            for (int q = 1; q <= N_Q; q++) {
                int res_idx = i * N_Q + (q - 1); 
                double resP = TMath::Sqrt((hRes->GetBinContent(1,res_idx+1)*hRes->GetBinContent(2,res_idx+1))/hRes->GetBinContent(3,res_idx+1));
                double resM = TMath::Sqrt((hRes->GetBinContent(1,res_idx+1)*hRes->GetBinContent(3,res_idx+1))/hRes->GetBinContent(2,res_idx+1));

                hn->GetAxis(1)->SetRange(p, p); 
                hn->GetAxis(2)->SetRange(min_cent[i]+1, max_cent[i]);
                hn->GetAxis(3)->SetRange(q, q);

                double local_sum_YiVi = 0, local_sum_Yi = 0;
                vector<double> loc_vi, loc_err;

                for (int sub = 0; sub < 2; sub++) {
                    hn->GetAxis(4)->SetRange(sub+1, sub+1);
                    hSP_Dist[i][p-1][q-1][sub] = new TH1D(Form("hSP_c%d_p%d_q%d_sub%d", i, p, q, sub), "Signal SP Distribution", N_SP_TOTAL, sp_edges);
                    double curRes = (sub == 0) ? resM : resP;

                    for (int b = 1; b <= N_SP_TOTAL; b++) {
                        hn->GetAxis(5)->SetRangeUser(sp_edges[b-1] + 0.0001, sp_edges[b] - 0.0001);

			TString uniqueProjName = Form("hMass_c%d_p%d_q%d_sub%d_sp%d", i, p, q, sub, b);
			TH1D *hSlice = (TH1D*)hn->Projection(0, uniqueProjName.Data());
			if (!hSlice) continue;
			hSlice->SetDirectory(0);
	
                        //double yield = 100, err = 10;
			double yield, err;
			ExtractYield(hSlice, yield, err, uniqueProjName, globalPars);

                        if (yield > 0) {
                            double vi_scaled = ((sp_edges[b-1]+sp_edges[b])/2.0) / curRes;
                            local_sum_YiVi += (yield * vi_scaled); local_sum_Yi += yield;
                            loc_vi.push_back(vi_scaled); loc_err.push_back(err);
                            global_sum_YiVi += (yield * vi_scaled); global_sum_Yi += yield;
                            glob_vi.push_back(vi_scaled); glob_err.push_back(err);

                            hSP_Dist[i][p-1][q-1][sub]->SetBinContent(b, yield);
                            hMassFits[i][p-1][q-1][sub][b-1] = hSlice;
                        }
			else delete hSlice;

                    } // -- SP loop ---
		    //cout << ">>> End SP loop <<<<" << endl;
		    hn->GetAxis(5)->SetRange(0, 0); 

                } // --- sub loop ---
		//cout << ">>> End Sub loop <<<<" << endl;

		    if (local_sum_Yi > 0) {
		      double vn = local_sum_YiVi / local_sum_Yi; double esq = 0;
		      for(size_t k=0; k<loc_vi.size(); k++) esq += TMath::Power(loc_vi[k]-vn, 2) * TMath::Power(loc_err[k], 2);
		      gVn_vs_Q[i][p-1]->SetPoint(q-1, q*10-5, vn);
		      gVn_vs_Q[i][p-1]->SetPointError(q-1, 0, TMath::Sqrt(esq)/local_sum_Yi);
		    }
            } // --- q loop ---

	    //cout << ">>> End Q loop <<<<" << endl;
	    
            if (global_sum_Yi > 0) {
                double vn = global_sum_YiVi / global_sum_Yi; double esq = 0;
                for(size_t k=0; k<glob_vi.size(); k++) esq += TMath::Power(glob_vi[k]-vn, 2) * TMath::Power(glob_err[k], 2);
                gVn_vs_Pt_Int[i]->SetPoint(p-1, (min_pt[p-1]+max_pt[p-1])/2.0, vn);
                gVn_vs_Pt_Int[i]->SetPointError(p-1, 0, TMath::Sqrt(esq)/global_sum_Yi);
            }
        }// -- loop before preFit ---
      } // --- pT loop ----
    }// -- cent loop ---


//==================================
// Store and save histograms
//==================================

TFile *fOut = new TFile(outputFilePath, "RECREATE");

// 1. Save Pre-Fit Checks (Only once per file)
 if (hPreFitResult) {
   fOut->mkdir("PreFit_Checks")->cd();
   hPreFitResult->Write();
 }
 
 // 2. Loop through centralities, but ONLY create/fill for the assigned one
 for (int i = 0; i < N_CENT; i++) {
   // CRITICAL: Only process the centrality assigned to this specific job
   if (specificCent != -1 && i != specificCent) continue;
   
   TString dirName = Form("Centrality_%d_%d", min_cent[i], max_cent[i]);
   TDirectory *dirCent = fOut->mkdir(dirName);
   
   // Create Subdirectories
   TDirectory *dirGr = dirCent->mkdir("PhysicsGraphs");
   TDirectory *dirSP = dirCent->mkdir("SP_Distributions");
   TDirectory *dirFits = dirCent->mkdir("MassFits");
   
   // Save Integrated vN vs Pt
   dirGr->cd();
   gVn_vs_Pt_Int[i]->SetName(Form("gV%d_vs_Pt_Integrated", harm));
   gVn_vs_Pt_Int[i]->Write();
   
   // Save vN vs Q for each pT bin
   for(int p=0; p < N_PT; p++) {
     if (gVn_vs_Q[i][p]) {
       gVn_vs_Q[i][p]->Write(Form("gV%d_vs_Q_p%d", harm, p+1));
     }
   }
   
   // Save Histograms (SP and Mass Fits)
   for(int p=0; p < N_PT; p++) {
     for(int q=0; q < N_Q; q++) {
       
       dirSP->cd();
       if(hSP_Dist[i][p][q][0]) hSP_Dist[i][p][q][0]->Write();
       if(hSP_Dist[i][p][q][1]) hSP_Dist[i][p][q][1]->Write();
       
       dirFits->cd();
       for(int s=0; s<2; s++) {
	 for(int b=0; b < N_SP_TOTAL; b++) {
	   if(hMassFits[i][p][q][s][b]) {
	     hMassFits[i][p][q][s][b]->Write();
	   }
	 }
       }
     }
   }
 }
 
 fOut->Close();
}
