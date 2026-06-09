//=======================================================
// OPTION 1: 10-Bin Comparison
// OPTION 2: Mean/RMS Ratios with Statistical Error Bars 
// OPTION 3: Multiplicity Proxy Ratio 
// OPTION 4: Resolution-Scaled q2 Comparison 
// OPTION 5: 1% Centrality Weight Comparison (0-10% Bins) 
// OPTION 6: q2 Weight ratio
//=======================================================

void compare_PbPb2018vs2023(int option = 6) {
    TFile *f18 = TFile::Open("hist_out_weight_PbPb2018.root");
    TFile *f23 = TFile::Open("hist_out_weight_PbPb2023.root");

    if (!f18 || !f23) {
        cout << "Check file paths! Ensure hist_out_2018.root and hist_out_2023.root are in this folder." << endl;
        return;
    }
    //====================================
    // --- OPTION 1: 10-Bin Comparison ---
    //====================================
    
    if (option == 1) {
        for (int i = 0; i < 10; i++) {
            TH1D *h18 = (TH1D*)f18->Get(Form("Slices/hRawQ2_c%d", i));
            TH1D *h23 = (TH1D*)f23->Get(Form("Slices/hRawQ2_c%d", i));
            if (!h18 || !h23) continue;

            h18->Scale(1.0/h18->Integral());
            h23->Scale(1.0/h23->Integral());

            TCanvas *c = new TCanvas(Form("c%d", i), "", 800, 600);
            h18->SetLineColor(kBlue); h23->SetLineColor(kRed);
            h18->GetXaxis()->SetRangeUser(0, 0.15);
            h18->Draw("HIST"); h23->Draw("SAME");
        }
    }

    //====================================
    // --- OPTION 2: Mean/RMS Ratios with Statistical Error Bars ---
    //====================================
    if (option == 2) {
        TH1D *hMeanRatio = new TH1D("hMeanRatio", "Ratio of Means (2023/2018);Centrality (%);Mean Ratio", 90, 0, 90);
        TH1D *hRMSRatio  = new TH1D("hRMSRatio",  "Ratio of RMS (2023/2018);Centrality (%);RMS Ratio", 90, 0, 90);

        for (int i = 0; i < 90; i++) {
            TH1D *h18 = (TH1D*)f18->Get(Form("Slices/hRawQ2_c%d", i));
            TH1D *h23 = (TH1D*)f23->Get(Form("Slices/hRawQ2_c%d", i));
            if (!h18 || !h23) continue;

            double m1 = h18->GetMean();
	    double m2 = h23->GetMean();
            double e1 = h18->GetMeanError();
	    double e2 = h23->GetMeanError();
            
            double r1 = h18->GetRMS();
	    double r2 = h23->GetRMS();
            double re1 = h18->GetRMSError();
	    double re2 = h23->GetRMSError();

            if (m1 > 0 && m2 > 0) {
                double ratio = m2 / m1;
                // Error propagation for Ratio R = A/B: sigma_R = R * sqrt((sigA/A)^2 + (sigB/B)^2)
                double err = ratio * sqrt(pow(e1/m1, 2) + pow(e2/m2, 2));
                hMeanRatio->SetBinContent(i+1, ratio);
                hMeanRatio->SetBinError(i+1, err);
            }
            if (r1 > 0 && r2 > 0) {
                double ratio = r2 / r1;
                double err = ratio * sqrt(pow(re1/r1, 2) + pow(re2/r2, 2));
                hRMSRatio->SetBinContent(i+1, ratio);
                hRMSRatio->SetBinError(i+1, err);
            }
        }

        TCanvas *c2 = new TCanvas("c2", "Trend Analysis with Errors", 1200, 500);
        c2->Divide(2, 1);
        c2->cd(1); gPad->SetGrid();
        hMeanRatio->SetMarkerStyle(20); hMeanRatio->GetYaxis()->SetRangeUser(0.8, 1.1);
        hMeanRatio->Draw("E1P");
        
        c2->cd(2); gPad->SetGrid();
        hRMSRatio->SetMarkerStyle(20); hRMSRatio->SetMarkerColor(kRed); hRMSRatio->GetYaxis()->SetRangeUser(0.8, 1.1);
        hRMSRatio->Draw("E1P");
    }

    //====================================
    // --- OPTION 3: Multiplicity Proxy Ratio ---
    //====================================
    if (option == 3) {
        TH1D *hEventRatio = new TH1D("hEventRatio", "Event Count Ratio (2023/2018);Centrality (%);Ratio", 90, 0, 90);
        for (int i = 0; i < 90; i++) {
            TH1D *h18 = (TH1D*)f18->Get(Form("Slices/hRawQ2_c%d", i));
            TH1D *h23 = (TH1D*)f23->Get(Form("Slices/hRawQ2_c%d", i));
            if (!h18 || !h23) continue;
            
            double n18 = h18->GetEntries();
            double n23 = h23->GetEntries();
            if (n18 > 0) hEventRatio->SetBinContent(i + 1, n23 / n18);
        }
        TCanvas *c3 = new TCanvas("c3", "Sample Size Ratio", 800, 600);
        hEventRatio->SetMarkerStyle(24);
        hEventRatio->Draw("P");
        cout << "Note: If your Ntuples have a 'mult' branch, plot that ratio instead for better physical insight." << endl;
    }

    //====================================
    // --- OPTION 4: Resolution-Scaled q2 Comparison ---
    //====================================
    if (option == 4) {
        // This checks if the ratio of means follows the 1/sqrt(N) scaling
        TH1D *hResScale = new TH1D("hResScale", "Scaling Check: (Mean Ratio) / sqrt(N18/N23);Cent;Double Ratio", 90, 0, 90);
        
        for (int i = 0; i < 90; i++) {
            TH1D *h18 = (TH1D*)f18->Get(Form("Slices/hRawQ2_c%d", i));
            TH1D *h23 = (TH1D*)f23->Get(Form("Slices/hRawQ2_c%d", i));
            if (!h18 || !h23) continue;

            double m18 = h18->GetMean();
            double m23 = h23->GetMean();
            double n18 = h18->GetEntries();
            double n23 = h23->GetEntries();

            if (m18 > 0 && n23 > 0) {
                double mean_ratio = m23 / m18;
                double stat_scaling = sqrt(n18 / n23);
                // If this value is 1.0, the shift is purely due to different event counts/stats
                hResScale->SetBinContent(i + 1, mean_ratio / stat_scaling);
            }
        }
        TCanvas *c4 = new TCanvas("c4", "Resolution Scaling Check", 800, 600);
        gPad->SetGrid();
        hResScale->GetYaxis()->SetRangeUser(0.8, 1.2);
        hResScale->SetMarkerStyle(21);
        hResScale->SetMarkerColor(kBlue+2);
        hResScale->Draw("P");
        
        TLine *l1 = new TLine(0, 1, 90, 1);
        l1->SetLineStyle(2); l1->Draw();
        
        cout << "Observation: If the blue points sit on the dashed line at 1.0, the shift is statistical." << endl;
    }

    //====================================
    // --- OPTION 5: 1% Centrality Weight Comparison (0-10% Bins) ---
    //====================================
    if (option == 5) {
        for (int centBin = 0; centBin < 10; centBin++) {
            // Access histograms from the internal 'Weights' directory
            TH1D *hw18 = (TH1D*)f18->Get(Form("Weights/hWeightQ2_c%d", centBin));
            TH1D *hw23 = (TH1D*)f23->Get(Form("Weights/hWeightQ2_c%d", centBin));

            if (!hw18 || !hw23) {
                cout << "Weight histograms not found for bin " << centBin << ". Skipping..." << endl;
                continue;
            }

            // Create a unique canvas for each bin
            TCanvas *cWeight = new TCanvas(Form("cWeight_cent%d", centBin), 
                                          Form("Weight Comparison Cent %d", centBin), 800, 600);
            gPad->SetGrid();
            gPad->SetLogy(); // Log scale is better for energy distributions

	    // Normalization to compare shapes
            hw18->Scale(1.0 / hw18->Integral());
            hw23->Scale(1.0 / hw23->Integral());

            // Formatting
            hw18->SetLineColor(kBlue);
            hw18->SetLineWidth(2);
            hw18->SetTitle(Form("HF Weight (E_{T}) Comparison (Cent %d-%d%%);E_{T} Sum;Probability", centBin, centBin+1));
            hw18->SetStats(0);

            hw23->SetLineColor(kRed);
            hw23->SetLineWidth(2);
            hw23->SetLineStyle(2); // Dashed line for 2023

	    // Drawing
            hw18->Draw("HIST");
            hw23->Draw("HIST SAME");

            // Legend
            TLegend *leg = new TLegend(0.55, 0.75, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->AddEntry(hw18, Form("2018 (Mean: %.1f)", hw18->GetMean()), "l");
            leg->AddEntry(hw23, Form("2023 (Mean: %.1f)", hw23->GetMean()), "l");
            leg->Draw();
        }
        cout << "Opened 10 canvases for weight comparison. Check the Mean values in the legends." << endl;
    }



    //====================================
    // --- OPTION 6: q2 Weight ratio 
    //====================================
    if (option == 6) {

      TH1D *hWeightRatio = new TH1D("hWeightRatio", "Ratio of Mean Weights (2023/2018);Centrality (%);Ratio <W_{2023}> / <W_{2018}>", 90, 0, 90);
        
        // Canvas for raw weight distributions of a few sample bins
        TCanvas *cRaw = new TCanvas("cRaw", "Sample Weight Distributions", 1200, 400);
        cRaw->Divide(3, 1);
        int sampleBins[3] = {10, 45, 80}; // Central, Mid, Peripheral

        for (int i = 0; i < 90; i++) {
            // Access histograms from the internal 'Weights' directory
            TH1D *hw18 = (TH1D*)f18->Get(Form("Weights/hWeightQ2_c%d", i));
            TH1D *hw23 = (TH1D*)f23->Get(Form("Weights/hWeightQ2_c%d", i));

	    if (!hw18 || !hw23) continue;

            double m18 = hw18->GetMean();
            double m23 = hw23->GetMean();
            double e18 = hw18->GetMeanError();
            double e23 = hw23->GetMeanError();

            if (m18 > 0) {
                double ratio = m23 / m18;
                double err = ratio * sqrt(pow(e18/m18, 2) + pow(e23/m23, 2));
                hWeightRatio->SetBinContent(i + 1, ratio);
                hWeightRatio->SetBinError(i + 1, err);
            }

	    // Plot sample distributions to see the shift visually
            for(int s=0; s<3; s++) {
                if(i == sampleBins[s]) {
                    cRaw->cd(s+1); gPad->SetLogy();
                    hw18->SetLineColor(kBlue); hw18->Scale(1.0/hw18->Integral());
                    hw23->SetLineColor(kRed);  hw23->Scale(1.0/hw23->Integral());
                    hw18->SetTitle(Form("Weight Dist Cent %d%%;E_{T} Sum;Prob", i));
                    hw18->Draw("HIST"); hw23->Draw("HIST SAME");
                }
            }
        }

	TCanvas *cRatio = new TCanvas("cRatio", "Weight Ratio Trend", 800, 600);
        gPad->SetGrid();
        hWeightRatio->SetMarkerStyle(20);
        hWeightRatio->SetMarkerColor(kGreen+2);
        hWeightRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
        hWeightRatio->Draw("E1P");

        TLine *l = new TLine(0, 1, 90, 1);
        l->SetLineStyle(2);
        l->Draw();

        cout << "--- Weight Analysis Complete ---" << endl;
        cout << "Observation: If Weight Ratio matches your q2 ratio trend, the shift is energy-scale driven." << endl;
    }



    
}// -- END ---
