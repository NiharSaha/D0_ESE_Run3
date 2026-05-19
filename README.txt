
## Prepare ingradients for the analysis!

1. Save q2, q3 and weights values in a ntuple  
   `/home/saha115/D0_ESE/CMSSW_13_2_11/src/Quantiles/make_q2_slices.C`

2. Extract quantiles and save in a header (one-time)  
   `/home/saha115/D0_ESE/CMSSW_13_2_11/src/Quantiles/Quantile/PbPb2023/ExtractQuantiles_sortFunc_new.C`

3. Calculate resolution factor (denominator of SP)  
   `/home/saha115/D0_ESE/CMSSW_13_2_11/src/Quantiles/Calculate_Resolution_latest.C`

4. Save all mass, v2, v3 etc in ntuple after applying optimized bdt cuts  
   `/home/saha115/D0_ESE/CMSSW_13_2_11/src/Quantiles/flow_Analysis_latest.C`

5. ploting macro (local)  
   `/Users/saha115/cernbox/Analysis/D0_ESE/compare_PbPb2018vs2023.C`  
   `/Users/saha115/cernbox/Analysis/D0_ESE/ESE_Diagnostics.C`

## Compare with publish results!

1. Calculate the SP resolution for coarser centrality bins (0-10%, 10-30%, 30-50%)  
   `/home/saha115/D0_ESE/CMSSW_13_2_11/src/Compare_PubResults/Calculate_Resolution_forPubResults.C`

2. Save ntuple with all flow ingradients for same centrality bins.  
   `/home/saha115/D0_ESE/CMSSW_13_2_11/src/Compare_PubResults/flow_Analysis_ingradients_forPubResults.C`

3. Make SP binning  
   `/home/saha115/D0_ESE/CMSSW_13_2_11/src/Compare_PubResults/save_vnbinning_to_txt.C`

4. Get the flow results  
   `/home/saha115/D0_ESE/CMSSW_13_2_11/src/Compare_PubResults/fit_mass_flow_forPubResults.C`

## Doing the analysis!

1. Get the vnbinning and save in a header file (one time)  
   `/home/saha115/D0_ESE/CMSSW_13_2_11/src/flow_extraction/save_vnbinning_to_txt.C`

2. Get 1% mass distribution and save in a output root file  
   `/home/saha115/D0_ESE/CMSSW_13_2_11/src/flow_extraction/save_mass_distributions.C`

3. Fit mass distributions to get flow results  
   `/home/saha115/D0_ESE/CMSSW_13_2_11/src/flow_extraction/fit_mass_and_flow.C`

4. Plotting/Showing macros/scripts

   - To merge all SP bin mass fit plots in pdf format  
     `/home/saha115/D0_ESE/CMSSW_13_2_11/src/flow_extraction/sbatch_combined_mass_pdf_v2v3.sh`

   - To draw SP distribution of candidates for different (cent, pT) bin  
     `/home/saha115/D0_ESE/CMSSW_13_2_11/src/flow_extraction/Draw_raw_SPv2v3_dist_CentPtBin.C`

   - To draw SP distribution of candidates for different (cent, pT, Qn) bin  
     `/home/saha115/D0_ESE/CMSSW_13_2_11/src/flow_extraction/Draw_raw_SPv2v3_dist_CentPtQuanBin.C`

5. Plotting macro (local)
