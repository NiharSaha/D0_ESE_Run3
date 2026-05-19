# D0 ESE Flow Analysis

## Table of Contents

- [Prepare Ingredients for the Analysis](#prepare-ingredients-for-the-analysis)
- [Compare with Published Results](#compare-with-published-results)
- [Doing the Analysis](#doing-the-analysis)

---

## Prepare Ingredients for the Analysis

1. **Save q2, q3 and weights values in a ntuple**
   ```
   Quantiles/make_q2_slices.C
   ```

2. **Extract quantiles and save in a header** *(one-time)*
   ```
   Quantiles/Quantile/PbPb2023/ExtractQuantiles_sortFunc_new.C
   ```

3. **Calculate resolution factor** *(denominator of SP)*
   ```
   Quantiles/Calculate_Resolution_latest.C
   ```

4. **Save all mass, v2, v3 etc. in ntuple** after applying optimized BDT cuts
   ```
   Quantiles/flow_Analysis_latest.C
   ```

5. **Plotting macros** *(local)*
   ```
   /Users/saha115/cernbox/Analysis/D0_ESE/compare_PbPb2018vs2023.C
   /Users/saha115/cernbox/Analysis/D0_ESE/ESE_Diagnostics.C
   ```

---

## Compare with Published Results

1. **Calculate the SP resolution** for coarser centrality bins (0–10%, 10–30%, 30–50%)
   ```
   Compare_PubResults/Calculate_Resolution_forPubResults.C
   ```

2. **Save ntuple with all flow ingredients** for same centrality bins
   ```
   Compare_PubResults/flow_Analysis_ingradients_forPubResults.C
   ```

3. **Make SP binning**
   ```
   Compare_PubResults/save_vnbinning_to_txt.C
   ```

4. **Get the flow results**
   ```
   Compare_PubResults/fit_mass_flow_forPubResults.C
   ```

---

## Doing the Analysis

1. **Get the vn binning and save in a header file** *(one-time)*
   ```
   flow_extraction/save_vnbinning_to_txt.C
   ```

2. **Get 1% mass distributions and save in an output ROOT file**
   ```
   flow_extraction/save_mass_distributions.C
   ```

3. **Fit mass distributions to get flow results**
   ```
   flow_extraction/fit_mass_and_flow.C
   ```

4. **Plotting / diagnostic macros**

   | Purpose | Script |
   |---|---|
   | Merge all SP bin mass fit plots into PDF | `flow_extraction/sbatch_combined_mass_pdf_v2v3.sh` |
   | Draw SP distribution of candidates per (cent, pT) bin | `flow_extraction/Draw_raw_SPv2v3_dist_CentPtBin.C` |
   | Draw SP distribution of candidates per (cent, pT, Qn) bin | `flow_extraction/Draw_raw_SPv2v3_dist_CentPtQuanBin.C` |

5. **Plotting macros** *(local)*
