# BDT-cut systematic (prompt-D0 v2 / v3)

Two approaches live in this directory. **Method B (data-driven cut variation)
is the current/primary one.** Method A (MC with/without-BDT closure) was tried
first, showed large bin-to-bin statistical variation (uncalibrated MC
Q-vectors + thin MC stats once split by the BDT-centrality basis), and is kept
for reference/cross-check.

---

## Method B (current): data-driven BDT cut variation

**Idea**: instead of an arbitrary flat shift of the BDT score, vary the
nominal cut in each (y, centrality, pT) analysis bin by an amount **justified
from real signal/background statistics**, then re-run the normal data
pipeline (mass fit + yield-weighted vn, exactly as for the nominal result)
with the varied cut, and take nominal-vs-varied as the systematic — the same
"vary an ingredient, re-run the chain, take the difference" pattern
`../centrality_sys/` uses for its centrality-calibration systematic.

### Step 1 — `scan_significance_vs_BDT.C`: size the variation
Reads `BDT_training/TMVA.root:dataset/{TestTree,TrainTree}` (the BDT training
sample), which carries `classID`, the evaluated `BDT` score, and — critically
— `npT`/`ncent`/`ny`, the **same kinematic variables `BDTHandler::get_bins`
bins on**. For each of the 104 (y,cent,pT) `BDTHandler` bins:
* builds a proxy figure of merit `FoM(cut) = S/√(S+B)` (unweighted training
  counts; the `weight` branch is 1 for every event in this sample);
* anchors at the **actual nominal cut** (from `bdt_cuts.csv`, loaded via the
  unmodified `BDTHandler`), and measures `δ_tight` = the score shift where
  `FoM` drops 1% (relative) on the well-posed, monotonically-declining
  *tighter* side — see the header comment in the macro for why only that
  side is well-posed here (the proxy's own maximum sits well away from the
  nominal cut, most likely because the analysis was actually tuned on a
  mass-fit significance `S/σ_S`, a different quantity from this raw-count
  proxy — confirmed by inspecting the curve shape). `δ_tight` is then
  **mirrored**: `cut_loose = nominal − δ_tight`, `cut_tight = nominal + δ_tight`.
* pT ≥ 20 GeV (BDTHandler bins 9-12) has zero training statistics (the
  training only covers `1≤pT<20 GeV`) — `δ_tight` from the last covered bin
  (15-20 GeV, same y/cent) is applied as an absolute score shift instead.
  Sentinel nominal cuts (`|v|≥1`, i.e. reject-all/accept-all placeholders)
  are never perturbed.

Interpretation: this is **not** a re-derivation of the original working
point, and the 1% figure is on a proxy metric, not the true `S/σ_S` the
analysis was tuned on. It should be read as a reproducible, bin-dependent
*local sensitivity scale* near the cut actually used, not as "the edge of
the optimal region." See the conversation/commit history for the full
reasoning; a stronger cross-check (building the true `S/σ_S` curve from real
data mass fits) is possible but not done here.

Run (local, seconds, no slurm):
```
g++ scan_significance_vs_BDT.C $(root-config --cflags --libs) -Wall -O2 -o scan_significance_vs_BDT.exe
./scan_significance_vs_BDT.exe
```
Writes `significance_scan.root`, `significance_plots/*.png` (FoM-vs-cut, one
panel per (y,cent), 9 pT curves each), a printed per-bin table, and
`bdt_cuts_loose.csv` / `bdt_cuts_tight.csv` (identical format to
`bdt_cuts.csv` — loadable by the unmodified `BDTHandler`, verified).

### Step 2 — re-run the data chain with the varied cut
`loose/` and `tight/` each hold a **self-contained copy** of the two pipeline
stages the BDT cut touches:

| file | copied from | what's different |
|---|---|---|
| `BDTHandler.h`/`.cc` | `Quantiles/for_20Qbin/` | ctor hardcodes `bdt_cuts_loose.csv` / `bdt_cuts_tight.csv` (`../bdt_cuts_<v>.csv`) instead of `bdt_cuts.csv` |
| `flow_Analysis_ingradients_BDT<v>.C` | `Quantiles/for_20Qbin/flow_Analysis_ingradients.C` | only the `BDTHandler` include path (→ the local copy above); quantile header + resolution-file path **unchanged/reused** (BDT doesn't affect either) |
| `fit_mass_and_flow_BDT<v>.C` | `flow_extraction/for_20Qbin/fit_mass_and_flow.C` | only `inputNtuple` (repointed at this variation's combined ntuple) and the output name tag; `Analysis_bin.h` include made absolute so it still resolves to the shared, unmodified nominal header (same `vnbinning`/quantile binning as data) |
| `scriptForSlurmJobs_optimized.py` | `Quantiles/for_20Qbin/` | `STEP`/`SOURCE_CODE_NAME`/`OUTPUT_BASE` renamed; **same nominal data file list**, same `N_FILES_PER_JOB=200` |
| `sbatch_mass_fit_flow_BDT<v>.sh` | `flow_extraction/for_20Qbin/sbatch_mass_fit_flow.sh` | `SRC`/`SCRATCH_BASE`/exec name renamed |

Both ingredients macros and both fit macros compile clean, and
`flow_Analysis_ingradients_BDT{loose,tight}.exe` were smoke-tested on 2 real
data files: loose keeps ~22% more candidates than tight there, in the
expected direction and a sensible magnitude for a 1%-significance-scale cut
shift.

**Scale warning**: the nominal data file list
(`Quantiles/for_20Qbin/inputFiles_PbPb2023_MB0to31_new.txt`) has ~212,000
files → ~1060 jobs at 200 files/job. This mirrors the nominal production
exactly (same file list, same job size) — it is **not** a lightweight step.
`IS_TEST_RUN = True` by default in both `scriptForSlurmJobs_optimized.py`
copies; flip to `False` deliberately.

Run, per variation (`loose` then `tight`, or in parallel):
```
cd loose   # or tight
python3 scriptForSlurmJobs_optimized.py        # stage 3: SP ingredients, BDT-varied
# when all ~1060 jobs finish:
hadd -f <OUT>/ROOT/flow_Analysis_out_combined.root <OUT>/ROOT/flow_Analysis_*_*.root
sbatch sbatch_mass_fit_flow_BDT<v>.sh           # stage 4: mass fit + yield-weighted vn (array 0-4)
```
`<OUT>` is `/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/BDT_sys_<v>_D0_SP_12NUQbin_diffq2q3_MB0to31_Sept10_v1`
(printed by the script). Stage 4's `SCRATCH_BASE` defaults to
`BDT_sys_<v>_D0_Flow_12Qbin_diffq2q3_output_Sept10_v1/output/`, producing
`Flow_output_BDT<v>_cen<0..4>.root` — same graph names/structure as the
nominal `Flow_output_cen<N>.root` from `flow_extraction/for_20Qbin/`.

### Step 3 — systematic (not automated; do per (cent,pT))
```
dvn     = vn_BDTloose/tight − vn_nominal
dvn/vn  = dvn / vn_nominal
```
using the matching graphs (`v2_graph_*`, `v3_graph_*`,
`v{2,3}_cen_inclusive*`, etc.) from the nominal `Flow_output_cen<N>.root`
(`flow_extraction/for_20Qbin/`) and from the two
`Flow_output_BDT{loose,tight}_cen<N>.root` files above — same envelope
convention as `../centrality_sys/`.

---

## Method A (earlier, exploratory): MC with/without-BDT closure

Measures prompt-D0 `v2(pT)`/`v3(pT)` on MC with the scalar-product method,
**with** and **without** the BDT cut, no mass fit — signal is simply
gen-matched, non-swap, prompt candidates
(`matchGEN==1 && isSwap==0 && Dgen_isPrompt==1`). Showed large statistical
fluctuation (see conversation), which is why Method B above was pursued
instead; kept here for reference.

### Files
| file | role |
|------|------|
| `Calculate_Resolution_MC.C` | stage 1 – SP 3-sub-event resolution ingredients `<Q_A·Q_B>`, q- and centrality-inclusive, one set of 6 weighted histograms. |
| `flow_MC_BDT_sys.C`         | stage 2 – reads the combined resolution file, one pass over the MC, fills `TProfile`s: `pv{2,3}_woBDT` (no cut) and `pv{2,3}_wBDT_<cen>` (per-class cut). No ntuple. |
| `make_BDT_systematic.C`     | stage 3 – reads the hadd-ed stage-2 output, builds `vn(pT)` graphs and the `dvn`/`dvn/vn` systematic graphs, prints a table. |
| `Analysis_bin_BDT.h`        | pT bin edges + centrality-class BDT-cut args + MC stitch-weight table + `USE_MC_WEIGHT` flag + default resolution-file path. |
| `BDTHandler.{h,cc}` / `bdt_cuts.csv` | local copy of the analysis BDT-cut lookup (nominal cuts; also read by `scan_significance_vs_BDT.C` above). |
| `scriptForSlurmJobs_optimized.py` | compiles + submits, switch `STAGE = "resolution"` / `"flow"`. |

### Method
Per event (weight `w = genWeight × stitchWeight`, `USE_MC_WEIGHT`):
* resolution (stage 1): `<Q2_HFm·Q2_HFp>`, `<Q2_HFm·Q2_Trk>`, `<Q2_HFp·Q2_Trk>`
  and the n=3 analogues; `TH1::StatOverflows(kTRUE)` for unbiased `GetMean()`.
* numerator (stage 2), per gen-matched prompt candidate, `pT>2`, `|y|<1`:
  `vn_obs = Re( u_n(D0) · Q*_n(HF_ref) )`, `HF_ref = HF-minus` for `y>0`,
  `HF-plus` for `y<0`.
```
res_plus  (y>0) = sqrt( <HFm·HFp> <HFm·Trk> / <HFp·Trk> )
res_minus (y<0) = sqrt( <HFm·HFp> <HFp·Trk> / <HFm·Trk> )
vn              = vn_obs / res(sign(y))
dvn     = vn_wBDT_<cen> - vn_woBDT
dvn/vn  = dvn / vn_woBDT
```
The BDT cut is applied the `../Prompt_frac/getDCA_fromMC.C` /
`../AccxEff/getEfficiency.C` way: candidate's true centrality ignored, cut
picked by analysis class (`cent_bdt_arg = {10,30,50,70,90}`) applied to the
full centrality-inclusive sample.

### Running
```
# STAGE = "resolution" in scriptForSlurmJobs_optimized.py
python3 scriptForSlurmJobs_optimized.py
hadd -f <OUT>/ROOT/Resolution_MC_combined.root <OUT>/ROOT/Resolution_MC_*_*.root
# STAGE = "flow"  (same DATE_tag)
python3 scriptForSlurmJobs_optimized.py
hadd -f <OUT2>/ROOT/flow_MC_BDT_out_combined.root <OUT2>/ROOT/flow_MC_BDT_sys_*_*.root
./make_BDT_systematic.exe <OUT2>/ROOT/flow_MC_BDT_out_combined.root BDT_sys_output.root
```
