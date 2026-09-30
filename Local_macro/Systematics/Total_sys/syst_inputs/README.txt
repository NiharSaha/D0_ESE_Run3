Put here:
  - the nominal pair: D0_Flow_*.root + Charge_flow_*.root (the same two files you'd pass
    to plot_ESE_scatter_and_pearson_combined.C as d0_file / chg_file for the nominal result).
  - for each of the 7 systematic sources, either ONE (D0, charged-particle) pair, or TWO
    such pairs if that source was evaluated with two variations (e.g. "Loose"/"Tight" or
    "Down"/"Up") — up to 15 pairs total.

Each file must have the exact same internal layout the nominal macro expects:
  D0 file:  directory "vn_vs_qbin" with graphs "v{2,3}_vs_q{2,3}bin_<cent>_pT<lo>to<hi>"
  Chg file: directory "vsQbin_TProfile" with profiles "hp_v{2,3}_vsq{2,3}_<pt_tag>_<cent>"
for <cent> in {cent0to10, cent10to20, cent20to30, cent30to40, cent40to50}.

Once your files are here, open ../plot_ESE_systematics.C and edit the "USER CONFIGURATION"
block right at the top of the file (just after the #include lines) to point nominal_d0_file /
nominal_chg_file / source_d0_file[][2] / source_chg_file[][2] / source_var_label[][2] /
source_label[] at your actual filenames and systematic-source names (e.g. "AccxEff", "BDT",
"PromptFrac", ...). The same block also has the 6 on/off switches, isRun3, the thresholds and
N_MC — everything you'd want to change lives there, so after editing it you just run:
  root -l -b -q 'plot_ESE_systematics.C()'
with no arguments needed.

For a source with only one variation, leave its [1] slot (both D0 and charged-particle) as
nullptr. For a source with two variations, fill both [0] and [1] with that variation's pair,
and set source_var_label[isrc] to e.g. {"Down","Up"} or {"Loose","Tight"} — at every
(harmonic, centrality, pT) bin the macro then compares BOTH variations to nominal and adopts
whichever deviates more (which side wins can differ from bin to bin) as that source's
systematic contribution.

This macro re-derives the Pearson r / slope / intercept for the nominal dataset and every
raw-file variation itself (same scatter-building, MC-resampling and linear-fit logic as
plot_ESE_scatter_and_pearson_combined.C) — it does not need any pre-computed
pearson_ESE_*.root file. Note: with N two-variation sources, the number of MC-resampling
computations grows accordingly (1 nominal + one-variation sources + 2x two-variation
sources), so runtime scales up correspondingly.
