#!/usr/bin/env python3
"""
Merge all partial counts_*.txt files produced by flow_Analysis_latest
and compute the final q2/q3 max/min ratio over the full sample.

Usage:
    python3 merge_counts.py /path/to/output_dir/LOG/
    python3 merge_counts.py /path/to/output_dir/LOG/ --root
"""

import sys, os, glob
import numpy as np

WRITE_ROOT = "--root" in sys.argv
try:
    import ROOT
    ROOT.gROOT.SetBatch(True)
except ImportError:
    WRITE_ROOT = False

# ── config ────────────────────────────────────────────────────────────────────
N_CENT  = 90
N_QBINS = 10
WARNING = 1.2

LOGDIR  = sys.argv[1] if len(sys.argv) > 1 else "."
OUTDIR  = LOGDIR

# ── accumulate counts from all partial files ──────────────────────────────────
total_q2 = np.zeros((N_CENT, N_QBINS), dtype=np.int64)
total_q3 = np.zeros((N_CENT, N_QBINS), dtype=np.int64)

files = sorted(glob.glob(os.path.join(LOGDIR, "counts_*.txt")))
if not files:
    print(f"ERROR: no counts_*.txt files found in {LOGDIR}")
    sys.exit(1)

print(f"Found {len(files)} partial count files.")
for f in files:
    with open(f) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            ic, iq, nq2, nq3 = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])
            total_q2[ic][iq] += nq2
            total_q3[ic][iq] += nq3

# ── compute per-cent ratios ────────────────────────────────────────────────────
SUMMARY_FILE = os.path.join(OUTDIR, "merged_ratio_summary.txt")

with open(SUMMARY_FILE, "w") as out:

    header = (f"\n{'cent':>6}  "
              + "  ".join(f"{'q2_iq'+str(iq):>10}" for iq in range(N_QBINS))
              + f"  {'q2_total':>12}  {'q2_max/min':>12}"
              + f"  {'q3_total':>12}  {'q3_max/min':>12}  {'warn':>8}\n")
    print(header)
    out.write(header)
    print("="*160)
    out.write("="*160 + "\n")

    # grand totals across all cent bins
    grand_q2 = np.zeros(N_QBINS, dtype=np.int64)
    grand_q3 = np.zeros(N_QBINS, dtype=np.int64)

    for ic in range(N_CENT):
        q2_row = total_q2[ic]
        q3_row = total_q3[ic]

        q2_sum = q2_row.sum()
        q3_sum = q3_row.sum()
        q2_max = q2_row.max(); q2_min = q2_row.min()
        q3_max = q3_row.max(); q3_min = q3_row.min()

        r2 = q2_max / q2_min if q2_min > 0 else float('nan')
        r3 = q3_max / q3_min if q3_min > 0 else float('nan')

        warn = ""
        if not np.isnan(r2) and r2 > WARNING: warn += "Q2! "
        if not np.isnan(r3) and r3 > WARNING: warn += "Q3!"

        q2_counts_str = "  ".join(f"{q2_row[iq]:>10}" for iq in range(N_QBINS))

        line = (f"{ic:>6}  {q2_counts_str}"
                f"  {q2_sum:>12}  {r2:>12.4f}"
                f"  {q3_sum:>12}  {r3:>12.4f}  {warn:>8}")
        print(line)
        out.write(line + "\n")

        grand_q2 += q2_row
        grand_q3 += q3_row

    # ── grand total row ───────────────────────────────────────────────────────
    print("="*160)
    out.write("="*160 + "\n")

    grand_q2_counts_str = "  ".join(f"{grand_q2[iq]:>10}" for iq in range(N_QBINS))
    grand_r2 = grand_q2.max() / grand_q2.min() if grand_q2.min() > 0 else float('nan')
    grand_r3 = grand_q3.max() / grand_q3.min() if grand_q3.min() > 0 else float('nan')

    grand_line = (f"{'ALL':>6}  {grand_q2_counts_str}"
                  f"  {grand_q2.sum():>12}  {grand_r2:>12.4f}"
                  f"  {grand_q3.sum():>12}  {grand_r3:>12.4f}")
    print(grand_line)
    out.write(grand_line + "\n\n")

    # ── extremes summary ──────────────────────────────────────────────────────
    def print_extremes(total_arr, label):
        ratios = np.array([
            total_arr[ic].max() / total_arr[ic].min()
            if total_arr[ic].min() > 0 else float('nan')
            for ic in range(N_CENT)
        ])
        valid  = ~np.isnan(ratios)
        imax   = np.nanargmax(ratios)
        imin   = np.nanargmin(ratios)
        warn_n = np.sum(ratios[valid] > WARNING)
        msg = (f"\n  [{label}]\n"
               f"    Overall MAX ratio = {ratios[imax]:.4f}  at cent1pct = {imax}\n"
               f"    Overall MIN ratio = {ratios[imin]:.4f}  at cent1pct = {imin}\n"
               f"    Bins with ratio > {WARNING}: {warn_n} / {N_CENT}\n"
               f"    Top-10 worst:\n")
        top10_idx = np.argsort(ratios[valid])[::-1][:10]
        valid_idx = np.where(valid)[0]
        for k, i in enumerate(top10_idx):
            ic = valid_idx[i]
            msg += (f"      #{k+1:2d}  cent={ic:3d}%  "
                    f"ratio={ratios[ic]:.4f}  "
                    f"max={total_arr[ic].max()}  "
                    f"min={total_arr[ic].min()}\n")
        print(msg)
        out.write(msg)

    print_extremes(total_q2, "q2 quantile balance (full sample)")
    print_extremes(total_q3, "q3 quantile balance (full sample)")

print(f"\n==> Summary written to: {SUMMARY_FILE}")

# ── optional ROOT output ───────────────────────────────────────────────────────
if WRITE_ROOT:
    ROOTFILE = os.path.join(OUTDIR, "merged_counts.root")
    rf = ROOT.TFile(ROOTFILE, "RECREATE")

    hq2 = ROOT.TH2D("h2_q2_counts",
                     "Q2 event counts (full sample);q-bin;Cent 1%",
                     N_QBINS, 0, N_QBINS, N_CENT, 0, N_CENT)
    hq3 = ROOT.TH2D("h2_q3_counts",
                     "Q3 event counts (full sample);q-bin;Cent 1%",
                     N_QBINS, 0, N_QBINS, N_CENT, 0, N_CENT)
    hr2 = ROOT.TH1D("h_q2_ratio", "Q2 max/min ratio (full sample);Cent 1%;max/min ratio",
                    N_CENT, 0, N_CENT)
    hr3 = ROOT.TH1D("h_q3_ratio", "Q3 max/min ratio (full sample);Cent 1%;max/min ratio",
                    N_CENT, 0, N_CENT)

    for ic in range(N_CENT):
        for iq in range(N_QBINS):
            hq2.SetBinContent(iq+1, ic+1, float(total_q2[ic][iq]))
            hq3.SetBinContent(iq+1, ic+1, float(total_q3[ic][iq]))
        r2 = total_q2[ic].max()/total_q2[ic].min() if total_q2[ic].min() > 0 else 0
        r3 = total_q3[ic].max()/total_q3[ic].min() if total_q3[ic].min() > 0 else 0
        hr2.SetBinContent(ic+1, r2)
        hr3.SetBinContent(ic+1, r3)

    rf.Write()
    rf.Close()
    print(f"==> ROOT file written to: {ROOTFILE}")
    print("    Open with:  root -l merged_counts.root  then  new TBrowser")
