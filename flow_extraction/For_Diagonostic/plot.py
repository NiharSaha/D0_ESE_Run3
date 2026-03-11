#!/usr/bin/env python3
"""
Parse diag_logs/log_cen{N}_pt{M}.txt produced by run_all_bins.sh
and save 2D/1D ratio plots into a single ROOT file.

Usage:
    python3 plot_ratio_maps.py                        # uses ./diag_logs/
    python3 plot_ratio_maps.py /path/to/diag_logs/   # custom log dir
"""

import sys
import os
import re
import numpy as np

# ---------------------------------------------------------------
# Try to import ROOT; fall back to matplotlib+uproot if unavailable
# ---------------------------------------------------------------
try:
    import ROOT
    ROOT.gROOT.SetBatch(True)
    HAS_ROOT = True
except ImportError:
    HAS_ROOT = False
    print("WARNING: ROOT not found, falling back to PNG output.")
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------
N_CENT  = 90
N_PT    = 10
PT_LABELS   = ["1-2","2-3","3-4","4-5","5-6","6-8","8-10","10-15","15-20","20-40"]
PT_LABELS_AXIS = ["1-2 GeV","2-3","3-4","4-5","5-6","6-8","8-10","10-15","15-20","20-40"]
WARNING_THRESHOLD = 1.2

LOGDIR  = sys.argv[1] if len(sys.argv) > 1 else "./diag_logs"
OUTDIR  = os.path.dirname(os.path.abspath(__file__))
OUTROOT = os.path.join(OUTDIR, "ratio_maps.root")

# ---------------------------------------------------------------
# Parse all log files -> fill 2D numpy arrays
# ---------------------------------------------------------------
q2_ratio = np.full((N_CENT, N_PT), -1.0)   # -1 = missing
q3_ratio = np.full((N_CENT, N_PT), -1.0)
q2_total = np.full((N_CENT, N_PT), -1.0)
q3_total = np.full((N_CENT, N_PT), -1.0)

parsed  = 0
missing = []

for icen in range(N_CENT):
    for ipt in range(N_PT):
        logfile = os.path.join(LOGDIR, f"log_cen{icen}_pt{ipt}.txt")
        if not os.path.isfile(logfile):
            missing.append((icen, ipt))
            continue
        with open(logfile) as f:
            text = f.read()

        m2 = re.search(r'q2 cand max/min ratio:\s*([\d.]+)', text)
        m3 = re.search(r'q3 cand max/min ratio:\s*([\d.]+)', text)
        mt = re.search(r'TOTAL\s+([\d.]+)\s+([\d.]+)', text)

        if m2: q2_ratio[icen, ipt] = float(m2.group(1))
        if m3: q3_ratio[icen, ipt] = float(m3.group(1))
        if mt:
            q2_total[icen, ipt] = float(mt.group(1))
            q3_total[icen, ipt] = float(mt.group(2))
        parsed += 1

print(f"Parsed {parsed} log files.  Missing: {len(missing)}")

# ---------------------------------------------------------------
# Helpers to print text summary
# ---------------------------------------------------------------
def print_extremes(data, label):
    valid = data[data >= 0]
    if len(valid) == 0:
        print(f"  No valid data for {label}"); return
    icens, ipts = np.where(data >= 0)
    vals = data[icens, ipts]
    idx_max = np.argmax(vals); idx_min = np.argmin(vals)
    warn_n  = np.sum(vals > WARNING_THRESHOLD)
    print(f"\n  [{label}]")
    print(f"    MAX = {vals[idx_max]:.4f}  cen1pct={icens[idx_max]}  pT={PT_LABELS[ipts[idx_max]]}")
    print(f"    MIN = {vals[idx_min]:.4f}  cen1pct={icens[idx_min]}  pT={PT_LABELS[ipts[idx_min]]}")
    print(f"    Bins with ratio > {WARNING_THRESHOLD}: {warn_n}")
    print(f"    Top-10 highest:")
    for k, i in enumerate(np.argsort(vals)[::-1][:10]):
        print(f"      #{k+1:2d}  ratio={vals[i]:.4f}  cen1pct={icens[i]:3d}  pT={PT_LABELS[ipts[i]]}")

print("\n=== RATIO EXTREMES ===")
print_extremes(q2_ratio, "q2 cand max/min ratio")
print_extremes(q3_ratio, "q3 cand max/min ratio")

# ---------------------------------------------------------------
# ==================  ROOT output  ==============================
# ---------------------------------------------------------------
if HAS_ROOT:
    outf = ROOT.TFile(OUTROOT, "RECREATE")

    # ------------------------------------------------------------------
    # Helper: fill a TH2D from a numpy array (x=pT, y=cent1pct)
    # ------------------------------------------------------------------
    def make_th2(arr, name, title):
        """
        arr shape: (N_CENT, N_PT)
        TH2D: x-axis = pT bin index (0..N_PT-1),  y-axis = cent 1% (0..N_CENT-1)
        Bins with arr==-1 are left at 0 (unfilled).
        """
        h = ROOT.TH2D(name, title,
                      N_PT,  0, N_PT,
                      N_CENT, 0, N_CENT)
        h.GetXaxis().SetTitle("p_{T} bin")
        h.GetYaxis().SetTitle("Centrality 1% bin")
        h.GetZaxis().SetTitle("max/min ratio")
        for ipt in range(N_PT):
            h.GetXaxis().SetBinLabel(ipt+1, PT_LABELS[ipt])
        for ic in range(N_CENT):
            h.GetYaxis().SetBinLabel(ic+1, str(ic))
        for ic in range(N_CENT):
            for ipt in range(N_PT):
                if arr[ic, ipt] >= 0:
                    h.SetBinContent(ipt+1, ic+1, arr[ic, ipt])
        return h

    # ------------------------------------------------------------------
    # Helper: fill a TH2D marking WARNING bins (ratio > threshold -> 1)
    # ------------------------------------------------------------------
    def make_th2_warn(arr, name, title, thr=WARNING_THRESHOLD):
        h = ROOT.TH2D(name, title,
                      N_PT,  0, N_PT,
                      N_CENT, 0, N_CENT)
        h.GetXaxis().SetTitle("p_{T} bin")
        h.GetYaxis().SetTitle("Centrality 1% bin")
        h.GetZaxis().SetTitle(f"ratio > {thr}")
        for ipt in range(N_PT):
            h.GetXaxis().SetBinLabel(ipt+1, PT_LABELS[ipt])
        for ic in range(N_CENT):
            h.GetYaxis().SetBinLabel(ic+1, str(ic))
        for ic in range(N_CENT):
            for ipt in range(N_PT):
                if arr[ic, ipt] > thr:
                    h.SetBinContent(ipt+1, ic+1, 1.0)
        return h

    # ------------------------------------------------------------------
    # Helper: 1D TH1D profile — ratio vs cent for one pT bin
    # ------------------------------------------------------------------
    def make_th1_vs_cent(arr, ipt, name, title):
        h = ROOT.TH1D(name, title, N_CENT, 0, N_CENT)
        h.GetXaxis().SetTitle("Centrality 1% bin")
        h.GetYaxis().SetTitle("max/min ratio")
        for ic in range(N_CENT):
            if arr[ic, ipt] >= 0:
                h.SetBinContent(ic+1, arr[ic, ipt])
        return h

    # ------------------------------------------------------------------
    # Helper: 1D TH1D profile — ratio vs pT for one cent 1% bin
    # ------------------------------------------------------------------
    def make_th1_vs_pt(arr, icen, name, title):
        h = ROOT.TH1D(name, title, N_PT, 0, N_PT)
        for ipt in range(N_PT):
            h.GetXaxis().SetBinLabel(ipt+1, PT_LABELS[ipt])
        h.GetXaxis().SetTitle("p_{T} bin")
        h.GetYaxis().SetTitle("max/min ratio")
        for ipt in range(N_PT):
            if arr[icen, ipt] >= 0:
                h.SetBinContent(ipt+1, arr[icen, ipt])
        return h

    # ------------------------------------------------------------------
    # Helper: TH1D profile — mean ratio vs pT for a centrality group
    # ------------------------------------------------------------------
    def make_th1_centgroup_vs_pt(arr, clo, chi, name, title):
        h = ROOT.TH1D(name, title, N_PT, 0, N_PT)
        for ipt in range(N_PT):
            h.GetXaxis().SetBinLabel(ipt+1, PT_LABELS[ipt])
        h.GetXaxis().SetTitle("p_{T} bin")
        h.GetYaxis().SetTitle("mean max/min ratio")
        for ipt in range(N_PT):
            slc = arr[clo:chi, ipt]
            slc = slc[slc >= 0]
            if len(slc): h.SetBinContent(ipt+1, float(np.mean(slc)))
        return h

    # ------------------------------------------------------------------
    # Write everything into directories inside the ROOT file
    # ------------------------------------------------------------------

    # ---- 1. 2D maps ----
    dir2d = outf.mkdir("2D_maps")
    dir2d.cd()

    hq2_map = make_th2(q2_ratio, "h2_q2_ratio_map",
                       "Q2 cand max/min ratio;p_{T} bin;Centrality 1%")
    hq3_map = make_th2(q3_ratio, "h2_q3_ratio_map",
                       "Q3 cand max/min ratio;p_{T} bin;Centrality 1%")
    hq2_warn = make_th2_warn(q2_ratio, "h2_q2_warning_map",
                             f"Q2 WARNING bins (ratio>{WARNING_THRESHOLD})")
    hq3_warn = make_th2_warn(q3_ratio, "h2_q3_warning_map",
                             f"Q3 WARNING bins (ratio>{WARNING_THRESHOLD})")

    hq2_map.Write()
    hq3_map.Write()
    hq2_warn.Write()
    hq3_warn.Write()

    # ---- 2. 1D: ratio vs cent, one histogram per pT bin ----
    dir_vs_cent_q2 = outf.mkdir("q2_ratio_vs_cent1pct")
    dir_vs_cent_q3 = outf.mkdir("q3_ratio_vs_cent1pct")

    for ipt in range(N_PT):
        dir_vs_cent_q2.cd()
        h = make_th1_vs_cent(q2_ratio, ipt,
                             f"h_q2_ratio_vs_cent_pt{ipt}",
                             f"Q2 ratio vs cent 1%  pT={PT_LABELS[ipt]};Cent 1%;max/min ratio")
        h.Write()

        dir_vs_cent_q3.cd()
        h = make_th1_vs_cent(q3_ratio, ipt,
                             f"h_q3_ratio_vs_cent_pt{ipt}",
                             f"Q3 ratio vs cent 1%  pT={PT_LABELS[ipt]};Cent 1%;max/min ratio")
        h.Write()

    # ---- 3. 1D: ratio vs pT, one histogram per cent 1% bin ----
    dir_vs_pt_q2 = outf.mkdir("q2_ratio_vs_pt_per_cent1pct")
    dir_vs_pt_q3 = outf.mkdir("q3_ratio_vs_pt_per_cent1pct")

    for icen in range(N_CENT):
        dir_vs_pt_q2.cd()
        h = make_th1_vs_pt(q2_ratio, icen,
                           f"h_q2_ratio_vs_pt_cen{icen}",
                           f"Q2 ratio vs pT  cent1pct={icen};p_{{T}} bin;max/min ratio")
        h.Write()

        dir_vs_pt_q3.cd()
        h = make_th1_vs_pt(q3_ratio, icen,
                           f"h_q3_ratio_vs_pt_cen{icen}",
                           f"Q3 ratio vs pT  cent1pct={icen};p_{{T}} bin;max/min ratio")
        h.Write()

    # ---- 4. 1D: mean ratio vs pT per cent group ----
    dir_grp = outf.mkdir("ratio_vs_pt_centgroups")
    dir_grp.cd()
    cent_groups = [(0,10,"0to10"), (10,30,"10to30"), (30,50,"30to50"), (50,90,"50to90")]
    for (clo, chi, clabel) in cent_groups:
        h = make_th1_centgroup_vs_pt(q2_ratio, clo, chi,
                                     f"h_q2_ratio_vs_pt_cent{clabel}",
                                     f"Q2 mean ratio vs pT  cent {clabel}%;p_{{T}};mean max/min ratio")
        h.Write()
        h = make_th1_centgroup_vs_pt(q3_ratio, clo, chi,
                                     f"h_q3_ratio_vs_pt_cent{clabel}",
                                     f"Q3 mean ratio vs pT  cent {clabel}%;p_{{T}};mean max/min ratio")
        h.Write()

    # ---- 5. Total candidate counts (2D maps) ----
    dir_cand = outf.mkdir("candidate_totals")
    dir_cand.cd()
    hq2_tot = make_th2(q2_total, "h2_q2_total_cands",
                       "Q2 total candidates;p_{T} bin;Centrality 1%")
    hq3_tot = make_th2(q3_total, "h2_q3_total_cands",
                       "Q3 total candidates;p_{T} bin;Centrality 1%")
    hq2_tot.Write()
    hq3_tot.Write()

    outf.Close()
    print(f"\n==> ROOT file saved: {OUTROOT}")
    print("    Directories inside:")
    print("      2D_maps/                    -- TH2D colour maps")
    print("      q2_ratio_vs_cent1pct/       -- TH1D per pT bin")
    print("      q3_ratio_vs_cent1pct/       -- TH1D per pT bin")
    print("      q2_ratio_vs_pt_per_cent1pct/-- TH1D per cent 1% bin")
    print("      q3_ratio_vs_pt_per_cent1pct/-- TH1D per cent 1% bin")
    print("      ratio_vs_pt_centgroups/     -- TH1D mean per cent group")
    print("      candidate_totals/           -- TH2D total candidate counts")
    print("\n    Open with:  root -l ratio_maps.root")
    print("    Then:  new TBrowser")
