#!/usr/bin/env python3
"""
analyze_bad_bins.py
Reads sigma_bad_bins_cent*.log and summarizes bad SP-bins
in the mid v-range (signal-contributing region).

Mid ranges (from Analysis_bin.h):
  v2 : i_v in [11, 30]  (tail: 0-10 and 31-41, N_SIDE_EDGES_V2-1=11 per side)
  v3 : i_v in [14, 33]  (tail: 0-13 and 34-47, N_SIDE_EDGES_V3-1=14 per side)

Log line format:
  type | i_cen | cen_name | i_q | i_pt | pt_name | i_v | reason

Usage:
  python3 analyze_bad_bins.py                        # reads logs in current dir
  python3 analyze_bad_bins.py --logdir /path/to/logs
  python3 analyze_bad_bins.py --all-reasons          # include low_stats too
  python3 analyze_bad_bins.py --show-tail            # also dump tail-bin summary
"""

import glob, sys, os, argparse
from collections import defaultdict

# ---- bin boundaries (must match Analysis_bin.h) ----
V2_MID_LO,  V2_MID_HI  = 11, 30   # inclusive
V3_MID_LO,  V3_MID_HI  = 14, 33
N_VBINS_V2, N_VBINS_V3 = 42, 48

PT_LABELS = [
    "1-2",  "2-3",  "3-4",  "4-5",
    "5-6",  "6-8",  "8-10", "10-15",
    "15-20","20-40","40-60","60-100"
]

def is_mid(vtype, iv):
    return (V2_MID_LO <= iv <= V2_MID_HI) if vtype == "v2" \
      else (V3_MID_LO <= iv <= V3_MID_HI)

# ---- parse ----
def parse_logs(log_dir):
    files = sorted(glob.glob(os.path.join(log_dir, "sigma_bad_bins_cent*.log")))
    if not files:
        sys.exit(f"[ERROR] No sigma_bad_bins_cent*.log files found in '{log_dir}'")
    entries = []
    for fpath in files:
        with open(fpath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = [p.strip() for p in line.split("|")]
                if len(parts) < 8:
                    continue
                try:
                    vtype    = parts[0]
                    i_cen    = int(parts[1])
                    cen_name = parts[2]
                    i_q      = int(parts[3])
                    i_pt     = int(parts[4])
                    pt_name  = parts[5]
                    i_v      = int(parts[6])
                    reason   = parts[7]
                except (ValueError, IndexError):
                    continue
                entries.append(dict(vtype=vtype, i_cen=i_cen, cen_name=cen_name,
                                    i_q=i_q, i_pt=i_pt, pt_name=pt_name,
                                    i_v=i_v, reason=reason,
                                    mid=is_mid(vtype, i_v)))
    print(f"Parsed {len(entries)} total bad-bin entries from {len(files)} file(s).\n")
    return entries

# ---- main report ----
def print_report(entries, skip_low_stats, show_tail):
    reasons_shown = "excluding low_stats" if skip_low_stats else "all reasons"

    for region, flag in [("MID (signal-contributing)", True),
                         ("TAIL (low-weight, informational)", False)]:
        if flag is False and not show_tail:
            continue

        subset = [e for e in entries if e["mid"] == flag]
        if skip_low_stats:
            subset = [e for e in subset if e["reason"] != "low_stats"]

        print("=" * 72)
        print(f"  {region}  [{reasons_shown}]  —  {len(subset)} bad bins")
        print("=" * 72)

        if not subset:
            print("  (none)\n")
            continue

        # --- group by (vtype, i_cen, i_q, i_pt) ---
        groups = defaultdict(list)
        for e in subset:
            key = (e["vtype"], e["i_cen"], e["cen_name"], e["i_q"], e["i_pt"], e["pt_name"])
            groups[key].append((e["i_v"], e["reason"]))

        sorted_keys = sorted(groups, key=lambda k: (k[1], k[0], k[3], k[4]))

        current_cen = None
        for key in sorted_keys:
            vtype, i_cen, cen_name, i_q, i_pt, pt_name = key
            if cen_name != current_cen:
                if current_cen is not None:
                    print()
                print(f"\n  --- {cen_name} ---")
                print(f"  {'TYPE':<4}  {'q-bin':>5}  {'pT (GeV)':>10}  bad i_v bins")
                print(f"  {'-'*4}  {'-'*5}  {'-'*10}  {'-'*40}")
                current_cen = cen_name

            bad = sorted(groups[key])
            pt_range = PT_LABELS[i_pt] if i_pt < len(PT_LABELS) else pt_name
            # compact: iv[reason_initial]
            iv_str = "  ".join(f"{iv}[{r[0]}]" for iv, r in bad)
            print(f"  {vtype:<4}  q={i_q:<3}  {pt_range:>10}  {iv_str}")

        # --- count summary table ---
        print(f"\n\n  Summary counts — {region}")
        print(f"  {'Centrality':<15}  {'v2 bad':>7}  {'v3 bad':>7}  {'total':>7}")
        print(f"  {'-'*15}  {'-'*7}  {'-'*7}  {'-'*7}")
        cen_order = sorted(set((e["i_cen"], e["cen_name"]) for e in subset))
        grand_v2 = grand_v3 = 0
        for i_cen, cen_name in cen_order:
            n_v2 = sum(len(v) for k, v in groups.items() if k[0] == "v2" and k[1] == i_cen)
            n_v3 = sum(len(v) for k, v in groups.items() if k[0] == "v3" and k[1] == i_cen)
            print(f"  {cen_name:<15}  {n_v2:>7}  {n_v3:>7}  {n_v2+n_v3:>7}")
            grand_v2 += n_v2
            grand_v3 += n_v3
        print(f"  {'TOTAL':<15}  {grand_v2:>7}  {grand_v3:>7}  {grand_v2+grand_v3:>7}")

        # --- most problematic pT bins ---
        pt_counts = defaultdict(int)
        for e in subset:
            pt_counts[(e["i_pt"], e["pt_name"])] += 1
        print(f"\n  Top pT bins by bad-bin count:")
        for (i_pt, pt_name), cnt in sorted(pt_counts.items(), key=lambda x: -x[1])[:6]:
            pt_range = PT_LABELS[i_pt] if i_pt < len(PT_LABELS) else pt_name
            print(f"    pT {pt_range:>10} GeV  ->  {cnt} bad bins")

        print()

# ---- reason key ----
def print_reason_key():
    print("\nReason codes:")
    print("  n = neg_yield  : fit converged but signal yield <= 0")
    print("  z = zero_error : fit error on yield is 0 (fit likely failed)")
    print("  l = low_stats  : histogram had < 10 entries (skipped by default)\n")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Summarize bad sigma bins from fit_mass_and_flow log files.")
    ap.add_argument("--logdir",      default=".", help="Directory containing sigma_bad_bins_cent*.log files (default: current dir)")
    ap.add_argument("--all-reasons", action="store_true", help="Include low_stats entries (excluded by default)")
    ap.add_argument("--show-tail",   action="store_true", help="Also print tail-bin summary (informational)")
    args = ap.parse_args()

    entries = parse_logs(args.logdir)
    print_report(entries, skip_low_stats=not args.all_reasons, show_tail=args.show_tail)
    print_reason_key()
