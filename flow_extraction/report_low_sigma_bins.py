#!/usr/bin/env python3
"""
report_low_sigma_bins.py
------------------------
Reads all sigma_bad_bins_*.log files from the analysis output directory and
produces a summary of SP bins rejected by the sigma > 1 criterion ("low_sigma").

Usage:
    python3 report_low_sigma_bins.py <output_dir>

Example:
    python3 report_low_sigma_bins.py \
        /scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Flow_output_FullStat_May19_v1/output
"""

import sys
import os
import glob
from collections import defaultdict

# ---------------------------------------------------------------------------
# Total SP bins per (flow_type, pt) — from Analysis_bin.h
# Used to compute the rejection fraction.
# ---------------------------------------------------------------------------
N_VBINS_V2 = 42
N_VBINS_V3 = 48
N_QBINS    = 10   # Q-bins per (cent, pT) group

# ---------------------------------------------------------------------------

def parse_logs(log_dir):
    """Return list of dicts, one per low_sigma line across all log files."""
    pattern = os.path.join(log_dir, "sigma_bad_bins_*.log")
    log_files = sorted(glob.glob(pattern))
    if not log_files:
        print(f"ERROR: No log files matching '{pattern}' found.")
        sys.exit(1)
    print(f"Found {len(log_files)} log file(s): {[os.path.basename(f) for f in log_files]}\n")

    records = []
    for fpath in log_files:
        with open(fpath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = [p.strip() for p in line.split("|")]
                if len(parts) < 8:
                    continue
                reason = parts[7]
                if reason != "low_sigma":
                    continue
                try:
                    rec = {
                        "type":     parts[0],          # v2 / v3
                        "i_cen":    int(parts[1]),
                        "cen_name": parts[2],
                        "i_q":      int(parts[3]),
                        "i_pt":     int(parts[4]),
                        "pt_name":  parts[5],
                        "i_v":      int(parts[6]),      # SP bin index
                        "mid_tail": parts[8] if len(parts) > 8 else "?",
                        "sigma":    float(parts[9].replace("sigma=", "")) if len(parts) > 9 else float("nan"),
                    }
                    records.append(rec)
                except (ValueError, IndexError):
                    continue
    return records


def make_report(records):
    # -----------------------------------------------------------------------
    # 1. Per (type, cen_name, pt_name) summary
    # -----------------------------------------------------------------------
    # key -> list of sigma values
    summary = defaultdict(list)
    # key -> set of SP bin indices that were rejected
    sp_bins = defaultdict(set)
    # key -> dict of {i_q: count}
    by_qbin = defaultdict(lambda: defaultdict(int))

    for r in records:
        key = (r["type"], r["i_cen"], r["cen_name"], r["i_pt"], r["pt_name"])
        summary[key].append(r["sigma"])
        sp_bins[key].add(r["i_v"])
        by_qbin[key][r["i_q"]] += 1

    # Sort keys: v2 before v3, then by cen index, then by pT index (low → high)
    sorted_keys = sorted(summary.keys(), key=lambda k: (k[0], k[1], k[3]))

    # -----------------------------------------------------------------------
    # 2. Print report
    # -----------------------------------------------------------------------
    SEP = "=" * 110

    print(SEP)
    print("  LOW-SIGMA REJECTION REPORT  (sigma <= 1 criterion)")
    print("  One entry per (flow_type, centrality, pT bin)")
    print("  'unique SP bins' = distinct SP bin indices rejected across ALL q-bins")
    print("  'fraction'       = unique SP bins / total SP bins in that flow type")
    print(SEP)

    current_type = None
    for key in sorted_keys:
        vtype, i_cen, cen, i_pt, pt = key
        sigmas = summary[key]
        n_total_vbins = N_VBINS_V2 if vtype == "v2" else N_VBINS_V3
        unique_sp = sp_bins[key]
        n_unique   = len(unique_sp)
        fraction   = n_unique / n_total_vbins

        sigma_min  = min(sigmas)
        sigma_max  = max(sigmas)
        sigma_mean = sum(sigmas) / len(sigmas)

        # Q-bin breakdown: which q-bins contribute most rejections
        qbin_counts = by_qbin[key]
        qbin_str = "  ".join(f"q{iq}:{cnt}" for iq, cnt in sorted(qbin_counts.items()))

        # Header when flow type changes
        if vtype != current_type:
            current_type = vtype
            print(f"\n{'─'*110}")
            print(f"  {'[[ ' + vtype.upper() + ' ]]':^106}")
            print(f"{'─'*110}")
            print(f"  {'Centrality':<14} {'pT bin':<14} {'unique SP bins':>14} {'/ total':>8} "
                  f"{'fraction':>9}   {'sigma min':>10} {'sigma mean':>11} {'sigma max':>10}")
            print(f"  {'-'*106}")

        print(f"  {cen:<14} {pt:<14} {n_unique:>14d} {n_total_vbins:>8d} "
              f"{fraction:>9.1%}   {sigma_min:>10.3f} {sigma_mean:>11.3f} {sigma_max:>10.3f}")
        print(f"    └─ q-bin breakdown: {qbin_str}")
        print(f"    └─ rejected SP bin indices: {sorted(unique_sp)}")

    print(f"\n{SEP}")

    # -----------------------------------------------------------------------
    # 3. Worst-offender summary: (cen, pt) pairs with > 50% rejection
    # -----------------------------------------------------------------------
    print("\n  FLAGGED: groups with > 50% of SP bins rejected (likely low statistics)\n")
    flagged = False
    for key in sorted_keys:
        vtype, i_cen, cen, i_pt, pt = key
        unique_sp     = sp_bins[key]
        n_total_vbins = N_VBINS_V2 if vtype == "v2" else N_VBINS_V3
        fraction      = len(unique_sp) / n_total_vbins
        if fraction > 0.5:
            sigmas = summary[key]
            print(f"  {vtype.upper():4s}  {cen:<14} {pt:<14}  {len(unique_sp):>3d}/{n_total_vbins} bins "
                  f"rejected ({fraction:.0%})  sigma range [{min(sigmas):.3f}, {max(sigmas):.3f}]")
            flagged = True
    if not flagged:
        print("  None — all (cent, pT) groups have <= 50% low-sigma rejection.")

    print(f"\n{SEP}\n")

    # -----------------------------------------------------------------------
    # 4. Overall statistics
    # -----------------------------------------------------------------------
    # Identify the last pT bin index for v2 and v3 (highest i_pt seen)
    last_ipt_v2 = max((k[3] for k in sorted_keys if k[0] == "v2"), default=None)
    last_ipt_v3 = max((k[3] for k in sorted_keys if k[0] == "v3"), default=None)
    # Retrieve the pt_name for the last bin (for display)
    last_pt_name_v2 = next((k[4] for k in sorted_keys if k[0] == "v2" and k[3] == last_ipt_v2), "?")
    last_pt_name_v3 = next((k[4] for k in sorted_keys if k[0] == "v3" and k[3] == last_ipt_v3), "?")

    def count_and_possible(vtype, exclude_last_ipt=None):
        """Return (n_rejections, n_possible_slots, n_groups) for the given filters."""
        n_vbins = N_VBINS_V2 if vtype == "v2" else N_VBINS_V3
        groups = set(
            k for k in sorted_keys
            if k[0] == vtype and (exclude_last_ipt is None or k[3] != exclude_last_ipt)
        )
        n_rej = sum(
            1 for r in records
            if r["type"] == vtype and (exclude_last_ipt is None or r["i_pt"] != exclude_last_ipt)
        )
        possible = len(groups) * N_QBINS * n_vbins
        return n_rej, possible, len(groups)

    tv2_all,  pv2_all,  gv2_all  = count_and_possible("v2")
    tv3_all,  pv3_all,  gv3_all  = count_and_possible("v3")
    tv2_excl, pv2_excl, gv2_excl = count_and_possible("v2", last_ipt_v2)
    tv3_excl, pv3_excl, gv3_excl = count_and_possible("v3", last_ipt_v3)

    def pct(n, d):
        return n / d * 100 if d else 0.0

    print(f"  Total low_sigma rejections  (incl. last pT bin):")
    print(f"    V2 : {tv2_all:>6d} / {pv2_all:>6d}  ({pct(tv2_all,pv2_all):.1f}%)   "
          f"[{gv2_all} cent-pT groups × {N_QBINS} q-bins × {N_VBINS_V2} SP bins]")
    print(f"    V3 : {tv3_all:>6d} / {pv3_all:>6d}  ({pct(tv3_all,pv3_all):.1f}%)   "
          f"[{gv3_all} cent-pT groups × {N_QBINS} q-bins × {N_VBINS_V3} SP bins]")
    print(f"    All: {tv2_all+tv3_all:>6d} / {pv2_all+pv3_all:>6d}  ({pct(tv2_all+tv3_all, pv2_all+pv3_all):.1f}%)")

    print(f"\n  Total low_sigma rejections  (excl. last pT bin — V2: {last_pt_name_v2}, V3: {last_pt_name_v3}):")
    print(f"    V2 : {tv2_excl:>6d} / {pv2_excl:>6d}  ({pct(tv2_excl,pv2_excl):.1f}%)   "
          f"[{gv2_excl} cent-pT groups × {N_QBINS} q-bins × {N_VBINS_V2} SP bins]")
    print(f"    V3 : {tv3_excl:>6d} / {pv3_excl:>6d}  ({pct(tv3_excl,pv3_excl):.1f}%)   "
          f"[{gv3_excl} cent-pT groups × {N_QBINS} q-bins × {N_VBINS_V3} SP bins]")
    print(f"    All: {tv2_excl+tv3_excl:>6d} / {pv2_excl+pv3_excl:>6d}  ({pct(tv2_excl+tv3_excl, pv2_excl+pv3_excl):.1f}%)")

    if records:
        all_sigma = [r["sigma"] for r in records]
        print(f"\n  Sigma value range of rejected bins: [{min(all_sigma):.4f}, {max(all_sigma):.4f}]")
        print(f"  Mean sigma of rejected bins:         {sum(all_sigma)/len(all_sigma):.4f}")
    print()


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    log_dir = sys.argv[1]
    if not os.path.isdir(log_dir):
        print(f"ERROR: Directory not found: {log_dir}")
        sys.exit(1)

    records = parse_logs(log_dir)
    if not records:
        print("No 'low_sigma' entries found in any log file.")
        sys.exit(0)
    print(f"Parsed {len(records)} low_sigma rejection entries.\n")
    make_report(records)


if __name__ == "__main__":
    main()
