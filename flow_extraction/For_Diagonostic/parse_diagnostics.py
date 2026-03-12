#!/usr/bin/env python3
"""
parse_diagnostics.py

Parses the stdout/log from save_mass_distributions and produces:
  - A full table of (cent_1pct, pT_bin, q2_ratio, q3_ratio, q2_candidates, q3_candidates)
  - The (cent, pT) where q2 and q3 max/min ratios are maximum and minimum.

Usage:
    python3 parse_diagnostics.py [log_file]
    # if log_file is omitted, reads from stdin
    ./save_mass_distributions | python3 parse_diagnostics.py

The script also prints per-iq candidate counts in a detailed section.
"""

import sys
import re
from collections import OrderedDict

# --------------------------------------------------------------------------
# Parsing helpers
# --------------------------------------------------------------------------
RE_CAND_HEADER = re.compile(
    r'\[cen1pct=(\d+),\s*pt=(\S+)\]\s+CANDIDATE entries in mass histograms'
)
RE_IQ_ROW = re.compile(
    r'^\s*(\d+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)'  # iq  q2_cand  q3_cand  (+ optional cols)
)
RE_TOTAL_ROW = re.compile(
    r'^\s*TOTAL\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)'
)
RE_Q2_RATIO = re.compile(r'q2 cand max/min ratio:\s*([\d.eE+\-]+)')
RE_Q3_RATIO = re.compile(r'q3 cand max/min ratio:\s*([\d.eE+\-]+)')

# --------------------------------------------------------------------------
# Main parser
# --------------------------------------------------------------------------
def parse_log(lines):
    """
    Returns a list of dicts, one per (cent_1pct, pt_name):
        {
          'cent': int,
          'pt':   str,
          'iq_q2': [float]*N_QBINS,
          'iq_q3': [float]*N_QBINS,
          'total_q2': float,
          'total_q3': float,
          'q2_ratio': float or None,
          'q3_ratio': float or None,
        }
    """
    records = []
    current = None
    in_cand_table = False

    for raw in lines:
        line = raw.rstrip()

        # Start of a new candidate block
        m = RE_CAND_HEADER.search(line)
        if m:
            current = {
                'cent': int(m.group(1)),
                'pt':   m.group(2),
                'iq_q2': [],
                'iq_q3': [],
                'total_q2': None,
                'total_q3': None,
                'q2_ratio': None,
                'q3_ratio': None,
            }
            records.append(current)
            in_cand_table = True
            continue

        if current is None:
            continue

        # Separator lines
        if re.match(r'^\s*[-=]+\s*$', line):
            continue

        # TOTAL row
        mt = RE_TOTAL_ROW.search(line)
        if mt:
            current['total_q2'] = float(mt.group(1))
            current['total_q3'] = float(mt.group(2))
            in_cand_table = False
            continue

        # iq data rows (only while we're inside the table)
        if in_cand_table:
            mr = RE_IQ_ROW.match(line)
            if mr:
                current['iq_q2'].append(float(mr.group(2)))
                current['iq_q3'].append(float(mr.group(3)))
                continue

        # Ratio lines
        mq2 = RE_Q2_RATIO.search(line)
        if mq2:
            current['q2_ratio'] = float(mq2.group(1))

        mq3 = RE_Q3_RATIO.search(line)
        if mq3:
            current['q3_ratio'] = float(mq3.group(1))

    return records


# --------------------------------------------------------------------------
# Report builders
# --------------------------------------------------------------------------
def ratio_summary_table(records):
    """Print a compact summary table of max/min ratios."""
    print("\n" + "="*90)
    print("  SUMMARY TABLE: q2/q3 candidate max/min ratios per (cent_1pct, pT_bin)")
    print("="*90)
    hdr = f"{'cent_1pct':>10}  {'pT_bin':<18}  {'q2_total':>12}  {'q3_total':>12}  {'q2_ratio':>10}  {'q3_ratio':>10}"
    print(hdr)
    print("-"*90)

    for r in records:
        q2r = f"{r['q2_ratio']:.4f}" if r['q2_ratio'] is not None else "  N/A   "
        q3r = f"{r['q3_ratio']:.4f}" if r['q3_ratio'] is not None else "  N/A   "
        tq2 = f"{r['total_q2']:.0f}" if r['total_q2'] is not None else "N/A"
        tq3 = f"{r['total_q3']:.0f}" if r['total_q3'] is not None else "N/A"
        print(f"{r['cent']:>10}  {r['pt']:<18}  {tq2:>12}  {tq3:>12}  {q2r:>10}  {q3r:>10}")

    print("="*90)


def extremes_report(records):
    """Print the (cent, pT) where q2/q3 ratios are highest and lowest."""
    valid_q2 = [(r['cent'], r['pt'], r['q2_ratio']) for r in records if r['q2_ratio'] is not None]
    valid_q3 = [(r['cent'], r['pt'], r['q3_ratio']) for r in records if r['q3_ratio'] is not None]

    print("\n" + "="*70)
    print("  EXTREMES: max and min max/min ratio across all (cent, pT) bins")
    print("="*70)

    if valid_q2:
        max_q2 = max(valid_q2, key=lambda x: x[2])
        min_q2 = min(valid_q2, key=lambda x: x[2])
        print(f"\n  Q2 candidate max/min ratio:")
        print(f"    Highest : cent_1pct={max_q2[0]:3d}, pT={max_q2[1]:<18s}  ratio={max_q2[2]:.4f}")
        print(f"    Lowest  : cent_1pct={min_q2[0]:3d}, pT={min_q2[1]:<18s}  ratio={min_q2[2]:.4f}")
    else:
        print("  No valid q2 ratios found.")

    if valid_q3:
        max_q3 = max(valid_q3, key=lambda x: x[2])
        min_q3 = min(valid_q3, key=lambda x: x[2])
        print(f"\n  Q3 candidate max/min ratio:")
        print(f"    Highest : cent_1pct={max_q3[0]:3d}, pT={max_q3[1]:<18s}  ratio={max_q3[2]:.4f}")
        print(f"    Lowest  : cent_1pct={min_q3[0]:3d}, pT={min_q3[1]:<18s}  ratio={min_q3[2]:.4f}")
    else:
        print("  No valid q3 ratios found.")

    print("="*70)


def per_iq_detail(records):
    """Print per-iq candidate counts for every (cent, pT)."""
    print("\n" + "="*90)
    print("  DETAILED PER-IQ CANDIDATE COUNTS")
    print("="*90)
    for r in records:
        print(f"\n  cent_1pct={r['cent']}  pT={r['pt']}")
        print(f"  {'iq':>4}  {'q2_cand':>14}  {'q3_cand':>14}")
        print("  " + "-"*40)
        n_iq = max(len(r['iq_q2']), len(r['iq_q3']))
        for iq in range(n_iq):
            q2c = r['iq_q2'][iq] if iq < len(r['iq_q2']) else float('nan')
            q3c = r['iq_q3'][iq] if iq < len(r['iq_q3']) else float('nan')
            print(f"  {iq:>4}  {q2c:>14.0f}  {q3c:>14.0f}")
        tq2 = r['total_q2'] if r['total_q2'] is not None else float('nan')
        tq3 = r['total_q3'] if r['total_q3'] is not None else float('nan')
        q2r = r['q2_ratio'] if r['q2_ratio'] is not None else float('nan')
        q3r = r['q3_ratio'] if r['q3_ratio'] is not None else float('nan')
        print("  " + "-"*40)
        print(f"  {'TOTAL':>4}  {tq2:>14.0f}  {tq3:>14.0f}")
        print(f"  q2 max/min ratio: {q2r:.4f}   q3 max/min ratio: {q3r:.4f}")


# --------------------------------------------------------------------------
# Imbalance flag table: which (cent, pT) bins have ratio > threshold
# --------------------------------------------------------------------------
def imbalance_flags(records, threshold=1.2):
    flagged = [(r['cent'], r['pt'], r['q2_ratio'], r['q3_ratio'])
               for r in records
               if (r['q2_ratio'] is not None and r['q2_ratio'] > threshold) or
                  (r['q3_ratio'] is not None and r['q3_ratio'] > threshold)]

    print(f"\n  IMBALANCED bins (ratio > {threshold}):")
    if not flagged:
        print(f"    None – all ratios <= {threshold}")
    else:
        print(f"  {'cent_1pct':>10}  {'pT_bin':<18}  {'q2_ratio':>10}  {'q3_ratio':>10}  flags")
        print("  " + "-"*70)
        for cent, pt, q2r, q3r in flagged:
            q2s = f"{q2r:.4f}" if q2r is not None else "  N/A "
            q3s = f"{q3r:.4f}" if q3r is not None else "  N/A "
            flags = []
            if q2r is not None and q2r > threshold: flags.append("Q2!")
            if q3r is not None and q3r > threshold: flags.append("Q3!")
            print(f"  {cent:>10}  {pt:<18}  {q2s:>10}  {q3s:>10}  {', '.join(flags)}")


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------
def main():
    if len(sys.argv) >= 2:
        fname = sys.argv[1]
        try:
            with open(fname) as fh:
                lines = fh.readlines()
        except FileNotFoundError:
            print(f"ERROR: cannot open file: {fname}", file=sys.stderr)
            sys.exit(1)
    else:
        lines = sys.stdin.readlines()

    records = parse_log(lines)

    if not records:
        print("No diagnostic records found in the input. "
              "Make sure the macro was compiled and run correctly.", file=sys.stderr)
        sys.exit(1)

    print(f"\nParsed {len(records)} (cent_1pct, pT) combinations.")

    # Print results
    ratio_summary_table(records)
    extremes_report(records)
    imbalance_flags(records, threshold=1.2)
    per_iq_detail(records)


if __name__ == "__main__":
    main()
