#!/usr/bin/env python3
"""
Parse all_bins_summary.log produced by run_all_bins.sh
and report:
  - q2 cand max/min ratio for every (cent, pT) bin
  - q3 cand max/min ratio for every (cent, pT) bin
  - which bin has the MAXIMUM and MINIMUM ratio
  - all bins flagged WARNING
"""

import sys
import re

def parse_log(logfile):
    results = []   # list of dicts: {cen, pt, q2_ratio, q3_ratio}
    current = {}

    with open(logfile) as f:
        for line in f:
            line = line.rstrip()

            # New block: --- cen=X  pt=Y ---
            m = re.match(r'---\s*cen=(\d+)\s+pt=(\d+)\s*---', line)
            if m:
                if current:
                    results.append(current)
                current = {
                    'cen': int(m.group(1)),
                    'pt':  int(m.group(2)),
                    'q2_ratio': None,
                    'q3_ratio': None,
                    'warnings': []
                }
                continue

            # q2 cand max/min ratio: X.XXXX
            m2 = re.search(r'q2 cand max/min ratio:\s*([\d.]+)', line)
            if m2 and current:
                current['q2_ratio'] = float(m2.group(1))
                continue

            # q3 cand max/min ratio: X.XXXX
            m3 = re.search(r'q3 cand max/min ratio:\s*([\d.]+)', line)
            if m3 and current:
                current['q3_ratio'] = float(m3.group(1))
                continue

            # WARNING lines
            if 'WARNING' in line and current:
                current['warnings'].append(line.strip())

    if current:
        results.append(current)

    return results


def main():
    logfile = sys.argv[1] if len(sys.argv) > 1 else "all_bins_summary.log"
    results = parse_log(logfile)

    if not results:
        print("No results parsed. Check the log file.")
        sys.exit(1)

    pt_names = ["pT1to2","pT2to3","pT3to4","pT4to5","pT5to6",
                "pT6to8","pT8to10","pT10to15","pT15to20","pT20to40"]

    def pt_label(i):
        return pt_names[i] if i < len(pt_names) else str(i)

    # ------------------------------------------------------------------
    # Full table
    # ------------------------------------------------------------------
    print("\n" + "="*72)
    print(f"  {'cen':>6}  {'pT':>12}  {'q2_ratio':>12}  {'q3_ratio':>12}")
    print("="*72)
    for r in results:
        q2s = f"{r['q2_ratio']:.4f}" if r['q2_ratio'] is not None else "  N/A  "
        q3s = f"{r['q3_ratio']:.4f}" if r['q3_ratio'] is not None else "  N/A  "
        print(f"  {r['cen']:>6}  {pt_label(r['pt']):>12}  {q2s:>12}  {q3s:>12}")
    print("="*72)

    # ------------------------------------------------------------------
    # Find extremes for q2
    # ------------------------------------------------------------------
    q2_valid = [r for r in results if r['q2_ratio'] is not None]
    q3_valid = [r for r in results if r['q3_ratio'] is not None]

    print("\n--- Q2 CANDIDATE RATIO EXTREMES ---")
    if q2_valid:
        rmax = max(q2_valid, key=lambda r: r['q2_ratio'])
        rmin = min(q2_valid, key=lambda r: r['q2_ratio'])
        print(f"  MAX q2 ratio: {rmax['q2_ratio']:.4f}  "
              f"at  cen1pct={rmax['cen']}  pt={pt_label(rmax['pt'])}")
        print(f"  MIN q2 ratio: {rmin['q2_ratio']:.4f}  "
              f"at  cen1pct={rmin['cen']}  pt={pt_label(rmin['pt'])}")

    print("\n--- Q3 CANDIDATE RATIO EXTREMES ---")
    if q3_valid:
        rmax = max(q3_valid, key=lambda r: r['q3_ratio'])
        rmin = min(q3_valid, key=lambda r: r['q3_ratio'])
        print(f"  MAX q3 ratio: {rmax['q3_ratio']:.4f}  "
              f"at  cen1pct={rmax['cen']}  pt={pt_label(rmax['pt'])}")
        print(f"  MIN q3 ratio: {rmin['q3_ratio']:.4f}  "
              f"at  cen1pct={rmin['cen']}  pt={pt_label(rmin['pt'])}")

    # ------------------------------------------------------------------
    # Top 10 highest ratios
    # ------------------------------------------------------------------
    print("\n--- TOP 10 HIGHEST q2 RATIOS ---")
    for r in sorted(q2_valid, key=lambda x: x['q2_ratio'], reverse=True)[:10]:
        print(f"  cen1pct={r['cen']:3d}  pt={pt_label(r['pt']):>12}  "
              f"q2_ratio={r['q2_ratio']:.4f}")

    print("\n--- TOP 10 HIGHEST q3 RATIOS ---")
    for r in sorted(q3_valid, key=lambda x: x['q3_ratio'], reverse=True)[:10]:
        print(f"  cen1pct={r['cen']:3d}  pt={pt_label(r['pt']):>12}  "
              f"q3_ratio={r['q3_ratio']:.4f}")

    # ------------------------------------------------------------------
    # All WARNING bins
    # ------------------------------------------------------------------
    warn_bins = [r for r in results if r['warnings']]
    print(f"\n--- BINS WITH WARNINGS ({len(warn_bins)} total) ---")
    if warn_bins:
        for r in warn_bins:
            print(f"  cen1pct={r['cen']:3d}  pt={pt_label(r['pt']):>12}  "
                  f"q2={r['q2_ratio']}  q3={r['q3_ratio']}")
            for w in r['warnings']:
                print(f"    {w}")
    else:
        print("  None found.")

    print()


if __name__ == "__main__":
    main()
