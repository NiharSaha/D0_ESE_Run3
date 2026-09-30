#!/usr/bin/env python3
"""
make_total_systematic_table.py

Reads text_output/systematic_uncertainty_v{2,3}.txt and writes the combined TOTAL
systematic uncertainty (quadratic sum over all active sources -- the TotalSyst_abs /
TotalSyst_rel% columns already in that file) in the same table* format used for the
per-source tables, for both harmonics and both abs/relative(%).

Usage:
    python3 make_total_systematic_table.py <outbase>
"""
import sys
import os
import importlib.util

spec = importlib.util.spec_from_file_location("srctab", os.path.join(os.path.dirname(__file__), "make_source_summary_tables.py"))
srctab = importlib.util.module_from_spec(spec)
spec.loader.exec_module(srctab)


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <outbase>", file=sys.stderr)
        sys.exit(1)
    outbase = sys.argv[1]

    blocks_by_ivn = {}
    for ivn in (2, 3):
        inpath = os.path.join(outbase, "text_output", f"systematic_uncertainty_v{ivn}.txt")
        blocks_by_ivn[ivn] = srctab.parse_file(inpath)

    outpath = os.path.join(outbase, "text_output", "total_systematic_tables.tex")
    with open(outpath, "w") as fout:
        for mode in ("abs", "rel"):
            for ivn in (2, 3):
                blocks = blocks_by_ivn[ivn]
                fout.write(f"% ===== Total : {mode} : v{ivn} =====\n")
                mode_word = "" if mode == "abs" else "relative "
                caption = (f"Summary of {mode_word}total systematic uncertainty (quadratic sum "
                           f"over all systematic sources) on the D0-vs-charged-particle Slope "
                           f"and Intercept for $v_{ivn}$")
                label_suffix = "" if mode == "abs" else "_rel"
                label = f"tab:sys_unc_q{ivn}_v{ivn}_Total{label_suffix}"
                srctab.write_table(fout, ivn, None, "Total", "", mode, blocks,
                                    caption_override=caption, label_override=label)
    print(f"Wrote {outpath}")


if __name__ == "__main__":
    main()
