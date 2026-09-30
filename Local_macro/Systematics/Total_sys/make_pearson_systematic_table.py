#!/usr/bin/env python3
"""
make_pearson_systematic_table.py

Reads text_output/pearson_systematic_v{2,3}.txt (STEP 8 output of plot_ESE_systematics.C --
Pearson r's systematic uncertainty from Gaussian MC resampling, NOT a per-source quadratic
sum; see that step's own header comment) and writes a table* in the same visual style as the
Slope/Intercept tables, but with a single value-block (Pearson r has no Slope/Intercept split).

SystLo == SystHi always in this file (USE_SYMMETRIC_INTERVAL = true in the macro), so a single
number per bin is used. Relative(%) = |Syst| / |Nominal| * 100, left as "--" wherever
|Nominal r| < REL_ABS_THRESHOLD_R (0.02) -- same omission convention as the Slope/Intercept
tables, applied here since the text dump itself has no precomputed relative column.

Usage:
    python3 make_pearson_systematic_table.py <outbase>
"""
import sys
import os

PT_ORDER = {
    2: ["2-3", "3-4", "4-5", "5-6", "6-8", "8-10", "10-15", "15-30", "30-100"],
    3: ["2-4", "4-6", "6-8", "8-10", "10-20", "20-50"],
}
CENT_ORDER = ["cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"]
CENT_DISPLAY = ["0--10\\%", "10--20\\%", "20--30\\%", "30--40\\%", "40--50\\%"]
REL_ABS_THRESHOLD_R = 0.02


def pt_display(pt):
    return pt.replace("-", "--", 1)


def parse_file(path):
    rows = []
    with open(path) as f:
        for line in f:
            tok = line.split()
            if len(tok) != 7:
                continue
            cent, pt, nominal, statlo, stathi, systlo, systhi = tok
            try:
                float(nominal)
            except ValueError:
                continue  # header row
            rows.append({"cent": cent, "pt": pt, "nominal": float(nominal), "syst": float(systhi)})
    return rows


def lookup(rows, cent, pt):
    for r in rows:
        if r["cent"] == cent and r["pt"] == pt:
            return r
    return None


def write_table(fout, ivn, mode, rows):
    pts = PT_ORDER[ivn]
    mode_word = "" if mode == "abs" else "relative "
    label_suffix = "" if mode == "abs" else "_rel"

    fout.write("\\begin{table*}[htbp]\n")
    fout.write("    \\centering\n")
    fout.write(f"    \\caption{{Summary of {mode_word}systematic uncertainty on the Pearson "
               f"correlation coefficient $r$ for $v_{ivn}$, from Gaussian MC resampling "
               f"(D$^0$ $v_{{{ivn}}}$ points smeared by their total per-source systematic "
               f"uncertainty, statistical weights kept fixed) -- not a per-source quadratic "
               f"sum, since $r$ is bounded to $[-1,1]$ by definition.}}\n")
    fout.write(f"    \\label{{tab:sys_unc_q{ivn}_v{ivn}_PearsonR{label_suffix}}}\n")
    fout.write("    \\resizebox{\\textwidth}{!}{\n")
    fout.write("    \\begin{tabular}{l|ccccc}\n")
    fout.write("        \\hline\n")
    cent_hdr = " & ".join(CENT_DISPLAY)
    fout.write(f"        Centrality & {cent_hdr} \\\\\n")
    fout.write("        \\hline\n")
    fout.write("        $p_{\\mathrm{T}}$ (GeV/$c$) & & & & & \\\\\n")
    fout.write("        \\hline\n")

    for pt in pts:
        cells = [pt_display(pt)]
        for cent in CENT_ORDER:
            row = lookup(rows, cent, pt)
            if row is None:
                cells.append("--")
                continue
            if mode == "abs":
                cells.append(f"{row['syst']:.3f}")
            else:
                if abs(row["nominal"]) < REL_ABS_THRESHOLD_R:
                    cells.append("--")
                else:
                    cells.append(f"{abs(row['syst'] / row['nominal']) * 100:.1f}\\%")
        fout.write("        " + " & ".join(cells) + " \\\\\n")
    fout.write("        \\hline\n")
    fout.write("    \\end{tabular}\n")
    fout.write("    }\n")
    fout.write("\\end{table*}\n\n")


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <outbase>", file=sys.stderr)
        sys.exit(1)
    outbase = sys.argv[1]

    rows_by_ivn = {}
    for ivn in (2, 3):
        inpath = os.path.join(outbase, "text_output", f"pearson_systematic_v{ivn}.txt")
        rows_by_ivn[ivn] = parse_file(inpath)

    outpath = os.path.join(outbase, "text_output", "pearson_systematic_tables.tex")
    with open(outpath, "w") as fout:
        for mode in ("abs", "rel"):
            for ivn in (2, 3):
                fout.write(f"% ===== PearsonR : {mode} : v{ivn} =====\n")
                write_table(fout, ivn, mode, rows_by_ivn[ivn])
    print(f"Wrote {outpath}")


if __name__ == "__main__":
    main()
