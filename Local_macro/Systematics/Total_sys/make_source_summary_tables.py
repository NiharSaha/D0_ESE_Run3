#!/usr/bin/env python3
"""
make_source_summary_tables.py

Reads text_output/systematic_uncertainty_v{2,3}.txt (STEP 6 output of
plot_ESE_systematics.C) and, for each requested systematic source, writes the "table*
with resizebox / 5+5 centrality columns for Slope+Intercept / one row per measured pT bin"
LaTeX table (absolute and relative-% versions) in the exact format given interactively for
SigPDF -- used here so the 20 remaining source/harmonic/abs-or-rel tables are generated
directly from the source-of-truth text dump instead of transcribed by hand.

Convention (confirmed with the user):
  - Cell value = |Adopted value| (unsigned magnitude): the source's adopted (whichever
    variation deviates more from nominal, for a two-variation source) deviation from nominal.
  - Absolute table: always filled when the source has a valid point in that bin.
  - Relative(%) table: filled with |Adopted_rel%|, EXCEPT left as "--" wherever
    |nominal Intercept| is below REL_ABS_THRESHOLD_INTER (marked "n/a*"/"n/a" in the text
    dump) -- same omission rule already used throughout plot_ESE_systematics.C's own
    reldiff/summary plots and text dump, so the two stay consistent.

Usage:
    python3 make_source_summary_tables.py <outbase>
e.g.
    python3 make_source_summary_tables.py ESE_systematics_20260916
"""
import sys
import os

SOURCE_LABELS = ["SigPDF", "BkgPDF", "BDT", "Centrality", "qn quantile", "PromptFraction", "Syst7"]
N_SOURCES = len(SOURCE_LABELS)
N_LEAD = 5
N_PER_SRC = 6
N_TRAIL = 2

PT_ORDER = {
    2: ["2-3", "3-4", "4-5", "5-6", "6-8", "8-10", "10-15", "15-30", "30-100"],
    3: ["2-4", "4-6", "6-8", "8-10", "10-20", "20-50"],
}
CENT_ORDER = ["cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"]
CENT_DISPLAY = ["0--10\\%", "10--20\\%", "20--30\\%", "30--40\\%", "40--50\\%"]

# One entry per source we're asked to tabulate here (SigPDF already done interactively).
SOURCES_TODO = [
    ("BkgPDF", "background shape variation"),
    ("BDT", "BDT selection (loose/tight) variation"),
    ("Centrality", "centrality bin border variation"),
    ("qn quantile", "ESE $q_n$-quantile binning variation"),
    ("PromptFraction", "non-prompt D$^0$ (prompt-fraction) correction"),
]


def pt_display(pt):
    return pt.replace("-", "--", 1)


def parse_quantity_block(lines):
    rows = []
    for line in lines:
        line = line.rstrip("\n")
        if not line.strip():
            break
        tok = line.split()
        if len(tok) != N_LEAD + N_PER_SRC * N_SOURCES + N_TRAIL:
            continue
        cent = tok[0]
        pt = tok[1]
        srcs = []
        pos = N_LEAD
        for _ in range(N_SOURCES):
            srcs.append((tok[pos + 4], tok[pos + 5]))  # (Adopted_abs, Adopted_rel)
            pos += N_PER_SRC
        total = (tok[pos], tok[pos + 1])  # (TotalSyst_abs, TotalSyst_rel%)
        rows.append({"cent": cent, "pt": pt, "sources": srcs, "total": total})
    return rows


def parse_file(path):
    with open(path) as f:
        lines = f.readlines()
    blocks = {}
    i = 0
    while i < len(lines):
        if lines[i].startswith("Quantity:"):
            qname = lines[i].split("Quantity:")[1].strip()
            i += 3
            blocks[qname] = parse_quantity_block(lines[i:])
        i += 1
    return blocks


def fmt_abs(tok):
    try:
        return f"{abs(float(tok)):.3f}"
    except ValueError:
        return "--"


def fmt_rel(tok):
    tok = tok.rstrip("*")
    try:
        return f"{abs(float(tok)):.1f}\\%"
    except ValueError:
        return "--"


def lookup(rows, cent, pt):
    for r in rows:
        if r["cent"] == cent and r["pt"] == pt:
            return r
    return None


def write_table(fout, ivn, src_idx, src_label, caption_desc, mode, blocks, caption_override=None,
                 label_override=None):
    """mode: 'abs' or 'rel'. src_idx=None means use the row's combined Total (quadratic sum
    over all sources) instead of one source's Adopted value."""
    pts = PT_ORDER[ivn]
    fmt = fmt_abs if mode == "abs" else fmt_rel
    mode_word = "" if mode == "abs" else "relative "
    label_suffix = "" if mode == "abs" else "_rel"
    safe_label = src_label.replace(" ", "")

    fout.write("\\begin{table*}[htbp]\n")
    fout.write("    \\centering\n")
    caption = caption_override or f"Summary of {mode_word}systematic uncertainties due to {caption_desc} for $v_{ivn}$"
    fout.write(f"    \\caption{{{caption}}}\n")
    label = label_override or f"tab:sys_unc_q{ivn}_v{ivn}_{safe_label}{label_suffix}"
    fout.write(f"    \\label{{{label}}}\n")
    fout.write("    \\resizebox{\\textwidth}{!}{\n")
    fout.write("    \\begin{tabular}{l|ccccc|ccccc}\n")
    fout.write("        \\hline\n")
    fout.write("        & \\multicolumn{5}{c|}{Uncertainty in slopes} & "
                "\\multicolumn{5}{c}{Uncertainty in intercepts} \\\\\n")
    fout.write("        \\cline{2-6} \\cline{7-11}\n")
    cent_hdr = " & ".join(CENT_DISPLAY)
    fout.write(f"        Centrality & {cent_hdr} & {cent_hdr} \\\\\n")
    fout.write("        \\hline\n")
    fout.write("        $p_{\\mathrm{T}}$ (GeV/$c$) & & & & & & & & & & \\\\\n")
    fout.write("        \\hline\n")

    for pt in pts:
        cells = [pt_display(pt)]
        for qname in ("Slope", "Intercept"):
            rows = blocks.get(qname, [])
            for cent in CENT_ORDER:
                row = lookup(rows, cent, pt)
                if row is None:
                    cells.append("--")
                    continue
                a, r = row["total"] if src_idx is None else row["sources"][src_idx]
                cells.append(fmt(a if mode == "abs" else r))
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

    blocks_by_ivn = {}
    for ivn in (2, 3):
        inpath = os.path.join(outbase, "text_output", f"systematic_uncertainty_v{ivn}.txt")
        blocks_by_ivn[ivn] = parse_file(inpath)

    # Ordered: grouped by source; within each source, abs-v2, abs-v3, rel-v2, rel-v3.
    outpath = os.path.join(outbase, "text_output", "source_summary_tables.tex")
    with open(outpath, "w") as fout:
        for src_label, caption_desc in SOURCES_TODO:
            src_idx = SOURCE_LABELS.index(src_label)
            for mode in ("abs", "rel"):
                for ivn in (2, 3):
                    blocks = blocks_by_ivn[ivn]
                    fout.write(f"% ===== {src_label} : {mode} : v{ivn} =====\n")
                    write_table(fout, ivn, src_idx, src_label, caption_desc, mode, blocks)
    print(f"Wrote {outpath}")


if __name__ == "__main__":
    main()
