#!/usr/bin/env python3
"""
tile_effW_comparison.py
------------------------
Tiles the individual with-eff-vs-without-eff comparison PDFs saved by
plot_ESE_effW_comparison.C into summary PDFs, one call per category
(handles both harmonics v2/v3 internally):

  d0      : vn_vs_qbin_D0comp_v{n}_{cent}_{pT}.pdf   -> summary_D0comp_v{n}.pdf
            (one page per centrality, pT bins tiled 3 cols)
  scatter : scatter_comp_v{n}_{cent}_{pT}.pdf        -> summary_scatterComp_v{n}.pdf
            (one page per centrality, pT bins tiled 3 cols)
  chg     : vn_vs_qbin_ChgComp_v{n}_{cent}.pdf       -> summary_ChgComp_v{n}.pdf
            (single page, all 6 centralities tiled 3 cols x 2 rows;
             charged vn is pT-integrated, so there is no pT dimension to tile)

Usage:
    python3 tile_effW_comparison.py d0      <input_dir> <output_prefix>
    python3 tile_effW_comparison.py scatter <input_dir> <output_prefix>
    python3 tile_effW_comparison.py chg     <input_dir> <output_prefix>
"""

import fitz
import os, sys
import math

DPI = 150
NCOLS = 3

CENT_ORDER = ["cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50", "cent0to50"]

PT_TAGS_V2 = ["pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8",
              "pT8to10", "pT10to15", "pT15to30", "pT30to100"]
PT_TAGS_V3 = ["pT2to4", "pT4to6", "pT6to8", "pT8to10", "pT10to20", "pT20to50"]


def tile_by_cent_pt(harmonic, input_dir, out_prefix, file_prefix):
    """One page per centrality, pT bins tiled left-to-right/top-to-bottom (NCOLS cols)."""
    pt_tags = PT_TAGS_V2 if harmonic == 2 else PT_TAGS_V3
    nrows = math.ceil(len(pt_tags) / NCOLS)

    ordered = []
    for cent in CENT_ORDER:
        for pt in pt_tags:
            fname = os.path.join(input_dir, f"{file_prefix}_v{harmonic}_{cent}_{pt}.pdf")
            if os.path.isfile(fname):
                ordered.append(fname)
            else:
                print(f"  MISSING: {fname}", file=sys.stderr)

    if not ordered:
        print(f"No PDFs found for {file_prefix} v{harmonic} in {input_dir}", file=sys.stderr)
        return

    src0 = fitz.open(ordered[0])
    pix0 = src0[0].get_pixmap(dpi=DPI)
    cell_w, cell_h = pix0.width, pix0.height
    src0.close()
    page_w, page_h = NCOLS * cell_w, nrows * cell_h

    out = fitz.open()
    plots_per_page = len(pt_tags)
    for i in range(0, len(ordered), plots_per_page):
        pg = out.new_page(width=page_w, height=page_h)
        batch = ordered[i:i + plots_per_page]
        for slot, fpath in enumerate(batch):
            col, row = slot % NCOLS, slot // NCOLS
            rect = fitz.Rect(col * cell_w, row * cell_h, (col + 1) * cell_w, (row + 1) * cell_h)
            src = fitz.open(fpath)
            pix = src[0].get_pixmap(dpi=DPI)
            pg.insert_image(rect, pixmap=pix)
            src.close()

    out_file = f"{out_prefix}_v{harmonic}.pdf"
    out.save(out_file)
    print(f"Saved {len(out)}-page summary to {out_file}")


def tile_chg(harmonic, input_dir, out_prefix):
    """Single page, all centralities tiled NCOLS cols (no pT dimension)."""
    ordered = []
    for cent in CENT_ORDER:
        fname = os.path.join(input_dir, f"vn_vs_qbin_ChgComp_v{harmonic}_{cent}.pdf")
        if os.path.isfile(fname):
            ordered.append(fname)
        else:
            print(f"  MISSING: {fname}", file=sys.stderr)

    if not ordered:
        print(f"No PDFs found for ChgComp v{harmonic} in {input_dir}", file=sys.stderr)
        return

    rows = math.ceil(len(ordered) / NCOLS)
    src0 = fitz.open(ordered[0])
    pix0 = src0[0].get_pixmap(dpi=DPI)
    cell_w, cell_h = pix0.width, pix0.height
    src0.close()

    out = fitz.open()
    pg = out.new_page(width=NCOLS * cell_w, height=rows * cell_h)
    for slot, fpath in enumerate(ordered):
        col, row = slot % NCOLS, slot // NCOLS
        rect = fitz.Rect(col * cell_w, row * cell_h, (col + 1) * cell_w, (row + 1) * cell_h)
        src = fitz.open(fpath)
        pix = src[0].get_pixmap(dpi=DPI)
        pg.insert_image(rect, pixmap=pix)
        src.close()

    out_file = f"{out_prefix}_v{harmonic}.pdf"
    out.save(out_file)
    print(f"Saved 1-page summary to {out_file}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 tile_effW_comparison.py <d0|scatter|chg> <input_dir> <output_prefix>")
        sys.exit(1)

    mode, input_dir, out_prefix = sys.argv[1], sys.argv[2], sys.argv[3]

    for harmonic in (2, 3):
        if mode == "d0":
            tile_by_cent_pt(harmonic, input_dir, out_prefix, "vn_vs_qbin_D0comp")
        elif mode == "scatter":
            tile_by_cent_pt(harmonic, input_dir, out_prefix, "scatter_comp")
        elif mode == "chg":
            tile_chg(harmonic, input_dir, out_prefix)
        else:
            print(f"Unknown mode: {mode}", file=sys.stderr)
            sys.exit(1)
