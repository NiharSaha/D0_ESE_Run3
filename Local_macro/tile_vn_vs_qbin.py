#!/usr/bin/env python3
"""
tile_vn_vs_qbin.py
Tiles the per-(cent x pT) v3/v2 PDFs from the specified input directory
into a single multi-page summary PDF.
Layout: 5 columns x 2 rows per page  (= 10 plots per page)
Cent order: 0-10, 10-20, 20-30, 30-40, 40-50
pT order (v3): 2-4, 4-6, 6-8, 8-10, 10-20, 20-50, 50-100
pT order (v2): 2-3, 3-4, 4-5, 5-6, 6-8, 8-10, 10-15, 15-30, 30-100
"""
import fitz
import os, sys
import math

# Require the 3 arguments passed from the C++ macro; an optional 4th
# restricts tiling to a single centrality tag (e.g. "cent0to50" for the
# inclusive-centrality summary), producing a 1-page PDF instead of one page
# per differential centrality.
if len(sys.argv) < 4:
    print("Usage: python3 tile_vn_vs_qbin.py <harmonic> <input_directory> <output_prefix> [cent_tag]")
    sys.exit(1)

harmonic = int(sys.argv[1])
input_dir = sys.argv[2]
out_prefix = sys.argv[3]
cent_filter = sys.argv[4] if len(sys.argv) > 4 else None

# Construct dynamic paths
TMP_DIR  = f"{input_dir}/v{harmonic}"
OUT_FILE = f"{out_prefix}.pdf"
DPI      = 150
NCOLS    = 3   # columns per page

# Automatically determine the file naming tag ("meanq" vs "qbin") based on the input directory
type_tag = "meanq" if "meanq" in input_dir else "qbin"

# Ordered list of centrality and pT tags matching the filename convention
cent_tags = ["cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"]
if harmonic == 2:
    pt_tags = ["pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8",
               "pT8to10", "pT10to15", "pT15to30", "pT30to100"]
else:
    pt_tags = ["pT2to4", "pT4to6", "pT6to8", "pT8to10", "pT10to20", "pT20to50"]

if cent_filter:
    cent_tags = [cent_filter]

NROWS = math.ceil(len(pt_tags) / NCOLS)   # rows needed per page

# Build ordered file list grouped by centrality (one group = one page)
ordered = []
for cent in cent_tags:
    for pt in pt_tags:
        # Dynamically inject the type_tag into the filename template
        fname = os.path.join(TMP_DIR, f"vn_vs_{type_tag}_v{harmonic}_{cent}_{pt}.pdf")
        if os.path.isfile(fname):
            ordered.append(fname)
        else:
            print(f"  MISSING: {fname}", file=sys.stderr)

if not ordered:
    print(f"ERROR: No PDFs found in {TMP_DIR}", file=sys.stderr)
    sys.exit(1)

print(f"Found {len(ordered)} PDFs to tile ({len(cent_tags)} pages, {len(pt_tags)} plots/page).")

# Determine cell size from first available PDF
src0 = fitz.open(ordered[0])
pix0 = src0[0].get_pixmap(dpi=DPI)
cell_w = pix0.width
cell_h = pix0.height
src0.close()

PAGE_W = NCOLS * cell_w
PAGE_H = NROWS * cell_h

out = fitz.open()
plots_per_page = len(pt_tags)   # one page per centrality

for page_idx in range(0, len(ordered), plots_per_page):
    pg = out.new_page(width=PAGE_W, height=PAGE_H)
    batch = ordered[page_idx : page_idx + plots_per_page]
    for slot, fpath in enumerate(batch):
        col = slot % NCOLS
        row = slot // NCOLS
        x0  = col * cell_w
        y0  = row * cell_h
        rect = fitz.Rect(x0, y0, x0 + cell_w, y0 + cell_h)
        src = fitz.open(fpath)
        pix = src[0].get_pixmap(dpi=DPI)
        pg.insert_image(rect, pixmap=pix)
        src.close()

out.save(OUT_FILE)
print(f"Saved {len(out)}-page summary to {OUT_FILE}")