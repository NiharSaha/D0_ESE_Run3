#!/usr/bin/env python3
"""
tile_SPdist_pdfs.py
--------------------
Tiles the individual single-pad SP-distribution PDFs (saved by
plot_SPdist_qbins.C) into multi-page PDFs, one page per (cent, pT)
combination with all 12 q bins arranged 4 cols x 3 rows.

  histogram_SP/v2_pdfs_tmp/sp_c{ic}_p{ip}_q{iq}.pdf   (5 x 9 x 12 = 540 files)
      -> histogram_SP/SPdist_v2_AllQ2Bins.pdf   (45 pages)

  histogram_SP/v3_pdfs_tmp/sp_c{ic}_p{ip}_q{iq}.pdf   (5 x 6 x 12 = 360 files)
      -> histogram_SP/SPdist_v3_AllQ3Bins.pdf   (30 pages)

Usage:
    python3 tile_SPdist_pdfs.py v2   # tile the 540 v2 plots
    python3 tile_SPdist_pdfs.py v3   # tile the 360 v3 plots
    python3 tile_SPdist_pdfs.py      # do both
"""

import fitz
import os, sys

DPI = 150  # render resolution for each individual plot

NCENT = 5
NQ    = 12
NPT   = {"v2": 9, "v3": 6}

OUT_DIR = "histogram_SP"


def get_cell_size(paths):
    """Grab width/height (post-rotation) from the first PDF that exists."""
    for p in paths:
        if os.path.exists(p):
            src = fitz.open(p)
            r = src[0].rect
            src.close()
            return r.width, r.height
    return None, None


def tile(pages_of_paths, out_file, ncols, nrows):
    """pages_of_paths: list of pages, each a list of source PDF paths in the
    order they should be placed left-to-right, top-to-bottom."""

    all_paths = [p for page in pages_of_paths for p in page]
    cell_w, cell_h = get_cell_size(all_paths)
    if not cell_w:
        print(f"ERROR: no individual PDFs found for {out_file}", file=sys.stderr)
        sys.exit(1)

    page_w, page_h = cell_w * ncols, cell_h * nrows
    print(f"[{out_file}] cell {cell_w:.0f}x{cell_h:.0f} pt -> "
          f"page {page_w:.0f}x{page_h:.0f} pt, {len(pages_of_paths)} pages")

    out_doc = fitz.open()
    for page_paths in pages_of_paths:
        page = out_doc.new_page(width=page_w, height=page_h)
        for i, fname in enumerate(page_paths):
            if not os.path.exists(fname):
                print(f"WARNING: {fname} not found, skipping", file=sys.stderr)
                continue
            col, row = i % ncols, i // ncols
            rect = fitz.Rect(col * cell_w, row * cell_h,
                              (col + 1) * cell_w, (row + 1) * cell_h)
            src = fitz.open(fname)
            pix = src[0].get_pixmap(dpi=DPI)
            page.insert_image(rect, pixmap=pix)
            src.close()

    out_doc.save(out_file, garbage=4, deflate=True)
    out_doc.close()
    print(f"Saved: {out_file}")


def tile_harmonic(tag):
    npt = NPT[tag]
    tmp_dir = os.path.join(OUT_DIR, f"{tag}_pdfs_tmp")
    pages = []
    for ic in range(NCENT):
        for ip in range(npt):
            pages.append([os.path.join(tmp_dir, f"sp_c{ic}_p{ip}_q{iq}.pdf")
                          for iq in range(NQ)])
    qtag = "2" if tag == "v2" else "3"
    out_file = os.path.join(OUT_DIR, f"SPdist_{tag}_AllQ{qtag}Bins.pdf")
    tile(pages, out_file, ncols=4, nrows=3)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which == "both":
        tile_harmonic("v2")
        tile_harmonic("v3")
    else:
        tile_harmonic(which)
