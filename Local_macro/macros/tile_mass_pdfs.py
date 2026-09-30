#!/usr/bin/env python3
"""
tile_mass_pdfs.py
------------------
Tiles the individual single-pad mass PDFs (saved by plot_mass_qbins.C) into
multi-page PDFs.

  data_pdfs_tmp/d_c{ic}_p{ip}_q{iq}.pdf   (5 x 9 x 12 = 540 files)
      -> DataMass_AllQ2Bins.pdf     (45 pages, 4 cols x 3 rows = 12/page,
                                      one page per centrality/pT combination)

  mc_pdfs_tmp/m_c{ic}_p{ip}.pdf           (5 x 9 = 45 files)
      -> MCMass_Sig_AllPtBins.pdf   (5 pages, 3 cols x 3 rows = 9/page,
                                      one page per centrality)

ROOT saves each canvas PDF with a rotation tag that PyMuPDF's show_pdf_page
does NOT apply, so we render each source page to a pixmap with get_pixmap()
(which does apply the rotation correctly) and embed that image instead.

DPI=150 gives good quality per cell while keeping file size reasonable.

Usage:
    python3 tile_mass_pdfs.py data   # tile the 540 data plots
    python3 tile_mass_pdfs.py mc     # tile the 45 MC plots
    python3 tile_mass_pdfs.py        # do both
"""

import fitz
import os, sys
import math

DPI = 150  # render resolution for each individual plot

NCENT = 5
NPT   = 9
NQ    = 12


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


def tile_data():
    tmp_dir = "data_pdfs_tmp"
    pages = []
    for ic in range(NCENT):
        for ip in range(NPT):
            pages.append([os.path.join(tmp_dir, f"d_c{ic}_p{ip}_q{iq}.pdf")
                          for iq in range(NQ)])
    tile(pages, "DataMass_AllQ2Bins.pdf", ncols=4, nrows=3)


def tile_mc():
    tmp_dir = "mc_pdfs_tmp"
    pages = []
    for ic in range(NCENT):
        pages.append([os.path.join(tmp_dir, f"m_c{ic}_p{ip}.pdf")
                      for ip in range(NPT)])
    tile(pages, "MCMass_Sig_AllPtBins.pdf", ncols=3, nrows=3)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which == "data":
        tile_data()
    elif which == "mc":
        tile_mc()
    else:
        tile_data()
        tile_mc()