#!/usr/bin/env python3
"""Arrange the per-bin Canvas_BestFit PDFs (saved by PlotBestFitPages.C) into a
5-page PDF, one page per centrality, each a 3x3 grid of the 9 pT bins.

Uses PyMuPDF's show_pdf_page(), which embeds each source PDF's vector content
directly into the output page (no rasterization), unlike the earlier
PNG-then-ImageMagick approach which visibly softened lines/text.
"""
import sys
import fitz  # PyMuPDF

tmp_dir, out_pdf = sys.argv[1], sys.argv[2]

cent_labels = ["0-10%", "10-20%", "20-30%", "30-40%", "40-50%"]
n_cent, n_pt = 5, 9
n_rows, n_cols = 3, 3

page_w, page_h = 1800, 1800
header_h = 60  # top strip reserved for the centrality label

grid_h = page_h - header_h
cell_w = page_w / n_cols
cell_h = grid_h / n_rows

out_doc = fitz.open()

for icent in range(n_cent):
    page = out_doc.new_page(width=page_w, height=page_h)
    label = f"Centrality {cent_labels[icent]}"
    fontsize = 28
    text_width = fitz.get_text_length(label, fontname="helv", fontsize=fontsize)
    page.insert_text(((page_w - text_width) / 2, header_h * 0.65), label,
                      fontsize=fontsize, fontname="helv")

    for ipt in range(n_pt):
        cell_path = f"{tmp_dir}/cell_{icent}_{ipt}.pdf"
        try:
            src = fitz.open(cell_path)
        except Exception:
            print(f"Missing cell PDF: {cell_path}")
            continue

        row, col = divmod(ipt, n_cols)
        x0 = col * cell_w
        y0 = header_h + row * cell_h
        rect = fitz.Rect(x0, y0, x0 + cell_w, y0 + cell_h)
        # ROOT's PDF driver stamps its pages with /Rotate 90. show_pdf_page() needs a
        # counter-rotation to display it correctly here, but its keep_proportion sizing
        # gets confused (clips the content) if the source page's own /Rotate flag is
        # still set while also passing a nonzero rotate= -- so strip that flag (metadata
        # only, the content stream itself is untouched) before embedding.
        orig_rotation = src[0].rotation
        src[0].set_rotation(0)
        page.show_pdf_page(rect, src, 0, rotate=-orig_rotation)
        src.close()

out_doc.save(out_pdf)
out_doc.close()
print(f"Wrote {out_pdf}")
