#!/usr/bin/env python3
"""
merge_yield_plots.py
--------------------
Combines the individual yield-vs-SP-bin PDF plots produced by fit_mass_and_flow
into two multi-page summary PDFs (one for v2, one for v3).

Page layout
-----------
  One page per (cent, pT) combination.
  Each page contains a 2 × 5 grid = 10 panels, one per q bin (0-9).

Outputs (written to OUTPUT_DIR)
-------
  yield_v2_summary.pdf
  yield_v3_summary.pdf

Usage
-----
  # use hard-coded paths below
  python3 merge_yield_plots.py

  # override input / output directory from the command line
  python3 merge_yield_plots.py <input_dir>
  python3 merge_yield_plots.py <input_dir> <output_dir>

PDF rendering back-ends (tried in order)
-----------------------------------------
  1. PyMuPDF  (fitz)      — pip install pymupdf
  2. pdf2image            — pip install pdf2image  (needs poppler)
  3. Ghostscript          — gs must be on PATH
"""

import io
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import warnings
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ── user-configurable paths ────────────────────────────────────────────────────
INPUT_DIR = (
    "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Flow_output_FullStat_May17_v1/output/yield_vs_spbin_plots"
)
OUTPUT_DIR = "."   # summary PDFs are written to the current working directory

# layout
N_QBINS   = 10   # panels per page
GRID_ROWS = 2
GRID_COLS = 5    # 2 × 5 = 10

# render resolution for each sub-panel image (higher → sharper but slower)
DPI = 150
# ──────────────────────────────────────────────────────────────────────────────


# ── PDF → PIL Image rendering ─────────────────────────────────────────────────

def _render_via_fitz(pdf_path, dpi):
    import fitz  # PyMuPDF
    from PIL import Image as PILImage
    doc = fitz.open(pdf_path)
    page = doc.load_page(0)
    mat = fitz.Matrix(dpi / 72.0, dpi / 72.0)
    pix = page.get_pixmap(matrix=mat, alpha=False)
    img = PILImage.open(io.BytesIO(pix.tobytes("png")))
    img.load()
    doc.close()
    return img


def _render_via_pdf2image(pdf_path, dpi):
    from pdf2image import convert_from_path
    return convert_from_path(pdf_path, dpi=dpi)[0]


def _render_via_ghostscript(pdf_path, dpi):
    import subprocess, tempfile
    from PIL import Image as PILImage
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        subprocess.run(
            ["gs", "-dNOPAUSE", "-dBATCH", "-dSAFER",
             "-sDEVICE=png16m", f"-r{dpi}",
             f"-sOutputFile={tmp_path}", pdf_path],
            check=True, capture_output=True,
        )
        img = PILImage.open(tmp_path)
        img.load()
    finally:
        os.unlink(tmp_path)
    return img


def render_pdf_page(pdf_path, dpi=DPI):
    """Return a PIL Image for the first page of *pdf_path*."""
    for renderer in (_render_via_fitz, _render_via_pdf2image, _render_via_ghostscript):
        try:
            return renderer(pdf_path, dpi)
        except Exception:
            continue
    raise RuntimeError(
        f"Cannot render {pdf_path}. Install PyMuPDF, pdf2image, or ghostscript."
    )


# ── filename parsing ──────────────────────────────────────────────────────────

# Filenames written by the sbatch script:
#   {TASK_ID}_yield_v2_{cen_name}_{pt_name}_q2bin{q}.pdf
#   {TASK_ID}_yield_v3_{cen_name}_{pt_name}_q3bin{q}.pdf
# The greedy (.+) before _q2bin/_q3bin captures "{cen_name}_{pt_name}".
_PAT_V2 = re.compile(r"^\d+_yield_v2_(.+)_q2bin(\d+)\.pdf$")
_PAT_V3 = re.compile(r"^\d+_yield_v3_(.+)_q3bin(\d+)\.pdf$")


def parse_files(input_dir):
    """
    Scan *input_dir* and return:
      v2_groups  dict[cen_pt_key] → sorted list of (q_bin, full_path)
      v3_groups  dict[cen_pt_key] → sorted list of (q_bin, full_path)
    """
    v2_groups = defaultdict(list)
    v3_groups = defaultdict(list)

    for fname in os.listdir(input_dir):
        m = _PAT_V2.match(fname)
        if m:
            v2_groups[m.group(1)].append(
                (int(m.group(2)), os.path.join(input_dir, fname))
            )
            continue
        m = _PAT_V3.match(fname)
        if m:
            v3_groups[m.group(1)].append(
                (int(m.group(2)), os.path.join(input_dir, fname))
            )

    for grp in (v2_groups, v3_groups):
        for k in grp:
            grp[k].sort(key=lambda x: x[0])

    return v2_groups, v3_groups


# ── natural sort (handles "cent0to10" < "cent10to20" correctly) ───────────────

def _natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


# ── page builder ──────────────────────────────────────────────────────────────

def build_summary_pdf(groups, output_path, vn_label, qbin_label, dpi):
    """
    Write a multi-page PDF to *output_path*.
    One page per key in *groups*; each page = GRID_ROWS × GRID_COLS sub-panels.
    """
    sorted_keys = sorted(groups.keys(), key=_natural_key)
    n_pages = len(sorted_keys)
    print(f"  → {n_pages} pages")

    with PdfPages(output_path) as pdf:
        for page_idx, key in enumerate(sorted_keys, 1):
            entries = groups[key]   # already sorted by q bin

            fig, axes = plt.subplots(
                GRID_ROWS, GRID_COLS,
                figsize=(GRID_COLS * 3.6, GRID_ROWS * 3.6),
            )
            axes_flat = axes.flatten()

            for panel in range(GRID_ROWS * GRID_COLS):
                ax = axes_flat[panel]
                if panel < len(entries):
                    q_bin, fpath = entries[panel]
                    if os.path.isfile(fpath):
                        try:
                            img = render_pdf_page(fpath, dpi=dpi)
                            ax.imshow(img, aspect="auto", interpolation="lanczos")
                            ax.set_title(f"{qbin_label} {q_bin}",
                                         fontsize=8, pad=2)
                        except Exception as exc:
                            ax.text(0.5, 0.5, f"render error\n{exc}",
                                    ha="center", va="center",
                                    transform=ax.transAxes,
                                    fontsize=6, color="red", wrap=True)
                    else:
                        ax.text(0.5, 0.5, "file missing",
                                ha="center", va="center",
                                transform=ax.transAxes,
                                fontsize=7, color="gray")
                else:
                    ax.set_visible(False)
                ax.axis("off")

            # Page title: "cent0to10  pt1to2" style
            title = f"{vn_label}    {key.replace('_', '   ')}"
            fig.suptitle(title, fontsize=10, y=1.005)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fig.tight_layout(rect=[0, 0, 1, 1])

            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            print(f"    page {page_idx}/{n_pages}: {key}")

    print(f"  Saved → {output_path}\n")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    input_dir  = sys.argv[1] if len(sys.argv) > 1 else INPUT_DIR
    output_dir = sys.argv[2] if len(sys.argv) > 2 else (
                 sys.argv[1] if len(sys.argv) > 1 else OUTPUT_DIR)

    input_dir  = str(Path(input_dir).expanduser().resolve())
    output_dir = str(Path(output_dir).expanduser().resolve())
    os.makedirs(output_dir, exist_ok=True)

    print(f"Input  : {input_dir}")
    print(f"Output : {output_dir}\n")

    v2_groups, v3_groups = parse_files(input_dir)
    print(f"v2 (cent, pT) combinations found : {len(v2_groups)}")
    print(f"v3 (cent, pT) combinations found : {len(v3_groups)}\n")

    if not v2_groups and not v3_groups:
        sys.exit("No matching PDF files found. Check INPUT_DIR and filename pattern.")

    if v2_groups:
        out_v2 = os.path.join(output_dir, "yield_v2_summary.pdf")
        print(f"Building v2 summary PDF ({len(v2_groups)} pages) ...")
        build_summary_pdf(v2_groups, out_v2,
                          vn_label="v2 yield vs SP bin index",
                          qbin_label="q2 bin", dpi=DPI)

    if v3_groups:
        out_v3 = os.path.join(output_dir, "yield_v3_summary.pdf")
        print(f"Building v3 summary PDF ({len(v3_groups)} pages) ...")
        build_summary_pdf(v3_groups, out_v3,
                          vn_label="v3 yield vs SP bin index",
                          qbin_label="q3 bin", dpi=DPI)


if __name__ == "__main__":
    main()
