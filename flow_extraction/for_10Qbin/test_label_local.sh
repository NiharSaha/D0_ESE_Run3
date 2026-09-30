#!/bin/bash
# Local test: verifies the eso-pic page label renders correctly in pdfjam output.
# Creates 4 tiny dummy PDFs, assembles them into a 2x2 grid page with a label,
# and opens the result for inspection.
# Usage: bash test_label_local.sh

set -e

TESTDIR="$(mktemp -d /tmp/pdfjam_label_test_XXXXXX)"
echo "Working in: $TESTDIR"
trap 'rm -rf "$TESTDIR"' EXIT

# --- 1. Check dependencies ---
for cmd in pdflatex pdfjam; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: '$cmd' not found. Load texlive first (e.g. module load texlive)."
        exit 1
    fi
done

# --- 2. Create 4 minimal single-page PDFs as stand-ins for mass-fit plots ---
echo "Creating dummy PDFs..."
for i in 1 2 3 4; do
    pdflatex -interaction=batchmode -output-directory="$TESTDIR" \
        <(printf '\\documentclass{minimal}\\begin{document}\\Huge Plot %d\\end{document}' "$i") \
        > /dev/null 2>&1 || true
done

# pdflatex doesn't accept process substitution well; use temp .tex files instead
for i in 1 2 3 4; do
    cat > "$TESTDIR/plot${i}.tex" <<EOF
\documentclass{minimal}
\begin{document}
\Huge Plot ${i}
\end{document}
EOF
    pdflatex -interaction=batchmode -output-directory="$TESTDIR" "$TESTDIR/plot${i}.tex" > /dev/null 2>&1
done

# Verify dummy PDFs were created
DUMMY_PDFS=( "$TESTDIR"/plot{1,2,3,4}.pdf )
for f in "${DUMMY_PDFS[@]}"; do
    [ -f "$f" ] || { echo "ERROR: Failed to create $f"; exit 1; }
done
echo "Dummy PDFs created: ${DUMMY_PDFS[*]}"

# --- 3. Run the exact pdfjam command from the script ---
LABEL="Cen: 0  |  pT: pT3to4  |  Q2bin: 5  (4 SP bins)"
OUTPUT="$TESTDIR/test_output.pdf"

echo "Running pdfjam with eso-pic label..."
TMPDIR="$TESTDIR" pdfjam "${DUMMY_PDFS[@]}" \
    --nup "2x2" \
    --papersize '{20in,20in}' \
    --scale 0.92 \
    --preamble "\\usepackage{eso-pic}" \
    --pagecommand "{\\AddToShipoutPictureFG*{\\AtPageUpperLeft{\\raisebox{-0.5in}{\\makebox[\\paperwidth][c]{\\Large\\bfseries ${LABEL}}}}}}" \
    --outfile "$OUTPUT" 2>"$TESTDIR/pdfjam_stderr.log"

if [ ! -f "$OUTPUT" ]; then
    echo "ERROR: pdfjam did not produce output. pdfjam stderr:"
    cat "$TESTDIR/pdfjam_stderr.log"
    exit 1
fi

echo ""
echo "SUCCESS: output PDF created."

# Copy to current dir so it survives the temp-dir cleanup
cp "$OUTPUT" "./test_label_output.pdf"
echo "Saved to: $(pwd)/test_label_output.pdf"
echo "Open it and confirm the label is visible at the top of the page."

# Try to open with a PDF viewer if available
for viewer in evince okular xdg-open; do
    if command -v "$viewer" &>/dev/null; then
        echo "Opening with $viewer ..."
        "$viewer" "./test_label_output.pdf" &
        break
    fi
done
