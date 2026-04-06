#!/bin/bash
#SBATCH --job-name=plot_stitch
#SBATCH --output=stitch_%j.out
#SBATCH --error=stitch_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH -A physics
#SBATCH -p cpu
#SBATCH -q standby

# Load necessary modules
module --force purge
module load texlive parallel

# --- CONFIGURATION ---
PLOT_SOURCE_DIR="/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Flow_output_Apr1_v0/output/prompt_mass_plot_withchi2_sigma"
GRID_V2="6x6"
GRID_V3="7x7"
FINAL_OUTPUT_V2="All_Plots_Combined_v2_Apr6_v0.pdf"
FINAL_OUTPUT_V3="All_Plots_Combined_v3_Apr6_v0.pdf"
TEMP_DIR_V2="temp_parallel_v2_$SLURM_JOB_ID"
TEMP_DIR_V3="temp_parallel_v3_$SLURM_JOB_ID"

# Fix Perl environment for GNU Parallel
unset PERL5LIB

# Setup
mkdir -p "$TEMP_DIR_V2" "$TEMP_DIR_V3"

# --- Binning definitions ---
CEN_IDS="0 1 2 3 4 5"
PT_NAMES="pT1to2 pT2to3 pT3to4 pT4to5 pT5to6 pT6to8 pT8to10 pT10to15 pT15to20 pT20to40 pT40to60 pT60to100"
N_QBINS=10
N_VBINS=62

echo "Job started on $(hostname)"

# --- Worker function: one page per (centrality, pT, qbin) combination ---
# Args: cen_id pt iq n_vbins vtag qtag temp_dir job_id source_dir grid seq_num
process_group() {
    unset PERL5LIB
    local cen_id=$1
    local pt=$2
    local iq=$3
    local n_vbins=$4
    local vtag=$5        # "v2" or "v3"
    local qtag=$6        # "q2" or "q3"
    local temp_dir=$7
    local job_id=$8
    local source_dir=$9
    local grid=${10}
    local seq_num=${11}

    local output_page
    output_page=$(printf "%s/page_%06d.pdf" "$temp_dir" "$seq_num")
    local file_list="list_${job_id}_${vtag}_${cen_id}_${pt}_${qtag}bin_${iq}.txt"

    # Collect all SP-bin (vbin) plots for this (cen, pT, qbin) combination, in order
    > "$file_list"
    for (( iv=0; iv<n_vbins; iv++ )); do
        local found
        found=$(ls "${source_dir}/${cen_id}_hmassfit_${pt}_"*"_${qtag}bin_${iq}_${vtag}bin_${iv}.pdf" 2>/dev/null | head -1)
        [ -n "$found" ] && echo "$found" >> "$file_list"
    done

    local count
    count=$(wc -l < "$file_list")

    if [ "$count" -eq 0 ]; then
        rm -f "$file_list"
        return 0
    fi

    # Build page label (flow harmonic number extracted from vtag last char)
    local qnum="${qtag: -1}"
    local label="Cen: ${cen_id}  |  pT: ${pt}  |  Q${qnum}bin: ${iq}  (${count} SP bins)"

    # Create the grid page with a centered title at the top of the 20in x 20in page.
    # TMPDIR is set to our writable scratch space so pdflatex temp files land there.
    # tikz overlay (remember picture + overlay) places the title as a zero-size node
    # anchored to page.north — it does not affect the pdfpages grid layout at all.
    # \usepackage{tikz} is loaded via --preamble; tikz is always present in texlive 2022.
    TMPDIR="$temp_dir" pdfjam $(cat "$file_list") \
        --nup "$grid" \
        --papersize '{20in,20in}' \
        --scale 0.92 \
        --preamble "\\usepackage{eso-pic}" \
        --pagecommand "{\\AddToShipoutPictureFG*{\\AtPageUpperLeft{\\raisebox{-0.5in}{\\makebox[\\paperwidth][c]{\\Large\\bfseries ${label}}}}}}" \
        --outfile "$output_page" 2>> "${temp_dir}/pdfjam_errors.log" > /dev/null

    # On failure: find and save the pdflatex log for diagnostics
    if [ ! -f "$output_page" ]; then
        local pdflatex_log
        pdflatex_log=$(find "$temp_dir" -name "a.log" 2>/dev/null | head -1)
        if [ -n "$pdflatex_log" ]; then
            echo "=== pdflatex log: ${label} ===" >> "${temp_dir}/pdflatex_full.log"
            cat "$pdflatex_log"              >> "${temp_dir}/pdflatex_full.log"
        fi
    fi

    rm -f "$file_list"
}
export -f process_group

# Split cores evenly between v2 and v3 streams
HALF_CORES=$(( SLURM_CPUS_PER_TASK / 2 ))

# --- Build ordered job lists (ordering: cen -> pT -> qbin) ---
JOBS_V2="jobs_v2_$SLURM_JOB_ID.txt"
JOBS_V3="jobs_v3_$SLURM_JOB_ID.txt"
> "$JOBS_V2"
> "$JOBS_V3"
for cen_id in $CEN_IDS; do
  for pt in $PT_NAMES; do
    for (( iq=0; iq<N_QBINS; iq++ )); do
      echo "$cen_id $pt $iq" >> "$JOBS_V2"
      echo "$cen_id $pt $iq" >> "$JOBS_V3"
    done
  done
done

TOTAL_GROUPS=$(wc -l < "$JOBS_V2")
echo "Groups per stream (cen x pT x Qbin): $TOTAL_GROUPS"
echo "Using $SLURM_CPUS_PER_TASK cores ($HALF_CORES per stream)."

# --- Process V2 in background ---
(
    echo "[V2] Processing $TOTAL_GROUPS groups (cen x pT x Q2bin), each page = all v2bins..."
    parallel --env process_group -j "$HALF_CORES" \
        --colsep ' ' \
        process_group {1} {2} {3} "$N_VBINS" "v2" "q2" \
            "$TEMP_DIR_V2" "$SLURM_JOB_ID" "$PLOT_SOURCE_DIR" "$GRID_V2" {#} \
        :::: "$JOBS_V2"

    if [ -s "${TEMP_DIR_V2}/pdflatex_full.log" ]; then
        echo "[V2] pdflatex error log (first 60 lines):"
        head -60 "${TEMP_DIR_V2}/pdflatex_full.log"
    elif [ -s "${TEMP_DIR_V2}/pdfjam_errors.log" ]; then
        echo "[V2] pdfjam stderr (first 20 lines):"
        head -20 "${TEMP_DIR_V2}/pdfjam_errors.log"
    fi

    PAGE_COUNT=$(ls "$TEMP_DIR_V2"/page_*.pdf 2>/dev/null | wc -l)
    if [ "$PAGE_COUNT" -eq 0 ]; then
        echo "[V2] No pages generated, skipping merge."
    else
        echo "[V2] Merging $PAGE_COUNT pages into $FINAL_OUTPUT_V2..."
        pdfjam "$TEMP_DIR_V2"/page_*.pdf --outfile "$FINAL_OUTPUT_V2" > /dev/null 2>&1
        echo "[V2] Done: $FINAL_OUTPUT_V2"
    fi
) &
PID_V2=$!

# --- Process V3 in background ---
(
    echo "[V3] Processing $TOTAL_GROUPS groups (cen x pT x Q3bin), each page = all v3bins..."
    parallel --env process_group -j "$HALF_CORES" \
        --colsep ' ' \
        process_group {1} {2} {3} "$N_VBINS" "v3" "q3" \
            "$TEMP_DIR_V3" "$SLURM_JOB_ID" "$PLOT_SOURCE_DIR" "$GRID_V3" {#} \
        :::: "$JOBS_V3"

    if [ -s "${TEMP_DIR_V3}/pdflatex_full.log" ]; then
        echo "[V3] pdflatex error log (first 60 lines):"
        head -60 "${TEMP_DIR_V3}/pdflatex_full.log"
    elif [ -s "${TEMP_DIR_V3}/pdfjam_errors.log" ]; then
        echo "[V3] pdfjam stderr (first 20 lines):"
        head -20 "${TEMP_DIR_V3}/pdfjam_errors.log"
    fi

    PAGE_COUNT=$(ls "$TEMP_DIR_V3"/page_*.pdf 2>/dev/null | wc -l)
    if [ "$PAGE_COUNT" -eq 0 ]; then
        echo "[V3] No pages generated, skipping merge."
    else
        echo "[V3] Merging $PAGE_COUNT pages into $FINAL_OUTPUT_V3..."
        pdfjam "$TEMP_DIR_V3"/page_*.pdf --outfile "$FINAL_OUTPUT_V3" > /dev/null 2>&1
        echo "[V3] Done: $FINAL_OUTPUT_V3"
    fi
) &
PID_V3=$!

# --- Wait for both streams to finish ---
wait $PID_V2
STATUS_V2=$?
wait $PID_V3
STATUS_V3=$?

[ $STATUS_V2 -ne 0 ] && echo "ERROR: V2 processing failed with status $STATUS_V2"
[ $STATUS_V3 -ne 0 ] && echo "ERROR: V3 processing failed with status $STATUS_V3"

# --- Cleanup ---
rm -rf "$TEMP_DIR_V2" "$TEMP_DIR_V3" "$JOBS_V2" "$JOBS_V3"

echo "============================="
echo "All done."
echo "V2 PDF : $FINAL_OUTPUT_V2  ($TOTAL_GROUPS pages)"
echo "V3 PDF : $FINAL_OUTPUT_V3  ($TOTAL_GROUPS pages)"
echo "============================="
