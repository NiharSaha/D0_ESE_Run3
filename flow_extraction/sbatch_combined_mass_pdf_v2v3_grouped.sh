#!/bin/bash
#SBATCH --job-name=plot_stitch
#SBATCH --output=stitch_%j.out
#SBATCH --error=stitch_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH -A physics
#SBATCH -p cpu
#SBATCH -q standby

# Load necessary modules
module --force purge
module load texlive parallel

# --- CONFIGURATION ---
PLOT_SOURCE_DIR="/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Flow_output_FullStat_May19_v0/output/prompt_mass_plot_withchi2_sigma"
GRID_V2="7x6"   # 7 cols x 6 rows = 42 = N_VBINS_V2 (one page per group)
GRID_V3="8x6"   # 8 cols x 6 rows = 48 = N_VBINS_V3 (one page per group)
FINAL_OUTPUT_V2="All_Plots_Combined_v2_FullStat_May19_v1.pdf"
FINAL_OUTPUT_V3="All_Plots_Combined_v3_FullStat_May19_v1.pdf"
TEMP_DIR_V2="temp_parallel_v2_$SLURM_JOB_ID"
TEMP_DIR_V3="temp_parallel_v3_$SLURM_JOB_ID"

# Fix Perl environment for GNU Parallel
unset PERL5LIB

# Setup
mkdir -p "$TEMP_DIR_V2" "$TEMP_DIR_V3"

# --- Binning definitions ---
CEN_NAMES="cent0to10 cent10to20 cent20to30 cent30to40 cent40to50 cent50to80"
# v2: 9 pT bins matching pt_name_v2[] in Analysis_bin.h
PT_NAMES_V2="pT2to3 pT3to4 pT4to5 pT5to6 pT6to8 pT8to10 pT10to15 pT15to30 pT30to100"
# v3: 7 pT bins matching pt_name_v3[] in Analysis_bin.h
PT_NAMES_V3="pT2to4 pT4to6 pT6to8 pT8to10 pT10to20 pT20to50 pT50to100"
N_QBINS=10
N_VBINS_V2=42   # N_VBINS_V2 in Analysis_bin.h
N_VBINS_V3=48   # N_VBINS_V3 in Analysis_bin.h

echo "Job started on $(hostname)"

# --- Worker function: one page per (centrality, pT, qbin) combination ---
# Args: cen_idx cen_name pt iq n_vbins vtag qtag temp_dir job_id source_dir grid seq_num
process_group() {
    unset PERL5LIB
    local cen_idx=$1
    local cen_name=$2
    local pt=$3
    local iq=$4
    local n_vbins=$5
    local vtag=$6        # "v2" or "v3"
    local qtag=$7        # "q2" or "q3"
    local temp_dir=$8
    local job_id=$9
    local source_dir=${10}
    local grid=${11}
    local seq_num=${12}

    local output_page
    output_page=$(printf "%s/page_%06d.pdf" "$temp_dir" "$seq_num")
    local file_list="list_${job_id}_${vtag}_${cen_name}_${pt}_${qtag}bin_${iq}.txt"

    # Collect all SP-bin (vbin) plots for this (cen, pT, qbin) combination, in order
    > "$file_list"
    for (( iv=0; iv<n_vbins; iv++ )); do
        local found
        found=$(ls "${source_dir}/${cen_idx}_hmassfit_${pt}_${cen_name}_${qtag}bin_${iq}_${vtag}bin_${iv}.pdf" 2>/dev/null | head -1)
        [ -n "$found" ] && echo "$found" >> "$file_list"
    done

    local count
    count=$(wc -l < "$file_list")

    if [ "$count" -eq 0 ]; then
        rm -f "$file_list"
        return 0
    fi

    # Build physics-style label with LaTeX math formatting.
    # Parse centrality from cen_name: "cent0to10" -> lo=0, hi=10
    local cent_stripped="${cen_name#cent}"
    local cent_lo="${cent_stripped%%to*}"
    local cent_hi="${cent_stripped##*to}"
    # Parse pT from pt name: "pT2to3" -> lo=2, hi=3
    local pt_stripped="${pt#pT}"
    local pt_lo="${pt_stripped%%to*}"
    local pt_hi="${pt_stripped##*to}"
    local qnum="${qtag: -1}"
    # \$ in bash double-quoted string gives literal $ for LaTeX math mode
    local label="\$${cent_lo} < Cent < ${cent_hi}\\%\$  |  \$${pt_lo} < p_{T} < ${pt_hi}\$ GeV/c  |  \$q_{${qnum}}\$ bin = ${iq}  (${count} SP bins)"

    # Paper size matched to grid aspect ratio (square ROOT plots) to minimise blank space.
    # v2 grid 7x6: width:height = 7:6; add 1in height for title -> {28in,25in}
    # v3 grid 8x6: width:height = 8:6; add 1in height for title -> {32in,25in}
    local papersize
    if [ "$vtag" = "v2" ]; then
        papersize='{28in,25in}'
    else
        papersize='{32in,25in}'
    fi

    TMPDIR="$temp_dir" pdfjam $(cat "$file_list") \
        --nup "$grid" \
        --papersize "$papersize" \
        --scale 0.95 \
        --preamble "\\usepackage{eso-pic}" \
        --pagecommand "{\\AddToShipoutPictureFG*{\\AtPageUpperLeft{\\raisebox{-0.85in}{\\makebox[\\paperwidth][c]{\\fontsize{36}{44}\\selectfont\\bfseries ${label}}}}}}" \
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
cen_idx=0
for cen_name in $CEN_NAMES; do
  for pt in $PT_NAMES_V2; do
    for (( iq=0; iq<N_QBINS; iq++ )); do
      echo "$cen_idx $cen_name $pt $iq" >> "$JOBS_V2"
    done
  done
  (( cen_idx++ ))
done
cen_idx=0
for cen_name in $CEN_NAMES; do
  for pt in $PT_NAMES_V3; do
    for (( iq=0; iq<N_QBINS; iq++ )); do
      echo "$cen_idx $cen_name $pt $iq" >> "$JOBS_V3"
    done
  done
  (( cen_idx++ ))
done

TOTAL_GROUPS_V2=$(wc -l < "$JOBS_V2")
TOTAL_GROUPS_V3=$(wc -l < "$JOBS_V3")
echo "V2 groups (cen x pT x Q2bin): $TOTAL_GROUPS_V2"
echo "V3 groups (cen x pT x Q3bin): $TOTAL_GROUPS_V3"
echo "Using $SLURM_CPUS_PER_TASK cores ($HALF_CORES per stream)."

# --- Process V2 in background ---
(
    echo "[V2] Processing $TOTAL_GROUPS_V2 groups (cen x pT x Q2bin), each page = all v2bins..."
    parallel --env process_group -j "$HALF_CORES" \
        --colsep ' ' \
        process_group {1} {2} {3} {4} "$N_VBINS_V2" "v2" "q2" \
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
    echo "[V3] Processing $TOTAL_GROUPS_V3 groups (cen x pT x Q3bin), each page = all v3bins..."
    parallel --env process_group -j "$HALF_CORES" \
        --colsep ' ' \
        process_group {1} {2} {3} {4} "$N_VBINS_V3" "v3" "q3" \
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
echo "V2 PDF : $FINAL_OUTPUT_V2  ($TOTAL_GROUPS_V2 pages)"
echo "V3 PDF : $FINAL_OUTPUT_V3  ($TOTAL_GROUPS_V3 pages)"
echo "============================="
