#!/bin/bash
#SBATCH --job-name=plot_stitch
#SBATCH --output=stitch_%j.out
#SBATCH --error=stitch_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -A physics
#SBATCH -p cpu
#SBATCH -q standby

# Load necessary modules
module --force purge
module load texlive parallel
unset PERL5LIB

# --- CONFIGURATION ---
PLOT_SOURCE_DIR="/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/D0_Flow_12Qbin_output_FullStat_Jun19_v2/output/prompt_mass_plot_withchi2_sigma"
GRID_V2="7x6"   # 7 cols x 6 rows = 42 = N_VBINS_V2 (one page per group)
GRID_V3="8x6"   # 8 cols x 6 rows = 48 = N_VBINS_V3 (one page per group)
FINAL_OUTPUT_V2="All_Plots_Combined_v2_FullStat_Jun19_12NUQbin.pdf"
FINAL_OUTPUT_V3="All_Plots_Combined_v3_FullStat_Jun19_12NUQbin.pdf"
TEMP_DIR="temp_stitch_$SLURM_JOB_ID"
mkdir -p "$TEMP_DIR"

# --- Binning definitions ---
CEN_NAMES="cent0to10 cent10to20 cent20to30 cent30to40 cent40to50"
PT_NAMES_V2="pT2to3 pT3to4 pT4to5 pT5to6 pT6to8 pT8to10 pT10to15 pT15to30 pT30to100"
PT_NAMES_V3="pT2to4 pT4to6 pT6to8 pT8to10 pT10to20 pT20to50 pT50to100"
N_QBINS=12
N_VBINS_V2=42
N_VBINS_V3=48

echo "Job started on $(hostname)"
echo "Using $SLURM_CPUS_PER_TASK cores (sequential: V2 then V3)."

# --- Worker function: one page per (centrality, pT, qbin) combination ---
# Args: cen_idx cen_name pt iq n_vbins vtag qtag temp_dir job_id source_dir grid seq_num
process_group() {
    unset PERL5LIB
    local cen_idx=$1 cen_name=$2 pt=$3 iq=$4 n_vbins=$5
    local vtag=$6 qtag=$7 temp_dir=$8 job_id=$9
    local source_dir=${10} grid=${11} seq_num=${12}

    # vtag-prefixed page name avoids any collision between v2/v3 page files
    local output_page
    output_page=$(printf "%s/${vtag}_page_%06d.pdf" "$temp_dir" "$seq_num")
    local file_list="${temp_dir}/list_${vtag}_${seq_num}.txt"

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

    # Build physics-style label
    local cent_stripped="${cen_name#cent}"
    local cent_lo="${cent_stripped%%to*}" cent_hi="${cent_stripped##*to}"
    local pt_stripped="${pt#pT}"
    local pt_lo="${pt_stripped%%to*}" pt_hi="${pt_stripped##*to}"
    local qnum="${qtag: -1}"
    local label="\$${cent_lo} < Cent < ${cent_hi}\\%\$  |  \$${pt_lo} < p_{T} < ${pt_hi}\$ GeV/c  |  \$q_{${qnum}}\$ bin = ${iq}  (${count} SP bins)"

    local papersize
    [ "$vtag" = "v2" ] && papersize='{28in,25in}' || papersize='{32in,25in}'

    TMPDIR="$temp_dir" pdfjam $(cat "$file_list") \
        --nup "$grid" \
        --papersize "$papersize" \
        --scale 0.95 \
        --preamble "\\usepackage{eso-pic}" \
        --pagecommand "{\\AddToShipoutPictureFG*{\\AtPageUpperLeft{\\raisebox{-0.85in}{\\makebox[\\paperwidth][c]{\\fontsize{36}{44}\\selectfont\\bfseries ${label}}}}}}" \
        --outfile "$output_page" 2>> "${temp_dir}/pdfjam_${vtag}.log" > /dev/null

    rm -f "$file_list"
}
export -f process_group

# --- Helper: merge pages using pdfunite (fast) or pdfjam (fallback) ---
merge_pages() {
    local vtag=$1 final_output=$2 temp_dir=$3
    local pages=( "${temp_dir}/${vtag}_page_"*.pdf )
    local count=${#pages[@]}
    if [ "$count" -eq 0 ] || [ ! -f "${pages[0]}" ]; then
        echo "[${vtag^^}] No pages generated, skipping merge."
        return 1
    fi
    echo "[${vtag^^}] Merging $count pages into $final_output..."
    if command -v pdfunite &>/dev/null; then
        pdfunite "${pages[@]}" "$final_output"
    else
        pdfjam "${pages[@]}" --outfile "$final_output" > /dev/null 2>&1
    fi
    echo "[${vtag^^}] Done: $final_output"
}

# --- Build ordered job lists (ordering: cen -> pT -> qbin) ---
JOBS_V2="$TEMP_DIR/jobs_v2.txt"
JOBS_V3="$TEMP_DIR/jobs_v3.txt"
> "$JOBS_V2"; > "$JOBS_V3"
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
echo "V2 groups: $TOTAL_GROUPS_V2 | V3 groups: $TOTAL_GROUPS_V3"

# =============================================================
# V2 — all cores, then V3 — all cores  (sequential, not split)
# This is more efficient: no idle cores while the other stream
# is still running.
# =============================================================

echo "[V2] Processing $TOTAL_GROUPS_V2 groups with $SLURM_CPUS_PER_TASK cores..."
parallel --env process_group -j "$SLURM_CPUS_PER_TASK" --colsep ' ' \
    process_group {1} {2} {3} {4} "$N_VBINS_V2" "v2" "q2" \
        "$TEMP_DIR" "$SLURM_JOB_ID" "$PLOT_SOURCE_DIR" "$GRID_V2" {#} \
    :::: "$JOBS_V2"
[ -s "$TEMP_DIR/pdfjam_v2.log" ] && { echo "[V2] pdfjam errors (first 20 lines):"; head -20 "$TEMP_DIR/pdfjam_v2.log"; }
merge_pages "v2" "$FINAL_OUTPUT_V2" "$TEMP_DIR"

echo "[V3] Processing $TOTAL_GROUPS_V3 groups with $SLURM_CPUS_PER_TASK cores..."
parallel --env process_group -j "$SLURM_CPUS_PER_TASK" --colsep ' ' \
    process_group {1} {2} {3} {4} "$N_VBINS_V3" "v3" "q3" \
        "$TEMP_DIR" "$SLURM_JOB_ID" "$PLOT_SOURCE_DIR" "$GRID_V3" {#} \
    :::: "$JOBS_V3"
[ -s "$TEMP_DIR/pdfjam_v3.log" ] && { echo "[V3] pdfjam errors (first 20 lines):"; head -20 "$TEMP_DIR/pdfjam_v3.log"; }
merge_pages "v3" "$FINAL_OUTPUT_V3" "$TEMP_DIR"

# --- Cleanup ---
rm -rf "$TEMP_DIR"

echo "============================="
echo "All done."
echo "V2 PDF : $FINAL_OUTPUT_V2  ($TOTAL_GROUPS_V2 pages)"
echo "V3 PDF : $FINAL_OUTPUT_V3  ($TOTAL_GROUPS_V3 pages)"
echo "============================="
