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
module purge
module load texlive parallel

# --- CONFIGURATION ---
PLOT_SOURCE_DIR="/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Flow_output_PubResults_Mar12_v1/output/prompt_mass_plot_withchi2_sigma"
PLOTS_PER_PAGE=100
GRID="10x10"
FINAL_OUTPUT_V2="All_Plots_Combined_v2_Mar12.pdf"
FINAL_OUTPUT_V3="All_Plots_Combined_v3_Mar12.pdf"
TEMP_DIR_V2="temp_parallel_v2_$SLURM_JOB_ID"
TEMP_DIR_V3="temp_parallel_v3_$SLURM_JOB_ID"

# Fix Perl environment for GNU Parallel
unset PERL5LIB

# Setup
mkdir -p "$TEMP_DIR_V2" "$TEMP_DIR_V3"

# --- Binning definitions ---
CEN_IDS="0 1 2"                          # numeric prefix in filename
CEN_NAMES="cent0to10 cent10to30 cent30to50"          # short cen name used in filename
PT_NAMES="pT1to2 pT2to3 pT3to4 pT4to5 pT5to6 pT6to8 pT8to10 pT10to15 pT15to20 pT20to40 pT40to60 pT60to100"
N_VBINS_V2=63
N_VBINS_V3=67

ORDERED_LIST_V2="ordered_files_v2_$SLURM_JOB_ID.txt"
ORDERED_LIST_V3="ordered_files_v3_$SLURM_JOB_ID.txt"
> "$ORDERED_LIST_V2"
> "$ORDERED_LIST_V3"

# --- Build ordered file list for V2: cent -> pt -> v2bin ---
# Pattern: {i_cen}_hmassfit_{pt_name}_{cen_short}_v2bin_{iv}.pdf
echo "Building v2 file list..."
for i_cen in 0 1 2; do
  cen_arr=($CEN_NAMES)
  cen="${cen_arr[$i_cen]}"
  for pt in $PT_NAMES; do
    for (( iv=0; iv<N_VBINS_V2; iv++ )); do
      found=$(ls "$PLOT_SOURCE_DIR"/${i_cen}_hmassfit_${pt}_${cen}_v2bin_${iv}.pdf 2>/dev/null | head -1)
      if [ -n "$found" ]; then
        echo "$found" >> "$ORDERED_LIST_V2"
      fi
    done
  done
done

# --- Build ordered file list for V3: cent -> pt -> v3bin ---
# Pattern: {i_cen}_hmassfit_{pt_name}_{cen_short}_v3bin_{iv}.pdf
echo "Building v3 file list..."
for i_cen in 0 1 2; do
  cen_arr=($CEN_NAMES)
  cen="${cen_arr[$i_cen]}"
  for pt in $PT_NAMES; do
    for (( iv=0; iv<N_VBINS_V3; iv++ )); do
      found=$(ls "$PLOT_SOURCE_DIR"/${i_cen}_hmassfit_${pt}_${cen}_v3bin_${iv}.pdf 2>/dev/null | head -1)
      if [ -n "$found" ]; then
        echo "$found" >> "$ORDERED_LIST_V3"
      fi
    done
  done
done

TOTAL_V2=$(wc -l < "$ORDERED_LIST_V2")
TOTAL_V3=$(wc -l < "$ORDERED_LIST_V3")
echo "Job started on $(hostname)"
echo "V2 plots found: $TOTAL_V2"
echo "V3 plots found: $TOTAL_V3"
echo "Using $SLURM_CPUS_PER_TASK cores ($((SLURM_CPUS_PER_TASK / 2)) per stream)."

# --- Guard: exit if no files found ---
if [ "$TOTAL_V2" -eq 0 ] && [ "$TOTAL_V3" -eq 0 ]; then
  echo "ERROR: No PDF files found in $PLOT_SOURCE_DIR. Exiting."
  exit 1
fi

# --- Worker function ---
process_batch() {
    unset PERL5LIB
    local start_line=$1
    local batch_size=$2
    local grid=$3
    local temp_dir=$4
    local job_id=$5
    local ordered_list=$6
    local tag=$7
    local batch_num=$((start_line / batch_size))
    local output_page
    output_page=$(printf "%s/page_%04d.pdf" "$temp_dir" "$batch_num")

    sed -n "$((start_line + 1)),$((start_line + batch_size))p" "$ordered_list" \
        > "list_${job_id}_${tag}_${batch_num}.txt"

    pdfjam $(cat "list_${job_id}_${tag}_${batch_num}.txt") \
        --nup "$grid" --papersize '{20in,20in}' \
        --outfile "$output_page" > /dev/null 2>&1

    rm -f "list_${job_id}_${tag}_${batch_num}.txt"
}
export -f process_batch

# Split cores evenly between v2 and v3
HALF_CORES=$(( SLURM_CPUS_PER_TASK / 2 ))

# --- Process V2 in background ---
(
  if [ "$TOTAL_V2" -gt 0 ]; then
    echo "[V2] Processing $TOTAL_V2 plots in batches of $PLOTS_PER_PAGE..."
    seq 0 $PLOTS_PER_PAGE $((TOTAL_V2 - 1)) | \
      parallel -j "$HALF_CORES" \
      process_batch {} $PLOTS_PER_PAGE $GRID "$TEMP_DIR_V2" "$SLURM_JOB_ID" "$ORDERED_LIST_V2" "v2"

    echo "[V2] Merging pages into $FINAL_OUTPUT_V2..."
    pdfjam "$TEMP_DIR_V2"/page_*.pdf --outfile "$FINAL_OUTPUT_V2" > /dev/null 2>&1
    echo "[V2] Done: $FINAL_OUTPUT_V2"
  else
    echo "[V2] No plots found, skipping."
  fi
) &
PID_V2=$!

# --- Process V3 in background ---
(
  if [ "$TOTAL_V3" -gt 0 ]; then
    echo "[V3] Processing $TOTAL_V3 plots in batches of $PLOTS_PER_PAGE..."
    seq 0 $PLOTS_PER_PAGE $((TOTAL_V3 - 1)) | \
      parallel -j "$HALF_CORES" \
      process_batch {} $PLOTS_PER_PAGE $GRID "$TEMP_DIR_V3" "$SLURM_JOB_ID" "$ORDERED_LIST_V3" "v3"

    echo "[V3] Merging pages into $FINAL_OUTPUT_V3..."
    pdfjam "$TEMP_DIR_V3"/page_*.pdf --outfile "$FINAL_OUTPUT_V3" > /dev/null 2>&1
    echo "[V3] Done: $FINAL_OUTPUT_V3"
  else
    echo "[V3] No plots found, skipping."
  fi
) &
PID_V3=$!

# --- Wait for both streams to finish ---
wait $PID_V2
STATUS_V2=$?
wait $PID_V3
STATUS_V3=$?

if [ $STATUS_V2 -ne 0 ]; then
  echo "ERROR: V2 processing failed with status $STATUS_V2"
fi
if [ $STATUS_V3 -ne 0 ]; then
  echo "ERROR: V3 processing failed with status $STATUS_V3"
fi

# --- Cleanup ---
rm -rf "$TEMP_DIR_V2" "$TEMP_DIR_V3" "$ORDERED_LIST_V2" "$ORDERED_LIST_V3"

echo "============================="
echo "All done."
echo "V2 PDF : $FINAL_OUTPUT_V2  ($TOTAL_V2 plots)"
echo "V3 PDF : $FINAL_OUTPUT_V3  ($TOTAL_V3 plots)"
echo "============================="
