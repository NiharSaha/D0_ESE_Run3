#!/bin/bash
#SBATCH --job-name=fit_v2_1pct_pt
#SBATCH --output=/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/fit_v2_outputs/logs/fit_v2_%A_%a.out
#SBATCH --error=/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/fit_v2_outputs/logs/fit_v2_%A_%a.err
#SBATCH --array=0-899%120
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=physics
#SBATCH --partition=cpu

set -euo pipefail

SRC="/home/saha115/D0_ESE/CMSSW_13_2_11/src/flow_extraction/fit_v2_prompt_dca_variable_v2_bin.C"
OUTDIR="${OUTDIR:-/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/fit_v2_outputs}"
BUILD_DIR="${OUTDIR}/build"
EXEC="${BUILD_DIR}/fit_v2_exec"
GEN_OUTNAME="prompt_vn_in_q2HF_total_graph_fullrange_withchi2_sigma_vs_v2v3.root"

module load root || true

mkdir -p "$BUILD_DIR" "$OUTDIR" "$OUTDIR/logs"
cd "$BUILD_DIR"

if [[ ! -x "$EXEC" || "$SRC" -nt "$EXEC" ]]; then
  echo "Compiling $SRC -> $EXEC"
  g++ -std=c++17 -O2 -pipe "$SRC" -o "$EXEC" $(root-config --cflags --libs) || { echo "Compile failed"; exit 1; }
fi

IDX=${SLURM_ARRAY_TASK_ID:-0}
N_PTBINS=10
CEN_IDX=$(( IDX / N_PTBINS ))
PT_IDX=$(( IDX % N_PTBINS ))

echo "Array idx=$IDX -> cent1pct=$CEN_IDX pt=$PT_IDX"
srun "$EXEC" "$CEN_IDX" "$PT_IDX"

if [[ -f "$BUILD_DIR/$GEN_OUTNAME" ]]; then
  mv -f "$BUILD_DIR/$GEN_OUTNAME" "${OUTDIR}/out_c${CEN_IDX}_pt${PT_IDX}.root"
elif [[ -f "$GEN_OUTNAME" ]]; then
  mv -f "$GEN_OUTNAME" "${OUTDIR}/out_c${CEN_IDX}_pt${PT_IDX}.root"
else
  echo "Warning: expected output $GEN_OUTNAME not found"
  exit 2
fi
