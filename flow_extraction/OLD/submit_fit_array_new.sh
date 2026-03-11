#!/bin/bash
#SBATCH --job-name=fit_v2_1pct
#SBATCH --output=logs/fit_v2_%A_%a.out
#SBATCH --error=logs/fit_v2_%A_%a.err
#SBATCH --array=0-89%60
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --account=physics
#SBATCH --partition=cpu

set -euo pipefail

SRC="/home/saha115/D0_ESE/CMSSW_13_2_11/src/flow_extraction/fit_v2_prompt_dca_variable_v2_bin.C"
OUTDIR="${OUTDIR:-/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/fit_v2_outputs}"
BUILD_DIR="${OUTDIR}/build"
EXEC="${BUILD_DIR}/fit_v2_exec"
GEN_OUTNAME="prompt_vn_in_q2HF_total_graph_fullrange_withchi2_sigma_vs_v2v3.root"

module load root || true

mkdir -p "$BUILD_DIR" "$OUTDIR"

# compile once using an atomic lock to avoid concurrent compiles
(
  exec 9>"${BUILD_DIR}/compile.lock"
  flock -x 9
  if [[ ! -x "$EXEC" || "$SRC" -nt "$EXEC" ]]; then
    echo "Compiling $SRC -> $EXEC"
    g++ -std=c++17 -O2 "$SRC" -o "$EXEC" $(root-config --cflags --libs) || { echo "Compile failed"; exit 1; }
  fi
) 

IDX=${SLURM_ARRAY_TASK_ID:-0}
echo "Running index $IDX"
srun "$EXEC" "$IDX"

if [[ -f "$GEN_OUTNAME" ]]; then
  mv -f "$GEN_OUTNAME" "${OUTDIR}/out_1pct_${IDX}.root"
  echo "Saved ${OUTDIR}/out_1pct_${IDX}.root"
else
  echo "Warning: expected output $GEN_OUTNAME not found"
  exit 2
fi
