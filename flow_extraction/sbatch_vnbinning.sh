#!/bin/bash
#SBATCH --job-name=save_vnbinning
#SBATCH --array=0-71%72
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=physics
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

SRC="/home/saha115/D0_ESE/CMSSW_13_2_11/src/flow_extraction/save_vnbinning_to_txt.C"
SCRATCH_BASE="${SCRATCH_BASE:-/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/save_vnbinning_outputs_Mar9_v0}"
LOG_DIR="${SCRATCH_BASE}/logs"
BUILD_DIR="${SCRATCH_BASE}/build"
#BUILD_DIR="${SCRATCH_BASE}/build_job${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-unknown}}_task${SLURM_ARRAY_TASK_ID:-0}"
OUT_DIR="${SCRATCH_BASE}/output"
EXEC="${BUILD_DIR}/save_vnbinning_to_txt"

mkdir -p "$LOG_DIR" "$BUILD_DIR" "$OUT_DIR"
cd "$BUILD_DIR"

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
JOB_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-unknown}}"
STDOUT_FILE="${LOG_DIR}/save_vnbinning_${JOB_ID}_${TASK_ID}.out"
STDERR_FILE="${LOG_DIR}/save_vnbinning_${JOB_ID}_${TASK_ID}.err"
exec 1>>"$STDOUT_FILE" 2>>"$STDERR_FILE"

# map array index (0..71) -> centrality (0..6) and pT (0..11)
CENT_IDX=$(( TASK_ID / 12 ))
PT_IDX=$(( TASK_ID % 12 ))

echo "JOB_ID=${JOB_ID} TASK_ID=${TASK_ID} CENT_IDX=${CENT_IDX} PT_IDX=${PT_IDX}"
echo "BUILD_DIR=${BUILD_DIR} OUT_DIR=${OUT_DIR} LOG_DIR=${LOG_DIR}"

# copy and compile
#cp "$SRC" .
#g++ -O2 -std=c++17 `root-config --cflags --libs` save_vnbinning_to_txt.C -o save_vnbinning_to_txt

# run program for the mapped (centrality, pT) pair
#./save_vnbinning_to_txt "${CENT_IDX}" "${PT_IDX}"

# compile once with a lock (replace the old cp+g++ block)
LOCKDIR="${BUILD_DIR}/compile.lock"
if [[ ! -x "$EXEC" || "$SRC" -nt "$EXEC" ]]; then
  if mkdir "$LOCKDIR" 2>/dev/null; then
    echo "Compiling $SRC -> $EXEC"
    g++ -O2 -std=c++17 $(root-config --cflags --libs) "$SRC" -o "$EXEC" || { echo "Compile failed"; rmdir "$LOCKDIR"; exit 1; }
    rmdir "$LOCKDIR"
  else
    echo "Waiting for compilation by another task..."
    while [[ -d "$LOCKDIR" ]]; do sleep 2; done
  fi
fi

# run program for the mapped (centrality, pT) pair
"$EXEC" "${CENT_IDX}" "${PT_IDX}"

#move any produced vnbinning_*.txt files to OUT_DIR (add job/task/timestamp to avoid collisions)
shopt -s nullglob
files=(vnbinning_*.txt)
if [ ${#files[@]} -gt 0 ]; then
  for f in "${files[@]}"; do
    # C++ already embeds job/task in the filename; just move as-is
    DEST="${OUT_DIR}/$(basename "${f}")"
    mv -- "${f}" "${DEST}"
    echo "Saved output to ${DEST}"
  done
else
  echo "Expected vnbinning_*.txt not found" >&2
  exit 1
fi
