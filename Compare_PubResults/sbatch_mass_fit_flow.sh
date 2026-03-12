#!/bin/bash
#SBATCH --job-name=fit_mass_and_flow
#SBATCH --array=0-2%3
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=physics
#SBATCH --partition=cpu
#SBATCH --time=02:00:00                                                                                                                                 
#SBATCH --qos=standby
#SBATCH --output=/dev/null                                                                                                                              
#SBATCH --error=/dev/null


set -euo pipefail
set -x

# SLURM ids (use these to make per-task build dir to avoid races)
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
JOB_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-unknown}}"

SRC="/home/saha115/D0_ESE/CMSSW_13_2_11/src/Compare_PubResults/fit_mass_flow_forPubResults.C"
SRC_DIR="$(dirname "$SRC")"

SCRATCH_BASE="${SCRATCH_BASE:-/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Flow_output_PubResults_Mar12_v1}"
LOG_DIR="${SCRATCH_BASE}/logs"
BUILD_DIR="${SCRATCH_BASE}/build/${JOB_ID}_${TASK_ID}"
OUT_DIR="${SCRATCH_BASE}/output"
PLOT_DIR="${OUT_DIR}/prompt_mass_plot_withchi2_sigma"
EXEC="${BUILD_DIR}/fit_mass_exec"

# create required folders (per-task build dir prevents concurrent cp/overwrite races)
mkdir -p "$LOG_DIR" "$BUILD_DIR" "$OUT_DIR" "$PLOT_DIR"
cd "$BUILD_DIR"

# redirect stdout/stderr to files under SCRATCH_BASE/logs
STDOUT_FILE="${LOG_DIR}/fit_mass_${JOB_ID}_${TASK_ID}.out"
STDERR_FILE="${LOG_DIR}/fit_mass_${JOB_ID}_${TASK_ID}.err"
exec 1>>"$STDOUT_FILE" 2>>"$STDERR_FILE"

# debug info
echo "JOB_ID=${JOB_ID} TASK_ID=${TASK_ID}"
echo "PWD=$(pwd)"
ls -la "$(pwd)"

# prepare CMSSW runtime
cd /home/saha115/D0_ESE/CMSSW_13_2_11/src || exit 1
eval "$(scramv1 runtime -sh)"

# compile from source path but run inside per-task build dir
cd "$BUILD_DIR"
g++ -O2 -std=c++17 "$SRC" $(root-config --cflags --libs) -I"$SRC_DIR" -I"$CMSSW_BASE/src" -o "$EXEC" || {
  echo "Compilation failed" >&2
  exit 2
}

# ensure we run in build dir so program writes outputs there (then we'll move them)
cd "$BUILD_DIR"
echo "Running in $(pwd)"
"$EXEC" "${TASK_ID}" || {
  echo "Execution failed for centrality index ${TASK_ID}" >&2
  exit 3
}

# move produced ROOT output and plots to OUT_DIR (keep per-task build dir clean)
mv "$BUILD_DIR"/*.root "$OUT_DIR"/ 2>/dev/null || true

# PDFs are inside the subdirectory created by gSystem->mkdir()
if [ -d "$BUILD_DIR/prompt_mass_plot_withchi2_sigma" ]; then
  find "$BUILD_DIR/prompt_mass_plot_withchi2_sigma" -name "*.pdf" | while read f; do
    fname="$(basename "$f")"
    cp "$f" "${PLOT_DIR}/${TASK_ID}_${fname}"   # prefix with task id to avoid collisions
  done
  rm -rf "$BUILD_DIR/prompt_mass_plot_withchi2_sigma"
fi

echo "Finished fit_mass_and_flow for centrality index ${TASK_ID}"
