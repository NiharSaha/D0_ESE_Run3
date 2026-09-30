#!/bin/bash
#SBATCH --job-name=save_mass_1pct_pt
#SBATCH --array=0-959%960
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=physics
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --qos=standby 
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

SRC="/home/saha115/D0_ESE/CMSSW_13_2_11/src/flow_extraction/save_mass_distributions.C"
SCRATCH_BASE="${SCRATCH_BASE:-/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/save_mass_outputs_Mar16_v0}"
LOG_DIR="${SCRATCH_BASE}/logs"
BUILD_DIR="${SCRATCH_BASE}/build"
OUT_DIR="${SCRATCH_BASE}/output"
EXEC="${BUILD_DIR}/save_mass_exec"

# create required folders
mkdir -p "$LOG_DIR" "$BUILD_DIR" "$OUT_DIR"
cd "$BUILD_DIR"

# redirect stdout/stderr to files under SCRATCH_BASE/logs (uses SLURM env vars)
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
JOB_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-unknown}}"
STDOUT_FILE="${LOG_DIR}/save_mass_${JOB_ID}_${TASK_ID}.out"
STDERR_FILE="${LOG_DIR}/save_mass_${JOB_ID}_${TASK_ID}.err"
# ensure log dir exists then redirect
exec 1>>"$STDOUT_FILE" 2>>"$STDERR_FILE"

#module load root
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}

# compile with lock to avoid races across array tasks
LOCKDIR="${BUILD_DIR}/compile.lock"
if [[ ! -x "$EXEC" || "$SRC" -nt "$EXEC" ]]; then
  if mkdir "$LOCKDIR" 2>/dev/null; then
    echo "Compiling $SRC -> $EXEC"
    g++ -std=c++17 -O2 -pipe "$SRC" -o "$EXEC" $(root-config --cflags --libs) || { echo "Compile failed"; rmdir "$LOCKDIR"; exit 1; }
    rmdir "$LOCKDIR"
  else
    echo "Waiting for compilation by another task..."
    while [[ -d "$LOCKDIR" ]]; do sleep 2; done
  fi
fi

# map array index to 1% centrality and pT bin
IDX_BASE=${JOB_OFFSET:-0}
IDX=$(( ${SLURM_ARRAY_TASK_ID:-0} + IDX_BASE ))
N_PTBINS=12
CEN_IDX=$(( IDX / N_PTBINS ))
PT_IDX=$(( IDX % N_PTBINS ))

echo "Array idx=$IDX -> cent1pct=$CEN_IDX pt=$PT_IDX"

srun --exclusive -n1 -c ${SLURM_CPUS_PER_TASK:-2} "$EXEC" "$CEN_IDX" "$PT_IDX"

# expected output (per-run) is mass_distributions[_1pct_X][_pt_Y].root
OUTNAME_BASE="mass_distributions"
if [[ -f "${BUILD_DIR}/${OUTNAME_BASE}_1pct_${CEN_IDX}_pt_${PT_IDX}.root" ]]; then
  mv -f "${BUILD_DIR}/${OUTNAME_BASE}_1pct_${CEN_IDX}_pt_${PT_IDX}.root" "${OUT_DIR}/out_c${CEN_IDX}_pt${PT_IDX}.root"
elif [[ -f "${BUILD_DIR}/${OUTNAME_BASE}.root" ]]; then
  mv -f "${BUILD_DIR}/${OUTNAME_BASE}.root" "${OUT_DIR}/out_c${CEN_IDX}_pt${PT_IDX}.root"
fi


echo "Saved output to ${OUT_DIR}/out_c${CEN_IDX}_pt${PT_IDX}.root and logs to ${LOG_DIR}"
