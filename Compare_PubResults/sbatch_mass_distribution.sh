#!/bin/bash
#SBATCH --job-name=save_mass_1pct_pt
#SBATCH --array=0-35%36
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=physics
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --qos=standby 
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

SRC="/home/saha115/D0_ESE/CMSSW_13_2_11/src/Compare_PubResults/save_mass_distributions_forPubResults.C"
SCRATCH_BASE="/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/save_mass_outputs_PubResults_Mar6_v0"
LOG_DIR="${SCRATCH_BASE}/logs"
BUILD_DIR="${SCRATCH_BASE}/build"
OUT_DIR="${SCRATCH_BASE}/output"
EXEC="${BUILD_DIR}/save_mass_exec"

mkdir -p "$LOG_DIR" "$BUILD_DIR" "$OUT_DIR"
cd "$BUILD_DIR"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}

# compile with lock to avoid races across array tasks
LOCKDIR="${BUILD_DIR}/compile.lock"
if [[ ! -x "$EXEC" || "$SRC" -nt "$EXEC" ]]; then
    if mkdir "$LOCKDIR" 2>/dev/null; then
        echo "Compiling $SRC -> $EXEC"
        g++ -std=c++17 -O2 -pipe "$SRC" -o "$EXEC" $(root-config --cflags --libs) \
            || { echo "Compile failed" >&2; rmdir "$LOCKDIR"; exit 1; }
        rmdir "$LOCKDIR"
    else
        echo "Waiting for compilation by another task..."
        while [[ -d "$LOCKDIR" ]]; do sleep 2; done
    fi
fi
# check exec exists after potential wait
if [[ ! -x "$EXEC" ]]; then
    echo "ERROR: exec missing after compilation step" >&2; exit 1
fi

# map array index → coarse centrality bin and pT bin
# 3 coarse cent bins × 12 pT bins = 36 total tasks
N_PTBINS=12
IDX=${SLURM_ARRAY_TASK_ID}
CEN_IDX=$(( IDX / N_PTBINS ))   # 0=0-10%, 1=10-30%, 2=30-50%
PT_IDX=$(( IDX % N_PTBINS ))    # 0-11

CEN_LABELS=("0to10" "10to30" "30to50")
CEN_LABEL="${CEN_LABELS[$CEN_IDX]}"

echo "Array idx=$IDX -> cen_coarse=$CEN_IDX ($(( CEN_IDX==0?0:CEN_IDX==1?10:30 ))-$(( CEN_IDX==0?10:CEN_IDX==1?30:50 ))%) pt_idx=$PT_IDX"

"$EXEC" "$CEN_IDX" "$PT_IDX"

# match output filename produced by macro:
# mass_distributions_pubcomp_cent<label>_pt<N>.root

EXPECTED="${BUILD_DIR}/mass_distributions_pubcomp_cent${CEN_LABEL}_pt${PT_IDX}.root"

if [[ -f "$EXPECTED" ]]; then
    mv -f "$EXPECTED" "${OUT_DIR}/out_c${CEN_IDX}_pt${PT_IDX}.root"
    echo "Saved: ${OUT_DIR}/out_c${CEN_IDX}_pt${PT_IDX}.root"
else
    echo "ERROR: expected output not found: $EXPECTED" >&2
    exit 1
fi
