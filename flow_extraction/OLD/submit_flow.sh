#!/bin/bash



#SBATCH --job-name=D0_ESE_Flow_pT
#SBATCH --account=physics
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --partition=cpu
#SBATCH --array=0-79
#SBATCH --output=OUTPUT/temp_logs/job_%A_%a.out
#SBATCH --error=OUTPUT/temp_logs/job_%A_%a.err

# --- Directory Setup ---
DATE_STR=$(date +%Y%m%d)
VERSION="v1"
BASE_DIR="OUTPUT/${DATE_STR}_${VERSION}"
LOG_DIR="${BASE_DIR}/logs"
ROOT_DIR="${BASE_DIR}/root_files"

mkdir -p $LOG_DIR $ROOT_DIR



INPUT_FILE="/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/flow_Jan19_MB15to20/ROOT/flow_Analysis_out_combined.root"

# --- Logic to Map Array ID ---
ID=${SLURM_ARRAY_TASK_ID}
if [ $ID -le 39 ]; then
    HARM=2
    TEMP_ID=$ID
else
    HARM=3
    TEMP_ID=$((ID - 40))
fi

CENT_IDX=$((TEMP_ID / 10))
PT_BIN=$(( (TEMP_ID % 10) + 1 ))

OUT_FILE="${ROOT_DIR}/Final_v${HARM}_c${CENT_IDX}_p${PT_BIN}.root"

# --- Execute ROOT ---
# We use 'tee' to capture the log into our LOG_DIR manually 
# since SBATCH header can't use shell variables easily.
root -l -b -q "Extract_Flow_Final_Directories.C(\"$INPUT_FILE\", $HARM, $CENT_IDX, $PT_BIN, \"$OUT_FILE\")" 2>&1 | tee ${LOG_DIR}/job_${ID}_v${HARM}_c${CENT_IDX}_p${PT_BIN}.log

mv "OUTPUT/temp_logs/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" "$LOG_DIR/"
mv "OUTPUT/temp_logs/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" "$LOG_DIR/"
