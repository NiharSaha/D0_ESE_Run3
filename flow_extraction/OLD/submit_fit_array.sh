#!/bin/bash
#SBATCH --job-name=D0_YieldAvg_Fit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=0-39            # 4 Cent Groups * 10 pT Bins = 40 tasks
#SBATCH --account=physics                                                                                                                                                
#SBATCH --partition=cpu            
#SBATCH --output=/dev/null      # Standard output handled inside script for custom naming

#SBATCH --job-name=D0_ESE_Flow_pT                                                                                                                                         

# --- 1. Directory Setup ---
# Define the base OUTPUT directory
BASE_OUTPUT="/scratch/negishi/saha115/D0_ESE_out/OUTPUT"

# Create a unique subdirectory for this specific run if it doesn't exist
# We use a date-based name: e.g., Run_20260122_2145
RUN_DIR="${BASE_OUTPUT}/Run_$(date +%Y%m%d_%H%M)"
LOG_DIR="${RUN_DIR}/logs"
RESULT_DIR="${RUN_DIR}/FitResults"

# Only the first task in the array needs to create the directories
if [ $SLURM_ARRAY_TASK_ID -eq 0 ]; then
    mkdir -p $LOG_DIR
    mkdir -p $RESULT_DIR
    echo "Created output directory: $RUN_DIR"
fi

# Wait a moment to ensure directory creation is propagated across nodes
sleep 5

# Redirect individual task output to the new log directory
exec > "${LOG_DIR}/fit_cent_pt_${SLURM_ARRAY_TASK_ID}.log" 2>&1

# --- 2. Map the Array ID to Indices ---
CENT_IDX=$((SLURM_ARRAY_TASK_ID / 10))
PT_IDX=$((SLURM_ARRAY_TASK_ID % 10))

# --- 3. Define Input ---
INPUT_FILE="/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/SP_MB11to21_Jan22/ROOT/flow_Analysis_out_combined.root"

echo "Processing Job Array ID: $SLURM_ARRAY_TASK_ID"
echo "Centrality Index: $CENT_IDX | pT Bin Index: $PT_IDX"
echo "Output will be saved in: $RESULT_DIR"

# --- 4. Execute the Fits ---
# Running both harmonics sequentially in the same job
echo ">>> Starting v2 Extraction..."
root -l -b -q "SP_new.C(\"$INPUT_FILE\", \"$RESULT_DIR\", 2, $CENT_IDX, $PT_IDX)"

echo ">>> Starting v3 Extraction..."
root -l -b -q "SP_new.C(\"$INPUT_FILE\", \"$RESULT_DIR\", 3, $CENT_IDX, $PT_IDX)"

echo ">>> Job $SLURM_ARRAY_TASK_ID completed successfully."
