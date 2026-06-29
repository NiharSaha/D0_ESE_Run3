#!/bin/bash
#SBATCH --job-name=ESE_Sort
#SBATCH --array=0-10
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --partition=cpu
#SBATCH --account=physics
#SBATCH --output=logs/job_%a.out
#SBATCH --error=logs/job_%a.err

# Create necessary directories
mkdir -p logs temp_cuts_MB0to1_charge temp_roots_MB0to1_charge

INPUT_DIR="/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Quantiles_MB0to1_Jun11_ForCheck/ROOT"

# NEW: Define your specific coarser centrality bins
BINS=(0 5 10 15 20 25 30 35 40 45 50)

# NEW: Map the Slurm array ID (0-10) to the corresponding bin number
CURRENT_BIN=${BINS[$SLURM_ARRAY_TASK_ID]}

# Execute the macro with centBin and input directory
root -l -b -q "ExtractQuantiles_sortFunc.C(${CURRENT_BIN}, ${CURRENT_BIN}+5, \"${INPUT_DIR}\")"
