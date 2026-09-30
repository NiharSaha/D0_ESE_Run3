#!/bin/bash
#SBATCH --job-name=ESE_quantiles
#SBATCH --array=0-49
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=cpu
#SBATCH --account=physics
#SBATCH --output=logs/job_%a.out
#SBATCH --error=logs/job_%a.err

# Create necessary directories
mkdir -p logs 

INPUT_DIR="/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/NEW_D0_Quantiles_MB0to31_FullStat_Aug23/ROOT"

OUT_NAME="HIST_NU_Aug23"

# when USE_SORT_METHOD = true
#mkdir -p temp_cuts_${OUT_NAME} temp_roots_${OUT_NAME}
#root -l -b -q "ExtractQuantiles.C(${SLURM_ARRAY_TASK_ID}, \"${INPUT_DIR}\", \"${OUT_NAME}\", true, true)"

# when USE_SORT_METHOD = false
mkdir -p temp_cuts_${OUT_NAME} temp_roots_${OUT_NAME}
root -l -b -q "ExtractQuantiles.C(${SLURM_ARRAY_TASK_ID}, \"${INPUT_DIR}\", \"${OUT_NAME}\", false, true)"
