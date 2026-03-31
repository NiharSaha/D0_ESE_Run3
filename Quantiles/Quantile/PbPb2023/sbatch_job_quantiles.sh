#!/bin/bash
#SBATCH --job-name=ESE_Sort
#SBATCH --array=0-89
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --partition=cpu
#SBATCH --account=physics
#SBATCH --output=logs/job_%a.out
#SBATCH --error=logs/job_%a.err

# Create necessary directories
mkdir -p logs temp_cuts_MB0to31 temp_roots_MB0to31

INPUT_DIR="/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/Quantiles_MB0to31_Mar30/ROOT"

# Execute the macro with centBin and input directory
root -l -b -q "ExtractQuantiles_sortFunc.C(${SLURM_ARRAY_TASK_ID}, \"${INPUT_DIR}\")"
