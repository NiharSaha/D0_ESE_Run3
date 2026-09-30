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

#INPUT_DIR="/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/NEW_Systematic_centralityUP_D0_Quantiles_MB0to31_Aug31_FullStat/ROOT"
INPUT_DIR="/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/NEW_Systematic_centralityDOWN_D0_Quantiles_MB0to31_Aug31_FullStat/ROOT"



# when USE_SORT_METHOD = true
#mkdir -p temp_cuts_SORT temp_roots_SORT
#root -l -b -q "ExtractQuantiles.C(${SLURM_ARRAY_TASK_ID}, \"${INPUT_DIR}\", \"SORT\", true, true)"

# when USE_SORT_METHOD = false
mkdir -p temp_cuts_HIST_NU_sysDOWN temp_roots_HIST_NU_sysDOWN
root -l -b -q "ExtractQuantiles.C(${SLURM_ARRAY_TASK_ID}, \"${INPUT_DIR}\", \"HIST_NU_sysDOWN\", false, true)"
