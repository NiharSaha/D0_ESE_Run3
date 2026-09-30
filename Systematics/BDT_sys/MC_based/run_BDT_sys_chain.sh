#!/bin/bash
# =============================================================================
#  run_BDT_sys_chain.sh
#
#  Full MC BDT-systematic chain for one DATE_tag, unattended:
#    resolution jobs -> wait -> hadd -> flow jobs -> wait -> hadd
#    -> make_BDT_systematic -> check_BDT_sys (vs a reference flow file)
#  Jobs whose output is missing after the queue drains (e.g. standby
#  pre-emption) are resubmitted once from their existing .sub file.
#
#  Usage (from MC_based/):
#     ./run_BDT_sys_chain.sh <DATE_tag> [<reference flow_MC_BDT_out_combined.root>]
#  Summary log: <flow output dir>/check_<DATE_tag>.log
# =============================================================================
set -u
TAG=${1:?usage: $0 DATE_tag [ref_flow_combined.root]}
REF=${2:-}
USER_NAME=saha115
OUT=/scratch/negishi/${USER_NAME}/D0_ESE_out/CMSSW_13_2_11/src
RES_DIR=${OUT}/BDT_sys_Resolution_MC_prompt_${TAG}
FLOW_DIR=${OUT}/BDT_sys_Flow_MC_prompt_${TAG}
HERE=$(cd "$(dirname "$0")" && pwd)

cd /home/${USER_NAME}/D0_ESE/CMSSW_13_2_11/src
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval $(scramv1 runtime -sh)
cd "${HERE}"

stamp() { echo "[$(date '+%F %T')] $*"; }

wait_queue() { # $1 = job-name prefix
  while [ "$(squeue -u ${USER_NAME} -h -o %j | grep -c "^$1")" -gt 0 ]; do
    sleep 120
  done
}

# $1 stage, $2 output dir, $3 per-job ROOT file prefix
run_stage() {
  local stage=$1 dir=$2 pfx=$3
  stamp "submitting ${stage} (${TAG})"
  python3 scriptForSlurmJobs_optimized.py "${stage}" "${TAG}" || { stamp "submission failed"; exit 1; }
  wait_queue "BDTsys_${stage}_"
  for attempt in 1 2; do
    local missing=()
    for sub in "${dir}"/job_files/BDTsys_${stage}_*.sub; do
      local rng=${sub##*BDTsys_${stage}_}; rng=${rng%.sub}
      [ -s "${dir}/ROOT/${pfx}_${rng}.root" ] || missing+=("${sub}")
    done
    [ ${#missing[@]} -eq 0 ] && break
    if [ ${attempt} -eq 1 ]; then
      stamp "${#missing[@]} ${stage} job(s) missing output, resubmitting once"
      for s in "${missing[@]}"; do sbatch "${s}"; done
      sleep 60; wait_queue "BDTsys_${stage}_"
    else
      stamp "WARNING: ${#missing[@]} ${stage} job(s) still missing output:"
      printf '    %s\n' "${missing[@]}"
    fi
  done
}

run_stage resolution "${RES_DIR}" Resolution_MC
hadd -f "${RES_DIR}/ROOT/Resolution_MC_combined.root" "${RES_DIR}"/ROOT/Resolution_MC_*_*.root > /dev/null \
  || { stamp "resolution hadd failed"; exit 1; }
stamp "resolution combined"

run_stage flow "${FLOW_DIR}" flow_MC_BDT_sys
FLOW=${FLOW_DIR}/ROOT/flow_MC_BDT_out_combined.root
hadd -f "${FLOW}" "${FLOW_DIR}"/ROOT/flow_MC_BDT_sys_*_*.root > /dev/null \
  || { stamp "flow hadd failed"; exit 1; }
stamp "flow combined"

./make_BDT_systematic.exe "${FLOW}" "${FLOW_DIR}/BDT_sys_output_${TAG}.root" > "${FLOW_DIR}/make_BDT_systematic_${TAG}.log"
g++ check_BDT_sys.C $(root-config --cflags --libs) -O2 -o check_BDT_sys.exe
./check_BDT_sys.exe "${FLOW}" ${REF} | tee "${FLOW_DIR}/check_${TAG}.log"
stamp "done: ${FLOW_DIR}/check_${TAG}.log"
