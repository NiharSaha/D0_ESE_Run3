import os
import subprocess
import time
import sys

# ==============================================================================
#  Slurm submission for the BDT systematic study (prompt-D0 MC).
#
#  Two stages, run one after the other (set STAGE below):
#
#   STAGE = "resolution"
#       Calculate_Resolution_MC.C  -- SP 3-sub resolution ingredients,
#       q-inclusive & centrality-inclusive over the whole 0-50% MC sample.
#       After all jobs finish:
#          hadd -f <OUT>/ROOT/Resolution_MC_combined.root <OUT>/ROOT/Resolution_MC_*_*.root
#       then point RES_FILE (below, and/or RES_FILE_DEFAULT in Analysis_bin_BDT.h)
#       at that combined file.
#
#   STAGE = "flow"
#       flow_MC_BDT_sys.C  -- reads the combined resolution file, does the whole
#       vn calculation in one pass, fills TProfiles of vn(pT) with / without the
#       BDT cut (no intermediate ntuple).
#       After all jobs finish:
#          hadd -f <OUT>/ROOT/flow_MC_BDT_out_combined.root <OUT>/ROOT/flow_MC_BDT_sys_*_*.root
#          ./make_BDT_systematic.exe <OUT>/ROOT/flow_MC_BDT_out_combined.root BDT_sys_output.root
# ==============================================================================

# ------------------------------- USER CONFIG ---------------------------------
#STAGE = "resolution"          # "resolution"  then  "flow"
STAGE = "flow"          # "resolution"  then  "flow"
DATE_tag = "Sept10_v3"        # bump on every re-run
DATASET  = "prompt"

IS_DRY_RUN  = False           # True: only write .sub files
IS_TEST_RUN = False            # True: first 10 files only

USERNAME    = "saha115"
CMSSW_DIR   = f"/home/{USERNAME}/D0_ESE/CMSSW_13_2_11/src"
INPUT_FLIST = f"{CMSSW_DIR}/Systematics/BDT_sys/MC_based/inputFiles_Prompt_MC_Sept5_noDpt8.txt"

# combined resolution file used by STAGE=="flow" (passed as 5th arg to the exe)
RES_FILE = ("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/"
            f"BDT_sys_Resolution_MC_{DATASET}_{DATE_tag}/ROOT/Resolution_MC_combined.root")

SLURM_ACCOUNT   = "physics"
SLURM_PARTITION = "cpu"
SLURM_TIME      = "04:00:00"
SLURM_QOS       = "standby"
N_FILES_PER_JOB = 10
MAX_JOBS_QUEUE  = 1000
PROXY_PATH      = f"/home/{USERNAME}/myproxy"
# --------------------------- END USER CONFIG --------------------------------

# optional command-line override (used by run_BDT_sys_chain.sh):
#   python3 scriptForSlurmJobs_optimized.py <resolution|flow> <DATE_tag>
if len(sys.argv) == 3:
    STAGE, DATE_tag = sys.argv[1], sys.argv[2]
    RES_FILE = ("/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/"
                f"BDT_sys_Resolution_MC_{DATASET}_{DATE_tag}/ROOT/Resolution_MC_combined.root")

if STAGE == "resolution":
    STEP             = "BDT_sys_Resolution_MC"
    SOURCE_CODE_NAME = "Calculate_Resolution_MC.C"
elif STAGE == "flow":
    STEP             = "BDT_sys_Flow_MC"
    SOURCE_CODE_NAME = "flow_MC_BDT_sys.C"
else:
    raise SystemExit(f"unknown STAGE '{STAGE}'")

OUTPUT_BASE = f"/scratch/negishi/{USERNAME}/D0_ESE_out/CMSSW_13_2_11/src/{STEP}_{DATASET}_{DATE_tag}"
EXECUTABLE  = f"{STEP}_out_{DATE_tag}.exe"
SOURCE_CODE = f"{SOURCE_CODE_NAME}"


def setup_environment():
    print("-" * 60)
    print(f">>> STAGE = {STAGE}   ({SOURCE_CODE})")
    try:
        timeleft = int(subprocess.check_output(
            f"voms-proxy-info -file {PROXY_PATH} -timeleft", shell=True).decode().strip())
    except Exception:
        timeleft = 0
    if timeleft < 3600:
        print(">>> Proxy missing/expiring. Initializing (192 hours)...")
        os.system(f"voms-proxy-init --voms cms --valid 192:00 --out {PROXY_PATH}")
    else:
        print(f">>> Proxy valid for {timeleft // 3600} hours.")

    print(f">>> Compiling {SOURCE_CODE}...")
    try:
        subprocess.run(
            f"g++ {SOURCE_CODE} $(root-config --cflags --libs) -Wall -O2 -o {EXECUTABLE}",
            shell=True, check=True)
        print(">>> Compilation Successful.")
    except subprocess.CalledProcessError:
        print("!!! Compilation FAILED.")
        exit(1)

    if STAGE == "flow":
        print(">>> Compiling make_BDT_systematic.C ...")
        subprocess.run(
            "g++ make_BDT_systematic.C $(root-config --cflags --libs) -Wall -O2 -o make_BDT_systematic.exe",
            shell=True, check=False)
    print("-" * 60)


class JobsSubmission:
    def __init__(self):
        self.dirs = {k: f"{OUTPUT_BASE}/{k}" for k in ("job_files", "log", "err", "ROOT")}
        for path in self.dirs.values():
            os.makedirs(path, exist_ok=True)

    def write_sub_file(self, job_name, istart, iend):
        sub_path = f"{self.dirs['job_files']}/{job_name}.sub"
        extra = f" {RES_FILE}" if STAGE == "flow" else ""
        with open(sub_path, "w") as f:
            f.write("#!/bin/bash\n#SBATCH --nodes=1\n")
            f.write(f"#SBATCH -o {self.dirs['log']}/{job_name}.out\n")
            f.write(f"#SBATCH -e {self.dirs['err']}/{job_name}.err\n")
            f.write(f"#SBATCH --job-name=\"{job_name}\"\n")
            f.write(f"#SBATCH --partition={SLURM_PARTITION}\n")
            f.write(f"#SBATCH -A {SLURM_ACCOUNT}\n")
            f.write(f"#SBATCH --time={SLURM_TIME}\n\n")
            f.write(f"#SBATCH --qos={SLURM_QOS}\n\n")
            f.write("export SCRAM_ARCH=el8_amd64_gcc12\n")
            f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
            f.write(f"cd {CMSSW_DIR}\n")
            f.write("eval `scramv1 runtime -sh` \n")
            f.write(f"cd {os.getcwd()}\n\n")
            f.write(f"export X509_USER_PROXY={PROXY_PATH}\n")
            f.write(f"./{EXECUTABLE} {INPUT_FLIST} {OUTPUT_BASE} {istart} {iend}{extra}\n")
        return sub_path

    def run(self):
        with open(INPUT_FLIST) as f:
            flist = [l.strip() for l in f if l.strip()]
        if IS_DRY_RUN or IS_TEST_RUN:
            flist = flist[:10]
            print(">>> MODE:", "DRY RUN" if IS_DRY_RUN else "TEST (10 files)")
        else:
            print(f">>> MODE: FULL PRODUCTION ({len(flist)} files)")

        njobs = (len(flist) + N_FILES_PER_JOB - 1) // N_FILES_PER_JOB
        print(f">>> {N_FILES_PER_JOB} files/job | Total Jobs: {njobs}")

        for i in range(njobs):
            if not IS_DRY_RUN and i > 0 and i % 100 == 0:
                while int(subprocess.check_output(
                        f"squeue -u {USERNAME} | wc -l", shell=True)) > MAX_JOBS_QUEUE:
                    print(">>> Queue full, waiting 60s..."); time.sleep(60)
            istart = i * N_FILES_PER_JOB
            iend = (i + 1) * N_FILES_PER_JOB
            job_name = f"BDTsys_{STAGE}_{istart}_{iend}"
            subfile = self.write_sub_file(job_name, istart, iend)
            if IS_DRY_RUN:
                print(f"Generated: {subfile}")
            else:
                os.system(f"sbatch {subfile}")
                if i % 10 == 0 or i == njobs - 1:
                    print(f"Submitted {i+1}/{njobs} jobs...")


if __name__ == "__main__":
    setup_environment()
    JobsSubmission().run()
    print("\n>>> Submission complete.")
    if STAGE == "resolution":
        print(f">>> When done:\n"
              f"    hadd -f {OUTPUT_BASE}/ROOT/Resolution_MC_combined.root {OUTPUT_BASE}/ROOT/Resolution_MC_*_*.root\n"
              f"    # then set STAGE='flow' (RES_FILE already points at that path)")
    else:
        print(f">>> When done:\n"
              f"    hadd -f {OUTPUT_BASE}/ROOT/flow_MC_BDT_out_combined.root {OUTPUT_BASE}/ROOT/flow_MC_BDT_sys_*_*.root\n"
              f"    ./make_BDT_systematic.exe {OUTPUT_BASE}/ROOT/flow_MC_BDT_out_combined.root BDT_sys_output.root")
