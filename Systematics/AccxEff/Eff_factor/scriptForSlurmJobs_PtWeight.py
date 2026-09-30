import os
import subprocess
import time

# ==============================================================================
# 1. USER CONFIGURATION (Everything you need to change is here)
# ==============================================================================
#--- pT-weighting ingredients (get_MassSpectra_PtWeight.C, macro 1) ---
STEP = "PtWeight"
DATE_tag = "Sept24_v2"  # include version here!!
SOURCE_CODE_NAME = "get_MassSpectra_PtWeight.C"

IS_DRY_RUN  = False  # True: Only generate .sub files. False: Ready to submit.
IS_TEST_RUN = False   # True: Only submit first 10 files (per sample). False: FULL RUN.

# Data and MC are independent chains (separate input lists, separate output files) --
# toggle these to (re)submit only the sample(s) you actually need. E.g. once the Data
# mass histograms are settled, set RUN_DATA = False and just resubmit MC after a code
# change, instead of re-submitting ~1000x more Data jobs for nothing.
RUN_DATA = False
RUN_MC   = True

USERNAME    = "saha115"
CMSSW_DIR   = f"/home/{USERNAME}/D0_ESE/CMSSW_13_2_11/src"
DATA_FLIST  = f"{CMSSW_DIR}/Systematics/AccxEff/inputFiles_PbPb2023_MB0to31_new.txt"
MC_FLIST    = f"{CMSSW_DIR}/Systematics/AccxEff/inputFiles_Prompt_MC_Sept5.txt"
OUTPUT_BASE = f"/scratch/negishi/{USERNAME}/D0_ESE_out/CMSSW_13_2_11/src/{STEP}_{DATE_tag}"

EXECUTABLE  = f"{STEP}_out_{DATE_tag}.exe"
SOURCE_CODE = f"{SOURCE_CODE_NAME}"

SLURM_ACCOUNT   = "physics"
SLURM_PARTITION = "cpu"
SLURM_TIME      = "04:00:00"
SLURM_QOS       = "standby"

# Data list is ~1000x larger than the MC list, so each sample gets its own
# files-per-job setting instead of sharing one N_FILES_PER_JOB.
N_FILES_PER_JOB_DATA = 200
N_FILES_PER_JOB_MC   = 10
MAX_JOBS_QUEUE  = 2000

# GRID PROXY
PROXY_PATH      = f"/home/{USERNAME}/myproxy"

# ================== USER CONFIGURATION ==============================================

def setup_environment():
    """Checks for a valid proxy and compiles the C++ code."""
    print("-" * 60)
    print(">>> STEP 1: Environment Setup")

    # Proxy Check
    try:
        timeleft = int(subprocess.check_output(f"voms-proxy-info -file {PROXY_PATH} -timeleft", shell=True).decode().strip())
    except:
        timeleft = 0

    if timeleft < 3600:
        print(">>> Proxy missing or expiring soon. Initializing (192 hours)...")
        os.system(f"voms-proxy-init --voms cms --valid 192:00 --out {PROXY_PATH}")
    else:
        print(f">>> Proxy is valid for {timeleft // 3600} hours.")

    # Compilation
    print(f">>> Compiling {SOURCE_CODE}...")
    compile_cmd = f"g++ {SOURCE_CODE} $(root-config --cflags --libs) -Wall -O2 -o {EXECUTABLE}"
    try:
        subprocess.run(compile_cmd, shell=True, check=True)
        print(">>> Compilation Successful.")
    except subprocess.CalledProcessError:
        print("!!! Compilation FAILED. Please check the errors above.")
        exit(1)
    print("-" * 60)

class JobsSubmission:
    def __init__(self):
        # Setup folder structure in scratch (shared by both the data and MC chains)
        self.dirs = {
            'sub':  f"{OUTPUT_BASE}/job_files",
            'log':  f"{OUTPUT_BASE}/log",
            'err':  f"{OUTPUT_BASE}/err",
            'root': f"{OUTPUT_BASE}/ROOT"
        }
        for path in self.dirs.values():
            os.makedirs(path, exist_ok=True)

    def write_sub_file(self, job_name, mode, input_flist, istart, iend):
        """Generates the .sub script with correct Negishi parameters.
        mode is 'data' or 'mc', matching get_MassSpectra_PtWeight.C's argv[1]."""
        sub_path = f"{self.dirs['sub']}/{job_name}.sub"
        with open(sub_path, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH --nodes=1\n")
            f.write(f"#SBATCH -o {self.dirs['log']}/{job_name}.out\n")
            f.write(f"#SBATCH -e {self.dirs['err']}/{job_name}.err\n")
            f.write(f"#SBATCH --job-name=\"{job_name}\"\n")
            f.write(f"#SBATCH --partition={SLURM_PARTITION}\n")
            f.write(f"#SBATCH -A {SLURM_ACCOUNT}\n")
            f.write(f"#SBATCH --time={SLURM_TIME}\n\n")
            f.write(f"#SBATCH --qos={SLURM_QOS}\n\n")

            # CMSSW / Environment setup
            f.write("export SCRAM_ARCH=el8_amd64_gcc12\n")
            f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
            f.write(f"cd {CMSSW_DIR}\n")
            f.write("eval `scramv1 runtime -sh` \n")
            f.write(f"cd {os.getcwd()}\n\n")

            # Execute
            f.write(f"export X509_USER_PROXY={PROXY_PATH}\n")
            f.write(f"./{EXECUTABLE} {mode} {input_flist} {OUTPUT_BASE} {istart} {iend}\n")
        return sub_path

    def submit_chain(self, mode, input_flist, n_files_per_job):
        """Submits one job array (all 'data' jobs, or all 'mc' jobs)."""
        with open(input_flist) as f:
            flist = [l.strip() for l in f if l.strip()]

        if IS_DRY_RUN or IS_TEST_RUN:
            label = "DRY RUN" if IS_DRY_RUN else "TEST SUBMISSION"
            print(f">>> MODE [{mode}]: {label} (10 files only)")
            flist = flist[:10]
        else:
            print(f">>> MODE [{mode}]: FULL PRODUCTION (Submitting all {len(flist)} files)")

        njobs = (len(flist) + n_files_per_job - 1) // n_files_per_job
        print(f">>> [{mode}] Configuration: {n_files_per_job} files/job | Total Jobs: {njobs}")

        for i in range(njobs):
            # Check queue every 100 submissions to prevent flooding
            if not IS_DRY_RUN and i > 0 and i % 100 == 0:
                while int(subprocess.check_output(f"squeue -u {USERNAME} | wc -l", shell=True)) > MAX_JOBS_QUEUE:
                    print(">>> Queue full, waiting 60s...")
                    time.sleep(60)

            istart = i * n_files_per_job
            iend = (i + 1) * n_files_per_job
            job_name = f"D0_ESE_{STEP}_{mode}_{istart}_{iend}"

            subfile = self.write_sub_file(job_name, mode, input_flist, istart, iend)

            if IS_DRY_RUN:
                print(f"Generated: {subfile}")
            else:
                os.system(f"sbatch {subfile}")
                if i % 10 == 0 or i == njobs - 1:
                    print(f">>> [{mode}] Submitted {i + 1}/{njobs} jobs...")

    def run(self):
        """Submits the data chain and/or the MC chain, per RUN_DATA / RUN_MC."""
        if RUN_DATA:
            self.submit_chain("data", DATA_FLIST, N_FILES_PER_JOB_DATA)
        else:
            print(">>> MODE [data]: SKIPPED (RUN_DATA = False)")

        if RUN_MC:
            self.submit_chain("mc", MC_FLIST, N_FILES_PER_JOB_MC)
        else:
            print(">>> MODE [mc]: SKIPPED (RUN_MC = False)")

if __name__ == "__main__":
    if not RUN_DATA and not RUN_MC:
        print(">>> RUN_DATA and RUN_MC are both False -- nothing to submit. Exiting.")
        exit(1)

    setup_environment()
    submitter = JobsSubmission()
    submitter.run()

    print("\n>>> Process Complete.")
    if IS_DRY_RUN:
        print(">>> Review scripts in scratch, then set IS_DRY_RUN = False.")
    elif IS_TEST_RUN:
        print(">>> Check 'squeue -u saha115'. If successful, set IS_TEST_RUN = False for full run.")
    else:
        print(f">>> Once all jobs finish, merge whichever sample(s) you (re)submitted, e.g.:")
        if RUN_DATA:
            print(f"    hadd DataMass_merged.root {OUTPUT_BASE}/ROOT/DataMass_*.root")
        if RUN_MC:
            print(f"    hadd MCSpectrum_merged.root {OUTPUT_BASE}/ROOT/MCSpectrum_*.root")
