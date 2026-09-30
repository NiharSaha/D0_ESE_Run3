import os
import subprocess
import time

# ==============================================================================
#  Slurm submission for the BDT-cut data systematic -- "tight" variation.
#  Stage 3 of the data pipeline: flow_Analysis_ingradients_BDTtight.C
#  (SP numerator ntuple, BDT cut = bdt_cuts_tight.csv instead of the nominal
#  bdt_cuts.csv; everything else -- input file list, quantile header,
#  resolution file -- reused UNCHANGED from the nominal for_20Qbin chain,
#  since the BDT cut does not affect the event-plane / quantile stages).
#
#  After all jobs finish:
#     hadd -f ${OUTPUT_BASE}/ROOT/flow_Analysis_out_combined.root ${OUTPUT_BASE}/ROOT/flow_Analysis_*_*.root
#  then run sbatch_mass_fit_flow_BDTtight.sh (stage 4) -- its SCRATCH_BASE/inputNtuple
#  already point at that combined file for the SAME DATE_tag.
# ==============================================================================

STEP        = "NEW_BDT_sys_tight_D0_SP_12NUQbin_diffq2q3"
DATE_tag    = "Sept11_v0"          # MUST match DATE_tag in fit_mass_and_flow_BDTtight.C's
                                    # inputNtuple string and in sbatch_mass_fit_flow_BDTtight.sh
DATASET     = "MB0to31"
SOURCE_CODE_NAME = "flow_Analysis_ingradients_BDTtight.C"

IS_DRY_RUN  = False   # True: only write .sub files, do not submit
IS_TEST_RUN = False    # True: only first 10 files (this is a ~212000-file, ~1060-job
                       # production when False -- confirm before flipping this)

USERNAME    = "saha115"
CMSSW_DIR   = f"/home/{USERNAME}/D0_ESE/CMSSW_13_2_11/src"
# same nominal data file list as Quantiles/for_20Qbin/ -- BDT cut doesn't change which
# files are read, only which candidates in them survive
INPUT_FLIST = f"{CMSSW_DIR}/Systematics/BDT_sys/Data_based/tight/inputFiles_PbPb2023_MB0to31_new.txt"
OUTPUT_BASE = f"/scratch/negishi/{USERNAME}/D0_ESE_out/CMSSW_13_2_11/src/{STEP}_{DATASET}_{DATE_tag}"

EXECUTABLE  = f"{STEP}_out_{DATE_tag}.exe"
SOURCE_CODE = f"{SOURCE_CODE_NAME}"

SLURM_ACCOUNT   = "physics"
SLURM_PARTITION = "cpu"
SLURM_TIME      = "04:00:00"
SLURM_QOS       = "normal"
#SLURM_QOS       = "standby"
N_FILES_PER_JOB = 200          # matches the nominal for_20Qbin submission
MAX_JOBS_QUEUE  = 2000

PROXY_PATH      = f"/home/{USERNAME}/myproxy"


def setup_environment():
    print("-" * 60)
    print(">>> STEP 1: Environment Setup")
    try:
        timeleft = int(subprocess.check_output(
            f"voms-proxy-info -file {PROXY_PATH} -timeleft", shell=True).decode().strip())
    except Exception:
        timeleft = 0
    if timeleft < 3600:
        print(">>> Proxy missing or expiring soon. Initializing (192 hours)...")
        os.system(f"voms-proxy-init --voms cms --valid 192:00 --out {PROXY_PATH}")
    else:
        print(f">>> Proxy is valid for {timeleft // 3600} hours.")

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
        self.dirs = {
            'sub':  f"{OUTPUT_BASE}/job_files",
            'log':  f"{OUTPUT_BASE}/log",
            'err':  f"{OUTPUT_BASE}/err",
            'root': f"{OUTPUT_BASE}/ROOT",
        }
        for path in self.dirs.values():
            os.makedirs(path, exist_ok=True)

    def write_sub_file(self, job_name, istart, iend):
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
            f.write("export SCRAM_ARCH=el8_amd64_gcc12\n")
            f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
            f.write(f"cd {CMSSW_DIR}\n")
            f.write("eval `scramv1 runtime -sh` \n")
            f.write(f"cd {os.getcwd()}\n\n")
            f.write(f"export X509_USER_PROXY={PROXY_PATH}\n")
            f.write(f"./{EXECUTABLE} {INPUT_FLIST} {OUTPUT_BASE} {istart} {iend}\n")
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
            job_name = f"BDTtight_{istart}_{iend}"
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
    print(f">>> When done:\n"
          f"    hadd -f {OUTPUT_BASE}/ROOT/flow_Analysis_out_combined.root {OUTPUT_BASE}/ROOT/flow_Analysis_*_*.root\n"
          f"    then run sbatch_mass_fit_flow_BDTtight.sh (stage 4)")
