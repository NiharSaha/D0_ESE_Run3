import os
import subprocess
import time

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

STEP = "NEW_Systematic_PromptFrac_MC_DCA_template"
DATE_tag = "Sept10_v0" 
DATASET = "MC2023"
SOURCE_CODE_NAME = "getDCA_fromMC.C"

IS_DRY_RUN = False
IS_TEST_RUN = False

USERNAME = "saha115"
CMSSW_DIR = f"/home/{USERNAME}/D0_ESE/CMSSW_13_2_11/src"

SAMPLES = {
    "Prompt":    (f"{CMSSW_DIR}/Systematics/Prompt_frac/inputFiles_Prompt_MC_Sept5.txt", 5),
    "NonPrompt": (f"{CMSSW_DIR}/Systematics/Prompt_frac/inputFiles_NonPrompt_MC_Sept5.txt", 5),
}



OUTPUT_BASE = f"/scratch/negishi/{USERNAME}/D0_ESE_out/CMSSW_13_2_11/src/{STEP}_{DATASET}_{DATE_tag}"

EXECUTABLE = f"{STEP}_out_{DATE_tag}.exe"
SOURCE_CODE = f"{SOURCE_CODE_NAME}"

SLURM_ACCOUNT   = "physics"
SLURM_PARTITION = "cpu"
SLURM_TIME      = "04:00:00"
SLURM_QOS       = "standby"
MAX_JOBS_QUEUE  = 2000

PROXY_PATH = f"/home/{USERNAME}/myproxy"

# ==============================================================================


def setup_environment():
    """Checks for a valid proxy and compiles the C++ code."""
    print("-" * 60)
    print(">>> STEP 1: Setup & Compilation")

    try:
        timeleft = int(subprocess.check_output(f"voms-proxy-info -file {PROXY_PATH} -timeleft", shell=True).decode().strip())
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
        print("!!! Compilation FAILED.")
        exit(1)
    print("-" * 60)


class JobsSubmission:
    def __init__(self):
        self.dirs = {
            'sub':            f"{OUTPUT_BASE}/job_files",
            'log':            f"{OUTPUT_BASE}/log",
            'err':            f"{OUTPUT_BASE}/err",
            'root_prompt':    f"{OUTPUT_BASE}/ROOT/Prompt",
            'root_nonprompt': f"{OUTPUT_BASE}/ROOT/NonPrompt"
        }
        for path in self.dirs.values():
            os.makedirs(path, exist_ok=True)

    def write_sub_file(self, job_name, input_flist, out_root_file, istart, iend, sample_type):
        """Generates the .sub script."""
        sub_path = f"{self.dirs['sub']}/{job_name}.sub"
        with open(sub_path, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH --nodes=1\n")
            f.write(f"#SBATCH -o {self.dirs['log']}/{job_name}.out\n")
            f.write(f"#SBATCH -e {self.dirs['err']}/{job_name}.err\n")
            f.write(f"#SBATCH --job-name=\"{job_name}\"\n")
            f.write(f"#SBATCH --partition={SLURM_PARTITION}\n")
            f.write(f"#SBATCH -A {SLURM_ACCOUNT}\n")
            f.write(f"#SBATCH --time={SLURM_TIME}\n")
            f.write(f"#SBATCH --qos={SLURM_QOS}\n\n")

            f.write("export SCRAM_ARCH=el8_amd64_gcc12\n")
            f.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
            f.write(f"cd {CMSSW_DIR}\n")
            f.write("eval `scramv1 runtime -sh`\n")
            f.write(f"cd {os.getcwd()}\n\n")

            f.write(f"export X509_USER_PROXY={PROXY_PATH}\n")
            f.write(f"./{EXECUTABLE} {input_flist} {out_root_file} {istart} {iend} {sample_type}\n")
        return sub_path

    def run(self):
        """Submits independent job sets for Prompt and NonPrompt."""
        total_submitted = 0

        for sample_type, (input_flist, n_files_per_job) in SAMPLES.items():
            with open(input_flist) as f:
                flist = [line.strip() for line in f if line.strip()]

            if IS_DRY_RUN:
                print(f">>> [{sample_type}] MODE: DRY RUN (10 files)")
                flist = flist[:10]
            elif IS_TEST_RUN:
                print(f">>> [{sample_type}] MODE: TEST SUBMISSION (10 files)")
                flist = flist[:10]
            else:
                print(f">>> [{sample_type}] MODE: FULL PRODUCTION ({len(flist)} files)")

            njobs = (len(flist) + n_files_per_job - 1) // n_files_per_job
            print(f">>> [{sample_type}] {n_files_per_job} files/job | Total Jobs: {njobs}")

            out_dir = self.dirs['root_prompt'] if sample_type == "Prompt" else self.dirs['root_nonprompt']

            for i in range(njobs):
                if not IS_DRY_RUN and total_submitted > 0 and total_submitted % 100 == 0:
                    while int(subprocess.check_output(f"squeue -u {USERNAME} | wc -l", shell=True)) > MAX_JOBS_QUEUE:
                        print(">>> Queue full, waiting 60s...")
                        time.sleep(60)

                istart = i * n_files_per_job
                iend = (i + 1) * n_files_per_job
                job_name = f"D0_DCA_{sample_type}_{istart}_{iend}"
                out_root_file = f"{out_dir}/Hist_DCA_{sample_type}_{istart}_{iend}.root"

                subfile = self.write_sub_file(job_name, input_flist, out_root_file, istart, iend, sample_type)

                if IS_DRY_RUN:
                    print(f"Generated: {subfile}")
                else:
                    os.system(f"sbatch {subfile}")
                    total_submitted += 1
                    if i % 10 == 0 or i == njobs - 1:
                        print(f"[{sample_type}] Submitted {i+1}/{njobs} jobs...")


if __name__ == "__main__":
    setup_environment()
    submitter = JobsSubmission()
    submitter.run()

    print("\n>>> Process Complete.")
