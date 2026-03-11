import os
import subprocess

# Step 1: Compile the executable
try:
    compile_command = [
        "g++",
        "../Ana_new.C",
        "$(root-config --cflags --libs)",
        "-Wall",
        "-O2",
        "-o",
        "D0_ESE_out.exe"
    ]
    subprocess.run(" ".join(compile_command), shell=True, check=True)
    print("Compilation successful!")
except subprocess.CalledProcessError:
    print("Compilation failed. Please check your source code.")
    exit(1)

# Step 2: Run the executable or proceed with other steps
print("Running subsequent operations...")

class jobs_submission():
    def __init__(self, flist_name):
        self.flist_name = flist_name
        self.max_njobs_per_round = 1500
        self.username = 'saha115'

    # Build job submission files 
    def build_job_files(self, job_name):
        dir_job_file = self.base_dir + '/job_files/'
        dir_job_out = self.base_dir + '/log/'
        dir_job_err = self.base_dir + '/err/'
        dir_job_output = self.base_dir + '/ROOT/'

        if not os.path.isdir(dir_job_file):
            os.makedirs(dir_job_file)
        if not os.path.isdir(dir_job_out):
            os.makedirs(dir_job_out)
        if not os.path.isdir(dir_job_err):
            os.makedirs(dir_job_err)
        if not os.path.isdir(dir_job_output):
            os.makedirs(dir_job_output)

        subfile = dir_job_file + job_name + '.sub'

        # SLURM job submission file
        with open(subfile, "w") as fsub:
            fsub.write('#!/bin/bash\n')
            fsub.write(f'#SBATCH --nodes=1\n')
            fsub.write(f'#SBATCH -o {dir_job_out}{job_name}.out\n')
            fsub.write(f'#SBATCH -e {dir_job_err}{job_name}.err\n')
            fsub.write(f'#SBATCH --job-name="{job_name}"\n')
            fsub.write('#SBATCH --no-requeue\n')
            fsub.write('#SBATCH --partition=bell-a\n')
            fsub.write('export SCRAM_ARCH=el8_amd64_gcc12\n')
            fsub.write('#SBATCH --export=ALL\n')
            fsub.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
            fsub.write('export X509_USER_PROXY=/home/saha115/myproxy\n')
            fsub.write(f'./D0_ESE_out.exe {self.flist_name} {self.base_dir} {self.istart} {self.iend} \n')
            fsub.close()

        return subfile

    # Submit jobs
    def submit_jobs(self, base_dir, njobs, nfile_per_job, dry_run=False):
        self.base_dir = base_dir

        dir_job_output = os.path.join(base_dir, 'ROOT/')


        count = 0
        for i in range(njobs):
            if count >= self.max_njobs_per_round: # hold until last round jobs are done
                self.check_running_jobs()
                count=0  # restart the counter for the next round
            else:
                count += 1

            print('job: ', i)
            self.istart = i*nfile_per_job
            self.iend = self.istart + nfile_per_job
            job_name = f"D0_ESE_out_{self.istart}_{self.iend}"
            subfile = self.build_job_files(job_name)
            if dry_run:
                print(f"Generated submission script: {subfile}")
            else:
                os.system(f"sbatch --time=04:00:00 -A standby {subfile}")

if __name__ == "__main__":

    # --- PROXY CHECK ---
    proxy_path = "/home/saha115/myproxy"
    min_proxy_hours = 24  # Ensure at least 24 hours of validity    
    print("Checking grid proxy...")
    proxy_check = f"voms-proxy-info -file {proxy_path} -timeleft"
    try:
        timeleft = int(subprocess.check_output(proxy_check, shell=True).decode('utf-8').strip())
    except (subprocess.CalledProcessError, ValueError):
        timeleft = 0
    if timeleft < (min_proxy_hours * 3600):
        print(f"Proxy is missing or expiring in less than {min_proxy_hours} hours.")
        os.system(f"voms-proxy-init --voms cms --valid 192:00 --out {proxy_path}")
    else:
        print(f"Proxy is valid for {timeleft // 3600} more hours. Proceeding...")
    # -------------------



    base_dir = "/scratch/negishi/saha115/D0_ESE_out/Jan16/"
    flist_name = "/home/saha115/D0_ESE/CMSSW_13_2_11/src/SlurmJob_submission/input_files.txt"
    #flist_name = "/home/saha115/Purdue/CMSSW_13_2_11/src/PbPbfiles_test.txt"

    sub = jobs_submission(flist_name)

    with open(flist_name) as f:
        flist = list(f)

    nfile_per_job = 5
    njobs = int(len(flist) / nfile_per_job) +1 
    print(f"Number of files in flist: {len(flist)}")
    print(f"Number of jobs: {njobs}")


    # Set dry_run to True to perform a dry run
    sub.submit_jobs(base_dir, njobs, nfile_per_job, dry_run=False) 

