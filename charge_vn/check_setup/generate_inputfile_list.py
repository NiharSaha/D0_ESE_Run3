import subprocess
import os

# --- Configuration ---
# The redirector for Purdue EOS
REDIRECTOR = "root://eos.cms.rcac.purdue.edu/"
# The base path on EOS (starting from /store/)
BASE_DIR = "/store/user/nsaha/2023PbPb_chargevn/HIPhysicsRawPrime1/"

# Output text file
OUTPUT_FILE = "inputFiles_charge_MB1.txt"

def get_files_recursive(directory):
    """Recursively lists files using the xrdfs ls command."""
    files_found = []
    
    # Run the xrdfs ls command to see contents
    cmd = ["xrdfs", REDIRECTOR, "ls", "-R", directory]
    
    try:
        print(f"Scanning {directory} ...")
        result = subprocess.check_output(cmd, stderr=subprocess.STDOUT).decode('utf-8')
        
        for line in result.splitlines():
            line = line.strip()
            # Only include .root files
            if line.endswith(".root"):
                # xrdfs returns the path without the redirector, so we prepend it
                files_found.append(f"{REDIRECTOR}/{line.lstrip('/')}")
                
    except subprocess.CalledProcessError as e:
        print(f"Error scanning directory: {e.output.decode('utf-8')}")
        
    return files_found

def main():
    # 1. Get the list of files
    all_files = get_files_recursive(BASE_DIR)
    
    # 2. Save to text file
    with open(OUTPUT_FILE, "w") as f:
        for file_path in all_files:
            f.write(file_path + "\n")
            
    print(f"\nSuccess! Found {len(all_files)} files.")
    print(f"File list saved to: {os.path.abspath(OUTPUT_FILE)}")

if __name__ == "__main__":
    main()
