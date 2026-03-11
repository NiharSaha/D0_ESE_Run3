import os

# --- SETTINGS ---
search_path = "/scratch/negishi/chand140/d0output_Milan/dataoutput_lowpt/"
output_filename = "root_files_Soumik.txt" 
# ----------------

def generate_file_list(input_dir, output_file):
    count = 0
    print(f"Scanning: {input_dir}...")
    
    try:
        with open(output_file, "w") as f:
            for root, dirs, files in os.walk(input_dir):
                for file in files:
                    if file.endswith(".root"):
                        # Construct the full absolute path
                        full_path = os.path.join(root, file)
                        f.write(full_path + "\n")
                        count += 1
        
        print(f"Success! {count} files found.")
        print(f"List saved to: {os.path.abspath(output_file)}")
        
    except PermissionError:
        print("Error: Permission denied. Check your access to the directory.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    generate_file_list(search_path, output_filename)
