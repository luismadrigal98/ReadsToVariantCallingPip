"""
Utilities for working with fastq files.

"""

import os

def decompress_and_split(input_dirs, output_dirs, job_dirs, lines=4000000):
    """
    Generates SLURM job scripts to decompress and split FASTQ files in given directories.
    
    Parameters:
    input_dirs (list of str): Directories containing the FASTQ files to be decompressed and split.
    output_dirs (list of str): Directories where the split files will be saved.
    job_dirs (list of str): Directories where the SLURM job scripts will be saved.
    lines (int): Number of lines to be included in each split file. Default is 4000000.
    
    Note: The number of lines should be a multiple of 4.
    """
    
    if not isinstance(input_dirs, list) or not isinstance(output_dirs, list) or not isinstance(job_dirs, list):
        raise ValueError("input_dirs, output_dirs, and job_dirs must be lists")
    
    if len(input_dirs) != len(output_dirs) or len(input_dirs) != len(job_dirs):
        raise ValueError("The number of input, output, and job directories must be the same")
    
    for dir_path in input_dirs + output_dirs + job_dirs:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
    
    for input_dir, output_dir, job_dir in zip(input_dirs, output_dirs, job_dirs):
        try:
            files = [f for f in os.listdir(input_dir) if f.endswith(".fastq.gz")]
            
            for file in files:
                name = file.replace(".fastq.gz", "_")
                job_script_path = os.path.join(job_dir, f"Split_{name}job.sh")
                
                fd = os.open(job_script_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o755)
                with os.fdopen(fd, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=split_{name}job\n")
                    script.write(f"#SBATCH --output=split_{name}output\n")
                    script.write("#SBATCH --partition=sixhour\n")
                    script.write("#SBATCH --nodes=1\n")
                    script.write("#SBATCH --ntasks=1\n")
                    script.write("#SBATCH --cpus-per-task=1\n")
                    script.write("#SBATCH --time=6:00:00\n")
                    script.write("#SBATCH --mail-user=l338m483@ku.edu\n")
                    script.write("#SBATCH --mail-type=END,FAIL\n")
                    script.write("#SBATCH --mem-per-cpu=5g\n")
                    script.write("\n")
                    script.write(f"cd {input_dir}\n")
                    script.write(f"zcat {file} | split -l {lines} - {output_dir}/{name}\n")
        
        except FileNotFoundError:
            print(f"Review the name format of the input files in {input_dir}")
        except Exception as e:
            print(f"An error occurred: {e}")