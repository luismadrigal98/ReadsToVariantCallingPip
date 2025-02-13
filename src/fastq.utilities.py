"""
Utilities for working with fastq files.

"""

import os

def decompress_and_split(input_dir, output_dir, job_dir, lines = 4000000):
    """
    Generates a SLURM job script to decompress and split in several smaller fils the input fastq files in a given directory.
    
    Parameters:
    
    input_dir (str): The directory or directories containing the fastq files to be decompressed and split.

    output_dir (str): The directory where the split files will be saved. Note that if multiple input directories are provided, the number of output directories should be the same as the number of input directories.

    job_dir (str): The directory where the SLURM job scripts will be saved.

    lines (int): The number of lines to be included in each split file. Default

    Note: The number of lines should be a multiple of 4.
    
    """
    
    try:
        
        # Does the input directory exist?
        assert os.path.exists(input_dir), "The input directory does not exist"

        # Does the output directory exist?
        try:
            for dir_ in range(len(output_dir)):
                assert os.path.exists(output_dir[dir_]), "The output directory does not exist"
        except:
            for dir_ in range(len(output_dir)):
                os.mkdir(output_dir[dir_])

        # Does the job directory exist?
        try:
            for dir_ in range(len(job_dir)):
                assert os.path.exists(job_dir[dir_]), "The job directory does not exist"
        except:
            for dir_ in range(len(job_dir)):
                os.mkdir(job_dir[dir_])

        assert len(input_dir) == len(output_dir), "The number of input directories should be the same as the number of output directories"
        
        try:
            for dir_ in range(len(input_dir)):
                files = []
                filenames = os.listdir(input_dir[dir_])

                for filename in filenames:
                    if filename.endswith(".fastq.gz"):
                        files.append(filename)

                for i in range(len(files)):
                    name = files[i].replace(".fastq.gz", "_")
                    out = f"{output_dir[dir_]}"

                    with open(f"{jobs[dir_]}/Split_{name}job.sh", "w") as scriptshell:
                        scriptshell.write("#!/bin/bash\n")
                        scriptshell.write(f"#SBATCH --job-name=split_{name}job\n") # Job name
                        scriptshell.write(f"#SBATCH --output=split_{name}output\n") # Output file
                        scriptshell.write("#SBATCH --partition=sixhour\n")
                        scriptshell.write("#SBATCH --nodes=1\n") # Number of nodes
                        scriptshell.write("#SBATCH --ntasks=1\n") # Number of tasks
                        scriptshell.write("#SBATCH --cpus-per-task=1\n") # Number of parallel processes (1 given that this is a serial job)
                        scriptshell.write("#SBATCH --time=6:00:00\n") # Time limit
                        scriptshell.write("#SBATCH --mail-user=l338m483@ku.edu\n") # Email for receiving notifications
                        scriptshell.write("#SBATCH --mail-type=END,FAIL\n") # Mail type
                        scriptshell.write("#SBATCH --mem-per-cpu=5g\n") # memory limit
                        scriptshell.write("\n")
                        scriptshell.write(f"cd {input_dir[dir_]}\n")
                        scriptshell.write(f"zcat {files[i]} | split -l {lines} - {out}/{name}\n")

                    os.system(f"chmod +x {job_dir[dir_]}/Split_{name}job.sh")

        except FileNotFoundError:
            print("Review the name format of the input files")
    except:
        AssertionError