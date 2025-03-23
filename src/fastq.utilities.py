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

def create_directory(directory):
    """Create directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)


def get_sample_structure(input_directories):
    """
    Analyze the input directories to determine sample structure.
    Returns a dictionary mapping sample IDs to their input directories.
    """
    samples = {}
    pattern = re.compile(r'(\d{4})([A-Z])?_sequences')
    
    for directory in input_directories:
        match = pattern.search(directory)
        if match:
            year = match.group(1)
            sample_id = match.group(2) if match.group(2) else "default"
            sample_key = f"{year}_{sample_id}"
            samples[sample_key] = directory
    
    return samples


def generate_split_jobs(input_dirs, output_dirs, job_dirs, lines=4000000, partition="sixhour", 
                        time="6:00:00", email="l338m483@ku.edu", mem_per_cpu="5g", cpus=10):
    """
    Generate SLURM job scripts to decompress and split FASTQ files.
    
    Parameters:
    input_dirs (list): List of directories containing FASTQ files
    output_dirs (list): List of directories for output split files
    job_dirs (list): List of directories for job scripts
    lines (int): Number of lines per split file (default 4000000)
    partition (str): SLURM partition to use
    time (str): Time limit for jobs
    email (str): Email for notifications
    mem_per_cpu (str): Memory per CPU
    cpus (int): Number of CPUs per task
    """
    
    for i, (input_dir, output_dir, job_dir) in enumerate(zip(input_dirs, output_dirs, job_dirs)):
        # Create directories if they don't exist
        create_directory(output_dir)
        create_directory(job_dir)
        
        try:
            files = [f for f in os.listdir(input_dir) if f.endswith(".fastq.gz")]
            
            for file in files:
                name = file.replace(".fastq.gz", "_")
                job_script_path = os.path.join(job_dir, f"Split_{name}job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=split_{name}job\n")
                    script.write(f"#SBATCH --output=split_{name}output\n")
                    script.write(f"#SBATCH --partition={partition}\n")
                    script.write("#SBATCH --nodes=1\n")
                    script.write("#SBATCH --ntasks=1\n")
                    script.write(f"#SBATCH --cpus-per-task={cpus}\n")
                    script.write(f"#SBATCH --time={time}\n")
                    script.write(f"#SBATCH --mail-user={email}\n")
                    script.write("#SBATCH --mail-type=END,FAIL\n")
                    script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                    script.write("\n")
                    script.write(f"cd {input_dir}\n")
                    script.write(f"zcat {file} | split -l {lines} - {output_dir}/{name}\n")
                
                os.chmod(job_script_path, 0o755)
                print(f"Created split job script: {job_script_path}")
        
        except FileNotFoundError:
            print(f"Review the name format of the input files in {input_dir}")
        except Exception as e:
            print(f"An error occurred: {e}")


def generate_compress_jobs(batch_dirs, job_dirs, partition="sixhour", 
                          time="6:00:00", email="l338m483@ku.edu", 
                          mem_per_cpu="5g", cpus=10):
    """
    Generate SLURM job scripts to compress split FASTQ files.
    
    Parameters:
    batch_dirs (list): List of directories containing split files
    job_dirs (list): List of directories for job scripts
    partition (str): SLURM partition to use
    time (str): Time limit for jobs
    email (str): Email for notifications
    mem_per_cpu (str): Memory per CPU
    cpus (int): Number of CPUs per task
    """
    
    for batch_dir, job_dir in zip(batch_dirs, job_dirs):
        # Create job directory if it doesn't exist
        create_directory(job_dir)
        
        try:
            filenames = os.listdir(batch_dir)
            
            # Get unique filename roots by removing the last part (usually aa, ab, ac, etc.)
            filenames_shortened = [name[:-3] for name in filenames if len(name) > 3]
            unique_filenames = list(set(filenames_shortened))
            
            for locator in unique_filenames:
                job_script_path = os.path.join(job_dir, f"gzip_in_batch_{locator}_compress_job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=gzip_in_batch_{locator}_job\n")
                    script.write(f"#SBATCH --output=gzip_in_batch_{locator}_output\n")
                    script.write(f"#SBATCH --partition={partition}\n")
                    script.write("#SBATCH --nodes=1\n")
                    script.write("#SBATCH --ntasks=1\n")
                    script.write(f"#SBATCH --cpus-per-task={cpus}\n")
                    script.write(f"#SBATCH --time={time}\n")
                    script.write(f"#SBATCH --mail-user={email}\n")
                    script.write("#SBATCH --mail-type=END,FAIL\n")
                    script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                    script.write("\n")
                    script.write(f"cd {batch_dir}\n")
                    script.write("\n")
                    script.write(f"gzip {locator}*\n")
                
                os.chmod(job_script_path, 0o755)
                print(f"Created compression job script: {job_script_path}")
        
        except FileNotFoundError:
            print(f"Directory not found: {batch_dir}")
        except Exception as e:
            print(f"An error occurred: {e}")


def generate_fastp_jobs(batch_dirs, output_dirs, job_dirs, fastp_path="/home/l338m483/fastp",
                        partition="sixhour", time="6:00:00", email="l338m483@ku.edu", 
                        mem_per_cpu="5g", cpus=10):
    """
    Generate SLURM job scripts to run fastp on paired FASTQ files.
    
    Parameters:
    batch_dirs (list): List of directories containing FASTQ files
    output_dirs (list): List of directories for fastp output
    job_dirs (list): List of directories for job scripts
    fastp_path (str): Path to fastp executable
    partition (str): SLURM partition to use
    time (str): Time limit for jobs
    email (str): Email for notifications
    mem_per_cpu (str): Memory per CPU
    cpus (int): Number of CPUs per task
    """
    
    for batch_dir, output_dir, job_dir in zip(batch_dirs, output_dirs, job_dirs):
        # Create directories if they don't exist
        create_directory(output_dir)
        create_directory(job_dir)
        
        try:
            filesR1 = []
            filesR2 = []
            
            filenames = os.listdir(batch_dir)
            
            for filename in filenames:
                if "R1" in filename:
                    filesR1.append(filename)
                elif "R2" in filename:
                    filesR2.append(filename)
            
            # Ensure filesR1 and filesR2 are sorted so paired files match up
            filesR1.sort()
            filesR2.sort()
            
            if len(filesR1) != len(filesR2):
                print(f"Warning: Number of R1 files ({len(filesR1)}) does not match number of R2 files ({len(filesR2)}) in {batch_dir}")
            
            for i in range(len(filesR1)):
                if i < len(filesR2):  # Make sure we have a matching R2 file
                    r1_file = filesR1[i]
                    r2_file = filesR2[i]
                    
                    out1 = r1_file.replace(".gz", "_preprocessed.fastq.gz")
                    out2 = r2_file.replace(".gz", "_preprocessed.fastq.gz")
                    name = r1_file.replace("_R1_001", "")
                    
                    fastp_command = (f"{fastp_path} -w {cpus} --correction -i {r1_file} -o {output_dir}/{out1} "
                                    f"-I {r2_file} -O {output_dir}/{out2} -3 --complexity_threshold=20 "
                                    f"--length_required=50 --cut_window_size=3 --cut_mean_quality=30")
                    
                    job_script_path = os.path.join(job_dir, f"fastp_in_batch_{name}_job.sh")
                    
                    with open(job_script_path, "w") as script:
                        script.write("#!/bin/bash\n")
                        script.write(f"#SBATCH --job-name=fastp_in_batch_{name}_job\n")
                        script.write(f"#SBATCH --output=fastp_in_batch_{name}_output\n")
                        script.write(f"#SBATCH --partition={partition}\n")
                        script.write("#SBATCH --nodes=1\n")
                        script.write("#SBATCH --ntasks=1\n")
                        script.write(f"#SBATCH --cpus-per-task={cpus}\n")
                        script.write(f"#SBATCH --time={time}\n")
                        script.write(f"#SBATCH --mail-user={email}\n")
                        script.write("#SBATCH --mail-type=FAIL\n")
                        script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                        script.write("\n")
                        script.write(f"cd {batch_dir}\n")
                        script.write("\n")
                        script.write(fastp_command)
                    
                    os.chmod(job_script_path, 0o755)
                    print(f"Created fastp job script: {job_script_path}")
        
        except FileNotFoundError:
            print(f"Directory not found: {batch_dir}")
        except Exception as e:
            print(f"An error occurred: {e}")