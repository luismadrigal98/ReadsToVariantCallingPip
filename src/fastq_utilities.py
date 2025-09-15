"""
Utilities for working with fastq files.

@author: Luis Javier Madrigal-Roca

"""

import os
import re
import logging

def create_directory(directory):
    """Create directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def detect_file_type(input_dir):
    """
    Analyze FASTQ files in a directory to determine if they are paired-end or single-end.
    
    Parameters:
    input_dir (str): Directory containing FASTQ files
    
    Returns:
    dict: Dictionary with lists of paired and single files
    """
    fastq_files = [f for f in os.listdir(input_dir) if (
                        f.endswith(('.fastq', '.fq')) or 
                        (f.endswith('.gz'))
                        )]
    
    # Initialize result structure
    result = {
        'paired_files': [],  # Will hold tuples of (R1, R2)
        'single_files': [],  # Will hold single-end files
    }
    
    # First identify potential pairs (Illumina naming convention)
    r1_files = [f for f in fastq_files if '_R1_' in f]
    r2_files = [f for f in fastq_files if '_R2_' in f]
    
    # Match pairs
    paired = set()
    for r1 in r1_files:
        r2_candidate = r1.replace('_R1_', '_R2_')
        if r2_candidate in r2_files:
            result['paired_files'].append((r1, r2_candidate))
            paired.add(r1)
            paired.add(r2_candidate)
    
    # All other files are considered single-end
    for f in fastq_files:
        if f not in paired:
            result['single_files'].append(f)
    
    logging.info(f"Found {len(result['paired_files'])} paired-end file sets and {len(result['single_files'])} single-end files")
    return result

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
    Handles both single-end and paired-end files.
    
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
            # Get all FASTQ files regardless of naming convention
            files = [f for f in os.listdir(input_dir) if f.endswith(('.fastq.gz', '.fq.gz'))]
            
            for file in files:
                # Create a clean name for the split files
                if file.endswith('.fastq.gz'):
                    name = file.replace(".fastq.gz", "_")
                elif file.endswith('.fq.gz'):
                    name = file.replace(".fq.gz", "_")
                
                job_script_path = os.path.join(job_dir, f"Split_{name}job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=split_{name}job\n")
                    script.write(f"#SBATCH --output={os.path.join(job_dir, f'split_{name}output')}\n")
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
                logging.info(f"Created split job script: {job_script_path}")
        
        except FileNotFoundError:
            logging.error(f"Review the name format of the input files in {input_dir}")
        except Exception as e:
            logging.error(f"An error occurred: {e}")

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

            if not filenames:
                logging.warning(f"No files found in {batch_dir} - this is expected when not using --submit")
            
            # Get unique filename roots by removing the last part (usually aa, ab, ac, etc.)
            filenames_shortened = [name[:-3] for name in filenames if len(name) > 3]
            unique_filenames = list(set(filenames_shortened))
            
            for locator in unique_filenames:
                job_script_path = os.path.join(job_dir, f"gzip_in_batch_{locator}_compress_job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=gzip_in_batch_{locator}_job\n")
                    script.write(f"#SBATCH --output={os.path.join(job_dir, f'gzip_in_batch_{locator}_output')}\n")
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
                logging.info(f"Created compression job script: {job_script_path}")
        
        except FileNotFoundError:
            logging.error(f"Directory not found: {batch_dir}")
        except Exception as e:
            logging.error(f"An error occurred: {e}")

def generate_fastp_jobs(batch_dirs, output_dirs, job_dirs, fastp_path="/home/l338m483/fastp",
                        fastp_control_param="",
                        partition="sixhour", time="6:00:00", email="l338m483@ku.edu", 
                        mem_per_cpu="5g", cpus=10):
    """
    Generate SLURM job scripts to run fastp on FASTQ files, handling both paired and single-end reads.
    
    Parameters:
    batch_dirs (list): List of directories containing FASTQ files
    output_dirs (list): List of directories for fastp output
    job_dirs (list): List of directories for job scripts
    fastp_path (str): Path to fastp executable
    fastp_options (str): Additional options to pass to fastp
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
            # Detect file types in the directory
            file_types = detect_file_type(batch_dir)
            
            # Process paired-end files
            for r1_file, r2_file in file_types['paired_files']:
                # Extract a meaningful name for the job
                # Try to get the sample ID from the filename (assumes format like SAMPLE_S1_L001_R1_001.fastq.gz)
                sample_id = r1_file.split('_')[0]
                
                out1 = r1_file.replace(".gz", "_preprocessed.fastq.gz")
                out2 = r2_file.replace(".gz", "_preprocessed.fastq.gz")
                
                # Generate paired-end fastp command
                fastp_command = (f"{fastp_path} -w {cpus} --correction "
                                    f"-i {os.path.join(batch_dir, r1_file)} "
                                    f"-o {output_dir}/{out1} "
                                    f"-I {os.path.join(batch_dir, r2_file)} "
                                    f"-O {output_dir}/{out2} {fastp_control_param} "
                                    f"-R {os.path.join(output_dir, f'{sample_id}_fastp_report')}")
                
                job_script_path = os.path.join(job_dir, f"fastp_paired_{sample_id}_job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=fastp_paired_{sample_id}_job\n")
                    script.write(f"#SBATCH --output={os.path.join(job_dir, f'fastp_paired_{sample_id}_output')}\n")
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
                logging.info(f"Created paired-end fastp job script: {job_script_path}")
            
            # Process single-end files
            for single_file in file_types['single_files']:
                # Extract a meaningful name for the job
                sample_id = single_file.split('.')[0]  # Assumes format like 1192c.fq.gz
                
                out_file = single_file.replace(".gz", "_preprocessed.fastq.gz")
                
                # Generate single-end fastp command (note: no -I or -O parameters)
                fastp_command = (f"{fastp_path} -w {cpus} "
                                    f"-i {os.path.join(batch_dir, single_file)} "
                                    f"-o {output_dir}/{out_file} {fastp_control_param} "
                                    f"-R {os.path.join(output_dir, f'{sample_id}_fastp_report')}")
                
                job_script_path = os.path.join(job_dir, f"fastp_single_{sample_id}_job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=fastp_single_{sample_id}_job\n")
                    script.write(f"#SBATCH --output={os.path.join(job_dir, f'fastp_single_{sample_id}_output')}\n")
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
                logging.info(f"Created single-end fastp job script: {job_script_path}")
                
        except FileNotFoundError:
            logging.error(f"Directory not found: {batch_dir}")
        except Exception as e:
            logging.error(f"An error occurred: {e}")