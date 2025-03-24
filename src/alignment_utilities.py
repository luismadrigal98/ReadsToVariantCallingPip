"""
Utilities for running alignment tools (BWA and Stampy) on FASTQ files.
"""

import os
import logging
import re

def create_directory(directory):
    """Create directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def detect_preprocessed_files(input_dir):
    """
    Analyze preprocessed FASTQ files in a directory to identify paired and single-end files.
    
    Parameters:
    input_dir (str): Directory containing preprocessed FASTQ files
    
    Returns:
    dict: Dictionary with lists of paired and single files
    """
    fastq_files = [f for f in os.listdir(input_dir) if f.endswith(('_preprocessed.fastq.gz'))]
    
    # Initialize result structure
    result = {
        'paired_files': [],  # Will hold tuples of (R1, R2)
        'single_files': [],  # Will hold single-end files
    }
    
    # Identify R1 files (for paired-end reads)
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

def index_bwa_reference(reference_genome, bwa_path):
    """
    Index the reference genome with BWA.
    
    Parameters:
    reference_genome (str): Path to the reference genome (FASTA)
    bwa_path (str): Path to BWA executable
    """
    try:
        # BWA index command
        index_command = f"{bwa_path} index {reference_genome}"
        
        # Run the command
        os.system(index_command)
        logging.info(f"Indexed reference genome with BWA: {reference_genome}")
        
    except FileNotFoundError:
        logging.error(f"File not found: {reference_genome}")
    except Exception as e:
        logging.error(f"An error occurred: {e}")

def detect_bam_files(input_dir):
    """
    Find sorted BAM files for Stampy processing.
    
    Parameters:
    input_dir (str): Directory containing BAM files
    
    Returns:
    list: List of sorted BAM files
    """
    bam_files = [f for f in os.listdir(input_dir) if f.endswith('_bwa-sorted.bam')]
    logging.info(f"Found {len(bam_files)} sorted BAM files for Stampy processing")
    return bam_files

def generate_bwa_jobs(input_dirs, output_dirs, job_dirs, reference_genome,
                    bwa_path="/kuhpc/sw/conda/latest/envs/bioconda/bin/bwa", 
                    partition="sixhour", time="6:00:00", 
                    email="l338m483@ku.edu", mem_per_cpu="5g", cpus=10):
    """
    Generate SLURM job scripts to run BWA on preprocessed FASTQ files.
    
    Parameters:
    input_dirs (list): List of directories containing preprocessed FASTQ files
    output_dirs (list): List of directories for BWA output (SAM files)
    job_dirs (list): List of directories for job scripts
    reference_genome (str): Path to the reference genome (indexed with BWA)
    bwa_path (str): Path to BWA executable (default: "bwa")
    partition (str): SLURM partition to use
    time (str): Time limit for jobs
    email (str): Email for notifications
    mem_per_cpu (str): Memory per CPU
    cpus (int): Number of CPUs per task
    """
    for input_dir, output_dir, job_dir in zip(input_dirs, output_dirs, job_dirs):
        # Create directories if they don't exist
        create_directory(output_dir)
        create_directory(job_dir)
        
        try:
            # Detect file types in the directory
            file_types = detect_preprocessed_files(input_dir)
            
            # Process paired-end files
            for r1_file, r2_file in file_types['paired_files']:
                # Extract a meaningful name for the job and output
                sample_name = re.sub(r'_R1.*preprocessed\.fastq\.gz', '', r1_file)
                
                # Output file name
                out_sam = f"{sample_name}_bwa-aligned.sam.gz"
                
                # BWA command for paired-end reads
                bwa_command = (f"{bwa_path} mem -t {cpus} {reference_genome} "
                                f"{os.path.join(input_dir, r1_file)} "
                                f"{os.path.join(input_dir, r2_file)} | "
                                f"gzip -3 > {os.path.join(output_dir, out_sam)}")
                
                job_script_path = os.path.join(job_dir, f"bwa_{sample_name}_alig_job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=bwa_{sample_name}_job\n")
                    script.write(f"#SBATCH --output={os.path.join(job_dir, f'bwa_{sample_name}_output')}\n")
                    script.write(f"#SBATCH --partition={partition}\n")
                    script.write("#SBATCH --nodes=1\n")
                    script.write("#SBATCH --ntasks=1\n")
                    script.write(f"#SBATCH --cpus-per-task={cpus}\n")
                    script.write(f"#SBATCH --time={time}\n")
                    script.write(f"#SBATCH --mail-user={email}\n")
                    script.write("#SBATCH --mail-type=FAIL\n")
                    script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                    script.write("\n")
                    script.write(f"cd {input_dir}\n")
                    script.write("\n")
                    script.write(bwa_command)
                
                os.chmod(job_script_path, 0o755)
                logging.info(f"Created BWA paired-end job script: {job_script_path}")
            
            # Process single-end files
            for single_file in file_types['single_files']:
                # Extract a meaningful name for the job and output
                sample_name = re.sub(r'_preprocessed\.fastq\.gz', '', single_file)
                
                # Output file name
                out_sam = f"{sample_name}_bwa-aligned.sam.gz"
                
                # BWA command for single-end reads
                bwa_command = (f"{bwa_path} mem -t {cpus} {reference_genome} "
                                f"{os.path.join(input_dir, single_file)} | "
                                f"gzip -3 > {os.path.join(output_dir, out_sam)}")
                
                job_script_path = os.path.join(job_dir, f"bwa_{sample_name}_alig_job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=bwa_{sample_name}_job\n")
                    script.write(f"#SBATCH --output={os.path.join(job_dir, f'bwa_{sample_name}_output')}\n")
                    script.write(f"#SBATCH --partition={partition}\n")
                    script.write("#SBATCH --nodes=1\n")
                    script.write("#SBATCH --ntasks=1\n")
                    script.write(f"#SBATCH --cpus-per-task={cpus}\n")
                    script.write(f"#SBATCH --time={time}\n")
                    script.write(f"#SBATCH --mail-user={email}\n")
                    script.write("#SBATCH --mail-type=FAIL\n")
                    script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                    script.write("\n")
                    script.write(f"cd {input_dir}\n")
                    script.write("\n")
                    script.write(bwa_command)
                
                os.chmod(job_script_path, 0o755)
                logging.info(f"Created BWA single-end job script: {job_script_path}")
                
        except FileNotFoundError:
            logging.error(f"Directory not found: {input_dir}")
        except Exception as e:
            logging.error(f"An error occurred: {e}")

def generate_sam_to_bam_jobs(input_dirs, output_dirs, job_dirs, 
                            samtools_path="samtools", partition="sixhour", 
                            time="6:00:00", email="l338m483@ku.edu", 
                            mem_per_cpu="5g", cpus=4):
    """
    Generate SLURM job scripts to convert SAM files to sorted BAM files.
    
    Parameters:
    input_dirs (list): List of directories containing SAM files (.sam.gz)
    output_dirs (list): List of directories for BAM output
    job_dirs (list): List of directories for job scripts
    samtools_path (str): Path to samtools executable (default: "samtools")
    partition (str): SLURM partition to use
    time (str): Time limit for jobs
    email (str): Email for notifications
    mem_per_cpu (str): Memory per CPU
    cpus (int): Number of CPUs per task
    """
    for input_dir, output_dir, job_dir in zip(input_dirs, output_dirs, job_dirs):
        # Create directories if they don't exist
        create_directory(output_dir)
        create_directory(job_dir)
        
        try:
            # Find all compressed SAM files
            sam_files = [f for f in os.listdir(input_dir) if f.endswith('_bwa-aligned.sam.gz')]
            
            for sam_file in sam_files:
                # Extract a meaningful name
                sample_name = sam_file.replace("_bwa-aligned.sam.gz", "")
                
                # Output BAM file name
                out_bam = f"{sample_name}_bwa-sorted.bam"
                
                # Command to convert SAM to sorted BAM
                bam_command = (f"zcat {os.path.join(input_dir, sam_file)} | "
                                f"{samtools_path} view -@ {cpus} -bS - | "
                                f"{samtools_path} sort -@ {cpus} -o {os.path.join(output_dir, out_bam)} -")
                
                job_script_path = os.path.join(job_dir, f"bam_{sample_name}_sort_job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=bam_{sample_name}_job\n")
                    script.write(f"#SBATCH --output={os.path.join(job_dir, f'bam_{sample_name}_output')}\n")
                    script.write(f"#SBATCH --partition={partition}\n")
                    script.write("#SBATCH --nodes=1\n")
                    script.write("#SBATCH --ntasks=1\n")
                    script.write(f"#SBATCH --cpus-per-task={cpus}\n")
                    script.write(f"#SBATCH --time={time}\n")
                    script.write(f"#SBATCH --mail-user={email}\n")
                    script.write("#SBATCH --mail-type=FAIL\n")
                    script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                    script.write("\n")
                    script.write(f"cd {input_dir}\n")
                    script.write("\n")
                    script.write(bam_command)
                
                os.chmod(job_script_path, 0o755)
                logging.info(f"Created SAM to BAM job script: {job_script_path}")
                
        except FileNotFoundError:
            logging.error(f"Directory not found: {input_dir}")
        except Exception as e:
            logging.error(f"An error occurred: {e}")

def index_stampy_reference(reference_genome, python_2_7_path, stampy_path):
    """
    Index the reference genome with Stampy.
    
    Parameters:
    reference_genome (str): Path to the reference genome (FASTA)
    python_2_7_path (str): Path to Python 2.7 executable
    stampy_path (str): Path to Stampy executable
    """
    try:
        # Reference name without extension (for Stampy -g -h options)
        ref_base = os.path.splitext(reference_genome)[0]
        
        # Stampy indexing command
        index_command = f"{python_2_7_path} {stampy_path} -G {ref_base} {reference_genome}"
        
        # Run the command
        os.system(index_command)
        logging.info(f"Indexed reference genome with Stampy: {reference_genome}")

        # Stampy hash command
        hash_command = f"{python_2_7_path} {stampy_path} -g {ref_base} -H {ref_base}"

        # Run the command
        os.system(hash_command)
        logging.info(f"Hashed reference genome with Stampy: {reference_genome}")
        
    except FileNotFoundError:
        logging.error(f"File not found: {reference_genome}")
    except Exception as e:
        logging.error(f"An error occurred: {e}")

def generate_stampy_jobs(input_dirs, output_dirs, job_dirs, reference_genome,
                        python_2_7_path,
                        stampy_path, 
                        samtools_path,
                        partition="sixhour", time="6:00:00", 
                        email="l338m483@ku.edu", mem_per_cpu="5g", cpus=3):
    """
    Generate SLURM job scripts to run Stampy on sorted BAM files.
    
    Parameters:
    input_dirs (list): List of directories containing sorted BAM files
    output_dirs (list): List of directories for Stampy output
    job_dirs (list): List of directories for job scripts
    reference_genome (str): Path to the reference genome (with .stidx and .sthash)
    python_27_path (str): Path to Python 2.7 executable
    stampy_path (str): Path to Stampy executable
    partition (str): SLURM partition to use
    time (str): Time limit for jobs
    email (str): Email for notifications
    mem_per_cpu (str): Memory per CPU
    cpus (int): Number of CPUs per task (note: Stampy may not use all cores effectively)
    """
    for input_dir, output_dir, job_dir in zip(input_dirs, output_dirs, job_dirs):
        # Create directories if they don't exist
        create_directory(output_dir)
        create_directory(job_dir)
        
        try:
            # Get all sorted BAM files
            bam_files = detect_bam_files(input_dir)
            
            # Reference name without extension (for Stampy -g -h options)
            ref_base = os.path.splitext(reference_genome)[0]
            
            for bam_file in bam_files:
                # Extract a meaningful name
                sample_name = bam_file.replace("_bwa-sorted.bam", "")
                
                # Output file name
                out_bam = f"{sample_name}_stampy.bam"
                
                # Stampy command
                stampy_command = (f"{python_2_7_path} {stampy_path} -t {cpus} --sensitive "
                                f"-g {ref_base} -h {ref_base} --bamkeepgoodreads "
                                f"-M {os.path.join(input_dir, bam_file)} | "
                                f"{samtools_path} view -Sb > {os.path.join(output_dir, out_bam)}")
                
                job_script_path = os.path.join(job_dir, f"stampy_{sample_name}_job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=stampy_{sample_name}_job\n")
                    script.write(f"#SBATCH --output={os.path.join(job_dir, f'stampy_{sample_name}_output')}\n")
                    script.write(f"#SBATCH --partition={partition}\n")
                    script.write("#SBATCH --nodes=1\n")
                    script.write("#SBATCH --ntasks=1\n")
                    script.write(f"#SBATCH --cpus-per-task={cpus}\n")
                    script.write(f"#SBATCH --time={time}\n")
                    script.write(f"#SBATCH --mail-user={email}\n")
                    script.write("#SBATCH --mail-type=FAIL\n")
                    script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                    script.write("\n")
                    script.write(f"cd {input_dir}\n")
                    script.write("\n")
                    script.write(stampy_command)
                
                os.chmod(job_script_path, 0o755)
                logging.info(f"Created Stampy job script: {job_script_path}")
                
        except FileNotFoundError:
            logging.error(f"Directory not found: {input_dir}")
        except Exception as e:
            logging.error(f"An error occurred: {e}")