"""
This module provides utility functions for handling BAM files and generating SLURM job scripts 
for various tasks such as merging, sorting, duplicate processing, and indexing.
Functions:
- create_directory(directory): Ensures a directory exists, creating it if necessary.
- detect_bam_files(input_dir, bam_type=None): Detects BAM files in a directory based on type.
- generate_merge_jobs(input_dirs, output_dirs, job_dirs, ...): Creates SLURM job scripts for merging and sorting BAM files.
- generate_duplicate_processing_jobs(input_dirs, output_dirs_base, job_dirs, ...): Creates SLURM job scripts for handling duplicates in BAM files.
- generate_indexing_jobs(input_dirs, job_dirs, ...): Creates SLURM job scripts for indexing BAM files.
Each function is designed to automate specific steps in BAM file processing pipelines, 
with support for SLURM job scheduling and customizable parameters.
"""

import os
import re
import logging
from pathlib import Path

def create_directory(directory):
    """Create directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)
        logging.info(f"Created directory: {directory}")

def detect_bam_files(input_dir, bam_type=None):
    """
    Find BAM files for merging.
    
    Parameters:
    input_dir (str): Directory containing BAM files
    bam_type (str): Type of BAM files to detect ('bwa', 'stampy', or None for both)
    
    Returns:
    list: List of BAM files
    """
    if bam_type == 'bwa':
        bam_files = [f for f in os.listdir(input_dir) if f.endswith('_bwa-sorted.bam')]
    elif bam_type == 'stampy':
        bam_files = [f for f in os.listdir(input_dir) if f.endswith('_stampy.bam')]
    else:
        bam_files = [f for f in os.listdir(input_dir) if f.endswith(('.bam'))]
        
    logging.info(f"Found {len(bam_files)} BAM files of type '{bam_type or 'any'}' for processing")
    return bam_files

def generate_merge_jobs(input_dirs, output_dirs, job_dirs, 
                       samtools_path="samtools", bam_type=None,
                       partition="sixhour", time="6:00:00", 
                       email="l338m483@ku.edu", mem_per_cpu="5g", cpus=14):
    """
    Generate SLURM job scripts to merge and sort BAM files.
    
    Parameters:
    input_dirs (list): List of directories containing BAM files
    output_dirs (list): List of directories for merged output
    job_dirs (list): List of directories for job scripts
    samtools_path (str): Path to samtools executable
    bam_type (str): Type of BAM files to merge ('bwa', 'stampy', or None for both)
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
            # Get all BAM files
            bam_files = detect_bam_files(input_dir, bam_type)
            
            if not bam_files:
                logging.warning(f"No {bam_type or 'BAM'} files found in {input_dir}")
                continue
                
            # Determine the output file name based on the first input file
            if bam_type == 'bwa' or (bam_type is None and 'bwa' in bam_files[0]):
                out = re.sub(r"_L[0-9]{3}_SET[0-9]_[a-z]{2}_bwa-sorted\.bam", "_bwa_merged.bam", bam_files[0])
                name = re.sub(r"_L[0-9]{3}_SET[0-9]_[a-z]{2}_bwa-sorted\.bam", "_bwa", bam_files[0])
            elif bam_type == 'stampy' or (bam_type is None and 'stampy' in bam_files[0]):
                out = re.sub(r"_L[0-9]{3}_SET[0-9]_[a-z]{2}_stampy\.bam", "_stampy_merged.bam", bam_files[0])
                name = re.sub(r"_L[0-9]{3}_SET[0-9]_[a-z]{2}_stampy\.bam", "_stampy", bam_files[0])
            else:
                # Generic fallback for other naming patterns
                out = f"{Path(input_dir).name}_merged.bam"
                name = f"{Path(input_dir).name}"
            
            # Full paths to files
            output_merged = os.path.join(output_dir, out)
            output_sorted = os.path.join(output_dir, out.replace('.bam', '_sorted.bam'))
            
            # Create commands
            bam_paths = ' '.join([os.path.join(input_dir, f) for f in bam_files])
            merging_command = f"{samtools_path} merge -@ {cpus} -o {output_merged} {bam_paths}"
            sorting_command = f"{samtools_path} sort -@ {cpus} -o {output_sorted} {output_merged}"
            
            # Create job script
            job_script_path = os.path.join(job_dir, f"{name}_merge_and_sort_bam_job.sh")
            
            with open(job_script_path, "w") as script:
                script.write("#!/bin/bash\n")
                script.write(f"#SBATCH --job-name=merge_{name}_job\n")
                script.write(f"#SBATCH --output={os.path.join(job_dir, f'merge_{name}_output')}\n")
                script.write(f"#SBATCH --partition={partition}\n")
                script.write("#SBATCH --nodes=1\n")
                script.write("#SBATCH --ntasks=1\n")
                script.write(f"#SBATCH --cpus-per-task={cpus}\n")
                script.write(f"#SBATCH --time={time}\n")
                script.write(f"#SBATCH --mail-user={email}\n")
                script.write("#SBATCH --mail-type=END,FAIL\n")
                script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                script.write("\n")
                script.write("# Ensure samtools is in PATH\n")
                script.write(f"export PATH=$(dirname {samtools_path}):$PATH\n")
                script.write(f"cd {input_dir}\n")
                script.write("\n")
                script.write(f"{merging_command}\n")
                script.write("\n")
                script.write(f"{sorting_command}\n")
            
            os.chmod(job_script_path, 0o755)
            logging.info(f"Created merge job script: {job_script_path}")
            
        except Exception as e:
            logging.error(f"Error processing directory {input_dir}: {e}")

def generate_duplicate_processing_jobs(input_dirs, output_dirs_base, job_dirs, 
                                     duplicate_mode="mark_remove", 
                                     samtools_path="samtools", picard_path="picard",
                                     partition="sixhour", time="16:00:00", 
                                     email="l338m483@ku.edu", mem_per_cpu="20g", cpus=5):
    """
    Generate SLURM job scripts for BAM duplicate processing and filtering.
    
    Parameters:
    input_dirs (list): List of directories containing merged BAM files
    output_dirs_base (list): List of base directories for output
    job_dirs (list): List of directories for job scripts
    duplicate_mode (str): How to handle duplicates ('preserve', 'mark', 'mark_remove')
    samtools_path (str): Path to samtools executable
    picard_path (str): Path to picard executable
    partition (str): SLURM partition to use
    time (str): Time limit for jobs
    email (str): Email for notifications
    mem_per_cpu (str): Memory per CPU
    cpus (int): Number of CPUs per task
    """
    # Map duplicate mode to directory names and processing commands
    mode_configs = {
        "preserve": {
            "dir_suffix": "Pipeline.with.duplicates",
            "name_suffix": "_pwd",
            "picard_cmd": None,  # No Picard command for this mode
        },
        "mark": {
            "dir_suffix": "Pipeline.with.marked.duplicates",
            "name_suffix": "_pwmd",
            "picard_cmd": f"{picard_path} MarkDuplicates TAGGING_POLICY=All {{input}} {{output}} {{metrics}} ASSUME_SORT_ORDER=coordinate",
        },
        "mark_remove": {
            "dir_suffix": "Pipeline.without.duplicates",
            "name_suffix": "_pwod",
            "picard_cmd": f"{picard_path} MarkDuplicates TAGGING_POLICY=All {{input}} {{output}} {{metrics}} ASSUME_SORT_ORDER=coordinate REMOVE_DUPLICATES=true",
        }
    }
    
    # Ensure the requested mode is valid
    if duplicate_mode not in mode_configs:
        logging.error(f"Invalid duplicate mode: {duplicate_mode}. Must be one of {list(mode_configs.keys())}")
        return
    
    # Get the configuration for the selected mode
    config = mode_configs[duplicate_mode]
    
    for input_dir, output_dir_base, job_dir in zip(input_dirs, output_dirs_base, job_dirs):
        create_directory(job_dir)
        
        try:
            # Find sorted merged BAM files
            merged_files = [f for f in os.listdir(input_dir) if f.endswith('_merged_sorted.bam')]
            
            if not merged_files:
                logging.warning(f"No merged sorted BAM files found in {input_dir}")
                continue
            
            for merged_file in merged_files:
                # Determine if this is BWA or Stampy
                if 'bwa' in merged_file:
                    output_dir = os.path.join(output_dir_base, config["dir_suffix"], "bwa")
                elif 'stampy' in merged_file:
                    output_dir = os.path.join(output_dir_base, config["dir_suffix"], "stampy")
                else:
                    output_dir = os.path.join(output_dir_base, config["dir_suffix"], "other")
                
                create_directory(output_dir)
                
                # Generate output file names
                name = merged_file.replace("_merged_sorted.bam", config["name_suffix"])
                out_1 = merged_file.replace("_merged_sorted.bam", f"{config['name_suffix']}_processed_merged_sorted.bam")
                out_2 = out_1.replace(f"{config['name_suffix']}_processed_merged_sorted.bam", 
                                     f"{config['name_suffix']}_filtered_merged_sorted.bam")
                
                # Create the commands
                if config["picard_cmd"] is None:
                    picard_command = None
                    filtering_command = f"{samtools_path} view -@ {cpus} -b -h -F 2308 -o {os.path.join(output_dir, out_2)} {os.path.join(input_dir, merged_file)}"
                else:
                    picard_command = config["picard_cmd"].format(
                        input=f"I={os.path.join(input_dir, merged_file)}",
                        output=f"O={os.path.join(output_dir, out_1)}",
                        metrics=f"M={os.path.join(output_dir, out_1)}.txt"
                    )
                    filtering_command = f"{samtools_path} view -@ {cpus} -b -h -F 2308 -o {os.path.join(output_dir, out_2)} {os.path.join(output_dir, out_1)}"
                
                # Create job script
                job_script_path = os.path.join(job_dir, f"{name}_duplicate_processing_job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=dup_proc_{name}_job\n")
                    script.write(f"#SBATCH --output={os.path.join(job_dir, f'dup_proc_{name}_output')}\n")
                    script.write(f"#SBATCH --partition={partition}\n")
                    script.write("#SBATCH --nodes=1\n")
                    script.write("#SBATCH --ntasks=1\n")
                    script.write(f"#SBATCH --cpus-per-task={cpus}\n")
                    script.write(f"#SBATCH --time={time}\n")
                    script.write(f"#SBATCH --mail-user={email}\n")
                    script.write("#SBATCH --mail-type=END,FAIL\n")
                    script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                    script.write("\n")
                    script.write("# Load required modules\n")
                    script.write(f"export PATH=$(dirname {samtools_path}):$PATH\n")
                    script.write("\n")
                    script.write(f"cd {input_dir}\n")
                    script.write("\n")
                    
                    if picard_command:
                        script.write(f"{picard_command}\n\n")
                    
                    script.write(f"{filtering_command}\n")
                
                os.chmod(job_script_path, 0o755)
                logging.info(f"Created duplicate processing job script: {job_script_path}")
                
        except Exception as e:
            logging.error(f"Error processing directory {input_dir}: {e}")

def generate_indexing_jobs(input_dirs, job_dirs, samtools_path="samtools", picard_path="picard",
                         partition="sixhour", time="16:00:00", 
                         email="l338m483@ku.edu", mem_per_cpu="5g", cpus=10):
    """
    Generate SLURM job scripts to create indexes for BAM files.
    
    Parameters:
    input_dirs (list): List of directories containing BAM files to index
    job_dirs (list): List of directories for job scripts
    samtools_path (str): Path to samtools executable
    picard_path (str): Path to picard executable
    partition (str): SLURM partition to use
    time (str): Time limit for jobs
    email (str): Email for notifications
    mem_per_cpu (str): Memory per CPU
    cpus (int): Number of CPUs per task
    """
    for input_dir, job_dir in zip(input_dirs, job_dirs):
        create_directory(job_dir)
        
        try:
            # Find BAM files to index
            bam_files = [f for f in os.listdir(input_dir) if f.endswith('.bam')]
            
            if not bam_files:
                logging.warning(f"No BAM files found in {input_dir}")
                continue
            
            for bam_file in bam_files:
                # Extract sample name for read groups
                if "12A" in bam_file:
                    year_ID = "2012A"
                    sample_name = re.sub(r'_.*$', '', bam_file)
                elif "12B" in bam_file:
                    year_ID = "2012B"
                    sample_name = re.sub(r'_.*$', '', bam_file)
                else:
                    year_ID = "unknown"
                    sample_name = re.sub(r'_.*$', '', bam_file)
                
                # Create commands
                picard_command = (f"{picard_path} AddOrReplaceReadGroups " 
                                f"I={os.path.join(input_dir, bam_file)} "
                                f"O={os.path.join(input_dir, bam_file.replace('.bam', '.wRG.bam'))} "
                                f"RGID={sample_name} RGLB={year_ID} RGPL=illumina RGPU=unit1 RGSM={sample_name}")
                
                mv_command = f"mv {os.path.join(input_dir, bam_file.replace('.bam', '.wRG.bam'))} {os.path.join(input_dir, bam_file)}"
                index_command = f"{samtools_path} index -b {os.path.join(input_dir, bam_file)}"
                
                # Create job script
                job_script_path = os.path.join(job_dir, f"indexing_{sample_name}_job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=index_{sample_name}_job\n")
                    script.write(f"#SBATCH --output={os.path.join(job_dir, f'index_{sample_name}_output')}\n")
                    script.write(f"#SBATCH --partition={partition}\n")
                    script.write("#SBATCH --nodes=1\n")
                    script.write("#SBATCH --ntasks=1\n")
                    script.write(f"#SBATCH --cpus-per-task={cpus}\n")
                    script.write(f"#SBATCH --time={time}\n")
                    script.write(f"#SBATCH --mail-user={email}\n")
                    script.write("#SBATCH --mail-type=END,FAIL\n")
                    script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                    script.write("\n")
                    script.write("# Load required modules\n")
                    script.write(f"export PATH=$(dirname {samtools_path}):$PATH\n")
                    script.write("\n")
                    script.write(f"cd {input_dir}\n")
                    script.write("\n")
                    script.write(f"{picard_command}\n")
                    script.write(f"rm {os.path.join(input_dir, bam_file)}\n")
                    script.write(f"{mv_command}\n")
                    script.write("\n")
                    script.write(f"{index_command}\n")
                
                os.chmod(job_script_path, 0o755)
                logging.info(f"Created indexing job script: {job_script_path}")
                
        except Exception as e:
            logging.error(f"Error processing directory {input_dir}: {e}")