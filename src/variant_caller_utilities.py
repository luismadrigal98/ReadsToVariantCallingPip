#!/usr/bin/env python3

import os
import logging
import re
from pathlib import Path

def create_directory(directory):
    """Create directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)
        logging.info(f"Created directory: {directory}")

def parse_reference_index(fai_path):
    """
    Parse a FASTA index file to get chromosome lengths.
    
    Parameters:
    fai_path (str): Path to the .fai file
    
    Returns:
    dict: Dictionary mapping chromosome names to their lengths
    """
    chrom_lengths = {}
    with open(fai_path, "r") as fai_file:
        for line in fai_file:
            cols = line.strip().split("\t")
            chrom_lengths[cols[0]] = int(cols[1])
    return chrom_lengths

def detect_bam_files(input_dir, bam_type=None, suffix="_filtered_merged_sorted.bam"):
    """
    Find BAM files in a directory, optionally filtering by type.
    
    Parameters:
    input_dir (str): Directory to search
    bam_type (str): Type of BAM files ('bwa', 'stampy', None for both)
    suffix (str): File suffix to match
    
    Returns:
    list: List of matching BAM files
    """
    if not os.path.exists(input_dir):
        logging.warning(f"Directory does not exist: {input_dir}")
        return []
    
    bam_files = []
    for filename in os.listdir(input_dir):
        # Check if file is a BAM file with the specified suffix
        if filename.endswith(suffix):
            # Filter by type if specified
            if bam_type is None:
                bam_files.append(filename)
            elif bam_type == 'bwa' and 'bwa' in filename:
                bam_files.append(filename)
            elif bam_type == 'stampy' and 'stampy' in filename:
                bam_files.append(filename)
    
    return sorted(bam_files)

def generate_variant_calling_jobs(input_dirs, output_dirs, job_dirs, 
                                regions, reference_fasta, fai_path,
                                window_size=1000000, variant_caller="freebayes",
                                caller_path='~/.conda/envs/Python2.7/bin/freebayes',
                                partition="sixhour", time="6:00:00",
                                email="l338m483@ku.edu", mem_per_cpu="30g",
                                caller_params=None):
    """
    Generate SLURM jobs for variant calling using a specified caller.
    
    Parameters:
    input_dirs (list): List of directories containing BAM files
    output_dirs (list): List of directories for variant calling output
    job_dirs (list): List of directories for job scripts
    regions (list): List of regions/chromosomes to process
    reference_fasta (str): Path to reference FASTA file
    fai_path (str): Path to FASTA index file
    window_size (int): Window size for chunking variant calling
    variant_caller (str): Variant caller to use ('freebayes', 'bcftools', etc.)
    caller_path (str): Path to variant caller executable
    partition (str): SLURM partition
    time (str): Job time limit
    email (str): Notification email
    mem_per_cpu (str): Memory per CPU
    caller_params (str): Additional parameters for the variant caller
    
    Returns:
    list: List of generated job scripts
    """
    # Default parameters for different callers
    default_params = {
        "freebayes": "-4 --limit-coverage=5000",
        "bcftools": "call -mv",
        "gatk": "HaplotypeCaller --emit-ref-confidence GVCF"
    }
    
    # Use default parameters if none specified
    if caller_params is None:
        caller_params = default_params.get(variant_caller, "")
    
    # Parse reference index
    chrom_lengths = parse_reference_index(fai_path)
    
    # Track generated job scripts
    job_scripts = []
    
    for input_dir, output_dir, job_dir in zip(input_dirs, output_dirs, job_dirs):
        # Create output and job directories if they don't exist
        create_directory(output_dir)
        create_directory(job_dir)
        
        # Get BAM files from input directory
        bam_files = detect_bam_files(input_dir)
        if not bam_files:
            logging.warning(f"No BAM files found in {input_dir}")
            continue
        
        # Process each region of interest
        for region in regions:
            if region not in chrom_lengths:
                logging.warning(f"Region {region} not found in reference index")
                continue
                
            region_length = chrom_lengths[region]
            position = 0
            
            # Generate jobs for windows across the region
            while position < region_length:
                start = position + 1
                end = position + window_size
                
                if end > region_length:
                    end = region_length
                
                interval = f"{region}:{start}-{end}"
                script_name = f"{variant_caller}_{interval.replace(':', '_')}"
                
                # Process each BAM file
                for bam_file in bam_files:
                    # Extract sample name from BAM file
                    sample_match = re.match(r'^(\d+[a-zA-Z]?)_', bam_file)
                    if sample_match:
                        sample_name = sample_match.group(1)
                    else:
                        sample_name = bam_file.split('_')[0]
                    
                    # Define output filename
                    output_vcf = f"{sample_name}_{interval.replace(':', '_')}_{variant_caller}_variants.vcf"
                    output_path = os.path.join(output_dir, output_vcf)
                    
                    # Create variant calling command
                    bam_path = os.path.join(input_dir, bam_file)
                    
                    if variant_caller == "freebayes":
                        call_cmd = f"{caller_path} {caller_params} -r {interval} -f {reference_fasta} {bam_path} > {output_path}"
                    elif variant_caller == "bcftools":
                        call_cmd = f"{caller_path} mpileup -Ou -r {interval} -f {reference_fasta} {bam_path} | bcftools {caller_params} -o {output_path}"
                    elif variant_caller == "gatk":
                        call_cmd = f"{caller_path} {caller_params} -R {reference_fasta} -I {bam_path} -L {interval} -O {output_path}"
                    else:
                        logging.error(f"Unsupported variant caller: {variant_caller}")
                        continue
                    
                    # Create job script
                    job_script_path = os.path.join(job_dir, f"{sample_name}_{script_name}_job.sh")
                    
                    with open(job_script_path, "w") as script:
                        script.write("#!/bin/bash\n")
                        script.write(f"#SBATCH --job-name={script_name}_{sample_name}_job\n")
                        script.write(f"#SBATCH --output={os.path.join(job_dir, f'{script_name}_{sample_name}_output')}\n")
                        script.write(f"#SBATCH --partition={partition}\n")
                        script.write("#SBATCH --nodes=1\n")
                        script.write("#SBATCH --ntasks=1\n")
                        script.write("#SBATCH --cpus-per-task=1\n")
                        script.write(f"#SBATCH --time={time}\n")
                        script.write(f"#SBATCH --mail-user={email}\n")
                        script.write("#SBATCH --mail-type=FAIL\n")
                        script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                        script.write("\n")
                        
                        # Load module based on variant caller
                        if variant_caller == "freebayes":
                            script.write("module load freebayes\n")
                        elif variant_caller == "bcftools":
                            script.write("module load bcftools\n")
                        elif variant_caller == "gatk":
                            script.write("module load gatk\n")
                        
                        script.write("\n")
                        script.write(f"{call_cmd}\n")
                    
                    # Make executable
                    os.chmod(job_script_path, 0o755)
                    job_scripts.append(job_script_path)
                    logging.info(f"Created variant calling job: {job_script_path}")
                
                # Update position for next window
                position = end
    
    return job_scripts