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

def create_fasta_index(reference_fasta, samtools_path):
    """
    Create a FASTA index file (.fai) for the given reference FASTA file.
    
    Parameters:
    reference_fasta (str): Path to the reference FASTA file
    """
    if not os.path.exists(reference_fasta):
        logging.error(f"Reference FASTA file does not exist: {reference_fasta}")
        return
    
    # Use samtools to create the index
    os.system(f"{samtools_path} faidx {reference_fasta}")
    logging.info(f"Created index for reference FASTA: {reference_fasta}.fai")

def generate_variant_calling_jobs(input_dirs, output_dirs, job_dirs, 
                               regions, reference_fasta, fai_path,
                               window_size=1000000, variant_caller="freebayes",
                               caller_path=None, partition="sixhour", time="6:00:00",
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
    list: List of generated job script paths
    """
    # Default paths for variant callers
    default_paths = {
        "freebayes": "~/.conda/envs/Python2.7/bin/freebayes",
        "bcftools": "bcftools",
        "gatk": "gatk"
    }
    
    # Default parameters for different callers
    default_params = {
        "freebayes": "-4 --limit-coverage=5000",
        "bcftools": "call -mv",
        "gatk": "HaplotypeCaller --emit-ref-confidence GVCF"
    }
    
    # Use default parameters if none specified
    if caller_params is None:
        caller_params = default_params.get(variant_caller, "")
    
    # Use default path if none specified
    if caller_path is None:
        caller_path = default_paths.get(variant_caller, variant_caller)
    
    # Parse reference index
    chrom_lengths = parse_reference_index(fai_path)
    
    # If "all" is in regions, use all chromosomes from the FAI file
    if "all" in regions:
        regions = list(chrom_lengths.keys())
        logging.info(f"Processing all {len(regions)} sequences from reference genome")
    
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
                        call_cmd = f"{caller_path} mpileup -Ou -r {interval} -f {reference_fasta} {bam_path} | {caller_path} call -mv -o {output_path}"
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
                        
                        # Add relevant module loading if needed and not using full paths
                        if not os.path.isabs(caller_path):
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

def extract_chromosome_position(filename):
    """
    Extract chromosome and position information from VCF filename.
    
    Parameters:
    filename (str): VCF filename to parse
    
    Returns:
    tuple: (chromosome_number, start_position, end_position)
    """
    parts = filename.split("_")
    # Find the part that contains chromosome info
    for i, part in enumerate(parts):
        if ":" in part and part.split(":")[0].replace("Chr_", "").isdigit():
            chrom_part = part
            # Try to convert chromosome number to int for proper sorting
            try:
                chrom = int(chrom_part.split(":")[0].replace("Chr_", ""))
            except ValueError:
                chrom = chrom_part.split(":")[0]  # Keep as string if not numeric
                
            pos = chrom_part.split(":")[1].split("-")
            start = int(pos[0])
            end = int(pos[1]) if len(pos) > 1 else start
            return (chrom, start, end)
    
    # Fallback if no position information found
    return (0, 0, 0)

def group_vcf_files_by_sample(vcf_files):
    """
    Group VCF files by sample name.
    
    Parameters:
    vcf_files (list): List of VCF filenames
    
    Returns:
    dict: Dictionary mapping sample names to lists of VCF files
    """
    sample_groups = {}
    
    for vcf_file in vcf_files:
        # Extract sample name - focus on the part before region information
        # Example: "1192c_pwmd_Chr_01:1-1000000_freebayes_variants.vcf"
        parts = vcf_file.split("_")
        
        # Extract the sample name (usually first part or a combination)
        sample_match = re.match(r'^(\d+[a-zA-Z]?)_', vcf_file)
        if sample_match:
            sample_name = sample_match.group(1)
        else:
            sample_name = parts[0]
        
        # Add aligner/processing info to create unique group key
        for i, part in enumerate(parts):
            if part in ["bwa", "stampy", "minimap2", "pwd", "pwmd", "pwod"]:
                sample_name = f"{sample_name}_{part}"
                break
        
        # Create variant caller signature
        for vc in ["freebayes", "bcftools", "gatk"]:
            if vc in vcf_file:
                sample_name = f"{sample_name}_{vc}"
                break
        
        # Add to proper group
        if sample_name not in sample_groups:
            sample_groups[sample_name] = []
        
        sample_groups[sample_name].append(vcf_file)
    
    return sample_groups

def generate_vcf_merge_jobs(input_dirs, output_dirs, job_dirs,
                         bcftools_path="bcftools", threads=10,
                         partition="sixhour", time="6:00:00",
                         email="l338m483@ku.edu", mem_per_cpu="5g"):
    """
    Generate SLURM jobs to merge VCF files by sample.
    
    Parameters:
    input_dirs (list): List of directories containing VCF files
    output_dirs (list): List of directories for merged output
    job_dirs (list): List of directories for job scripts
    bcftools_path (str): Path to bcftools executable
    threads (int): Number of threads for bcftools concat
    partition (str): SLURM partition
    time (str): Job time limit
    email (str): Notification email
    mem_per_cpu (str): Memory per CPU
    
    Returns:
    list: List of generated job scripts
    """
    job_scripts = []
    
    for input_dir, output_dir, job_dir in zip(input_dirs, output_dirs, job_dirs):
        # Create job directory if it doesn't exist
        create_directory(job_dir)
        create_directory(output_dir)
        
        # Get all VCF files from directory
        vcf_files = []
        try:
            for filename in os.listdir(input_dir):
                if filename.endswith(".vcf") or "variants" in filename:
                    vcf_files.append(filename)
        except FileNotFoundError:
            logging.warning(f"Directory not found: {input_dir}")
            continue
        
        if not vcf_files:
            logging.warning(f"No VCF files found in {input_dir}")
            continue
        
        # Group VCF files by sample
        sample_groups = group_vcf_files_by_sample(vcf_files)
        logging.info(f"Found {len(sample_groups)} sample groups to merge in {input_dir}")
        
        # Process each sample group
        for sample_name, files in sample_groups.items():
            # Sort files by chromosome and position
            sorted_files = sorted(files, key=extract_chromosome_position)
            
            # Define output filename
            merged_vcf = f"{sample_name}_merged.vcf"
            merged_vcf_path = os.path.join(output_dir, merged_vcf)
            
            # Create merge command
            file_paths = ' '.join([os.path.join(input_dir, f) for f in sorted_files])
            merge_command = f"{bcftools_path} concat --threads={threads} -Ov -o {merged_vcf_path} {file_paths}"
            
            # Create job script path
            job_script_path = os.path.join(job_dir, f"{sample_name}_vcf_merge_job.sh")
            
            # Write job script
            with open(job_script_path, "w") as script:
                script.write("#!/bin/bash\n")
                script.write(f"#SBATCH --job-name=merge_vcfs_{sample_name}_job\n")
                script.write(f"#SBATCH --output={os.path.join(job_dir, f'merge_vcfs_{sample_name}_output')}\n")
                script.write(f"#SBATCH --partition={partition}\n")
                script.write("#SBATCH --nodes=1\n")
                script.write("#SBATCH --ntasks=1\n")
                script.write(f"#SBATCH --cpus-per-task={threads}\n")
                script.write(f"#SBATCH --time={time}\n")
                script.write(f"#SBATCH --mail-user={email}\n")
                script.write("#SBATCH --mail-type=FAIL\n")
                script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                script.write("\n")
                
                # Add module loading if bcftools path is not absolute
                if not os.path.isabs(bcftools_path):
                    script.write("module load bcftools\n")
                
                script.write("\n")
                script.write(f"cd {input_dir}\n")
                script.write("\n")
                script.write(f"{merge_command}\n")
            
            # Make script executable
            os.chmod(job_script_path, 0o755)
            job_scripts.append(job_script_path)
            logging.info(f"Created VCF merge job script for {sample_name}: {job_script_path}")
    
    return job_scripts