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
        "freebayes": "/home/l338m483/.conda/envs/Python2.7/bin/freebayes",
        "bcftools": "bcftools",
        "gatk": "gatk"
    }
    
    # Default parameters for different callers
    default_params = {
        "freebayes": "-4 --max-coverage=5000",
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
    
    # If "all" in regions, use all chromosomes
    if "all" in regions:
        regions = list(chrom_lengths.keys())
        logging.info(f"Processing all {len(regions)} sequences from reference genome")
    
    # Track generated job scripts
    job_scripts = []
    
    for input_dir, output_dir, job_dir in zip(input_dirs, output_dirs, job_dirs):
        # Create output and job directories if they don't exist
        create_directory(output_dir)
        create_directory(job_dir)
        
        # Get ALL BAM files from input directory
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
                
                # Define output filename for the region
                output_vcf = f"{interval.replace(':', '_')}_{variant_caller}_variants.vcf"
                output_path = os.path.join(output_dir, output_vcf)
                
                # Create variant calling command with ALL BAM files
                bam_paths = ' '.join([os.path.join(input_dir, f) for f in bam_files])
                
                if variant_caller == "freebayes":
                    call_cmd = f"{caller_path} {caller_params} -r {interval} -f {reference_fasta} {bam_paths} > {output_path}"
                elif variant_caller == "bcftools":
                    call_cmd = f"{caller_path} mpileup -Ou -r {interval} -f {reference_fasta} {bam_paths} | {caller_path} call -mv -o {output_path}"
                elif variant_caller == "gatk":
                    # For GATK, we need to format the input differently
                    gatk_inputs = ' '.join([f"-I {os.path.join(input_dir, f)}" for f in bam_files])
                    call_cmd = f"{caller_path} {caller_params} -R {reference_fasta} {gatk_inputs} -L {interval} -O {output_path}"
                else:
                    logging.error(f"Unsupported variant caller: {variant_caller}")
                    continue
                
                # Create job script
                job_script_path = os.path.join(job_dir, f"{script_name}_job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name={script_name}_job\n")
                    script.write(f"#SBATCH --output={os.path.join(job_dir, f'{script_name}_output')}\n")
                    script.write(f"#SBATCH --partition={partition}\n")
                    script.write("#SBATCH --nodes=1\n")
                    script.write("#SBATCH --ntasks=1\n")
                    script.write("#SBATCH --cpus-per-task=1\n")
                    script.write(f"#SBATCH --time={time}\n")
                    script.write(f"#SBATCH --mail-user={email}\n")
                    script.write("#SBATCH --mail-type=FAIL\n")
                    script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                    script.write("\n")
                    
                    # Add module loading if needed
                    if not os.path.isabs(caller_path) and caller_path != "bcftools" and caller_path != "gatk":
                        if variant_caller == "freebayes":
                            script.write("module load freebayes\n")
                        elif variant_caller == "bcftools":
                            script.write("module load bcftools\n")
                        elif variant_caller == "gatk":
                            script.write("module load gatk\n")
                        script.write("\n")
                    
                    # Write the command
                    script.write(f"{call_cmd}\n")
                
                # Make executable
                os.chmod(job_script_path, 0o755)
                job_scripts.append(job_script_path)
                logging.info(f"Created joint variant calling job for {len(bam_files)} samples in region {interval}: {job_script_path}")
                
                # Update position for next window
                position = end
    
    return job_scripts

def sort_vcf_files_by_region(vcf_files):
    """
    Sort VCF files by chromosome and position.
    
    Parameters:
    vcf_files (list): List of VCF file paths
    
    Returns:
    list: Sorted list of VCF file paths
    """
    def extract_region_info(filepath):
        # Extract chromosome and position from filename
        filename = os.path.basename(filepath)
        parts = filename.split('_')
        
        # Try to find parts with Chr_ or : in them
        chrom_parts = [p for p in parts if p.startswith('Chr_') or ':' in p]
        if chrom_parts:
            region_part = chrom_parts[0]
            # Handle formats like "Chr_10:8000001-9000000"
            if ':' in region_part:
                chrom_str, pos_range = region_part.split(':', 1)
                chrom = chrom_str.replace('Chr_', '')
                start, end = map(int, pos_range.split('-'))
                return (chrom, start, end)
        
        # Fallback - just return filename for lexicographical sorting
        return (filename, 0, 0)
    
    return sorted(vcf_files, key=extract_region_info)

def merge_vcf_files_jobs_generator(input_dirs, output_dirs, job_dirs,
                                  bcftools_path='/kuhpc/sw/conda/latest/envs/bioconda/bin/bcftools',
                                  partition="sixhour", time="6:00:00",
                                  email="l338m483@ku.edu", mem_per_cpu="30g",
                                  threads=1, merge_command='concat',
                                  sample_names=None, merge_mode='by_chromosome'):
    """
    Generate jobs to merge VCF files globally or by chromosome.
    
    Parameters:
    input_dirs (list): List of directories containing VCF files
    output_dirs (list): List of directories for merged output
    job_dirs (list): List of directories for job scripts
    bcftools_path (str): Path to bcftools executable
    partition (str): SLURM partition
    time (str): Job time limit
    email (str): Notification email
    mem_per_cpu (str): Memory per CPU
    threads (int): Number of threads to use
    merge_command (str): Merge command to use (default: 'concat')
    sample_names (list): List of sample names to use for output files
    merge_mode (str): 'global' for one output file, 'by_chromosome' for one per chromosome
    
    Returns:
    list: List of generated job script paths
    """
    job_scripts = []
    
    # For each input/output/job directory set
    for idx, (input_dir, output_dir, job_dir) in enumerate(zip(input_dirs, output_dirs, job_dirs)):
        # Create output and job directories if they don't exist
        create_directory(output_dir)
        create_directory(job_dir)
        
        # Get sample name if provided, otherwise use directory name
        sample_name = None
        if sample_names and idx < len(sample_names):
            sample_name = sample_names[idx]
        
        # Find all VCF files in the input directory
        vcf_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.vcf')]
        if not vcf_files:
            logging.warning(f"No VCF files found in {input_dir}")
            continue
        
        if merge_mode == 'global':
            # GLOBAL MODE: Merge all files into one output file
            
            # Sort all VCF files by chromosome and position
            sorted_files = sort_vcf_files_by_region(vcf_files)
            
            # Create output filename
            if sample_name:
                output_filename = f"{sample_name}_merged_all.vcf"
            else:
                output_filename = "merged_all.vcf"
                
            merged_vcf = os.path.join(output_dir, output_filename)
            
            # Create job script path
            job_script_path = os.path.join(job_dir, "merge_all_job.sh")
            
            # Create job script
            with open(job_script_path, "w") as script:
                script.write("#!/bin/bash\n")
                script.write(f"#SBATCH --job-name=merge_{sample_name}_all_job\n")
                script.write(f"#SBATCH --output={os.path.join(job_dir, 'merge_all_output')}\n")
                script.write(f"#SBATCH --partition={partition}\n")
                script.write("#SBATCH --nodes=1\n")
                script.write("#SBATCH --ntasks=1\n")
                script.write(f"#SBATCH --cpus-per-task={threads}\n")
                script.write(f"#SBATCH --time={time}\n")
                script.write(f"#SBATCH --mail-user={email}\n")
                script.write("#SBATCH --mail-type=FAIL\n")
                script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                script.write("\n")
                
                # Add module loading if needed
                if not os.path.isabs(bcftools_path):
                    script.write("module load bcftools\n")
                    script.write("\n")
                
                # Write the command to merge VCF files
                vcf_files_str = ' '.join(sorted_files)
                if merge_command == 'concat':
                    merge_cmd = f"{bcftools_path} concat --threads {threads} -Ov -o {merged_vcf} {vcf_files_str}"
                else:
                    merge_cmd = f"{bcftools_path} merge --threads {threads} -Ov -o {merged_vcf} {vcf_files_str}"
                
                script.write(f"{merge_cmd}\n")
            
            # Make executable
            os.chmod(job_script_path, 0o755)
            
            # Add to job scripts
            job_scripts.append(job_script_path)
            
            logging.info(f"Created job to merge all VCF files into {merged_vcf}")
            
        else:
            # BY_CHROMOSOME MODE: One merged file per chromosome
            
            # Group VCF files by chromosome
            chrom_groups = {}
            for vcf_file in vcf_files:
                # Parse the chromosome from the filename
                # Format: Chr_14_13000001-14000000_freebayes_variants.vcf
                filename = os.path.basename(vcf_file)
                parts = filename.split('_')
                if len(parts) >= 2 and parts[0] == 'Chr':
                    chromosome = parts[1]  # Extract chromosome number
                    if chromosome not in chrom_groups:
                        chrom_groups[chromosome] = []
                    chrom_groups[chromosome].append(vcf_file)
            
            # Sort chromosomes numerically if possible
            chromosomes = sorted(chrom_groups.keys(), key=lambda x: int(x) if x.isdigit() else x)
            
            # Process each chromosome
            for chromosome in chromosomes:
                # Sort files by position
                def extract_position(filepath):
                    filename = os.path.basename(filepath)
                    parts = filename.split('_')
                    if len(parts) >= 3:
                        pos_part = parts[2]
                        if '-' in pos_part:
                            start, _ = pos_part.split('-', 1)
                            return int(start)
                    return 0
                
                sorted_files = sorted(chrom_groups[chromosome], key=extract_position)
                
                # Create output filename
                if sample_name:
                    output_filename = f"{sample_name}_Chr_{chromosome}_merged.vcf"
                else:
                    output_filename = f"Chr_{chromosome}_merged.vcf"
                    
                merged_vcf = os.path.join(output_dir, output_filename)
                
                # Create job script path
                job_script_path = os.path.join(job_dir, f"merge_Chr_{chromosome}_job.sh")
                
                # Create job script
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name=merge_{sample_name}_{chromosome}_job\n")
                    script.write(f"#SBATCH --output={os.path.join(job_dir, f'merge_Chr_{chromosome}_output')}\n")
                    script.write(f"#SBATCH --partition={partition}\n")
                    script.write("#SBATCH --nodes=1\n")
                    script.write("#SBATCH --ntasks=1\n")
                    script.write(f"#SBATCH --cpus-per-task={threads}\n")
                    script.write(f"#SBATCH --time={time}\n")
                    script.write(f"#SBATCH --mail-user={email}\n")
                    script.write("#SBATCH --mail-type=FAIL\n")
                    script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                    script.write("\n")
                    
                    # Add module loading if needed
                    if not os.path.isabs(bcftools_path):
                        script.write("module load bcftools\n")
                        script.write("\n")
                    
                    # Write the command to merge VCF files
                    vcf_files_str = ' '.join(sorted_files)
                    if merge_command == 'concat':
                        merge_cmd = f"{bcftools_path} concat --threads {threads} -Ov -o {merged_vcf} {vcf_files_str}"
                    else:
                        merge_cmd = f"{bcftools_path} merge --threads {threads} -Ov -o {merged_vcf} {vcf_files_str}"
                    
                    script.write(f"{merge_cmd}\n")
                
                # Make executable
                os.chmod(job_script_path, 0o755)
                
                # Add to job scripts
                job_scripts.append(job_script_path)
                
                if sample_name:
                    logging.info(f"Created merge job for chromosome {chromosome} in sample {sample_name}: {job_script_path}")
                else:
                    logging.info(f"Created merge job for chromosome {chromosome}: {job_script_path}")
    
    return job_scripts