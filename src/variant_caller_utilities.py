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

def validate_bam_file(bam_path):
    """
    Validate that a BAM file exists and is not empty.
    
    Parameters:
    bam_path (str): Path to BAM file
    
    Returns:
    bool: True if valid, False otherwise
    """
    if not os.path.exists(bam_path):
        logging.error(f"BAM file does not exist: {bam_path}")
        return False
    
    if os.path.getsize(bam_path) == 0:
        logging.error(f"BAM file is empty: {bam_path}")
        return False
    
    # Check for BAM index
    bai_path = bam_path + ".bai"
    if not os.path.exists(bai_path):
        logging.warning(f"BAM index not found: {bai_path}. Variant calling may be slower.")
    
    return True

def validate_vcf_files(vcf_files):
    """
    Validate that VCF files exist and are not empty.
    
    Parameters:
    vcf_files (list): List of VCF file paths
    
    Returns:
    tuple: (valid_files, invalid_files)
    """
    valid_files = []
    invalid_files = []
    
    for vcf_file in vcf_files:
        if not os.path.exists(vcf_file):
            logging.warning(f"VCF file does not exist: {vcf_file}")
            invalid_files.append(vcf_file)
        elif os.path.getsize(vcf_file) == 0:
            logging.warning(f"VCF file is empty: {vcf_file}")
            invalid_files.append(vcf_file)
        else:
            valid_files.append(vcf_file)
    
    return valid_files, invalid_files

def check_reference_compatibility(reference_fasta, bam_file):
    """
    Check if BAM file contigs are compatible with reference FASTA.
    Requires samtools to be available.
    
    Parameters:
    reference_fasta (str): Path to reference FASTA
    bam_file (str): Path to BAM file
    
    Returns:
    bool: True if compatible or cannot verify, False if incompatible
    """
    try:
        import subprocess
        
        # Get BAM header contigs
        result = subprocess.run(
            ['samtools', 'view', '-H', bam_file],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode != 0:
            logging.warning(f"Could not read BAM header for {bam_file}")
            return True  # Assume compatible if we can't check
        
        # Parse reference names from @SQ lines
        bam_contigs = set()
        for line in result.stdout.split('\n'):
            if line.startswith('@SQ'):
                for part in line.split('\t'):
                    if part.startswith('SN:'):
                        bam_contigs.add(part[3:])
        
        # Get reference contigs from .fai file
        fai_path = reference_fasta + ".fai"
        if os.path.exists(fai_path):
            ref_contigs = set()
            with open(fai_path, 'r') as fai:
                for line in fai:
                    ref_contigs.add(line.split('\t')[0])
            
            # Check if there's overlap
            common_contigs = bam_contigs.intersection(ref_contigs)
            if not common_contigs:
                logging.error(f"No common contigs between BAM {bam_file} and reference {reference_fasta}")
                logging.error(f"BAM contigs: {sorted(list(bam_contigs))[:5]}...")
                logging.error(f"Reference contigs: {sorted(list(ref_contigs))[:5]}...")
                return False
            
            logging.debug(f"Found {len(common_contigs)} common contigs between BAM and reference")
        
        return True
        
    except Exception as e:
        logging.warning(f"Could not verify reference compatibility: {e}")
        return True  # Assume compatible if we can't check

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

def extract_sample_names_from_bams(bam_files, suffix="_filtered_merged_sorted.bam"):
    """
    Extract sample names from BAM filenames, preserving library IDs.
    
    This function extracts sample names while preserving year-based library identifiers
    (e.g., "2007A", "2012B") that are used as sample names in VCF files.
    
    Parameters:
    bam_files (list): List of BAM filenames
    suffix (str): Suffix to remove from BAM filenames (or None to try multiple patterns)
    
    Returns:
    list: List of cleaned sample names with library IDs preserved
    """
    sample_names = []
    
    # Common BAM file suffixes to try
    suffixes_to_try = [
        suffix,
        "_filtered_merged_sorted.bam",
        "_dedup.bam",
        "_marked_duplicates.bam",
        "_rg.bam",
        ".bam"
    ]
    
    for bam_file in bam_files:
        base_name = bam_file
        
        # Try to remove suffix
        for suf in suffixes_to_try:
            if base_name.endswith(suf):
                base_name = base_name.replace(suf, "")
                break
        
        # Extract the primary identifier (preserve library IDs like "2007A", "2012B")
        # Split by underscore and take the first part which typically contains the year+ID
        parts = base_name.split('_')
        if parts:
            sample_name = parts[0]
        else:
            sample_name = base_name
        
        # Remove any leading/trailing whitespace
        sample_name = sample_name.strip()
        sample_names.append(sample_name)
    
    return sample_names

def detect_bam_files(input_dir, bam_type=None, suffix=None):
    """
    Find BAM files in a directory, optionally filtering by type.
    
    This function can detect BAM files with various suffixes including:
    - Filtered/merged: _filtered_merged_sorted.bam
    - Deduplicated: _dedup.bam, _marked_duplicates.bam
    - With read groups: _rg.bam
    - Any .bam file
    
    Parameters:
    input_dir (str): Directory to search
    bam_type (str): Type of BAM files ('bwa', 'stampy', None for both)
    suffix (str or list): File suffix(es) to match. If None, tries common patterns.
    
    Returns:
    list: List of matching BAM files
    """
    if not os.path.exists(input_dir):
        logging.warning(f"Directory does not exist: {input_dir}")
        return []
    
    # Define common BAM file suffixes in order of preference
    if suffix is None:
        suffixes = [
            "_dedup.bam",                      # Deduplicated
            "_marked_duplicates.bam",           # Marked duplicates
            "_filtered_merged_sorted.bam",      # Filtered/merged
            "_rg.bam",                          # With read groups
            "_sorted.bam",                      # Sorted
            ".bam"                              # Any BAM file
        ]
    elif isinstance(suffix, list):
        suffixes = suffix
    else:
        suffixes = [suffix]
    
    bam_files = []
    all_files = set(os.listdir(input_dir))
    
    # Try each suffix in order
    for suf in suffixes:
        for filename in all_files:
            # Check if file is a BAM file with the specified suffix
            if filename.endswith(suf) and filename not in bam_files:
                # Filter by type if specified
                if bam_type is None:
                    bam_files.append(filename)
                elif bam_type == 'bwa' and 'bwa' in filename:
                    bam_files.append(filename)
                elif bam_type == 'stampy' and 'stampy' in filename:
                    bam_files.append(filename)
        
        # If we found files with this suffix, stop searching
        if bam_files:
            logging.info(f"Found {len(bam_files)} BAM files with suffix '{suf}' in {input_dir}")
            break
    
    if not bam_files:
        logging.warning(f"No BAM files found in {input_dir} with any recognized suffix")
    
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
                                caller_params=None, paired_input_dirs=None, threads=1, constraint=None):
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
    paired_input_dirs (list): List of paired directories for joint calling (optional)
    threads (int): Number of threads to use
    
    Returns:
    list: List of generated job script paths
    """
    # Default paths for variant callers
    default_paths = {
        "freebayes": "/home/l338m483/.conda/envs/Python2.7/bin/freebayes",
        "bcftools": "/kuhpc/sw/conda/latest/envs/bioconda/bin/bcftools",
        "gatk": "gatk"
    }
    
    # Default parameters for different callers
    default_params = {
        "freebayes": "-4 --max-coverage=5000",
        "bcftools": "-Ou -a FORMAT/AD --max-depth 5000",
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
    
    # Check if we're doing paired analysis
    is_paired = paired_input_dirs is not None and len(paired_input_dirs) == len(input_dirs)
    
    # Process directories
    for i, (input_dir, output_dir, job_dir) in enumerate(zip(input_dirs, output_dirs, job_dirs)):
        # Create output and job directories if they don't exist
        create_directory(output_dir)
        create_directory(job_dir)
        
        # Get BAM files from input directory
        bam_files = detect_bam_files(input_dir)
        if not bam_files:
            logging.warning(f"No BAM files found in {input_dir}")
            continue
        
        # Validate BAM files
        valid_bam_files = []
        for bam_file in bam_files:
            bam_path = os.path.join(input_dir, bam_file)
            if validate_bam_file(bam_path):
                # Check compatibility with reference
                if check_reference_compatibility(reference_fasta, bam_path):
                    valid_bam_files.append(bam_file)
                else:
                    logging.error(f"BAM file {bam_file} is incompatible with reference {reference_fasta}")
            else:
                logging.error(f"BAM file validation failed: {bam_file}")
        
        if not valid_bam_files:
            logging.error(f"No valid BAM files found in {input_dir}")
            continue
        
        bam_files = valid_bam_files
        logging.info(f"Validated {len(bam_files)} BAM files in {input_dir}")
        
        # Get paired BAM files if applicable
        paired_bam_files = []
        if is_paired:
            paired_input_dir = paired_input_dirs[i]
            paired_bam_files = detect_bam_files(paired_input_dir)
            if not paired_bam_files:
                logging.warning(f"No BAM files found in paired directory {paired_input_dir}")
                continue
            
            # Validate paired BAM files
            valid_paired_bam_files = []
            for bam_file in paired_bam_files:
                bam_path = os.path.join(paired_input_dir, bam_file)
                if validate_bam_file(bam_path):
                    if check_reference_compatibility(reference_fasta, bam_path):
                        valid_paired_bam_files.append(bam_file)
                    else:
                        logging.error(f"Paired BAM file {bam_file} is incompatible with reference")
                else:
                    logging.error(f"Paired BAM file validation failed: {bam_file}")
            
            if not valid_paired_bam_files:
                logging.error(f"No valid paired BAM files found in {paired_input_dir}")
                continue
            
            paired_bam_files = valid_paired_bam_files
            logging.info(f"Validated {len(paired_bam_files)} paired BAM files in {paired_input_dir}")
        
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
                if len(bam_files) == 1:
                    # Extract sample name from BAM filename
                    sample_name = bam_files[0].replace("_filtered_merged_sorted.bam", "")
                    output_vcf = f"{sample_name}_{interval.replace(':', '_')}_{variant_caller}_variants.vcf"
                else:
                    output_vcf = f"{interval.replace(':', '_')}_{variant_caller}_variants.vcf"
                
                output_path = os.path.join(output_dir, output_vcf)
                
                # Create variant calling command
                if variant_caller == "freebayes":
                    # For Freebayes
                    bam_paths = ' '.join([os.path.join(input_dir, f) for f in bam_files])
                    if is_paired:
                        # Add paired BAM files
                        paired_bam_paths = ' '.join([os.path.join(paired_input_dirs[i], f) for f in paired_bam_files])
                        bam_paths = f"{bam_paths} {paired_bam_paths}"
                    
                    call_cmd = f"{caller_path} {caller_params} -r {interval} -f {reference_fasta} {bam_paths} > {output_path}"
                
                elif variant_caller == "bcftools":
                    # For bcftools with mpileup+call pipeline
                    if len(bam_files) == 1 and is_paired and len(paired_bam_files) == 1:
                        # Special case for paired calling with one sample per directory
                        bam_path = os.path.join(input_dir, bam_files[0])
                        paired_bam_path = os.path.join(paired_input_dirs[i], paired_bam_files[0])
                        
                        # Extract sample name (preserve library IDs like "2007A")
                        sample_names = extract_sample_names_from_bams([bam_files[0]])
                        sample_name = sample_names[0] if sample_names else "sample"
                        
                        output_vcf = f"{sample_name}_{interval.replace(':', '_')}_bcftools_variants.vcf"
                        output_path = os.path.join(output_dir, output_vcf)
                        
                        # Construct command with mpileup | call pipeline
                        # Note: bcftools will use read group sample names from BAM files
                        call_cmd = (f"{caller_path} mpileup --threads={threads} {caller_params} "
                                f"-r {interval} -f {reference_fasta} {bam_path} {paired_bam_path} | "
                                f"{caller_path} call -mv -Ov --threads={threads} -o {output_path}")
                    else:
                        # Standard case (single directory or multiple samples)
                        bam_paths = ' '.join([os.path.join(input_dir, f) for f in bam_files])
                        all_bam_files = bam_files.copy()
                        
                        if is_paired:
                            # Add paired BAM files
                            paired_bam_paths = ' '.join([os.path.join(paired_input_dirs[i], f) for f in paired_bam_files])
                            bam_paths = f"{bam_paths} {paired_bam_paths}"
                            all_bam_files.extend(paired_bam_files)
                        
                        # Construct command with mpileup | call pipeline
                        # bcftools will automatically use read group sample names from BAM @RG headers
                        call_cmd = (f"{caller_path} mpileup --threads={threads} {caller_params} "
                                    f"-r {interval} -f {reference_fasta} {bam_paths} | "
                                    f"{caller_path} call -mv -Ov --threads={threads} -o {output_path}")
                
                elif variant_caller == "gatk":
                    # For GATK, we need to format the input differently
                    gatk_inputs = ' '.join([f"-I {os.path.join(input_dir, f)}" for f in bam_files])
                    if is_paired:
                        # Add paired BAM files
                        paired_gatk_inputs = ' '.join([f"-I {os.path.join(paired_input_dirs[i], f)}" for f in paired_bam_files])
                        gatk_inputs = f"{gatk_inputs} {paired_gatk_inputs}"
                    
                    # Calculate Java heap memory (80% of SLURM allocation, minimum 8GB)
                    mem_value = float(mem_per_cpu.rstrip('gGmM'))
                    total_mem_gb = mem_value * threads
                    java_mem_gb = max(8, int(total_mem_gb * 0.8))
                    
                    call_cmd = f"{caller_path} --java-options '-Xmx{java_mem_gb}g' {caller_params} -R {reference_fasta} {gatk_inputs} -L {interval} -O {output_path}"
                else:
                    logging.error(f"Unsupported variant caller: {variant_caller}")
                    continue
                
                # Create job script
                job_script_path = os.path.join(job_dir, f"{script_name}_{i}_job.sh")
                
                with open(job_script_path, "w") as script:
                    script.write("#!/bin/bash\n")
                    script.write(f"#SBATCH --job-name={script_name}_{i}_job\n")
                    script.write(f"#SBATCH --output={os.path.join(job_dir, f'{script_name}_{i}_output')}\n")
                    script.write(f"#SBATCH --partition={partition}\n")
                    if constraint:
                        script.write(f"#SBATCH --constraint={constraint}\n")
                    script.write("#SBATCH --nodes=1\n")
                    script.write("#SBATCH --ntasks=1\n")
                    script.write(f"#SBATCH --cpus-per-task={threads}\n")
                    script.write(f"#SBATCH --time={time}\n")
                    script.write(f"#SBATCH --mail-user={email}\n")
                    script.write("#SBATCH --mail-type=FAIL\n")
                    script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                    script.write("\n")
                    
                    # Write the command
                    script.write(f"{call_cmd}\n")
                
                # Make executable
                os.chmod(job_script_path, 0o755)
                job_scripts.append(job_script_path)
                
                if is_paired:
                    logging.info(f"Created joint variant calling job for region {interval} using paired directories: {job_script_path}")
                else:
                    logging.info(f"Created variant calling job for {len(bam_files)} samples in region {interval}: {job_script_path}")
                
                # Update position for next window
                position = end
    
    return job_scripts

def generate_global_variant_calling_jobs(input_dirs, output_dir, job_dir,
                                        regions, reference_fasta, fai_path,
                                        window_size=1000000, variant_caller="bcftools",
                                        caller_path=None, partition="sixhour", time="6:00:00",
                                        email="l338m483@ku.edu", mem_per_cpu="30g",
                                        caller_params=None, threads=1, constraint=None,
                                        output_prefix="global_variants", bam_suffix="_filtered_merged_sorted.bam"):
    """
    Generate SLURM jobs for GLOBAL variant calling - all BAM files from all directories called together.
    
    This creates one VCF per genomic window containing ALL samples from ALL input directories.
    
    Parameters:
    input_dirs (list): List of directories containing BAM files (all will be combined)
    output_dir (str): Single directory for variant calling output
    job_dir (str): Single directory for job scripts
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
    threads (int): Number of threads to use
    constraint (str): SLURM node constraints
    output_prefix (str): Prefix for output VCF files
    bam_suffix (str): BAM file suffix to select (e.g., '_filtered_merged_sorted.bam')
    
    Returns:
    list: List of generated job script paths
    """
    # Default paths for variant callers
    default_paths = {
        "freebayes": "/home/l338m483/.conda/envs/Python2.7/bin/freebayes",
        "bcftools": "/kuhpc/sw/conda/latest/envs/bioconda/bin/bcftools",
        "gatk": "gatk"
    }
    
    # Default parameters for different callers
    default_params = {
        "freebayes": "-4 --max-coverage=5000",
        "bcftools": "-Ou -a FORMAT/AD --max-depth 5000",
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
    
    # Create output and job directories if they don't exist
    create_directory(output_dir)
    create_directory(job_dir)
    
    # Collect ALL BAM files from ALL input directories
    all_bam_files = []
    all_bam_paths = []
    
    logging.info(f"Collecting BAM files from {len(input_dirs)} directories for global variant calling")
    logging.info(f"Filtering for BAM files ending with: {bam_suffix}")
    
    for input_dir in input_dirs:
        bam_files = detect_bam_files(input_dir)
        if not bam_files:
            logging.warning(f"No BAM files found in {input_dir}")
            continue
        
        # Filter BAM files by suffix
        filtered_bam_files = [f for f in bam_files if f.endswith(bam_suffix)]
        
        if not filtered_bam_files:
            logging.warning(f"No BAM files with suffix '{bam_suffix}' found in {input_dir}")
            continue
        
        logging.info(f"Found {len(filtered_bam_files)} BAM files with suffix '{bam_suffix}' in {input_dir}")
        
        # Validate BAM files
        for bam_file in filtered_bam_files:
            bam_path = os.path.join(input_dir, bam_file)
            if validate_bam_file(bam_path):
                # Check compatibility with reference
                if check_reference_compatibility(reference_fasta, bam_path):
                    all_bam_files.append(bam_file)
                    all_bam_paths.append(bam_path)
                else:
                    logging.error(f"BAM file {bam_file} is incompatible with reference {reference_fasta}")
            else:
                logging.error(f"BAM file validation failed: {bam_file}")
    
    if not all_bam_paths:
        logging.error("No valid BAM files found in any input directory")
        return []
    
    logging.info(f"Found {len(all_bam_paths)} valid BAM files across all directories for global calling")
    
    # Track generated job scripts
    job_scripts = []
    
    # Create a single BAM paths string for all variant calling commands
    bam_paths_str = ' '.join(all_bam_paths)
    
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
            output_vcf = f"{output_prefix}_{interval.replace(':', '_')}_{variant_caller}_variants.vcf"
            output_path = os.path.join(output_dir, output_vcf)
            
            # Create variant calling command
            if variant_caller == "freebayes":
                # For Freebayes
                call_cmd = f"{caller_path} {caller_params} -r {interval} -f {reference_fasta} {bam_paths_str} > {output_path}"
            
            elif variant_caller == "bcftools":
                # For bcftools with mpileup+call pipeline
                # bcftools will automatically use read group sample names from BAM @RG headers
                call_cmd = (f"{caller_path} mpileup --threads={threads} {caller_params} "
                            f"-r {interval} -f {reference_fasta} {bam_paths_str} | "
                            f"{caller_path} call -mv -Ov --threads={threads} -o {output_path}")
            
            elif variant_caller == "gatk":
                # For GATK, we need to format the input differently
                gatk_inputs = ' '.join([f"-I {bam_path}" for bam_path in all_bam_paths])
                
                # Calculate Java heap memory (80% of SLURM allocation, minimum 8GB)
                mem_value = float(mem_per_cpu.rstrip('gGmM'))
                total_mem_gb = mem_value * threads
                java_mem_gb = max(8, int(total_mem_gb * 0.8))
                
                call_cmd = f"{caller_path} --java-options '-Xmx{java_mem_gb}g' {caller_params} -R {reference_fasta} {gatk_inputs} -L {interval} -O {output_path}"
            else:
                logging.error(f"Unsupported variant caller: {variant_caller}")
                continue
            
            # Create job script
            job_script_path = os.path.join(job_dir, f"{script_name}_global_job.sh")
            
            with open(job_script_path, "w") as script:
                script.write("#!/bin/bash\n")
                script.write(f"#SBATCH --job-name={script_name}_global_job\n")
                script.write(f"#SBATCH --output={os.path.join(job_dir, f'{script_name}_global_output')}\n")
                script.write(f"#SBATCH --partition={partition}\n")
                if constraint:
                    script.write(f"#SBATCH --constraint={constraint}\n")
                script.write("#SBATCH --nodes=1\n")
                script.write("#SBATCH --ntasks=1\n")
                script.write(f"#SBATCH --cpus-per-task={threads}\n")
                script.write(f"#SBATCH --time={time}\n")
                script.write(f"#SBATCH --mail-user={email}\n")
                script.write("#SBATCH --mail-type=FAIL\n")
                script.write(f"#SBATCH --mem-per-cpu={mem_per_cpu}\n")
                script.write("\n")
                
                # Write the command
                script.write(f"{call_cmd}\n")
            
            # Make executable
            os.chmod(job_script_path, 0o755)
            job_scripts.append(job_script_path)
            
            logging.info(f"Created global variant calling job for {len(all_bam_paths)} samples in region {interval}: {job_script_path}")
            
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
                                  sample_names=None, merge_mode='by_chromosome', constraint=None):
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
        
        # Validate VCF files
        valid_vcf_files, invalid_vcf_files = validate_vcf_files(vcf_files)
        if not valid_vcf_files:
            logging.error(f"No valid VCF files found in {input_dir}")
            continue
        
        if invalid_vcf_files:
            logging.warning(f"Skipping {len(invalid_vcf_files)} invalid VCF files in {input_dir}")
        
        vcf_files = valid_vcf_files
        logging.info(f"Validated {len(vcf_files)} VCF files in {input_dir}")
        
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
                if constraint:
                    script.write(f"#SBATCH --constraint={constraint}\n")
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
                    if constraint:
                        script.write(f"#SBATCH --constraint={constraint}\n")
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