#!/usr/bin/env python3

"""
This program runs the alignment workflow using BWA and Stampy on a high performance cluster.

The program integrates all required steps from BWA alignment through Stampy processing.
It supports both single-end and paired-end reads and various file formats.

Author: Luis Javier Madrigal-Roca
Date: 2025-03-23
"""

import os
import sys
import argparse
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Include the source directory in the path
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))

# Import utilities
from alignment_utilities import *
from slurm_utilities import *

# Python 2.7 path for Stampy
default_python=os.path.expanduser("~/.conda/envs/Python2.7/bin/python")

def main():
    parser = argparse.ArgumentParser(description="Alignment pipeline using BWA and Stampy (optional)")
    
    # Common arguments
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument("--partition", type=str, default="sixhour",
                            help="SLURM partition to use")
    common_parser.add_argument("--time", type=str, default="6:00:00",
                            help="Time limit for jobs")
    common_parser.add_argument("--email", type=str, default="l338m483@ku.edu",
                            help="Email for notifications")
    common_parser.add_argument("--mem-per-cpu", type=str, default="5g",
                            help="Memory per CPU")
    common_parser.add_argument("--cpus", type=int, default=10,
                            help="Number of CPUs per task")
    common_parser.add_argument("--submit", action="store_true",
                            help="Automatically submit jobs after generation")
    common_parser.add_argument("--max-jobs", type=int, default=5000,
                            help="Maximum number of jobs to have in queue at once")
    common_parser.add_argument("--check-interval", type=int, default=4000,
                            help="Time between job status checks in seconds")
    common_parser.add_argument("--max-wait-time", type=int, default=86400,
                            help="Maximum time to wait for jobs in seconds (default: 24 hours)")
    common_parser.add_argument("--max-retries", type=int, default=1,
                            help="Maximum number of retry attempts for failed jobs")
    common_parser.add_argument("--abort-on-failure", action="store_true",
                            help="Abort workflow if jobs fail after retries")
    
    # Subparsers for different commands
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")
    
    # BWA command
    bwa_parser = subparsers.add_parser("bwa", parents=[common_parser], help="Run BWA alignment")
    bwa_parser.add_argument("--input-dirs", type=str, nargs="+", required=True,
                            help="Directories containing preprocessed FASTQ files")
    bwa_parser.add_argument("--output-dirs", type=str, nargs="+", required=True,
                            help="Directories for BWA output (SAM files)")
    bwa_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    bwa_parser.add_argument("--reference", type=str, required=True,
                            help="Path to the reference genome (indexed with BWA)")
    bwa_parser.add_argument("--bwa-path", type=str, default="/kuhpc/sw/conda/latest/envs/bioconda/bin/bwa",
                            help="Path to BWA executable")
    
    # SAM to BAM command
    sambam_parser = subparsers.add_parser("sambam", parents=[common_parser], help="Convert SAM to sorted BAM")
    sambam_parser.add_argument("--input-dirs", type=str, nargs="+", required=True,
                            help="Directories containing BWA output (SAM files)")
    sambam_parser.add_argument("--output-dirs", type=str, nargs="+", required=True,
                            help="Directories for sorted BAM output")
    sambam_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    sambam_parser.add_argument("--samtools-path", type=str, 
                            default="/kuhpc/sw/conda/latest/envs/bioconda/bin/samtools",
                            help="Path to samtools executable")
    
    # Stampy command
    stampy_parser = subparsers.add_parser("stampy", parents=[common_parser], help="Run Stampy on sorted BAM files")
    stampy_parser.add_argument("--input-dirs", type=str, nargs="+", required=True,
                            help="Directories containing sorted BAM files")
    stampy_parser.add_argument("--output-dirs", type=str, nargs="+", required=True,
                            help="Directories for Stampy output")
    stampy_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    stampy_parser.add_argument("--reference", type=str, required=True,
                            help="Path to the reference genome (with .stidx and .sthash)")
    stampy_parser.add_argument("--python_2_7_path", type=str, default=str(default_python),
                            help="Path to Python 2.7 executable. This is required for Stampy")
    stampy_parser.add_argument("--stampy-path", default="/home/l338m483/stampy/stampy.py",
                            type=str,
                            help="Path to Stampy executable")
    stampy_parser.add_argument("--samtools-path", type=str,
                            default="/kuhpc/sw/conda/latest/envs/bioconda/bin/samtools",
                            help="Path to samtools executable")
    stampy_parser.add_argument("--stampy-cpus", type=int, default=3,
                            help="Number of CPUs for Stampy (default: 3)")
    
    # Workflow command (runs all steps)
    workflow_parser = subparsers.add_parser("workflow", parents=[common_parser], 
                                        help="Run complete workflow (BWA, SAM to BAM, Stampy). Notice that if you were to run the pipeline without Stampy, you can just run BWA and then convert SAM to BAM.")
    workflow_parser.add_argument("--input-dirs", type=str, nargs="+", required=True,
                            help="Directories containing preprocessed FASTQ files")
    workflow_parser.add_argument("--bwa-dirs", type=str, nargs="+", required=True,
                            help="Directories for BWA output (SAM files)")
    workflow_parser.add_argument("--bam-dirs", type=str, nargs="+", required=True,
                            help="Directories for sorted BAM output")
    workflow_parser.add_argument("--stampy-dirs", type=str, nargs="+", required=True,
                            help="Directories for Stampy output")
    workflow_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    workflow_parser.add_argument("--reference", type=str, required=True,
                            help="Path to the reference genome")
    workflow_parser.add_argument("--bwa-path", type=str, default="/kuhpc/sw/conda/latest/envs/bioconda/bin/bwa",
                            help="Path to BWA executable")
    workflow_parser.add_argument("--samtools-path", type=str, default="/kuhpc/sw/conda/latest/envs/bioconda/bin/samtools",
                            help="Path to samtools executable")
    # Make sure this is correct in the workflow parser section
    workflow_parser.add_argument("--python_2_7_path", type=str, default=str(default_python),
                            help="Path to Python 2.7 executable. This is required for Stampy")
    workflow_parser.add_argument("--stampy-path", type=str, default="/home/l338m483/stampy/stampy.py", required=False,
                            help="Path to Stampy executable")
    workflow_parser.add_argument("--stampy-cpus", type=int, default=3,
                            help="Number of CPUs for Stampy (default: 3)")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Check if bwa-generated index exists
    if args.command == "bwa" or args.command == "workflow":
        if not os.path.exists(args.reference + ".bwt"):
            logging.warning("WARNING: BWA index not found. The reference is going to be indexed.")
            index_bwa_reference(args.reference, args.bwa_path)

    # Check if indexed reference files exist
    if args.command == "stampy" or args.command == "workflow":
        # Use normpath to ensure consistent path format
        ref_base = os.path.normpath(os.path.splitext(args.reference)[0])
        stidx_file = f"{ref_base}.stidx"
        sthash_file = f"{ref_base}.sthash"
        
        # More robust check with debug info
        if not os.path.exists(stidx_file) or not os.path.exists(sthash_file):
            logging.warning(f"WARNING: Stampy reference files not found at {stidx_file} or {sthash_file}")
            logging.warning("The reference is going to be indexed.")
            index_stampy_reference(args.reference, args.python_2_7_path, args.stampy_path)
        else:
            logging.info(f"Using existing Stampy index files: {stidx_file} and {sthash_file}")

    # Execute appropriate command
    if args.command == "bwa":
        if len(args.input_dirs) != len(args.output_dirs) or len(args.input_dirs) != len(args.job_dirs):
            logging.error("Error: Number of input, output, and job directories must match")
            sys.exit(1)
        
        generate_bwa_jobs(args.input_dirs, args.output_dirs, args.job_dirs, 
                        args.reference, bwa_path=args.bwa_path,
                        partition=args.partition, time=args.time, email=args.email,
                        mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Find and submit all BWA jobs from all directories with retry capability
            logging.info("=== Submitting BWA jobs from all directories ===")
            all_bwa_jobs = []
            
            # Collect all BWA jobs from all directories
            for job_dir in args.job_dirs:
                bwa_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir)
                            if f.startswith("bwa_") and f.endswith("_alig_job.sh")]
                all_bwa_jobs.extend(bwa_jobs)
                logging.info(f"Found {len(bwa_jobs)} BWA jobs in {job_dir}")
            
            # Submit BWA jobs with retry
            logging.info(f"Submitting {len(all_bwa_jobs)} total BWA jobs with retry capability")
            bwa_result = submit_jobs_with_retry(
                job_scripts=all_bwa_jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            if not bwa_result['all_successful']:
                error_msg = f"BWA jobs failed: {len(bwa_result['final_failures'])} jobs failed after {args.max_retries} retry attempts"
                logging.error(error_msg)
                if args.abort_on_failure:
                    logging.error("Aborting due to BWA failures")
                    sys.exit(1)
            else:
                logging.info("All BWA jobs completed successfully")
    
    elif args.command == "sambam":
        if len(args.input_dirs) != len(args.output_dirs) or len(args.input_dirs) != len(args.job_dirs):
            logging.error("Error: Number of input, output, and job directories must match")
            sys.exit(1)
        
        generate_sam_to_bam_jobs(args.input_dirs, args.output_dirs, args.job_dirs,
                                samtools_path=args.samtools_path,
                                partition=args.partition, time=args.time, email=args.email,
                                mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Find and submit all SAM to BAM jobs from all directories with retry capability
            logging.info("=== Submitting SAM to BAM jobs from all directories ===")
            all_bam_jobs = []
            
            # Collect all SAM to BAM jobs from all directories
            for job_dir in args.job_dirs:
                bam_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir)
                            if f.startswith("bam_") and f.endswith("_sort_job.sh")]
                all_bam_jobs.extend(bam_jobs)
                logging.info(f"Found {len(bam_jobs)} SAM to BAM jobs in {job_dir}")
            
            # Submit SAM to BAM jobs with retry
            logging.info(f"Submitting {len(all_bam_jobs)} total SAM to BAM jobs with retry capability")
            bam_result = submit_jobs_with_retry(
                job_scripts=all_bam_jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            if not bam_result['all_successful']:
                error_msg = f"SAM to BAM jobs failed: {len(bam_result['final_failures'])} jobs failed after {args.max_retries} retry attempts"
                logging.error(error_msg)
                if args.abort_on_failure:
                    logging.error("Aborting due to SAM to BAM failures")
                    sys.exit(1)
            else:
                logging.info("All SAM to BAM jobs completed successfully")
    
    elif args.command == "stampy":
        if len(args.input_dirs) != len(args.output_dirs) or len(args.input_dirs) != len(args.job_dirs):
            logging.error("Error: Number of input, output, and job directories must match")
            sys.exit(1)
        
        generate_stampy_jobs(args.input_dirs, args.output_dirs, args.job_dirs,
                    args.reference, args.python_2_7_path, args.stampy_path,
                    args.samtools_path,
                    partition=args.partition, time=args.time, email=args.email,
                    mem_per_cpu=args.mem_per_cpu, cpus=args.stampy_cpus)
        
        if args.submit:
            # Find and submit all Stampy jobs from all directories with retry capability
            logging.info("=== Submitting Stampy jobs from all directories ===")
            all_stampy_jobs = []
            
            # Collect all Stampy jobs from all directories
            for job_dir in args.job_dirs:
                stampy_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir)
                            if f.startswith("stampy_") and f.endswith("_job.sh")]
                all_stampy_jobs.extend(stampy_jobs)
                logging.info(f"Found {len(stampy_jobs)} Stampy jobs in {job_dir}")
            
            # Submit Stampy jobs with retry
            logging.info(f"Submitting {len(all_stampy_jobs)} total Stampy jobs with retry capability")
            stampy_result = submit_jobs_with_retry(
                job_scripts=all_stampy_jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            if not stampy_result['all_successful']:
                error_msg = f"Stampy jobs failed: {len(stampy_result['final_failures'])} jobs failed after {args.max_retries} retry attempts"
                logging.error(error_msg)
                if args.abort_on_failure:
                    logging.error("Aborting due to Stampy failures")
                    sys.exit(1)
            else:
                logging.info("All Stampy jobs completed successfully")
    
    elif args.command == "workflow":
        # Verify directory counts match
        if (len(args.input_dirs) != len(args.bwa_dirs) or
            len(args.input_dirs) != len(args.bam_dirs) or
            len(args.input_dirs) != len(args.stampy_dirs) or
            len(args.input_dirs) != len(args.job_dirs)):
            logging.error("Error: Number of directories must match across all arguments")
            sys.exit(1)
        
        # STEP 1: BWA alignment
        logging.info("=== STEP 1: Generating BWA alignment jobs ===")
        generate_bwa_jobs(args.input_dirs, args.bwa_dirs, args.job_dirs, 
                        args.reference, bwa_path=args.bwa_path,
                        partition=args.partition, time=args.time, email=args.email,
                        mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Submit all BWA jobs from all directories with retry capability
            logging.info("\n=== Submitting BWA jobs from all directories ===")
            all_bwa_jobs = []
            
            # Collect all BWA jobs from all directories
            for job_dir in args.job_dirs:
                bwa_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir)
                        if f.startswith("bwa_") and f.endswith("_alig_job.sh")]
                all_bwa_jobs.extend(bwa_jobs)
                logging.info(f"Found {len(bwa_jobs)} BWA jobs in {job_dir}")
            
            # Submit BWA jobs with retry
            logging.info(f"Submitting {len(all_bwa_jobs)} total BWA jobs with retry capability")
            bwa_result = submit_jobs_with_retry(
                job_scripts=all_bwa_jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            if not bwa_result['all_successful']:
                error_msg = f"BWA step failed: {len(bwa_result['final_failures'])} jobs failed after {args.max_retries} retry attempts"
                logging.error(error_msg)
                if args.abort_on_failure:
                    logging.error("Aborting workflow due to BWA failures")
                    sys.exit(1)
                else:
                    logging.warning("Continuing workflow despite BWA failures - data may be incomplete")
            else:
                logging.info("All BWA jobs completed successfully")
        
        # STEP 2: SAM to BAM conversion
        logging.info("\n=== STEP 2: Generating SAM to BAM conversion jobs ===")
        generate_sam_to_bam_jobs(args.bwa_dirs, args.bam_dirs, args.job_dirs,
                                samtools_path=args.samtools_path,
                                partition=args.partition, time=args.time, email=args.email,
                                mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Submit all SAM to BAM jobs from all directories with retry capability
            logging.info("\n=== Submitting SAM to BAM conversion jobs from all directories ===")
            all_bam_jobs = []
            
            # Collect all SAM to BAM jobs from all directories
            for job_dir in args.job_dirs:
                bam_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir)
                        if f.startswith("bam_") and f.endswith("_sort_job.sh")]
                all_bam_jobs.extend(bam_jobs)
                logging.info(f"Found {len(bam_jobs)} SAM to BAM jobs in {job_dir}")
            
            # Submit SAM to BAM jobs with retry
            logging.info(f"Submitting {len(all_bam_jobs)} total SAM to BAM jobs with retry capability")
            bam_result = submit_jobs_with_retry(
                job_scripts=all_bam_jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            if not bam_result['all_successful']:
                error_msg = f"SAM to BAM step failed: {len(bam_result['final_failures'])} jobs failed after {args.max_retries} retry attempts"
                logging.error(error_msg)
                if args.abort_on_failure:
                    logging.error("Aborting workflow due to SAM to BAM failures")
                    sys.exit(1)
                else:
                    logging.warning("Continuing workflow despite SAM to BAM failures - data may be incomplete")
            else:
                logging.info("All SAM to BAM jobs completed successfully")
        
        # STEP 3: Stampy processing
        logging.info("\n=== STEP 3: Generating Stampy jobs ===")
        generate_stampy_jobs(args.bam_dirs, args.stampy_dirs, args.job_dirs,
                        args.reference, args.python_2_7_path, args.stampy_path,
                        args.samtools_path,
                        partition=args.partition, time=args.time, email=args.email,
                        mem_per_cpu=args.mem_per_cpu, cpus=args.stampy_cpus)
        
        if args.submit:
            # Submit all Stampy jobs from all directories with retry capability
            logging.info("\n=== Submitting Stampy jobs from all directories ===")
            all_stampy_jobs = []
            
            # Collect all Stampy jobs from all directories
            for job_dir in args.job_dirs:
                stampy_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir)
                            if f.startswith("stampy_") and f.endswith("_job.sh")]
                all_stampy_jobs.extend(stampy_jobs)
                logging.info(f"Found {len(stampy_jobs)} Stampy jobs in {job_dir}")
            
            # Submit Stampy jobs with retry
            logging.info(f"Submitting {len(all_stampy_jobs)} total Stampy jobs with retry capability")
            stampy_result = submit_jobs_with_retry(
                job_scripts=all_stampy_jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            if not stampy_result['all_successful']:
                error_msg = f"Stampy step failed: {len(stampy_result['final_failures'])} jobs failed after {args.max_retries} retry attempts"
                logging.error(error_msg)
                if args.abort_on_failure:
                    logging.error("Aborting workflow due to Stampy failures")
                    sys.exit(1)
                else:
                    logging.warning("Workflow completed despite Stampy failures - data may be incomplete")
            else:
                logging.info("All Stampy jobs completed successfully")
            
            logging.info("\nAlignment workflow completed successfully!")
        else:
            logging.info("\nWorkflow job generation complete. Submit jobs in order:")
            logging.info("1. Run BWA jobs")
            logging.info("2. Run SAM to BAM conversion jobs")
            logging.info("3. Run Stampy jobs")
    
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()