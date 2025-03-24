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
    parser = argparse.ArgumentParser(description="Alignment pipeline using BWA and Stampy")
    
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
    common_parser.add_argument("--check-interval", type=int, default=300,
                            help="Time between job status checks in seconds")
    common_parser.add_argument("--max-wait-time", type=int, default=86400,
                            help="Maximum time to wait for jobs in seconds (default: 24 hours)")
    
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
    stampy_parser.add_argument("--python_2.7_path", type=str, default=str(default_python),
                            help="Path to Python 2.7 executable. This is required for Stampy")
    stampy_parser.add_argument("--stampy-path", default="/home/l338m483/stampy/stampy.py",
                            type=str,
                            help="Path to Stampy executable")
    stampy_parser.add_argument("--stampy-cpus", type=int, default=3,
                            help="Number of CPUs for Stampy (default: 3)")
    
    # Workflow command (runs all steps)
    workflow_parser = subparsers.add_parser("workflow", parents=[common_parser], 
                                        help="Run complete workflow (BWA, SAM to BAM, Stampy)")
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
    stampy_parser.add_argument("--python_2.7_path", type=str, default=str(default_python),
                            help="Path to Python 2.7 executable. This is required for Stampy")
    workflow_parser.add_argument("--stampy-path", type=str, default="/home/l338m483/stampy/stampy.py", required=False,
                            help="Path to Stampy executable")
    workflow_parser.add_argument("--stampy-cpus", type=int, default=3,
                            help="Number of CPUs for Stampy (default: 3)")
    
    # Parse arguments
    args = parser.parse_args()
    
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
            # Submit BWA jobs
            for job_dir in args.job_dirs:
                bwa_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir)
                            if f.startswith("bwa_") and f.endswith("_alig_job.sh")]
                
                bwa_job_ids = submit_jobs_with_limit(bwa_jobs, args.max_jobs)
                logging.info(f"Submitted {len(bwa_job_ids)} BWA jobs from {job_dir}")
                
                # Wait for BWA jobs to complete
                logging.info("Waiting for BWA jobs to complete...")
                wait_for_jobs_to_complete(job_ids=bwa_job_ids, 
                                        check_interval=args.check_interval,
                                        max_wait_time=args.max_wait_time)
    
    elif args.command == "sambam":
        if len(args.input_dirs) != len(args.output_dirs) or len(args.input_dirs) != len(args.job_dirs):
            logging.error("Error: Number of input, output, and job directories must match")
            sys.exit(1)
        
        generate_sam_to_bam_jobs(args.input_dirs, args.output_dirs, args.job_dirs,
                                samtools_path=args.samtools_path,
                                partition=args.partition, time=args.time, email=args.email,
                                mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Submit SAM to BAM jobs
            for job_dir in args.job_dirs:
                bam_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir)
                            if f.startswith("bam_") and f.endswith("_sort_job.sh")]
                
                bam_job_ids = submit_jobs_with_limit(bam_jobs, args.max_jobs)
                logging.info(f"Submitted {len(bam_job_ids)} SAM to BAM jobs from {job_dir}")
                
                # Wait for SAM to BAM jobs to complete
                logging.info("Waiting for SAM to BAM jobs to complete...")
                wait_for_jobs_to_complete(job_ids=bam_job_ids,
                                        check_interval=args.check_interval,
                                        max_wait_time=args.max_wait_time)
    
    elif args.command == "stampy":
        if len(args.input_dirs) != len(args.output_dirs) or len(args.input_dirs) != len(args.job_dirs):
            logging.error("Error: Number of input, output, and job directories must match")
            sys.exit(1)
        
        generate_stampy_jobs(args.input_dirs, args.output_dirs, args.job_dirs,
                    args.reference, args.stampy_path, args.python_2_7_path,
                    partition=args.partition, time=args.time, email=args.email,
                    mem_per_cpu=args.mem_per_cpu, cpus=args.stampy_cpus)
        
        if args.submit:
            # Submit Stampy jobs
            for job_dir in args.job_dirs:
                stampy_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir)
                            if f.startswith("stampy_") and f.endswith("_job.sh")]
                
                stampy_job_ids = submit_jobs_with_limit(stampy_jobs, args.max_jobs)
                logging.info(f"Submitted {len(stampy_job_ids)} Stampy jobs from {job_dir}")
                
                # Wait for Stampy jobs to complete
                logging.info("Waiting for Stampy jobs to complete...")
                wait_for_jobs_to_complete(job_ids=stampy_job_ids,
                                    check_interval=args.check_interval,
                                    max_wait_time=args.max_wait_time)
    
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
            # Submit BWA jobs and wait for completion
            logging.info("\n=== Submitting BWA jobs ===")
            for job_dir in args.job_dirs:
                bwa_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir)
                        if f.startswith("bwa_") and f.endswith("_alig_job.sh")]
                
                bwa_job_ids = submit_jobs_with_limit(bwa_jobs, args.max_jobs)
                logging.info(f"Submitted {len(bwa_job_ids)} BWA jobs from {job_dir}")
                
                # Wait for BWA jobs to complete
                logging.info("Waiting for BWA jobs to complete...")
                wait_for_jobs_to_complete(job_ids=bwa_job_ids,
                                    check_interval=args.check_interval,
                                    max_wait_time=args.max_wait_time)
        
        # STEP 2: SAM to BAM conversion
        logging.info("\n=== STEP 2: Generating SAM to BAM conversion jobs ===")
        generate_sam_to_bam_jobs(args.bwa_dirs, args.bam_dirs, args.job_dirs,
                                samtools_path=args.samtools_path,
                                partition=args.partition, time=args.time, email=args.email,
                                mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Submit SAM to BAM jobs and wait for completion
            logging.info("\n=== Submitting SAM to BAM conversion jobs ===")
            for job_dir in args.job_dirs:
                bam_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir)
                        if f.startswith("bam_") and f.endswith("_sort_job.sh")]
                
                bam_job_ids = submit_jobs_with_limit(bam_jobs, args.max_jobs)
                logging.info(f"Submitted {len(bam_job_ids)} SAM to BAM jobs from {job_dir}")
                
                # Wait for SAM to BAM jobs to complete
                logging.info("Waiting for SAM to BAM jobs to complete...")
                wait_for_jobs_to_complete(job_ids=bam_job_ids,
                                    check_interval=args.check_interval,
                                    max_wait_time=args.max_wait_time)
        
        # STEP 3: Stampy processing
        logging.info("\n=== STEP 3: Generating Stampy jobs ===")
        generate_stampy_jobs(args.bam_dirs, args.stampy_dirs, args.job_dirs,
                        args.reference, args.python_2_7_path, args.stampy_path,
                        partition=args.partition, time=args.time, email=args.email,
                        mem_per_cpu=args.mem_per_cpu, cpus=args.stampy_cpus)
        
        if args.submit:
            # Submit Stampy jobs and wait for completion
            logging.info("\n=== Submitting Stampy jobs ===")
            for job_dir in args.job_dirs:
                stampy_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir)
                            if f.startswith("stampy_") and f.endswith("_job.sh")]
                
                stampy_job_ids = submit_jobs_with_limit(stampy_jobs, args.max_jobs)
                logging.info(f"Submitted {len(stampy_job_ids)} Stampy jobs from {job_dir}")
                
                # Wait for Stampy jobs to complete
                logging.info("Waiting for Stampy jobs to complete...")
                wait_for_jobs_to_complete(job_ids=stampy_job_ids,
                                    check_interval=args.check_interval,
                                    max_wait_time=args.max_wait_time)
            
            logging.info("\nAll alignment jobs have been submitted and completed!")
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