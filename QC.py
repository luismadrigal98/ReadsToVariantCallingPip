#!/usr/bin/env python3

"""
This program will be used to run fastp over a set of fastq files on a high performance cluster.

This essentially will be a wrapper script that will integrate all required steps from preprocessing of fastq files prior to
the mapping step using fastp.

Author: Luis Javier Madrigal-Roca

Date: 2025-03-11

"""

import os
import argparse
import sys
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Include the source directory in the path
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))

# Import the fastq utilities module
from fastq_utilities import *
from slurm_utilities import *

def main():
    parser = argparse.ArgumentParser(description="FASTQ preprocessing for variant calling pipeline")
    
    # Required subparsers for different operations
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")
    
    # Common arguments for all commands
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument("--email", type=str, default="l338m483@ku.edu", help="Email for job notifications")
    common_parser.add_argument("--partition", type=str, default="sixhour", help="SLURM partition to use")
    common_parser.add_argument("--time", type=str, default="6:00:00", help="Time limit for jobs")
    common_parser.add_argument("--mem-per-cpu", type=str, default="5g", help="Memory per CPU")
    common_parser.add_argument("--cpus", type=int, default=10, help="Number of CPUs per task")
    
    # Split command
    split_parser = subparsers.add_parser("split", parents=[common_parser], help="Split FASTQ files")
    split_parser.add_argument("--input-dirs", type=str, nargs="+", required=True, 
                            help="Directories containing FASTQ files")
    split_parser.add_argument("--output-dirs", type=str, nargs="+", required=True,
                            help="Directories for split output files")
    split_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    split_parser.add_argument("--lines", type=int, default=4000000,
                            help="Number of lines per split file (default: 4000000)")
    
    # Compress command
    compress_parser = subparsers.add_parser("compress", parents=[common_parser], help="Compress split FASTQ files")
    compress_parser.add_argument("--batch-dirs", type=str, nargs="+", required=True,
                            help="Directories containing split files")
    compress_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    
    # Fastp command
    fastp_parser = subparsers.add_parser("fastp", parents=[common_parser], help="Run fastp on FASTQ files")
    fastp_parser.add_argument("--batch-dirs", type=str, nargs="+", required=True,
                            help="Directories containing FASTQ files")
    fastp_parser.add_argument("--output-dirs", type=str, nargs="+", required=True,
                            help="Directories for fastp output")
    fastp_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    fastp_parser.add_argument("--fastp-path", type=str, default="/home/l338m483/fastp",
                            help="Path to fastp executable")
    fastp_parser.add_argument("--fastp-control-param", type=str,  # Change underscore to hyphen
                        default="-3 --complexity_threshold=20 --length_required=50 --cut_window_size=3 --cut_mean_quality=30",
                        help="Control parameters common to single and double end reads. Notice that you must provide an unique string separated by spaces.")

    # Workflow command (runs all steps)
    workflow_parser = subparsers.add_parser("workflow", parents=[common_parser], 
                                        help="Run complete workflow (split, compress, fastp). Notice that this is meant to be run with the submit flag enabled.")
    workflow_parser.add_argument("--input-dirs", type=str, nargs="+", required=True, 
                            help="Directories containing original FASTQ files")
    workflow_parser.add_argument("--split-dirs", type=str, nargs="+", required=True,
                            help="Directories for split output files (intermediate)")
    workflow_parser.add_argument("--fastp-dirs", type=str, nargs="+", required=True,
                            help="Directories for final fastp output")
    workflow_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    workflow_parser.add_argument("--lines", type=int, default=4000000,
                            help="Number of lines per split file (default: 4000000)")
    workflow_parser.add_argument("--fastp-path", type=str, default="/home/l338m483/fastp",
                            help="Path to fastp executable")
    workflow_parser.add_argument("--fastp-control-param", type=str, 
                            default="-3 --complexity_threshold=20 --length_required=50 --cut_window_size=3 --cut_mean_quality=30",
                            help="Control parameters common to single and double end reads")
    workflow_parser.add_argument("--submit", action="store_true",
                            help="Automatically submit jobs after generation")
    workflow_parser.add_argument("--max-jobs", type=int, default=5000,
                            help="Maximum number of jobs to have in queue at once (default: 5000)")
    workflow_parser.add_argument("--check-interval", type=int, default=300,
                            help="Time between job status checks in seconds (default: 300)")
    workflow_parser.add_argument("--max-wait-time", type=int, default=86400,
                            help="Maximum time to wait for jobs in seconds (default: 24 hours)")

    # Parse arguments
    args = parser.parse_args()
    
    # Execute appropriate command
    if args.command == "split":
        if len(args.input_dirs) != len(args.output_dirs) or len(args.input_dirs) != len(args.job_dirs):
            logging.error("Error: Number of input directories, output directories, and job directories must match")
            sys.exit(1)
        generate_split_jobs(args.input_dirs, args.output_dirs, args.job_dirs, lines=args.lines,
                        partition=args.partition, time=args.time, email=args.email, 
                        mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
    
    elif args.command == "compress":
        if len(args.batch_dirs) != len(args.job_dirs):
            logging.error("Error: Number of batch directories and job directories must match")
            sys.exit(1)
        generate_compress_jobs(args.batch_dirs, args.job_dirs, partition=args.partition,
                            time=args.time, email=args.email, mem_per_cpu=args.mem_per_cpu,
                            cpus=args.cpus)
    
    elif args.command == "fastp":
        if len(args.batch_dirs) != len(args.output_dirs) or len(args.batch_dirs) != len(args.job_dirs):
            logging.error("Error: Number of batch directories, output directories, and job directories must match")
            sys.exit(1)
        generate_fastp_jobs(args.batch_dirs, args.output_dirs, args.job_dirs, fastp_path=args.fastp_path,
                        fastp_control_param=args.fastp_control_param,
                        partition=args.partition, time=args.time, email=args.email,
                        mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
    
    elif args.command == "workflow":
        if len(args.input_dirs) != len(args.split_dirs) or len(args.input_dirs) != len(args.fastp_dirs) or len(args.input_dirs) != len(args.job_dirs):
            logging.error("Error: Number of directories must match across all arguments")
            sys.exit(1)
        
        logging.info("=== STEP 1: Generating split jobs ===")
        generate_split_jobs(args.input_dirs, args.split_dirs, args.job_dirs, lines=args.lines,
                        partition=args.partition, time=args.time, email=args.email, 
                        mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        logging.info("\n=== STEP 2: Generating compression jobs ===")
        generate_compress_jobs(args.split_dirs, args.job_dirs, partition=args.partition,
                            time=args.time, email=args.email, mem_per_cpu=args.mem_per_cpu,
                            cpus=args.cpus)
        
        logging.info("\n=== STEP 3: Generating fastp jobs ===")
        generate_fastp_jobs(args.split_dirs, args.fastp_dirs, args.job_dirs, fastp_path=args.fastp_path,
                        fastp_control_param=args.fastp_control_param,
                        partition=args.partition, time=args.time, email=args.email,
                        mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        logging.info("\nWorkflow job generation complete. Submit jobs in order:")
        logging.info("1. Run split jobs")
        logging.info("2. Run compression jobs")
        logging.info("3. Run fastp jobs")

        if args.submit:
            logging.info("\nSubmitting jobs automatically...")
            
            # Loop through each job directory and submit jobs in sequence
            for job_dir in args.job_dirs:
                # Step 1: Submit split jobs
                logging.info(f"\n=== Submitting split jobs from {job_dir} ===")
                split_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir) 
                            if f.startswith("Split_") and f.endswith(".sh")]
                
                if not split_jobs:
                    logging.warning(f"No split jobs found in {job_dir}")
                    continue
                
                split_job_ids = submit_jobs_with_limit(split_jobs, args.max_jobs)
                logging.info(f"Submitted {len(split_job_ids)} split jobs from {job_dir}")
                
                # Wait for split jobs to complete
                logging.info(f"Waiting for split jobs to complete...")
                wait_success = wait_for_jobs_to_complete(
                    job_ids=split_job_ids,
                    check_interval=args.check_interval,
                    max_wait_time=args.max_wait_time
                )
                
                if not wait_success:
                    logging.error(f"Timed out waiting for split jobs in {job_dir}. Skipping next steps.")
                    continue
                
                # Step 2: Submit compression jobs
                logging.info(f"\n=== Submitting compression jobs from {job_dir} ===")
                compress_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir) 
                                if f.startswith("gzip_in_batch_") and f.endswith("_compress_job.sh")]
                
                if not compress_jobs:
                    logging.warning(f"No compression jobs found in {job_dir}")
                    continue
                
                compress_job_ids = submit_jobs_with_limit(compress_jobs, args.max_jobs)
                logging.info(f"Submitted {len(compress_job_ids)} compression jobs from {job_dir}")
                
                # Wait for compression jobs to complete
                logging.info(f"Waiting for compression jobs to complete...")
                wait_success = wait_for_jobs_to_complete(
                    job_ids=compress_job_ids,
                    check_interval=args.check_interval,
                    max_wait_time=args.max_wait_time
                )
                
                if not wait_success:
                    logging.error(f"Timed out waiting for compression jobs in {job_dir}. Skipping next steps.")
                    continue
                
                # Step 3: Submit fastp jobs
                logging.info(f"\n=== Submitting fastp jobs from {job_dir} ===")
                # Change to:
                fastp_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir) 
                                if (f.startswith("fastp_paired_") or f.startswith("fastp_single_")) and f.endswith(".sh")]
                
                if not fastp_jobs:
                    logging.warning(f"No fastp jobs found in {job_dir}")
                    continue
                
                fastp_job_ids = submit_jobs_with_limit(fastp_jobs, args.max_jobs)
                logging.info(f"Submitted {len(fastp_job_ids)} fastp jobs from {job_dir}")
                
                # Wait for fastp jobs to complete (optional)
                logging.info(f"Waiting for fastp jobs to complete...")
                wait_for_jobs_to_complete(
                    job_ids=fastp_job_ids,
                    check_interval=args.check_interval,
                    max_wait_time=args.max_wait_time
                )
            
            logging.info("\nAll jobs have been submitted and completed!")

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()