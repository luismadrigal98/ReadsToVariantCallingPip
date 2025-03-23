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
import subprocess
import time

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Include the source directory in the path
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))

# Import the fastq utilities module
from fastq_utilities import *

def get_job_count(user=None):
    """Get the number of jobs in the SLURM queue for a specific user."""
    if user is None:
        # Get current username if not specified
        user = subprocess.check_output("whoami", shell=True).decode().strip()
    
    cmd = f"squeue -u {user} -h | wc -l"
    try:
        count = int(subprocess.check_output(cmd, shell=True).decode().strip())
        return count
    except subprocess.CalledProcessError:
        logging.error(f"Error checking job count for user {user}")
        return 0

def submit_jobs_with_limit(job_files, max_jobs=4900, sleep_time=60):
    """
    Submit jobs while respecting a maximum job limit.
    
    Parameters:
    job_files (list): List of job script paths to submit
    max_jobs (int): Maximum number of jobs allowed in queue
    sleep_time (int): Time to wait between submission batches in seconds
    
    Returns:
    list: List of job IDs that were submitted
    """
    submitted_job_ids = []
    
    for job_file in job_files:
        # Check current job count
        while get_job_count() >= max_jobs:
            logging.info(f"Job limit reached ({max_jobs}). Waiting for jobs to complete...")
            time.sleep(sleep_time)
        
        # Submit the job
        try:
            output = subprocess.check_output(f"sbatch {job_file}", shell=True).decode().strip()
            job_id = output.split()[-1]  # Extract job ID from sbatch output
            submitted_job_ids.append(job_id)
            logging.info(f"Submitted job {job_id} from {job_file}")
            
            # Small delay to avoid overwhelming the scheduler
            time.sleep(0.5)
            
        except subprocess.CalledProcessError as e:
            logging.error(f"Error submitting job {job_file}: {e}")
    
    return submitted_job_ids

def wait_for_jobs_to_complete(job_ids=None, job_name_pattern=None, check_interval=60, max_wait_time=86400):
    """
    Wait for SLURM jobs to complete.
    
    Parameters:
    job_ids (list): List of job IDs to wait for
    job_name_pattern (str): Pattern to match job names (alternative to job_ids)
    check_interval (int): Time between job status checks in seconds
    max_wait_time (int): Maximum time to wait in seconds (default: 24 hours)
    
    Returns:
    bool: True if all jobs completed successfully, False otherwise
    """
    if job_ids is None and job_name_pattern is None:
        logging.error("Either job_ids or job_name_pattern must be provided")
        return False
    
    user = subprocess.check_output("whoami", shell=True).decode().strip()
    start_time = time.time()
    
    while time.time() - start_time < max_wait_time:
        if job_ids:
            # Check if any of the specific job IDs are still running
            running_jobs = []
            for job_id in job_ids:
                cmd = f"squeue -j {job_id} -u {user} -h"
                try:
                    output = subprocess.check_output(cmd, shell=True).decode().strip()
                    if output:  # If there's output, the job is still running
                        running_jobs.append(job_id)
                except subprocess.CalledProcessError:
                    # Job not found in queue, which means it completed (or errored)
                    pass
            
            if not running_jobs:
                logging.info("All jobs completed")
                return True
            else:
                logging.info(f"Waiting for {len(running_jobs)} jobs to complete...")
        
        elif job_name_pattern:
            # Check if any jobs matching the pattern are still running
            cmd = f"squeue -u {user} -h -n {job_name_pattern} | wc -l"
            try:
                count = int(subprocess.check_output(cmd, shell=True).decode().strip())
                if count == 0:
                    logging.info(f"All jobs matching '{job_name_pattern}' completed")
                    return True
                else:
                    logging.info(f"Waiting for {count} jobs matching '{job_name_pattern}' to complete...")
            except subprocess.CalledProcessError:
                logging.error(f"Error checking job status for pattern '{job_name_pattern}'")
        
        time.sleep(check_interval)
    
    logging.error(f"Timed out waiting for jobs to complete after {max_wait_time/3600:.1f} hours")
    return False

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
    
    # Add this code to your main function

    # Workflow command (runs all steps)
    workflow_parser = subparsers.add_parser("workflow", parents=[common_parser], 
                                        help="Run complete workflow (split, compress, fastp)")
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
                                if f.startswith("gzip_in_batch_") and f.endswith(".sh")]
                
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
                fastp_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir) 
                            if f.startswith("fastp_in_batch_") and f.endswith(".sh")]
                
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