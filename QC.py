#!/usr/bin/env python3

"""
This program will be used to run fastp over a set of fastq files on a high perform    compress_parser.add_argument("--max-    fastp_parser.add_argument("--max-jobs", type=int, default=500,
                            help="Maximum number of jobs to submit at once (reduced for stability)")
    fastp_parser.add_argument("--check-interval", type=int, default=120,
                            help="Time between job status checks in seconds (reduced for faster detection)")
    fastp_parser.add_argument("--max-wait-time", type=int, default=86400,
                            help="Maximum time to wait for jobs in seconds")
    fastp_parser.add_argument("--max-retries", type=int, default=2,
                            help="Maximum number of retry attempts for failed jobs")ype=int, default=500,
                               help="Maximum number of jobs to submit at once (reduced for stability)")
    compress_parser.add_argument("--check-interval", type=int, default=120,
                               help="Time between job status checks in seconds (reduced for faster detection)")
    compress_parser.add_argument("--max-wait-time", type=int, default=86400,
                               help="Maximum time to wait for jobs in seconds")
    compress_parser.add_argument("--max-retries", type=int, default=2,
                               help="Maximum number of retry attempts for failed jobs")ster.

This essentially will be a wrapper script that will integrate all required steps for the preprocessing of fastq files prior to
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

def print_job_summary(step_name, result, total_jobs):
    """
    Print a comprehensive summary of job execution results.
    
    Parameters:
    step_name (str): Name of the pipeline step (e.g., "Split", "Compression", "FastP")
    result (dict): Result dictionary from submit_jobs_with_retry
    total_jobs (int): Total number of jobs submitted
    """
    logging.info(f"\n{'='*60}")
    logging.info(f"{step_name.upper()} JOBS SUMMARY")
    logging.info(f"{'='*60}")
    
    completed = len(result['completed'])
    failed = len(result['final_failures'])
    success_rate = (completed / total_jobs * 100) if total_jobs > 0 else 0
    
    logging.info(f"Total jobs submitted:    {total_jobs}")
    logging.info(f"Successfully completed:  {completed}")
    logging.info(f"Final failures:          {failed}")
    logging.info(f"Success rate:            {success_rate:.1f}%")
    
    if result['all_successful']:
        logging.info(f"✅ ALL {step_name.upper()} JOBS COMPLETED SUCCESSFULLY!")
    else:
        logging.warning(f"⚠️  {failed} {step_name.upper()} JOBS FAILED AFTER RETRIES")
        if result['final_failures']:
            logging.warning("Failed job scripts:")
            for job in result['final_failures'][:5]:  # Show first 5 failures
                logging.warning(f"  - {job}")
            if len(result['final_failures']) > 5:
                logging.warning(f"  ... and {len(result['final_failures']) - 5} more")
    
    logging.info(f"{'='*60}\n")
    return result['all_successful']

def main():
    parser = argparse.ArgumentParser(description="FASTQ preprocessing for variant calling pipeline")
    
    # Required subparsers for different operations
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")
    
    # Common arguments for all commands
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument("--email", type=str, default=None, help="Email for job notifications (optional)")
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
    split_parser.add_argument("--submit", action="store_true",
                            help="Submit jobs to SLURM")
    split_parser.add_argument("--max-jobs", type=int, default=500,
                            help="Maximum number of jobs to submit at once (reduced for stability)")
    split_parser.add_argument("--check-interval", type=int, default=120,
                            help="Time between job status checks in seconds (reduced for faster detection)")
    split_parser.add_argument("--max-wait-time", type=int, default=86400,
                            help="Maximum time to wait for jobs in seconds")
    split_parser.add_argument("--max-retries", type=int, default=2,
                            help="Maximum number of retry attempts for failed jobs")
    split_parser.add_argument("--abort-on-failure", action="store_true",
                            help="Abort if jobs fail after retries")
    
    # Compress command
    compress_parser = subparsers.add_parser("compress", parents=[common_parser], help="Compress split FASTQ files")
    compress_parser.add_argument("--batch-dirs", type=str, nargs="+", required=True,
                            help="Directories containing split files")
    compress_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    compress_parser.add_argument("--submit", action="store_true",
                            help="Submit jobs to SLURM")
    compress_parser.add_argument("--max-jobs", type=int, default=5000,
                            help="Maximum number of jobs to submit at once")
    compress_parser.add_argument("--check-interval", type=int, default=300,
                            help="Time between job status checks in seconds")
    compress_parser.add_argument("--max-wait-time", type=int, default=86400,
                            help="Maximum time to wait for jobs in seconds")
    compress_parser.add_argument("--max-retries", type=int, default=1,
                            help="Maximum number of retry attempts for failed jobs")
    compress_parser.add_argument("--abort-on-failure", action="store_true",
                            help="Abort if jobs fail after retries")
    
    # Fastp command
    fastp_parser = subparsers.add_parser("fastp", parents=[common_parser], help="Run fastp on FASTQ files")
    fastp_parser.add_argument("--batch-dirs", type=str, nargs="+", required=True,
                            help="Directories containing FASTQ files")
    fastp_parser.add_argument("--output-dirs", type=str, nargs="+", required=True,
                            help="Directories for fastp output")
    fastp_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    fastp_parser.add_argument("--fastp-path", type=str, default="fastp",
                            help="Path to fastp executable (default: 'fastp' - assumes it's in PATH)")
    fastp_parser.add_argument("--fastp-control-param", type=str,  # Change underscore to hyphen
                        default="-3 --complexity_threshold=20 --length_required=50 --cut_window_size=3 --cut_mean_quality=30",
                        help="Control parameters common to single and double end reads. Notice that you must provide an unique string separated by spaces.")
    fastp_parser.add_argument("--r1-pattern", type=str, default="_R1_",
                            help="Pattern to identify R1 files (default: '_R1_')")
    fastp_parser.add_argument("--r2-pattern", type=str, default="_R2_",
                            help="Pattern to identify R2 files (default: '_R2_')")
    fastp_parser.add_argument("--sample-id-method", type=str, default="auto",
                            choices=["auto", "prefix", "remove_extensions", "full"],
                            help="Method for sample ID extraction (default: 'auto')")
    fastp_parser.add_argument("--submit", action="store_true",
                            help="Submit jobs to SLURM")
    fastp_parser.add_argument("--max-jobs", type=int, default=5000,
                            help="Maximum number of jobs to submit at once")
    fastp_parser.add_argument("--check-interval", type=int, default=300,
                            help="Time between job status checks in seconds (default: 300)")
    fastp_parser.add_argument("--max-wait-time", type=int, default=86400,
                            help="Maximum time to wait for jobs in seconds (default: 24 hours)")
    fastp_parser.add_argument("--max-retries", type=int, default=2,
                            help="Maximum number of retry attempts for failed jobs (default: 2)")
    fastp_parser.add_argument("--abort-on-failure", action="store_true",
                            help="Abort if jobs fail after retries")
    fastp_parser.add_argument("--retry-delay", type=int, default=60,
                            help="Delay in seconds before retrying failed jobs (default: 60)")

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
    workflow_parser.add_argument("--fastp-path", type=str, default="fastp",
                            help="Path to fastp executable (default: 'fastp' - assumes it's in PATH)")
    workflow_parser.add_argument("--fastp-control-param", type=str, 
                            default="-3 --complexity_threshold=20 --length_required=50 --cut_window_size=3 --cut_mean_quality=30",
                            help="Control parameters common to single and double end reads")
    workflow_parser.add_argument("--r1-pattern", type=str, default="_R1_",
                            help="Pattern to identify R1 files (default: '_R1_')")
    workflow_parser.add_argument("--r2-pattern", type=str, default="_R2_",
                            help="Pattern to identify R2 files (default: '_R2_')")
    workflow_parser.add_argument("--sample-id-method", type=str, default="auto",
                            choices=["auto", "prefix", "remove_extensions", "full"],
                            help="Method for sample ID extraction (default: 'auto')")
    workflow_parser.add_argument("--submit", action="store_true",
                            help="Automatically submit jobs after generation")
    workflow_parser.add_argument("--max-jobs", type=int, default=5000,
                            help="Maximum number of jobs to have in queue at once (default: 5000)")
    workflow_parser.add_argument("--check-interval", type=int, default=300,
                            help="Time between job status checks in seconds (default: 300)")
    workflow_parser.add_argument("--max-wait-time", type=int, default=86400,
                            help="Maximum time to wait for jobs in seconds (default: 24 hours)")
    workflow_parser.add_argument("--max-retries", type=int, default=1,
                            help="Maximum number of retry attempts for failed jobs")
    workflow_parser.add_argument("--abort-on-failure", action="store_true",
                            help="Abort workflow if jobs fail after retries")

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
        
        if args.submit:
            # Submit split jobs with retry
            logging.info("=== Submitting split jobs from all directories ===")
            all_split_jobs = []
            
            # Collect all split jobs from all directories
            for job_dir in args.job_dirs:
                split_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir) 
                            if f.startswith("Split_") and f.endswith(".sh")]
                all_split_jobs.extend(split_jobs)
                logging.info(f"Found {len(split_jobs)} split jobs in {job_dir}")
            
            # Submit split jobs with retry
            logging.info(f"Submitting {len(all_split_jobs)} total split jobs with retry capability")
            result = submit_jobs_with_retry(
                job_scripts=all_split_jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            # Print comprehensive summary
            success = print_job_summary("Split", result, len(all_split_jobs))
            
            if not success and args.abort_on_failure:
                logging.error("Aborting due to split failures")
                sys.exit(1)
    
    elif args.command == "compress":
        if len(args.batch_dirs) != len(args.job_dirs):
            logging.error("Error: Number of batch directories and job directories must match")
            sys.exit(1)
        generate_compress_jobs(args.batch_dirs, args.job_dirs, partition=args.partition,
                            time=args.time, email=args.email, mem_per_cpu=args.mem_per_cpu,
                            cpus=args.cpus)
        
        if args.submit:
            # Submit compression jobs with retry capability
            logging.info("=== Submitting compression jobs from all directories ===")
            all_compress_jobs = []
            
            # Collect all compression jobs from all directories
            for job_dir in args.job_dirs:
                compress_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir) 
                                if f.startswith("gzip_in_batch_") and f.endswith("_compress_job.sh")]
                all_compress_jobs.extend(compress_jobs)
                logging.info(f"Found {len(compress_jobs)} compression jobs in {job_dir}")
            
            # Submit compression jobs with retry
            logging.info(f"Submitting {len(all_compress_jobs)} total compression jobs with retry capability")
            result = submit_jobs_with_retry(
                job_scripts=all_compress_jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            # Print comprehensive summary
            success = print_job_summary("Compression", result, len(all_compress_jobs))
            
            if not success and args.abort_on_failure:
                logging.error("Aborting due to compression failures")
                sys.exit(1)
    
    elif args.command == "fastp":
        if len(args.batch_dirs) != len(args.output_dirs) or len(args.batch_dirs) != len(args.job_dirs):
            logging.error("Error: Number of batch directories, output directories, and job directories must match")
            sys.exit(1)
        generate_fastp_jobs(args.batch_dirs, args.output_dirs, args.job_dirs, fastp_path=args.fastp_path,
                        fastp_control_param=args.fastp_control_param,
                        partition=args.partition, time=args.time, email=args.email,
                        mem_per_cpu=args.mem_per_cpu, cpus=args.cpus,
                        r1_pattern=args.r1_pattern, r2_pattern=args.r2_pattern,
                        sample_id_method=args.sample_id_method)
        
        if args.submit:
            # Submit fastp jobs with retry capability
            logging.info("=== Submitting fastp jobs from all directories ===")
            all_fastp_jobs = []
            
            # Collect all fastp jobs from all directories
            for job_dir in args.job_dirs:
                fastp_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir) 
                            if (f.startswith("fastp_paired_") or f.startswith("fastp_single_")) and f.endswith(".sh")]
                all_fastp_jobs.extend(fastp_jobs)
                logging.info(f"Found {len(fastp_jobs)} fastp jobs in {job_dir}")
            
            # Submit fastp jobs with retry
            logging.info(f"Submitting {len(all_fastp_jobs)} total fastp jobs with retry capability")
            result = submit_jobs_with_retry(
                job_scripts=all_fastp_jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            # Print comprehensive summary
            success = print_job_summary("FastP", result, len(all_fastp_jobs))
            
            if not success and args.abort_on_failure:
                logging.error("Aborting due to fastp failures")
                sys.exit(1)
            
            logging.info("FastP processing completed!")
    
    elif args.command == "workflow":
    # Verify directory counts match
        if len(args.input_dirs) != len(args.split_dirs) or len(args.input_dirs) != len(args.fastp_dirs) or len(args.input_dirs) != len(args.job_dirs):
            logging.error("Error: Number of directories must match across all arguments")
            sys.exit(1)
        
        # STEP 1: Split files
        logging.info("=== STEP 1: Generating split jobs ===")
        generate_split_jobs(args.input_dirs, args.split_dirs, args.job_dirs, lines=args.lines,
                        partition=args.partition, time=args.time, email=args.email, 
                        mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Submit split jobs from all directories with retry capability
            logging.info("\n=== Submitting split jobs from all directories ===")
            all_split_jobs = []
            
            # Collect all split jobs from all directories
            for job_dir in args.job_dirs:
                split_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir) 
                            if f.startswith("Split_") and f.endswith(".sh")]
                all_split_jobs.extend(split_jobs)
                logging.info(f"Found {len(split_jobs)} split jobs in {job_dir}")
            
            # Submit split jobs with retry
            logging.info(f"Submitting {len(all_split_jobs)} total split jobs with retry capability")
            split_result = submit_jobs_with_retry(
                job_scripts=all_split_jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            # Print comprehensive summary
            split_success = print_job_summary("Split", split_result, len(all_split_jobs))
            
            if not split_success:
                if args.abort_on_failure:
                    logging.error("Aborting workflow due to split failures")
                    sys.exit(1)
                else:
                    logging.warning("Continuing workflow despite split failures - data may be incomplete")
        
        # STEP 2: Generate compression jobs AFTER split is done
        logging.info("\n=== STEP 2: Generating compression jobs ===")
        generate_compress_jobs(args.split_dirs, args.job_dirs, partition=args.partition,
                            time=args.time, email=args.email, mem_per_cpu=args.mem_per_cpu,
                            cpus=args.cpus)
        
        if args.submit:
            # Submit compression jobs from all directories with retry capability
            logging.info("\n=== Submitting compression jobs from all directories ===")
            all_compress_jobs = []
            
            # Collect all compression jobs from all directories
            for job_dir in args.job_dirs:
                compress_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir) 
                                if f.startswith("gzip_in_batch_") and f.endswith("_compress_job.sh")]
                all_compress_jobs.extend(compress_jobs)
                logging.info(f"Found {len(compress_jobs)} compression jobs in {job_dir}")
            
            # Submit compression jobs with retry
            logging.info(f"Submitting {len(all_compress_jobs)} total compression jobs with retry capability")
            compress_result = submit_jobs_with_retry(
                job_scripts=all_compress_jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            # Print comprehensive summary
            compress_success = print_job_summary("Compression", compress_result, len(all_compress_jobs))
            
            if not compress_success:
                if args.abort_on_failure:
                    logging.error("Aborting workflow due to compression failures")
                    sys.exit(1)
                else:
                    logging.warning("Continuing workflow despite compression failures - data may be incomplete")
        
        # STEP 3: Generate fastp jobs AFTER compression is done
        logging.info("\n=== STEP 3: Generating fastp jobs ===")
        generate_fastp_jobs(args.split_dirs, args.fastp_dirs, args.job_dirs, fastp_path=args.fastp_path,
                        fastp_control_param=args.fastp_control_param,
                        partition=args.partition, time=args.time, email=args.email,
                        mem_per_cpu=args.mem_per_cpu, cpus=args.cpus,
                        r1_pattern=args.r1_pattern, r2_pattern=args.r2_pattern,
                        sample_id_method=args.sample_id_method)
        
        if args.submit:
            # Submit fastp jobs from all directories with retry capability
            logging.info("\n=== Submitting fastp jobs from all directories ===")
            all_fastp_jobs = []
            
            # Collect all fastp jobs from all directories
            for job_dir in args.job_dirs:
                fastp_jobs = [os.path.join(job_dir, f) for f in os.listdir(job_dir) 
                            if (f.startswith("fastp_paired_") or f.startswith("fastp_single_")) and f.endswith(".sh")]
                all_fastp_jobs.extend(fastp_jobs)
                logging.info(f"Found {len(fastp_jobs)} fastp jobs in {job_dir}")
            
            # Submit fastp jobs with retry
            logging.info(f"Submitting {len(all_fastp_jobs)} total fastp jobs with retry capability")
            fastp_result = submit_jobs_with_retry(
                job_scripts=all_fastp_jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            # Print comprehensive summary
            fastp_success = print_job_summary("Fastp", fastp_result, len(all_fastp_jobs))
            
            if not fastp_success:
                if args.abort_on_failure:
                    logging.error("Aborting workflow due to fastp failures")
                    sys.exit(1)
                else:
                    logging.warning("Workflow completed despite fastp failures - data may be incomplete")
            
            logging.info("\nQC workflow completed successfully!")
        else:
            logging.info("\nWorkflow job generation complete. Submit jobs in order:")
            logging.info("1. Run split jobs")
            logging.info("2. Run compression jobs")
            logging.info("3. Run fastp jobs")

    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()