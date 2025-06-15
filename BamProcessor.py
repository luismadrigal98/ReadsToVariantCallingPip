#!/usr/bin/env python3
"""
BamProcessor.py - A tool for processing BAM files after alignment

This script     workflow_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                                help="Directories for job scripts")
    workflow_parser.add_argument("--bam-type", type=str, choices=["bwa", "stampy", "all"],
                                default="all", help="Type of BAM files to merge")
    workflow_parser.add_argument("--duplicate-mode", type=str,rates three processing steps:
1. Merging individual BAM files and sorting the merged file
2. Optional duplicate processing (preserve, mark, or remove)
3. Indexing of BAM files for downstream analysis

This tool is part of the StampyToFreeBayesPip pipeline.
"""

import os
import argparse
import logging
import sys
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Include the source directory in the path
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))

# Import utilities
from bam_utilities import *
from slurm_utilities import *

def main():
    """Main function to parse arguments and execute the workflow."""
    parser = argparse.ArgumentParser(description='Process BAM files (merge, filter, index)')
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    
    # Common arguments for all commands
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument("--partition", type=str, default="sixhour",
                            help="SLURM partition to use")
    common_parser.add_argument("--time", type=str, default="6:00:00",
                            help="Time limit for jobs (format: HH:MM:SS)")
    common_parser.add_argument("--email", type=str, default="l338m483@ku.edu",
                            help="Email for job notifications")
    common_parser.add_argument("--mem-per-cpu", type=str, default="5g",
                            help="Memory per CPU")
    common_parser.add_argument("--cpus", type=int, default=10,
                            help="Number of CPUs per task")
    common_parser.add_argument("--submit", action="store_true",
                            help="Submit jobs to SLURM")
    common_parser.add_argument("--max-jobs", type=int, default=5000,
                            help="Maximum number of jobs to submit at once")
    common_parser.add_argument("--check-interval", type=int, default=60,
                            help="Seconds between job status checks")
    common_parser.add_argument("--max-wait-time", type=int, default=86400,
                            help="Maximum seconds to wait for jobs")
    common_parser.add_argument("--samtools-path", type=str, 
                            default="/kuhpc/sw/conda/latest/envs/bioconda/bin/samtools",
                            help="Path to samtools executable")
    common_parser.add_argument("--picard-path", type=str, 
                            default="~/.conda/envs/Python2.7/bin/picard",
                            help="Path to picard executable")
    
    # Command: merge
    merge_parser = subparsers.add_parser("merge", parents=[common_parser],
                                        help="Merge and sort BAM files")
    merge_parser.add_argument("--input-dirs", type=str, nargs="+", required=True,
                            help="Directories containing BAM files")
    merge_parser.add_argument("--output-dirs", type=str, nargs="+", required=True,
                            help="Directories for merged BAM output")
    merge_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    merge_parser.add_argument("--bam-type", type=str, choices=["bwa", "stampy", "all"],
                            default="all", help="Type of BAM files to merge")
    
    # Command: dedup
    dedup_parser = subparsers.add_parser("dedup", parents=[common_parser],
                                        help="Process duplicates in merged BAM files")
    dedup_parser.add_argument("--input-dirs", type=str, nargs="+", required=True,
                            help="Directories containing merged BAM files")
    dedup_parser.add_argument("--output-dirs-base", type=str, nargs="+", required=True,
                            help="Base directories for processed BAM output")
    dedup_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    dedup_parser.add_argument("--duplicate-mode", type=str, 
                            choices=["preserve", "mark", "mark_remove"],
                            default="mark",
                            help="How to handle duplicate reads")
    
    # Command: index
    index_parser = subparsers.add_parser("index", parents=[common_parser],
                                        help="Index BAM files")
    index_parser.add_argument("--input-dirs", type=str, nargs="+", required=True,
                            help="Directories containing BAM files to index")
    index_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    
    # Command: workflow
    workflow_parser = subparsers.add_parser("workflow", parents=[common_parser],
                                            help="Run complete workflow (merge, dedup, index)")
    workflow_parser.add_argument("--input-dirs", type=str, nargs="+", required=True,
                                help="Directories containing BAM files")
    workflow_parser.add_argument("--merged-dirs", type=str, nargs="+", required=True,
                                help="Directories for merged BAM output")
    workflow_parser.add_argument("--processed-dirs-base", type=str, nargs="+", required=True,
                                help="Base directories for processed BAM output")
    workflow_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                                help="Directories for job scripts")
    workflow_parser.add_argument("--bam-type", type=str, choices=["bwa", "stampy", "all"],
                                default="all", help="Type of BAM files to merge")
    workflow_parser.add_argument("--duplicate-mode", type=str, 
                                choices=["preserve", "mark", "mark_remove"],
                                default="mark",
                                help="How to handle duplicate reads")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Execute the appropriate command
    if args.command == "merge":
        if len(args.input_dirs) != len(args.output_dirs) or len(args.input_dirs) != len(args.job_dirs):
            logging.error("Error: Number of input, output, and job directories must match")
            sys.exit(1)
        
        bam_type = None if args.bam_type == "all" else args.bam_type
        
        generate_merge_jobs(args.input_dirs, args.output_dirs, args.job_dirs,
                            samtools_path=args.samtools_path, bam_type=bam_type,
                            partition=args.partition, time=args.time, email=args.email,
                            mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Find and submit merge jobs
            merge_jobs = []
            for job_dir in args.job_dirs:
                merge_jobs.extend([
                    os.path.join(job_dir, f) for f in os.listdir(job_dir)
                    if f.endswith("_merge_and_sort_bam_job.sh")
                ])
            
            merge_job_ids = submit_jobs_with_limit(merge_jobs, args.max_jobs)
            
            if merge_job_ids:
                logging.info("Waiting for merge jobs to complete...")
                wait_for_jobs_to_complete(merge_job_ids, 
                                        check_interval=args.check_interval,
                                        max_wait_time=args.max_wait_time)
    
    elif args.command == "dedup":
        if len(args.input_dirs) != len(args.output_dirs_base) or len(args.input_dirs) != len(args.job_dirs):
            logging.error("Error: Number of input, output, and job directories must match")
            sys.exit(1)
        
        generate_duplicate_processing_jobs(args.input_dirs, args.output_dirs_base, args.job_dirs,
                                            duplicate_mode=args.duplicate_mode,
                                            samtools_path=args.samtools_path, picard_path=args.picard_path,
                                            partition=args.partition, time=args.time, email=args.email,
                                            mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Find and submit duplicate processing jobs
            dedup_jobs = []
            for job_dir in args.job_dirs:
                dedup_jobs.extend([
                    os.path.join(job_dir, f) for f in os.listdir(job_dir)
                    if f.endswith("_duplicate_processing_job.sh")
                ])
            
            dedup_job_ids = submit_jobs_with_limit(dedup_jobs, args.max_jobs)
            
            if dedup_job_ids:
                logging.info("Waiting for duplicate processing jobs to complete...")
                wait_for_jobs_to_complete(dedup_job_ids, 
                                        check_interval=args.check_interval,
                                        max_wait_time=args.max_wait_time)
    
    elif args.command == "index":
        if len(args.input_dirs) != len(args.job_dirs):
            logging.error("Error: Number of input and job directories must match")
            sys.exit(1)
        
        generate_indexing_jobs(args.input_dirs, args.job_dirs,
                                samtools_path=args.samtools_path, picard_path=args.picard_path,
                                partition=args.partition, time=args.time, email=args.email,
                                mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Find and submit indexing jobs
            index_jobs = []
            for job_dir in args.job_dirs:
                index_jobs.extend([
                    os.path.join(job_dir, f) for f in os.listdir(job_dir)
                    if f.endswith("_job.sh") and "indexing" in f
                ])
            
            index_job_ids = submit_jobs_with_limit(index_jobs, args.max_jobs)
            
            if index_job_ids:
                logging.info("Waiting for indexing jobs to complete...")
                wait_for_jobs_to_complete(index_job_ids, 
                                        check_interval=args.check_interval,
                                        max_wait_time=args.max_wait_time)
    
    elif args.command == "workflow":
        if (len(args.input_dirs) != len(args.merged_dirs) or 
            len(args.input_dirs) != len(args.processed_dirs_base) or 
            len(args.input_dirs) != len(args.job_dirs)):
            logging.error("Error: Number of directories must match across all parameters")
            sys.exit(1)
        
        bam_type = None if args.bam_type == "all" else args.bam_type
        
        logging.info("=== STEP 1: Merging BAM files ===")
        generate_merge_jobs(args.input_dirs, args.merged_dirs, args.job_dirs,
                            samtools_path=args.samtools_path, bam_type=bam_type,
                            partition=args.partition, time=args.time, email=args.email,
                            mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Find and submit merge jobs
            merge_jobs = []
            for job_dir in args.job_dirs:
                merge_jobs.extend([
                    os.path.join(job_dir, f) for f in os.listdir(job_dir)
                    if f.endswith("_merge_and_sort_bam_job.sh")
                ])
            
            merge_job_ids = submit_jobs_with_limit(merge_jobs, args.max_jobs)
            
            if merge_job_ids:
                logging.info("Waiting for merge jobs to complete...")
                wait_for_jobs_to_complete(merge_job_ids, 
                                        check_interval=args.check_interval,
                                        max_wait_time=args.max_wait_time)
        
        logging.info("=== STEP 2: Indexing merged BAM files ===")
        generate_indexing_jobs(args.merged_dirs, args.job_dirs,
                                samtools_path=args.samtools_path, picard_path=args.picard_path,
                                partition=args.partition, time=args.time, email=args.email,
                                mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Find and submit indexing jobs for merged files
            index_jobs = []
            for job_dir in args.job_dirs:
                index_jobs.extend([
                    os.path.join(job_dir, f) for f in os.listdir(job_dir)
                    if f.endswith("_job.sh") and "indexing" in f
                ])
            
            index_job_ids = submit_jobs_with_limit(merge_jobs, args.max_jobs)
            
            if index_job_ids:
                logging.info("Waiting for indexing jobs to complete...")
                wait_for_jobs_to_complete(index_job_ids, 
                                        check_interval=args.check_interval,
                                        max_wait_time=args.max_wait_time)
        
        logging.info("=== STEP 3: Processing duplicates ===")
        generate_duplicate_processing_jobs(args.merged_dirs, args.processed_dirs_base, args.job_dirs,
                                            duplicate_mode=args.duplicate_mode,
                                            samtools_path=args.samtools_path, picard_path=args.picard_path,
                                            partition=args.partition, time=args.time, email=args.email,
                                            mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Find and submit duplicate processing jobs
            dedup_jobs = []
            for job_dir in args.job_dirs:
                dedup_jobs.extend([
                    os.path.join(job_dir, f) for f in os.listdir(job_dir)
                    if f.endswith("_duplicate_processing_job.sh")
                ])
            
            dedup_job_ids = submit_jobs_with_limit(dedup_jobs, args.max_jobs)
            
            if dedup_job_ids:
                logging.info("Waiting for duplicate processing jobs to complete...")
                wait_for_jobs_to_complete(dedup_job_ids, 
                                        check_interval=args.check_interval,
                                        max_wait_time=args.max_wait_time)
        
        logging.info("=== STEP 4: Indexing processed BAM files ===")
        # Collect all processed output directories
        logging.info("=== STEP 4: Indexing processed BAM files ===")
        # Collect all processed output directories that actually exist
        processed_dirs = []
        for base_dir in args.processed_dirs_base:
            if args.duplicate_mode == "preserve":
                suffix = "Pipeline.with.duplicates"
            elif args.duplicate_mode == "mark":
                suffix = "Pipeline.with.marked.duplicates"
            else:
                suffix = "Pipeline.without.duplicates"
            
            # Only add directories that actually exist
            if bam_type is None or bam_type == "bwa":
                bwa_dir = os.path.join(base_dir, suffix, "bwa")
                if os.path.exists(bwa_dir):
                    processed_dirs.append(bwa_dir)
                else:
                    logging.warning(f"Skipping non-existent directory: {bwa_dir}")
            
            if bam_type is None or bam_type == "stampy":
                stampy_dir = os.path.join(base_dir, suffix, "stampy")
                if os.path.exists(stampy_dir):
                    processed_dirs.append(stampy_dir)
                else:
                    logging.warning(f"Skipping non-existent directory: {stampy_dir}")
        
        generate_indexing_jobs(processed_dirs, args.job_dirs,
                                samtools_path=args.samtools_path, picard_path=args.picard_path,
                                partition=args.partition, time=args.time, email=args.email,
                                mem_per_cpu=args.mem_per_cpu, cpus=args.cpus)
        
        if args.submit:
            # Find and submit indexing jobs for processed files
            index_jobs = []
            for job_dir in args.job_dirs:
                new_index_jobs = [
                    os.path.join(job_dir, f) for f in os.listdir(job_dir)
                    if f.endswith("_job.sh") and "indexing" in f and not f in index_jobs
                ]
                index_jobs.extend(new_index_jobs)
            
            index_job_ids = submit_jobs_with_limit(index_jobs, args.max_jobs)
            
            if index_job_ids:
                logging.info("Waiting for indexing jobs to complete...")
                wait_for_jobs_to_complete(index_job_ids, 
                                        check_interval=args.check_interval,
                                        max_wait_time=args.max_wait_time)
        
        logging.info("\nAll BAM processing jobs have been submitted and completed!")
    
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()