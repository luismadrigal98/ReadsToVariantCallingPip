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
logging.basicConfig(level=logging.INFO)

# Include the source directory in the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "src")))

# Import the fastq utilities module
from fastq.utilities import *

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

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()