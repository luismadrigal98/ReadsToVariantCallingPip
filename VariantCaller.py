#!/usr/bin/env python3

"""
VariantCaller.py - A tool for running variant calling on processed BAM files

This script enables variant calling using multiple callers (FreeBayes, bcftools, GATK)
with support for:
- Region-based chunking for parallel execution
- Joint variant calling across samples
- Paired input directories for cohort analysis
- VCF file merging across genomic regions

This tool is part of the StampyToFreeBayesPip pipeline.
"""

import os
import argparse
import logging
import sys

# Set up logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Include the source directory in the path
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))

# Import utilities
from src.variant_caller_utilities import *
from src.slurm_utilities import *
from src.miscellaneous_utilities import atomize_vcf_file

def print_job_summary(step_name, result, total_jobs):
    """
    Print a comprehensive summary of job execution results.
    
    Parameters:
    step_name (str): Name of the pipeline step (e.g., "FreeBayes", "BCFtools", "GATK")
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
    """Main function to parse arguments and execute variant calling."""
    parser = argparse.ArgumentParser(
        description='Run variant calling on processed BAM files')
    
    # Common arguments for all commands
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument("--partition", type=str, default="sixhour",
                            help="SLURM partition to use")
    common_parser.add_argument("--constraint", type=str, default="avx2&noib",
                            help="SLURM constraint for node selection (default: CPU-only nodes)")
    common_parser.add_argument("--time", type=str, default="6:00:00",
                            help="Time limit for jobs (format: HH:MM:SS)")
    common_parser.add_argument("--email", type=str, default="l338m483@ku.edu",
                            help="Email for job notifications")
    common_parser.add_argument("--mem-per-cpu", type=str, default="20g",
                            help="Memory per CPU")
    common_parser.add_argument("--threads", type=int, default=2,
                            help="Number of threads per task")
    common_parser.add_argument("--submit", action="store_true",
                            help="Submit jobs to SLURM")
    common_parser.add_argument("--max-jobs", type=int, default=500,
                            help="Maximum number of jobs to submit at once")
    common_parser.add_argument("--check-interval", type=int, default=120,
                            help="Seconds between job status checks")
    common_parser.add_argument("--max-wait-time", type=int, default=86400,
                            help="Maximum seconds to wait for jobs")
    common_parser.add_argument("--max-retries", type=int, default=2,
                            help="Maximum number of retry attempts for failed jobs")
    common_parser.add_argument("--abort-on-failure", action="store_true",
                            help="Abort workflow if jobs fail after retries")
    common_parser.add_argument("--regions", type=str, nargs="+", default=["all"],
                            help="Genomic regions to process (e.g., Chr_01 Chr_02)")
    common_parser.add_argument("--window-size", type=int, default=1000000,
                            help="Window size for variant calling")
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    
    # Call command
    call_parser = subparsers.add_parser("call", parents=[common_parser],
                                        help="Run variant calling")
    call_parser.add_argument("--reference", type=str, required=True,
                            help="Path to reference FASTA file")
    call_parser.add_argument("--input-dirs", type=str, nargs="+", required=True,
                            help="Directories containing BAM files")
    call_parser.add_argument("--output-dirs", type=str, nargs="+", required=True,
                            help="Directories for variant calling output")
    call_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    call_parser.add_argument("--variant-caller", type=str, default="freebayes",
                            choices=["freebayes", "bcftools", "gatk"],
                            help="Variant caller to use")
    call_parser.add_argument("--caller-path", type=str, default=None,
                            help="Path to variant caller executable")
    call_parser.add_argument("--caller-params", type=str, default=None,
                            help="Additional parameters for variant caller")
    call_parser.add_argument("--samtools-path", type=str,
                            default="/kuhpc/sw/conda/latest/envs/bioconda/bin/samtools",
                            help="Path to samtools executable")
    
    # Joint call command
    joint_parser = subparsers.add_parser("joint-call", parents=[common_parser],
                                        help="Run joint variant calling with paired directories")
    joint_parser.add_argument("--reference", type=str, required=True,
                            help="Path to reference FASTA file")
    joint_parser.add_argument("--input-dirs-1", type=str, nargs="+", required=True,
                            help="First set of directories containing BAM files")
    joint_parser.add_argument("--input-dirs-2", type=str, nargs="+", required=True,
                            help="Second set of directories containing BAM files")
    joint_parser.add_argument("--output-dirs", type=str, nargs="+", required=True,
                            help="Directories for variant calling output")
    joint_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    joint_parser.add_argument("--variant-caller", type=str, default="bcftools",
                            choices=["freebayes", "bcftools", "gatk"],
                            help="Variant caller to use")
    joint_parser.add_argument("--caller-path", type=str, default=None,
                            help="Path to variant caller executable")
    joint_parser.add_argument("--caller-params", type=str, default=None,
                            help="Additional parameters for variant caller")
    
    # Merge command
    merge_parser = subparsers.add_parser("merge", parents=[common_parser],
                                       help="Merge VCF files")
    merge_parser.add_argument("--input-dirs", type=str, nargs="+", required=True,
                            help="Directories containing VCF files")
    merge_parser.add_argument("--output-dirs", type=str, nargs="+", required=True,
                            help="Directories for merged output")
    merge_parser.add_argument("--job-dirs", type=str, nargs="+", required=True,
                            help="Directories for job scripts")
    merge_parser.add_argument("--bcftools-path", type=str,
                            default="/kuhpc/sw/conda/latest/envs/bioconda/bin/bcftools",
                            help="Path to bcftools executable")
    merge_parser.add_argument("--merge-command", type=str, default="concat",
                            choices=["concat", "merge"],
                            help="bcftools command to use for merging")
    merge_parser.add_argument("--merge-mode", type=str, default="by_chromosome",
                            choices=["global", "by_chromosome"],
                            help="How to merge VCF files")
    merge_parser.add_argument("--sample-names", type=str, nargs="+", default=None,
                            help="Sample names for output files")
    
    # Atomize command
    atomize_parser = subparsers.add_parser("atomize",
                                            help="Run atomization on input files")
    atomize_parser.add_argument("--input-file", type=str, required=True,
                                help="Input vcf file to atomize")
    atomize_parser.add_argument("--output-file", type=str, required=True,
                                help="Output vcf file with atomized variants")
    atomize_parser.add_argument('--bcftools_path', type=str, default='bcftools',
                                help="Path to bcftools executable")
    atomize_parser.add_argument('--verbose', action='store_true',
                                help="Enable verbose logging")
    atomize_parser.add_argument('--debug-mode', action='store_true',
                                help="Enable debug mode. Command is not executed, but the function print the command.")

    # Parse arguments
    args = parser.parse_args()
    
    # Execute appropriate command
    if args.command == "call":
        if len(args.input_dirs) != len(args.output_dirs) or len(args.input_dirs) != len(args.job_dirs):
            logging.error("Error: Number of input, output, and job directories must match")
            sys.exit(1)
        
        # Check if the reference exists and create an index if needed
        if not os.path.exists(args.reference):
            logging.error(f"Reference file not found: {args.reference}")
            sys.exit(1)
        
        fai_path = args.reference + ".fai"
        if not os.path.exists(fai_path):
            logging.info(f"Reference index not found, creating one...")
            create_fasta_index(args.reference, args.samtools_path)

        # Generate variant calling jobs
        jobs = generate_variant_calling_jobs(
            input_dirs=args.input_dirs,
            output_dirs=args.output_dirs,
            job_dirs=args.job_dirs,
            regions=args.regions,
            reference_fasta=args.reference,
            fai_path=fai_path,
            window_size=args.window_size,
            variant_caller=args.variant_caller,
            caller_path=args.caller_path,
            partition=args.partition,
            time=args.time,
            email=args.email,
            mem_per_cpu=args.mem_per_cpu,
            caller_params=args.caller_params,
            threads=args.threads,
            constraint=args.constraint
        )
        
        if args.submit and jobs:
            # Submit variant calling jobs with retry capability
            logging.info(f"Submitting {len(jobs)} variant calling jobs with retry capability")
            result = submit_jobs_with_retry(
                job_scripts=jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            vc_successful = print_job_summary("Variant Calling", result, len(jobs))
            if not vc_successful:
                if args.abort_on_failure:
                    logging.error("Aborting due to variant calling failures")
                    sys.exit(1)
                else:
                    logging.warning("Variant calling completed with some failures - data may be incomplete")
    
    elif args.command == "joint-call":
        if (len(args.input_dirs_1) != len(args.input_dirs_2) or 
            len(args.input_dirs_1) != len(args.output_dirs) or 
            len(args.input_dirs_1) != len(args.job_dirs)):
            logging.error("Error: Number of directories must match across all parameters")
            sys.exit(1)
        
        # Check if the reference exists and create an index if needed
        if not os.path.exists(args.reference):
            logging.error(f"Reference file not found: {args.reference}")
            sys.exit(1)
        
        fai_path = args.reference + ".fai"
        if not os.path.exists(fai_path):
            logging.info(f"Reference index not found, creating one...")
            create_fasta_index(args.reference, args.samtools_path)
            
        # Generate joint variant calling jobs
        jobs = generate_variant_calling_jobs(
            input_dirs=args.input_dirs_1,
            output_dirs=args.output_dirs,
            job_dirs=args.job_dirs,
            regions=args.regions,
            reference_fasta=args.reference,
            fai_path=fai_path,
            window_size=args.window_size,
            variant_caller=args.variant_caller,
            caller_path=args.caller_path,
            partition=args.partition,
            time=args.time,
            email=args.email,
            mem_per_cpu=args.mem_per_cpu,
            caller_params=args.caller_params,
            paired_input_dirs=args.input_dirs_2,
            threads=args.threads,
            constraint=args.constraint
        )
        
        if args.submit and jobs:
            # Submit joint variant calling jobs with retry capability
            logging.info(f"Submitting {len(jobs)} joint variant calling jobs with retry capability")
            result = submit_jobs_with_retry(
                job_scripts=jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            joint_successful = print_job_summary("Joint Variant Calling", result, len(jobs))
            if not joint_successful:
                if args.abort_on_failure:
                    logging.error("Aborting due to joint variant calling failures")
                    sys.exit(1)
                else:
                    logging.warning("Joint variant calling completed with some failures - data may be incomplete")
    
    elif args.command == "merge":
        if len(args.input_dirs) != len(args.output_dirs) or len(args.input_dirs) != len(args.job_dirs):
            logging.error("Error: Number of input, output, and job directories must match")
            sys.exit(1)
        
        # Generate merge jobs
        jobs = merge_vcf_files_jobs_generator(
            input_dirs=args.input_dirs,
            output_dirs=args.output_dirs,
            job_dirs=args.job_dirs,
            bcftools_path=args.bcftools_path,
            partition=args.partition,
            time=args.time,
            email=args.email,
            mem_per_cpu=args.mem_per_cpu,
            threads=args.threads,
            merge_command=args.merge_command,
            sample_names=args.sample_names,
            merge_mode=args.merge_mode
        )
        
        if args.submit and jobs:
            # Submit VCF merge jobs with retry capability
            logging.info(f"Submitting {len(jobs)} VCF merge jobs with retry capability")
            result = submit_jobs_with_retry(
                job_scripts=jobs,
                max_jobs=args.max_jobs,
                max_retries=args.max_retries,
                check_interval=args.check_interval,
                max_wait_time=args.max_wait_time
            )
            
            merge_successful = print_job_summary("VCF Merge", result, len(jobs))
            if not merge_successful:
                if args.abort_on_failure:
                    logging.error("Aborting due to VCF merge failures")
                    sys.exit(1)
                else:
                    logging.warning("VCF merge completed with some failures - data may be incomplete")

    elif args.command == "atomize":
        # Run atomization
        atomize_vcf_file(
            input_file=args.input_file,
            output_file=args.output_file,
            bcftools_path=args.bcftools_path,
            verbose=args.verbose,
            debug_mode=args.debug_mode
        )

    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
