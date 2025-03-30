#!/usr/bin/env python3

import argparse
import logging
import sys
from src.variant_caller_utilities import *
from src.slurm_utilities import *

def main():
    """Main function to parse arguments and execute variant calling workflow."""
    parser = argparse.ArgumentParser(description='Call variants on processed BAM files')
    
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
    common_parser.add_argument("--mem-per-cpu", type=str, default="30g",
                            help="Memory per CPU")
    common_parser.add_argument("--submit", action="store_true",
                            help="Submit jobs to SLURM")
    common_parser.add_argument("--max-jobs", type=int, default=5000,
                            help="Maximum number of jobs to submit at once")
    common_parser.add_argument("--check-interval", type=int, default=60,
                            help="Seconds between job status checks")
    common_parser.add_argument("--max-wait-time", type=int, default=86400,
                            help="Maximum seconds to wait for jobs")
    common_parser.add_argument("--samtools_path", type=str, default='/kuhpc/sw/conda/latest/envs/bioconda/bin/samtools',
                            required=False, help="Path to samtools executable")
    
    # Call variants command
    call_parser = subparsers.add_parser('call', parents=[common_parser],
                                    help='Call variants on BAM files')
    call_parser.add_argument("--input-dirs", type=str, nargs='+', required=True,
                            help="Directories containing processed BAM files")
    call_parser.add_argument("--output-dirs", type=str, nargs='+', required=True,
                            help="Directories for variant calling output")
    call_parser.add_argument("--job-dirs", type=str, nargs='+', required=True,
                            help="Directories for job scripts")
    call_parser.add_argument("--regions", type=str, nargs='+', default=["all"],
                        help="Chromosomes or regions to process ('all' for entire genome)")
    call_parser.add_argument("--reference", type=str, required=True,
                            help="Path to reference FASTA file")
    call_parser.add_argument("--fai", type=str, default=None,
                            help="Path to FASTA index file (default: [reference].fai)")
    call_parser.add_argument("--window-size", type=int, default=1000000,
                            help="Window size for chunking variant calling")
    call_parser.add_argument("--variant-caller", type=str, default="freebayes",
                            choices=["freebayes", "bcftools", "gatk"],
                            help="Variant caller to use")
    call_parser.add_argument("--caller-params", type=str, default=None,
                            help="Additional parameters for the variant caller")
    call_parser.add_argument("--variant-caller-path", type=str, default='~/.conda/envs/Python2.7/bin/freebayes',
                            help="Path to variant caller executable")
    
    # Process arguments
    args = parser.parse_args()
    
    # Configure logging
    logging.basicConfig(level=logging.INFO, 
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # Check if command is specified
    if args.command is None:
        parser.print_help()
        sys.exit(1)
    
    # Handle variant calling
    if args.command == 'call':
        # Validate arguments
        if len(args.input_dirs) != len(args.output_dirs) or len(args.input_dirs) != len(args.job_dirs):
            logging.error("Number of input, output, and job directories must match")
            sys.exit(1)
        
        # Set default FAI path if not provided
        if args.fai is None:
            args.fai = f"{args.reference}.fai"

        # Check if index file exist, if not, create it
        if not os.path.exists(args.fai):
            logging.info(f"Creating index file for reference: {args.reference}")
            create_fasta_index(args.reference, args.samtools_path)
            args.fai = f"{args.reference}.fai"
            logging.info(f"Index file created: {args.fai}")
        else:
            logging.info(f"Using existing index file: {args.fai}")
        # Check if input directories exist
        for input_dir in args.input_dirs:
            if not os.path.exists(input_dir):
                logging.error(f"Input directory does not exist: {input_dir}")
                sys.exit(1)
        # Check if output directories exist, if not, create them
        for output_dir in args.output_dirs:
            if not os.path.exists(output_dir):
                logging.info(f"Creating output directory: {output_dir}")
                os.makedirs(output_dir, exist_ok=True)
            else:
                logging.info(f"Output directory already exists: {output_dir}")
        # Check if job directories exist, if not, create them
        for job_dir in args.job_dirs:
            if not os.path.exists(job_dir):
                logging.info(f"Creating job directory: {job_dir}")
                os.makedirs(job_dir, exist_ok=True)
            else:
                logging.info(f"Job directory already exists: {job_dir}")
        # Check if variant caller executable exists
        if not os.path.exists(args.variant_caller_path):
            logging.error(f"Variant caller executable does not exist: {args.variant_caller_path}")
            sys.exit(1)
        
        # Generate variant calling jobs
        try:
            logging.info(f"Generating variant calling jobs using {args.variant_caller}")
            job_scripts = generate_variant_calling_jobs(
                input_dirs=args.input_dirs,
                output_dirs=args.output_dirs,
                job_dirs=args.job_dirs,
                regions=args.regions,
                reference_fasta=args.reference,
                fai_path=args.fai,
                window_size=args.window_size,
                variant_caller=args.variant_caller,
                caller_path=args.variant_caller_path,
                partition=args.partition,
                time=args.time,
                email=args.email,
                mem_per_cpu=args.mem_per_cpu,
                caller_params=args.caller_params
            )
            
            logging.info(f"Generated {len(job_scripts)} variant calling jobs")
            
            # Submit jobs if requested
            if args.submit:
                logging.info("Submitting variant calling jobs to SLURM")
                job_ids = submit_jobs_with_limit(job_scripts, args.max_jobs)
                
                # Wait for jobs to complete
                logging.info(f"Waiting for {len(job_ids)} variant calling jobs to complete")
                wait_for_jobs(job_ids, args.check_interval, args.max_wait_time)
                logging.info("All variant calling jobs have completed")
            else:
                logging.info("Jobs created but not submitted. Use --submit to submit jobs.")
                
        except Exception as e:
            logging.error(f"Error generating variant calling jobs: {e}")
            sys.exit(1)

if __name__ == "__main__":
    main()