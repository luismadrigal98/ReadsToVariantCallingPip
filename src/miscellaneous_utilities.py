"""
Miscellaneous utilities for variant calling

@author: Luis Javier Madrigal Roca

"""

import subprocess
import sys
import logging


def atomize_vcf_file(input_file, output_file, bcftools_path="bcftools", verbose=False, debug_mode=False):
    """Atomize a VCF file into individual variant calls."""
    
    # Option 1: Use shell=True with string command (most similar to terminal)
    command_str = f"{bcftools_path} norm -a {input_file} > {output_file}"
    
    if verbose:
        logging.info(f"Running command: {command_str}")
    
    if not debug_mode:
        # Use shell=True to handle redirection properly
        # Capture stderr to control what's shown
        result = subprocess.run(command_str, shell=True, check=True, 
                              capture_output=False, text=True)