"""
Miscellaneous utilities for variant calling

@author: Luis Javier Madrigal Roca

"""

import subprocess
import sys
import logging


def atomize_vcf_file(input_file, output_file, bcftools_path="bcftools", verbose = False):
    """Atomize a VCF file into individual variant calls."""
    command = [bcftools_path, "norm", "-a", input_file, ">", output_file]
    
    if verbose:
        logging.info(f"Running command: {' '.join(command)}")
    
    subprocess.run(command, check=True)