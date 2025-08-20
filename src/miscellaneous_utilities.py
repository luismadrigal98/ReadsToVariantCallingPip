"""
Miscellaneous utilities for variant calling

@author: Luis Javier Madrigal Roca

"""

import subprocess
import sys

def atomize_vcf_file(input_file, output_file, bcftools_path="bcftools"):
    """Atomize a VCF file into individual variant calls."""
    command = [bcftools_path, "norm", "-a", input_file, ">", output_file]
    subprocess.run(command, check=True)