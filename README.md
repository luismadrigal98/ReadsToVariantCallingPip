## StampyToFreeBayesPip: Genomic Alignment Pipeline

# Overview

StampyToFreeBayesPip is a comprehensive genomic alignment pipeline designed for high-performance computing environments. It automates the process from raw FASTQ files to analysis-ready BAM files, leveraging BWA and Stampy aligners for improved accuracy. The pipeline is optimized for variant detection and is particularly useful for non-model organisms or divergent reference genomes.

# Features

- FASTQ preprocessing with quality control via fastp
- Efficient alignment using BWA with Stampy refinement
- Flexible duplicate handling strategies
- Parallel processing via SLURM job submission
- Modular design allowing execution of specific pipeline components
- Support for both paired-end and single-end sequencing data

# Installation

# Prerequisites

Python 3.6+
SLURM workload manager
Python 2.7 (for Stampy)
The following bioinformatics tools:
fastp
BWA
Stampy
samtools
Picard

# Setup

Clone the repository:

```
git clone https://github.com/luismadrigal98/StampyToFreeBayesPip.git
cd StampyToFreeBayesPip
```

Ensure all dependencies are installed and in your PATH.

Make the pipeline scripts executable:

```
chmod +x QC.py StampyRunner.py BamProcessor.py
```

# Pipeline Components

1. QC.py - FASTQ Preprocessing
Handles preprocessing of FASTQ files using fastp, with options for splitting large files.

2. StampyRunner.py - Sequence Alignment
Performs alignment using BWA followed by Stampy for improved mapping quality.

3. BamProcessor.py - BAM File Processing
Processes aligned reads, including merging, duplicate handling, and indexing.

# Usage

Quality Control and Preprocessing (QC.py)

```
# Run complete preprocessing workflow
./QC.py workflow \
  --input-dirs /path/to/fastq/dir1 /path/to/fastq/dir2 \
  --split-dirs /path/to/split/dir1 /path/to/split/dir2 \
  --fastp-dirs /path/to/output/dir1 /path/to/output/dir2 \
  --job-dirs /path/to/job/dir1 /path/to/job/dir2 \
  --submit \
  --fastp-control-param "-3 --complexity_threshold=20 --length_required=50"
```

Alignment (StampyRunner.py)
```
# Run complete alignment workflow
./StampyRunner.py workflow \
  --input-dirs /path/to/preprocessed/dir1 /path/to/preprocessed/dir2 \
  --bwa-dirs /path/to/bwa/dir1 /path/to/bwa/dir2 \
  --bam-dirs /path/to/bam/dir1 /path/to/bam/dir2 \
  --stampy-dirs /path/to/stampy/dir1 /path/to/stampy/dir2 \
  --job-dirs /path/to/job/dir1 /path/to/job/dir2 \
  --reference /path/to/reference.fasta \
  --submit
```

BAM Processing (BamProcessor.py)