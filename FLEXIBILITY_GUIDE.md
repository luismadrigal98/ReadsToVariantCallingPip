# Pipeline Flexibility and Portability Guide

## Overview

The StampyToFreeBayesPip has been designed with flexibility in mind to accommodate different:
- FASTQ naming conventions
- Computing environments
- User preferences
- Institutional setups

## Key Flexible Features

### 1. File Pattern Recognition

The pipeline can now handle various FASTQ naming conventions through customizable patterns:

#### Default Illumina Convention
```bash
# Files like: Sample1_S1_L001_R1_001.fastq.gz, Sample1_S1_L001_R2_001.fastq.gz
python QC.py fastp --batch-dirs /data --output-dirs /output --job-dirs /jobs
```

#### Custom Patterns
```bash
# For files like: Sample1.read1.fastq.gz, Sample1.read2.fastq.gz
python QC.py fastp ... --r1-pattern ".read1." --r2-pattern ".read2."

# For files like: Sample1_1.fastq.gz, Sample1_2.fastq.gz  
python QC.py fastp ... --r1-pattern "_1." --r2-pattern "_2."

# For files like: Sample1.R1.fastq.gz, Sample1.R2.fastq.gz
python QC.py fastp ... --r1-pattern ".R1." --r2-pattern ".R2."
```

### 2. Sample ID Extraction Methods

Control how unique job identifiers are generated:

#### Auto (Default) - Smart Detection
```bash
# Automatically detects split files vs. original files
# Split files: Uses full name (JKK-16A_S3_L004_SET3_001_fu)
# Original files: Uses prefix (JKK-16A)
python QC.py fastp ... --sample-id-method auto
```

#### Prefix - First Part Only
```bash
# Always uses first part before underscore
# JKK-16A_S3_L004_SET3_R1_001_fu.gz → JKK-16A
python QC.py fastp ... --sample-id-method prefix
```

#### Remove Extensions - Full Base Name
```bash
# Uses complete filename without extensions
# JKK-16A_S3_L004_SET3_R1_001_fu.gz → JKK-16A_S3_L004_SET3_R1_001_fu
python QC.py fastp ... --sample-id-method remove_extensions
```

#### Full - Complete Filename
```bash
# Uses original filename as-is
python QC.py fastp ... --sample-id-method full
```

### 3. Environment Configuration

#### Tool Paths
```bash
# Default: assumes tools are in PATH
python QC.py fastp --fastp-path fastp

# Custom installation
python QC.py fastp --fastp-path /opt/bioinformatics/fastp/bin/fastp
python QC.py fastp --fastp-path /home/user/anaconda3/bin/fastp
```

#### Email Notifications
```bash
# With email notifications
python QC.py fastp ... --email researcher@university.edu

# Without email notifications (default)
python QC.py fastp ...
```

#### SLURM Configuration
```bash
# Default settings (suitable for many clusters)
python QC.py fastp ... --partition sixhour --time 6:00:00 --mem-per-cpu 5g --cpus 10

# Custom for high-memory nodes
python QC.py fastp ... --partition highmem --time 12:00:00 --mem-per-cpu 16g --cpus 8

# Custom for long-running jobs
python QC.py fastp ... --partition long --time 48:00:00 --mem-per-cpu 8g --cpus 16
```

## Institution-Specific Setup Examples

### Example 1: University of Kansas
```bash
# Your original setup (now default-free)
python QC.py fastp \
  --batch-dirs /home/user/scratch/data \
  --output-dirs /home/user/scratch/output \
  --job-dirs /home/user/scratch/jobs \
  --email user@ku.edu \
  --partition "sixhour,eeb,kelly,kucg" \
  --fastp-path /home/user/fastp \
  --submit
```

### Example 2: Generic HPC Cluster
```bash
python QC.py fastp \
  --batch-dirs /scratch/project/data \
  --output-dirs /scratch/project/output \
  --job-dirs /home/user/jobs \
  --partition compute \
  --time 24:00:00 \
  --mem-per-cpu 8g \
  --cpus 12 \
  --fastp-path fastp \
  --submit
```

### Example 3: Different Sequencing Center
```bash
# Files named: SP001_1.fq.gz, SP001_2.fq.gz
python QC.py fastp \
  --batch-dirs /data/sequencing/batch1 \
  --output-dirs /data/processed/batch1 \
  --job-dirs /jobs/batch1 \
  --r1-pattern "_1." \
  --r2-pattern "_2." \
  --sample-id-method prefix \
  --partition normal \
  --submit
```

## Migration Guide

### For Existing Users
1. **No changes required**: Your existing commands will work with sensible defaults
2. **Optional improvements**: Add `--email your@email.edu` for notifications
3. **Tool paths**: If fastp is in PATH, remove hardcoded paths

### For New Users
1. **Start simple**: Use defaults and add customizations as needed
2. **Test patterns**: Use a small dataset to verify file pattern matching
3. **Check job names**: Ensure sample IDs are unique with your naming convention

### Converting from Old Version
```bash
# Old hardcoded version
python QC.py fastp --batch-dirs /data --output-dirs /output --job-dirs /jobs --submit

# New flexible version (same result)
python QC.py fastp \
  --batch-dirs /data \
  --output-dirs /output \
  --job-dirs /jobs \
  --email your@email.edu \
  --submit
```

## Best Practices

### File Organization
- Use consistent naming conventions within projects
- Test with a small subset before processing large datasets
- Keep job directories organized by project/date

### Resource Management
- Start with default resources and adjust based on performance
- Monitor job completion rates and adjust `--max-retries`
- Use appropriate `--max-jobs` for your cluster policies

### Troubleshooting
1. **Check file detection**: Run without `--submit` first to see generated jobs
2. **Verify patterns**: List files and confirm R1/R2 patterns match
3. **Test sample IDs**: Ensure unique job names are generated
4. **Check resources**: Verify SLURM settings work for your cluster

## Future Extensions

The flexible design allows for easy future enhancements:
- Additional file format support
- More sophisticated pattern matching
- Configuration file support
- Integration with other workflow managers

## Support

For questions about flexibility features:
1. Check the pattern matching with a test dataset
2. Review generated job scripts before submitting
3. Contact the developer with specific use cases for further customization