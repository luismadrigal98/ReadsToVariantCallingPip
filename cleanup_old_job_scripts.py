#!/usr/bin/env python3

"""
cleanup_old_job_scripts.py - Clean up legacy job scripts referencing the 'eeb' partition

This script safely removes job scripts (.sh files) that contain references to the 'eeb' partition.
It provides a dry-run mode for safety and detailed logging of what will be/was removed.

Usage:
    python cleanup_old_job_scripts.py [--dry-run] [--job-dirs DIR1 DIR2 ...]
    
Examples:
    # Check what would be removed (safe preview)
    python cleanup_old_job_scripts.py --dry-run
    
    # Actually remove the scripts
    python cleanup_old_job_scripts.py
    
    # Clean specific directories only
    python cleanup_old_job_scripts.py --job-dirs /path/to/jobs1 /path/to/jobs2
"""

import os
import argparse
import logging
import glob
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def find_job_directories():
    """Find common job directory patterns in the current workspace."""
    common_patterns = [
        "*job*", "*Job*", "*jobs*", "*Jobs*",
        "*/job*", "*/Job*", "*/jobs*", "*/Jobs*"
    ]
    
    job_dirs = set()
    for pattern in common_patterns:
        matches = glob.glob(pattern)
        for match in matches:
            if os.path.isdir(match):
                job_dirs.add(os.path.abspath(match))
    
    return sorted(list(job_dirs))

def find_eeb_job_scripts(job_dirs):
    """Find all .sh job scripts that reference the 'eeb' partition."""
    eeb_scripts = []
    
    for job_dir in job_dirs:
        if not os.path.exists(job_dir):
            logging.warning(f"Job directory does not exist: {job_dir}")
            continue
            
        logging.info(f"Searching in: {job_dir}")
        sh_files = glob.glob(os.path.join(job_dir, "*.sh"))
        
        for sh_file in sh_files:
            try:
                with open(sh_file, 'r') as f:
                    content = f.read()
                    
                # Check for 'eeb' partition references in SLURM directives
                if 'eeb' in content.lower():
                    lines = content.split('\n')
                    for line in lines:
                        line = line.strip()
                        if (line.startswith('#SBATCH') and 
                            ('--partition' in line or '-p' in line) and 
                            'eeb' in line.lower()):
                            eeb_scripts.append(sh_file)
                            break
                        elif '--partition=eeb' in line.lower() or '-p eeb' in line.lower():
                            eeb_scripts.append(sh_file)
                            break
                            
            except Exception as e:
                logging.warning(f"Error reading {sh_file}: {e}")
                
    return eeb_scripts

def cleanup_scripts(scripts, dry_run=True):
    """Remove the specified job scripts."""
    if not scripts:
        logging.info("✅ No legacy 'eeb' job scripts found - workspace is clean!")
        return
        
    logging.info(f"\n{'='*60}")
    if dry_run:
        logging.info("DRY RUN - Preview of files that WOULD be removed:")
    else:
        logging.info("REMOVING legacy 'eeb' job scripts:")
    logging.info(f"{'='*60}")
    
    removed_count = 0
    total_size = 0
    
    for script in scripts:
        try:
            file_size = os.path.getsize(script)
            total_size += file_size
            
            if dry_run:
                logging.info(f"  WOULD REMOVE: {script} ({file_size} bytes)")
            else:
                os.remove(script)
                logging.info(f"  ✓ REMOVED: {script} ({file_size} bytes)")
                removed_count += 1
                
        except Exception as e:
            logging.error(f"  ✗ ERROR with {script}: {e}")
    
    if dry_run:
        logging.info(f"\n📊 SUMMARY (DRY RUN):")
        logging.info(f"  Would remove: {len(scripts)} files")
        logging.info(f"  Total size: {total_size:,} bytes ({total_size/1024:.1f} KB)")
        logging.info(f"\n💡 To actually remove these files, run without --dry-run")
    else:
        logging.info(f"\n📊 CLEANUP SUMMARY:")
        logging.info(f"  Successfully removed: {removed_count}/{len(scripts)} files")
        logging.info(f"  Total space freed: {total_size:,} bytes ({total_size/1024:.1f} KB)")
        if removed_count < len(scripts):
            logging.warning(f"  Failed to remove: {len(scripts) - removed_count} files")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Clean up legacy job scripts referencing the 'eeb' partition"
    )
    
    parser.add_argument(
        "--dry-run", 
        action="store_true",
        help="Preview files that would be removed without actually deleting them"
    )
    
    parser.add_argument(
        "--job-dirs",
        nargs="+",
        help="Specific job directories to clean (default: auto-detect)"
    )
    
    args = parser.parse_args()
    
    # Determine job directories to search
    if args.job_dirs:
        job_dirs = [os.path.abspath(d) for d in args.job_dirs]
        logging.info(f"Using specified job directories: {job_dirs}")
    else:
        job_dirs = find_job_directories()
        if job_dirs:
            logging.info(f"Auto-detected job directories: {job_dirs}")
        else:
            logging.error("No job directories found. Specify directories with --job-dirs")
            sys.exit(1)
    
    # Find legacy scripts
    logging.info("Searching for legacy 'eeb' job scripts...")
    eeb_scripts = find_eeb_job_scripts(job_dirs)
    
    # Clean up
    cleanup_scripts(eeb_scripts, dry_run=args.dry_run)
    
    logging.info("\n🎉 Cleanup operation completed!")

if __name__ == "__main__":
    main()