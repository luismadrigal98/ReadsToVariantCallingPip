"""
This module provides utility functions for interacting with the SLURM workload manager. 
It includes functions to check the number of jobs in the queue, submit jobs while respecting 
a maximum job limit, and wait for jobs to complete.
Functions:
- get_job_count(user=None): Retrieves the number of jobs in the SLURM queue for a specific user.
- submit_jobs_with_limit(job_files, max_jobs=4900, sleep_time=60): Submits job scripts to SLURM 
    while ensuring the total number of jobs in the queue does not exceed a specified limit.
- wait_for_jobs_to_complete(job_ids=None, job_name_pattern=None, check_interval=60, max_wait_time=86400): 
    Waits for specified SLURM jobs to complete, either by job IDs or by matching a job name pattern.
Dependencies:
- subprocess: Used for executing shell commands and retrieving their output.
- logging: Used for logging informational and error messages.
- time: Used for implementing delays and measuring elapsed time.
Usage:
This module is intended for use in environments where SLURM is available and configured. 
It is particularly useful for managing job submissions and monitoring job statuses in 
high-performance computing (HPC) clusters.

"""

import subprocess
import time
import logging

def get_job_count(user=None):
    """Get the number of jobs in the SLURM queue for a specific user."""
    if user is None:
        # Get current username if not specified
        user = subprocess.check_output("whoami", shell=True).decode().strip()
    
    cmd = f"squeue -u {user} -h | wc -l"
    try:
        count = int(subprocess.check_output(cmd, shell=True).decode().strip())
        return count
    except subprocess.CalledProcessError:
        logging.error(f"Error checking job count for user {user}")
        return 0

def submit_jobs_with_limit(job_files, max_jobs=4900, sleep_time=60):
    """
    Submit jobs while respecting a maximum job limit.
    
    Parameters:
    job_files (list): List of job script paths to submit
    max_jobs (int): Maximum number of jobs allowed in queue
    sleep_time (int): Time to wait between submission batches in seconds
    
    Returns:
    list: List of job IDs that were submitted
    """
    submitted_job_ids = []
    
    for job_file in job_files:
        # Check current job count
        while get_job_count() >= max_jobs:
            logging.info(f"Job limit reached ({max_jobs}). Waiting for jobs to complete...")
            time.sleep(sleep_time)
        
        # Submit the job
        try:
            output = subprocess.check_output(f"sbatch {job_file}", shell=True).decode().strip()
            # More robust job ID extraction
            if "Submitted batch job" in output:
                job_id = output.split()[-1]
                # Verify it's a number
                if job_id.isdigit():
                    submitted_job_ids.append(job_id)
                    logging.info(f"Submitted job {job_id} from {job_file}")
                else:
                    logging.warning(f"Could not extract job ID from output: {output}")
            else:
                logging.warning(f"Unexpected sbatch output: {output}")
            
            # Small delay to avoid overwhelming the scheduler
            time.sleep(0.5)
            
        except subprocess.CalledProcessError as e:
            logging.error(f"Error submitting job {job_file}: {e}")
    
    return submitted_job_ids

def wait_for_jobs_to_complete(job_ids=None, job_name_pattern=None, check_interval=60, max_wait_time=86400):
    """
    Wait for SLURM jobs to complete.
    
    Parameters:
    job_ids (list): List of job IDs to wait for
    job_name_pattern (str): Pattern to match job names (alternative to job_ids)
    check_interval (int): Time between job status checks in seconds
    max_wait_time (int): Maximum time to wait in seconds (default: 24 hours)
    
    Returns:
    bool: True if all jobs completed successfully, False otherwise
    """
    
    # Check if there are no jobs to wait for
    if job_ids is not None and len(job_ids) == 0:
        logging.info("No jobs to wait for - continuing immediately")
        return True
    
    # Check if job IDs or job name pattern is provided
    if job_ids is None and job_name_pattern is None:
        logging.error("Either job_ids or job_name_pattern must be provided")
        return False
    
    user = subprocess.check_output("whoami", shell=True).decode().strip()
    start_time = time.time()
    
    while time.time() - start_time < max_wait_time:
        if job_ids:
            # Check if any of the specific job IDs are still running
            running_jobs = []
            for job_id in job_ids:
                cmd = f"squeue -j {job_id} -u {user} -h 2>/dev/null"
                try:
                    output = subprocess.check_output(cmd, shell=True).decode().strip()
                    if output:  # If there's output, the job is still running
                        running_jobs.append(job_id)
                except subprocess.CalledProcessError:
                    # Job not found in queue, which means it completed (or errored)
                    pass
            
            if not running_jobs:
                logging.info("All jobs completed")
                return True
            else:
                logging.info(f"Waiting for {len(running_jobs)} jobs to complete...")
        
        elif job_name_pattern:
            # Check if any jobs matching the pattern are still running
            cmd = f"squeue -u {user} -h -n {job_name_pattern} | wc -l"
            try:
                count = int(subprocess.check_output(cmd, shell=True).decode().strip())
                if count == 0:
                    logging.info(f"All jobs matching '{job_name_pattern}' completed")
                    return True
                else:
                    logging.info(f"Waiting for {count} jobs matching '{job_name_pattern}' to complete...")
            except subprocess.CalledProcessError:
                logging.error(f"Error checking job status for pattern '{job_name_pattern}'")
        
        time.sleep(check_interval)
    
    logging.error(f"Timed out waiting for jobs to complete after {max_wait_time/3600:.1f} hours")
    return False