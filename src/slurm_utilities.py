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
    Wait for SLURM jobs to complete and detect failed jobs
    
    Parameters:
    job_ids (list): List of job IDs to wait for
    job_name_pattern (str): Pattern to match job names (alternative to job_ids)
    check_interval (int): Time between job status checks in seconds
    max_wait_time (int): Maximum time to wait in seconds (default: 24 hours)
    
    Returns:
    dict: {'completed': [job_ids], 'failed': [job_ids], 'timeout': [job_ids]}
    """
    
    # Check if there are no jobs to wait for
    if job_ids is not None and len(job_ids) == 0:
        logging.info("No jobs to wait for - continuing immediately")
        return {'completed': [], 'failed': [], 'timeout': []}
    
    # Check if job IDs or job name pattern is provided
    if job_ids is None and job_name_pattern is None:
        logging.error("Either job_ids or job_name_pattern must be provided")
        return {'completed': [], 'failed': [], 'timeout': []}
    
    user = subprocess.check_output("whoami", shell=True).decode().strip()
    start_time = time.time()
    
    completed_jobs = []
    failed_jobs = []
    pending_jobs = set(job_ids) if job_ids else set()
    
    while time.time() - start_time < max_wait_time:
        if job_ids:
            # Check status of each job
            still_running = []
            
            for job_id in list(pending_jobs):
                # Check if job is still in queue
                cmd = f"squeue -j {job_id} -u {user} -h 2>/dev/null"
                try:
                    output = subprocess.check_output(cmd, shell=True).decode().strip()
                    if output:  # Job is still running/pending
                        still_running.append(job_id)
                    else:
                        # Job completed - check if it was successful
                        # Use sacct to check job exit status
                        sacct_cmd = f"sacct -j {job_id} --format=State,ExitCode --parsable2 --noheader 2>/dev/null"
                        try:
                            sacct_output = subprocess.check_output(sacct_cmd, shell=True).decode().strip()
                            lines = sacct_output.split('\n')
                            
                            # Look for the main job status (not sub-tasks)
                            job_failed = False
                            for line in lines:
                                if line and not '.batch' in line and not '.extern' in line:
                                    parts = line.split('|')
                                    if len(parts) >= 2:
                                        state = parts[0].strip()
                                        exit_code = parts[1].strip()
                                        
                                        # Check for failure states
                                        if state in ['FAILED', 'CANCELLED', 'TIMEOUT', 'NODE_FAIL', 'PREEMPTED']:
                                            job_failed = True
                                            break
                                        elif exit_code != '0:0' and exit_code != '':
                                            job_failed = True
                                            break
                            
                            if job_failed:
                                failed_jobs.append(job_id)
                                logging.warning(f"Job {job_id} failed with state: {state}, exit code: {exit_code}")
                            else:
                                completed_jobs.append(job_id)
                                logging.info(f"Job {job_id} completed successfully")
                                
                        except subprocess.CalledProcessError:
                            # If sacct fails, assume job completed (might be too old)
                            completed_jobs.append(job_id)
                            logging.warning(f"Could not check status of job {job_id} - assuming completed")
                        
                        pending_jobs.remove(job_id)
                        
                except subprocess.CalledProcessError:
                    # Job not found in queue - check with sacct
                    # This handles the case where job completed but we missed it
                    sacct_cmd = f"sacct -j {job_id} --format=State,ExitCode --parsable2 --noheader 2>/dev/null"
                    try:
                        sacct_output = subprocess.check_output(sacct_cmd, shell=True).decode().strip()
                        if sacct_output:
                            lines = sacct_output.split('\n')
                            job_failed = False
                            for line in lines:
                                if line and not '.batch' in line and not '.extern' in line:
                                    parts = line.split('|')
                                    if len(parts) >= 2:
                                        state = parts[0].strip()
                                        exit_code = parts[1].strip()
                                        
                                        if state in ['FAILED', 'CANCELLED', 'TIMEOUT', 'NODE_FAIL', 'PREEMPTED']:
                                            job_failed = True
                                            break
                                        elif exit_code != '0:0' and exit_code != '':
                                            job_failed = True
                                            break
                            
                            if job_failed:
                                failed_jobs.append(job_id)
                                logging.warning(f"Job {job_id} failed")
                            else:
                                completed_jobs.append(job_id)
                                logging.info(f"Job {job_id} completed successfully")
                        else:
                            # No record found, assume completed
                            completed_jobs.append(job_id)
                            logging.warning(f"No record found for job {job_id} - assuming completed")
                    except subprocess.CalledProcessError:
                        completed_jobs.append(job_id)
                        logging.warning(f"Could not check status of job {job_id} - assuming completed")
                    
                    if job_id in pending_jobs:
                        pending_jobs.remove(job_id)
            
            if not pending_jobs:
                logging.info(f"All jobs completed. Successful: {len(completed_jobs)}, Failed: {len(failed_jobs)}")
                return {'completed': completed_jobs, 'failed': failed_jobs, 'timeout': []}
            else:
                logging.info(f"Waiting for {len(pending_jobs)} jobs to complete... (Completed: {len(completed_jobs)}, Failed: {len(failed_jobs)})")
        
        elif job_name_pattern:
            # For job name pattern, just check if any are running (legacy behavior)
            cmd = f"squeue -u {user} -h -n {job_name_pattern} | wc -l"
            try:
                count = int(subprocess.check_output(cmd, shell=True).decode().strip())
                if count == 0:
                    logging.info(f"All jobs matching '{job_name_pattern}' completed")
                    return {'completed': [], 'failed': [], 'timeout': []}
                else:
                    logging.info(f"Waiting for {count} jobs matching '{job_name_pattern}' to complete...")
            except subprocess.CalledProcessError:
                logging.error("Error checking job status")
                return {'completed': [], 'failed': [], 'timeout': []}
        
        time.sleep(check_interval)
    
    # Timeout reached
    timeout_jobs = list(pending_jobs) if job_ids else []
    logging.error(f"Timeout reached after {max_wait_time} seconds. Jobs still pending: {timeout_jobs}")
    return {'completed': completed_jobs, 'failed': failed_jobs, 'timeout': timeout_jobs}


def retry_failed_jobs(failed_job_scripts, max_jobs, max_retries=1):
    """
    Retry failed jobs up to max_retries times
    
    Parameters:
    failed_job_scripts (list): List of job script paths that failed
    max_jobs (int): Maximum number of concurrent jobs
    max_retries (int): Maximum number of retry attempts (default: 1)
    
    Returns:
    dict: {'success': [job_scripts], 'final_failures': [job_scripts]}
    """
    
    if not failed_job_scripts:
        return {'success': [], 'final_failures': []}
    
    retry_count = 0
    jobs_to_retry = failed_job_scripts.copy()
    final_failures = []
    successful_retries = []
    
    while retry_count < max_retries and jobs_to_retry:
        retry_count += 1
        logging.info(f"=== RETRY ATTEMPT {retry_count}/{max_retries} for {len(jobs_to_retry)} failed jobs ===")
        
        # Submit the failed jobs again
        retry_job_ids = submit_jobs_with_limit(jobs_to_retry, max_jobs)
        logging.info(f"Resubmitted {len(retry_job_ids)} jobs for retry attempt {retry_count}")
        
        if not retry_job_ids:
            logging.error("Failed to resubmit any jobs")
            final_failures.extend(jobs_to_retry)
            break
        
        # Wait for retry jobs to complete
        retry_result = wait_for_jobs_to_complete(
            job_ids=retry_job_ids,
            check_interval=60,
            max_wait_time=86400
        )
        
        # Map job IDs back to job scripts for successful retries
        successful_job_ids = retry_result['completed']
        failed_job_ids = retry_result['failed'] + retry_result['timeout']
        
        # Find which job scripts succeeded
        jobs_to_retry_next = []
        for i, job_id in enumerate(retry_job_ids):
            job_script = jobs_to_retry[i] if i < len(jobs_to_retry) else None
            if job_script:
                if job_id in successful_job_ids:
                    successful_retries.append(job_script)
                    logging.info(f"Retry successful for: {job_script}")
                elif job_id in failed_job_ids:
                    jobs_to_retry_next.append(job_script)
                    logging.warning(f"Retry failed for: {job_script}")
        
        jobs_to_retry = jobs_to_retry_next
        
        if not jobs_to_retry:
            logging.info("All retry attempts successful!")
            break
    
    # Any remaining jobs are final failures
    final_failures.extend(jobs_to_retry)
    
    if final_failures:
        logging.error(f"FINAL FAILURES after {retry_count} retry attempts: {len(final_failures)} jobs")
        for job in final_failures:
            logging.error(f"  - {job}")
    
    return {'success': successful_retries, 'final_failures': final_failures}


def submit_jobs_with_retry(job_scripts, max_jobs, max_retries=1, check_interval=60, max_wait_time=86400):
    """
    Submit jobs with automatic retry on failure
    
    Parameters:
    job_scripts (list): List of job script paths to submit
    max_jobs (int): Maximum number of concurrent jobs
    max_retries (int): Maximum number of retry attempts for failed jobs
    check_interval (int): Time between job status checks in seconds
    max_wait_time (int): Maximum time to wait for jobs in seconds
    
    Returns:
    dict: {'all_successful': bool, 'completed': [scripts], 'final_failures': [scripts]}
    """
    
    if not job_scripts:
        return {'all_successful': True, 'completed': [], 'final_failures': []}
    
    logging.info(f"Submitting {len(job_scripts)} jobs with retry capability (max retries: {max_retries})")
    
    # Initial submission
    job_ids = submit_jobs_with_limit(job_scripts, max_jobs)
    if not job_ids:
        logging.error("Failed to submit any jobs")
        return {'all_successful': False, 'completed': [], 'final_failures': job_scripts}
    
    logging.info(f"Initial submission: {len(job_ids)} jobs")
    
    # Wait for initial jobs to complete
    result = wait_for_jobs_to_complete(
        job_ids=job_ids,
        check_interval=check_interval,
        max_wait_time=max_wait_time
    )
    
    completed_jobs = result['completed']
    failed_job_ids = result['failed'] + result['timeout']
    
    # Map failed job IDs back to job scripts
    failed_job_scripts = []
    successful_job_scripts = []
    
    for i, job_id in enumerate(job_ids):
        job_script = job_scripts[i] if i < len(job_scripts) else None
        if job_script:
            if job_id in failed_job_ids:
                failed_job_scripts.append(job_script)
            elif job_id in completed_jobs:
                successful_job_scripts.append(job_script)
    
    logging.info(f"Initial results: {len(successful_job_scripts)} successful, {len(failed_job_scripts)} failed")
    
    # Retry failed jobs if any
    if failed_job_scripts and max_retries > 0:
        retry_result = retry_failed_jobs(failed_job_scripts, max_jobs, max_retries)
        successful_job_scripts.extend(retry_result['success'])
        final_failures = retry_result['final_failures']
    else:
        final_failures = failed_job_scripts
    
    all_successful = len(final_failures) == 0
    
    logging.info(f"FINAL RESULTS: {len(successful_job_scripts)} successful, {len(final_failures)} failed")
    
    return {
        'all_successful': all_successful,
        'completed': successful_job_scripts,
        'final_failures': final_failures
    }