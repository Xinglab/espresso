import argparse
import datetime
import os.path
import subprocess
import sys
import time

import cluster_commands


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'job_id', help='the cluster id of the job to check the status of')
    parser.add_argument(
        '--resource-usage-dir',
        help='a directory for storing the file paths where the resource usage'
        ' of each job should be logged')
    parser.add_argument(
        '--resource-usage-min-interval',
        type=float,
        default=120,
        help='only log the resource usage if it has been at least this many'
        ' seconds since the last log')
    args = parser.parse_args()

    return args


def run_command(command):
    completed_process = subprocess.run(command,
                                       check=False,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
    decoded_stdout = completed_process.stdout.decode()
    decoded_stderr = completed_process.stderr.decode()
    if completed_process.returncode != 0:
        print('stdout:\n{}\n\nstderr:\n{}'.format(decoded_stdout,
                                                  decoded_stderr),
              file=sys.stderr)
        sys.exit(1)

    return decoded_stdout


def extract_job_info(stdout, job_id):
    info, error = cluster_commands.try_extract_job_info_from_status_output(
        stdout, job_id)
    if error:
        print('error: {}\n{}'.format(error, stdout), file=sys.stderr)
        sys.exit(1)

    now = datetime.datetime.now()
    formatted_now = now.isoformat()
    info['resource_usage'] = 'current_time: {}, {}'.format(
        formatted_now, info['resource_usage'])

    return info


def update_resource_log(status, resource_usage, resource_dir,
                        min_interval_seconds, job_id):
    if not (resource_dir and os.path.isdir(resource_dir)):
        return

    resource_dir_job_file = os.path.join(resource_dir, '{}.txt'.format(job_id))
    if not os.path.exists(resource_dir_job_file):
        return

    with open(resource_dir_job_file, 'rt') as f_handle:
        resource_log_file = f_handle.read().strip()

    is_final_update = status != 'running'
    update_resource_log_with_file(resource_usage, min_interval_seconds,
                                  resource_log_file, is_final_update)

    if is_final_update:
        os.remove(resource_dir_job_file)


def update_resource_log_with_file(resource_usage, min_interval_seconds,
                                  resource_log_file, is_final_update):
    if not resource_usage:
        return

    # resource_log_file is created and written to by this function.
    # log_dir should have been created by snakemake when the job was submitted.
    log_dir = os.path.dirname(resource_log_file)
    if not (resource_log_file.endswith('.cluster.usage')
            and os.path.isdir(log_dir)):
        return

    # The first write creates the file and does not check min_interval_seconds
    if ((min_interval_seconds and os.path.exists(resource_log_file)
         and not is_final_update)):
        mod_time_seconds = os.stat(resource_log_file).st_mtime
        current_seconds = time.time()
        diff_seconds = current_seconds - mod_time_seconds
        if diff_seconds < min_interval_seconds:
            return

    with open(resource_log_file, 'at') as f_handle:
        f_handle.write('{}\n'.format(resource_usage))


def main():
    args = parse_args()
    job_id = args.job_id
    command = cluster_commands.status_command(job_id)
    stdout = run_command(command)
    job_info = extract_job_info(stdout, job_id)
    status = job_info['status']
    update_resource_log(status, job_info['resource_usage'],
                        args.resource_usage_dir,
                        args.resource_usage_min_interval, job_id)
    print(status)


if __name__ == '__main__':
    main()
