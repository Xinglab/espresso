import argparse
import datetime
import os.path
import sys
import time

import cluster_commands
import try_command


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'job_id', help='the cluster id of the job to check the status of')
    parser.add_argument(
        '--retry-status-interval-seconds',
        default='',
        help='a "," separated list of integers representing'
        ' the number of seconds to wait after sequential failed'
        ' job status commands before retrying')
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
    parser.add_argument(
        '--max-job-days',
        type=int,
        default=60,
        help=('cluster job IDs can start being reused after reaching the'
              ' max id. If a job is older than --max-job-days assume that the'
              ' time is actually from an old job ID'))
    args = parser.parse_args()

    retry_status_interval_seconds = list()
    for int_str in args.retry_status_interval_seconds.split(','):
        retry_status_interval_seconds.append(int(int_str))

    return {
        'job_id': args.job_id,
        'retry_status_interval_seconds': retry_status_interval_seconds,
        'resource_usage_dir': args.resource_usage_dir,
        'resource_usage_min_interval': args.resource_usage_min_interval,
        'max_job_days': args.max_job_days,
    }


def extract_job_info(stdout, job_id, max_job_days):
    info, error = cluster_commands.try_extract_job_info_from_status_output(
        stdout, job_id, max_job_days)
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


def run_status_command(command, retry_status_interval_seconds):
    stdout, error = try_command.try_command(command,
                                            retry_status_interval_seconds)
    if error:
        sys.exit(error)

    return stdout


def main():
    args = parse_args()
    job_id = args['job_id']
    command = cluster_commands.status_command(job_id)
    stdout = run_status_command(command, args['retry_status_interval_seconds'])
    job_info = extract_job_info(stdout, job_id, args['max_job_days'])
    status = job_info['status']
    update_resource_log(status, job_info['resource_usage'],
                        args['resource_usage_dir'],
                        args['resource_usage_min_interval'], job_id)
    print(status)


if __name__ == '__main__':
    main()
