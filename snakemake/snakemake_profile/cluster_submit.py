import argparse
import math
import os
import os.path
import sys

from snakemake.utils import read_job_properties

import cluster_commands
import try_command


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('jobscript',
                        help='the script to be executed on the cluster')
    parser.add_argument(
        '--retry-submit-interval-seconds',
        default='',
        help='a "," separated list of integers representing'
        ' the number of seconds to wait after sequential failed'
        ' job submission commands before retrying')
    parser.add_argument(
        '--resource-usage-dir',
        help='a directory for storing the file paths where cluster_status.py'
        ' should log the resource usage of each job')
    args = parser.parse_args()
    jobscript = args.jobscript

    job_properties = read_job_properties(jobscript)

    retry_submit_interval_seconds = list()
    for int_str in args.retry_submit_interval_seconds.split(','):
        retry_submit_interval_seconds.append(int(int_str))

    resource_usage_dir = args.resource_usage_dir
    if resource_usage_dir:
        os.makedirs(resource_usage_dir, exist_ok=True)

    return {
        'jobscript': jobscript,
        'job_properties': job_properties,
        'retry_submit_interval_seconds': retry_submit_interval_seconds,
        'resource_usage_dir': resource_usage_dir,
    }


def get_base_path_from_jobscript(jobscript):
    with open(jobscript) as f_handle:
        for line in f_handle:
            tokens = line.split()
            if len(tokens) == 4:
                if ((tokens[0] == 'cd' and tokens[2] == '&&'
                     and tokens[3] == '\\')):
                    base_path = tokens[1]
                    if os.path.isdir(base_path):
                        return os.path.abspath(base_path)

    return None


def get_cluster_log_paths(jobscript, job_properties):
    cluster_logs = {'out': os.devnull, 'err': os.devnull, 'usage': os.devnull}
    base_path = get_base_path_from_jobscript(jobscript)
    orig_logs = job_properties.get('log')
    if (not base_path) or (not orig_logs):
        return cluster_logs

    if len(orig_logs) == 1:
        orig_log_out = orig_logs[0]
        orig_log_err = orig_log_out
    else:
        orig_log_out = orig_logs[0]
        orig_log_err = orig_logs[1]

    cluster_logs['out'] = os.path.join(base_path,
                                       '{}.cluster.out'.format(orig_log_out))
    cluster_logs['err'] = os.path.join(base_path,
                                       '{}.cluster.err'.format(orig_log_err))
    cluster_logs['usage'] = os.path.join(
        base_path, '{}.cluster.usage'.format(orig_log_out))
    return cluster_logs


def build_submit_command(jobscript, job_properties, cluster_log_out,
                         cluster_log_err):
    threads = job_properties.get('threads')
    resources = job_properties.get('resources')
    time_hours = None
    mem_mb = None
    mem_mb_per_thread = None
    if resources:
        time_hours = resources.get('time_hours')
        gpus = resources.get('gpus')
        gpu_name = resources.get('gpu_name')
        mem_mb = resources.get('mem_mb')
        mem_mb_per_thread = mem_mb
        if mem_mb and threads:
            mem_mb_per_thread /= float(threads)
            mem_mb_per_thread = math.ceil(mem_mb_per_thread)

    return cluster_commands.submit_command(cluster_log_out, cluster_log_err,
                                           threads, time_hours, mem_mb,
                                           mem_mb_per_thread, gpus, gpu_name,
                                           jobscript)


def run_submit_command(command, retry_submit_interval_seconds):
    stdout, error = try_command.try_command(command,
                                            retry_submit_interval_seconds)
    if error:
        sys.exit(error)

    return stdout


def extract_job_id(stdout):
    job_id, error = cluster_commands.try_extract_job_id_from_submit_output(
        stdout)
    if error:
        print('error: {}\n{}'.format(error, stdout), file=sys.stderr)
        sys.exit(1)

    return job_id


def record_usage_file(job_id, cluster_log_usage, resource_usage_dir):
    job_file_path = os.path.join(resource_usage_dir, '{}.txt'.format(job_id))
    with open(job_file_path, 'wt') as f_handle:
        f_handle.write('{}\n'.format(cluster_log_usage))


def main():
    parsed_args = parse_args()
    jobscript = parsed_args['jobscript']
    job_properties = parsed_args['job_properties']
    cluster_logs = get_cluster_log_paths(jobscript, job_properties)
    command = build_submit_command(jobscript, job_properties,
                                   cluster_logs['out'], cluster_logs['err'])
    stdout = run_submit_command(command,
                                parsed_args['retry_submit_interval_seconds'])
    job_id = extract_job_id(stdout)
    record_usage_file(job_id, cluster_logs['usage'],
                      parsed_args['resource_usage_dir'])
    print(job_id)


if __name__ == '__main__':
    main()
