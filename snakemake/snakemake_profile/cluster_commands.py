import os
import os.path


def submit_command(log_out, log_err, threads, time_hours, mem_mb,
                   mem_mb_per_thread, jobscript):
    # sbatch requires that the directories for the log files already exist
    for log_path in [log_out, log_err]:
        log_dir = os.path.dirname(log_path)
        if log_dir != '':
            os.makedirs(log_dir, exist_ok=True)

    command = ['sbatch', '-o', log_out, '-e', log_err]
    if threads:
        command.extend(['-c', str(threads)])

    if time_hours:
        days, hours_float = divmod(time_hours, 24)
        hours_whole, hours_part = divmod(hours_float, 1)
        minutes = hours_part * 60
        time_str = '{}-{}:{}'.format(int(days), int(hours_whole), int(minutes))
        command.extend(['--time', time_str])

    if mem_mb:
        command.append('--mem={}M'.format(mem_mb))

    command.append(jobscript)
    return command


def try_extract_job_id_from_submit_output(stdout):
    tokens = stdout.split()
    if len(tokens) < 4:
        return None, 'expected at least 4 tokens'

    if (((tokens[0] != 'Submitted') or (tokens[1] != 'batch')
         or (tokens[2] != 'job'))):
        return None, 'expected output to look like "Submitted batch job ..."'

    try:
        job_id = int(tokens[3])
    except ValueError as e:
        return None, 'could not parse {} as an int: {}'.format(tokens[3], e)

    return job_id, None


def status_command(job_id):
    status_fields = [
        'ElapsedRaw',
        'End',
        'ExitCode',
        'JobIDRaw',
        'MaxDiskRead',
        'MaxDiskWrite',
        'MaxRss',
        'MaxVMSize',
        'Start',
        'State',
        'Submit',
        'TotalCPU',
    ]
    return [
        'sacct', '--parsable', '-j', job_id,
        '--format={}'.format(','.join(status_fields))
    ]


def try_extract_job_info_from_status_output(stdout, job_id):
    rows, error = _parse_rows(stdout)
    if error:
        return None, error

    if not rows:
        # the output may be empty if the job was submitted very recently
        return {'status': 'running', 'resource_usage': None}, None

    parent_rows = list()
    batch_rows = list()
    other_rows = list()
    for row in rows:
        row_job_id, row_job_step = _get_job_id_and_step(row.get('JobIDRaw'))
        if row_job_id != job_id:
            continue

        if row_job_step is None:
            parent_rows.append(row)
        elif row_job_step == 'batch':
            batch_rows.append(row)
        else:
            other_rows.append(row)

    usage = {
        'cpu': None,
        'end_time': None,
        'exit_code': None,
        'exit_signal': None,
        'max_disk_read': None,
        'max_disk_write': None,
        'max_rss': None,
        'max_vmem': None,
        'start_time': None,
        'state': None,
        'submit_time': None,
        'wallclock': None,
    }
    status, error = _update_from_parent_rows(parent_rows, usage)
    if error:
        return None, error

    status, error = _update_from_batch_rows(batch_rows, status, usage)
    if error:
        return None, error

    status, error = _update_from_other_rows(other_rows, status, usage)
    if error:
        return None, error

    if status is None:
        return None, 'no status found'

    resource_usage = ('cpu: {cpu},'
                      ' end_time: {end_time},'
                      ' exit_code: {exit_code},'
                      ' exit_signal: {exit_signal},'
                      ' max_disk_read: {max_disk_read},'
                      ' max_disk_write: {max_disk_write},'
                      ' max_rss: {max_rss},'
                      ' max_vmem: {max_vmem}'
                      ' start_time: {start_time},'
                      ' state: {state},'
                      ' submit_time: {submit_time},'
                      ' wallclock: {wallclock}'.format(**usage))
    return {'status': status, 'resource_usage': resource_usage}, None


def _parse_rows(stdout):
    lines = stdout.splitlines()
    if not lines:
        return list(), None

    header = lines[0]
    header_cols = header.split('|')
    rows = list()
    for i, line in enumerate(lines[1:]):
        row_cols = line.split('|')
        if len(header_cols) != len(row_cols):
            return None, 'row {} had {} columns but expected {}'.format(
                i, len(row_cols), len(header_cols))

        row = dict(zip(header_cols, row_cols))
        rows.append(row)

    return rows, None


def _get_job_id_and_step(job_id_raw):
    job_id_sep_index = job_id_raw.find('.')
    if job_id_sep_index <= 0:
        return job_id_raw, None

    job_id_base = job_id_raw[:job_id_sep_index]
    job_id_step = job_id_raw[job_id_sep_index + 1:]
    return job_id_base, job_id_step


def _update_from_parent_rows(rows, usage):
    if not rows:
        return None, None

    if len(rows) > 1:
        return None, 'expected at most 1 parent row'

    row = rows[0]
    parsed_values = _parse_values(row)
    # parent row is handled first.
    # Add starting values which may be overwritten later
    _add_if_not_none_keys(parsed_values, usage, [
        'cpu', 'exit_code', 'exit_signal', 'wallclock', 'submit_time',
        'start_time', 'end_time', 'state', 'max_disk_read', 'max_disk_write',
        'max_rss', 'max_vmem'
    ])

    parsed_status = parsed_values.get('state_for_snakemake')
    return parsed_status, None


def _update_from_batch_rows(rows, status, usage):
    if not rows:
        return status, None

    if len(rows) > 1:
        return None, 'expected at most 1 batch row'

    row = rows[0]
    parsed_values = _parse_values(row)
    # the batch row seems to have more details for these fields
    _overwrite_if_not_none_keys(parsed_values, usage, [
        'cpu', 'exit_code', 'exit_signal', 'max_disk_read', 'max_disk_write',
        'max_rss', 'max_vmem'
    ])
    # prefer the info from the parent row for these fields
    _add_if_not_none_keys(
        parsed_values, usage,
        ['wallclock', 'submit_time', 'start_time', 'end_time', 'state'])

    # prefer the parent status
    if status is None:
        status = parsed_values.get('state_for_snakemake')

    return status, None


def _update_from_other_rows(rows, status, usage):
    for row in rows:
        parsed_values = _parse_values(row)
        # use the "other" rows to fill in missing information
        _add_if_not_none_keys(parsed_values, usage, [
            'cpu', 'exit_code', 'exit_signal', 'wallclock', 'submit_time',
            'start_time', 'end_time', 'state', 'max_disk_read',
            'max_disk_write', 'max_rss', 'max_vmem'
        ])

        if status is None:
            status = parsed_values.get('state_for_snakemake')

    return status, None


def _parse_values(row):
    values = dict()
    values['cpu'] = _parse_cpu_time(row)
    exit_code, exit_signal = _parse_exit_code(row)
    values['exit_code'] = exit_code
    values['exit_signal'] = exit_signal
    values['wallclock'] = _parse_wallclock(row)
    values['submit_time'] = _parse_submit_time(row)
    values['start_time'] = _parse_start_time(row)
    values['end_time'] = _parse_end_time(row)
    raw_state, state_for_snakemake = _parse_state(row)
    values['state'] = raw_state
    values['state_for_snakemake'] = state_for_snakemake
    values['max_disk_read'] = _parse_max_disk_read(row)
    values['max_disk_write'] = _parse_max_disk_write(row)
    values['max_rss'] = _parse_max_rss(row)
    values['max_vmem'] = _parse_max_vmem(row)
    return values


def _parse_state(row):
    raw = row.get('State')
    if not raw:
        return None, None

    # translate the slurm state into snakemake terms
    for_snakemake = None
    if ((raw.startswith('RUNNING') or raw.startswith('PENDING')
         or raw.startswith('REQUEUED') or raw.startswith('RESIZING')
         or raw.startswith('SUSPENDED'))):
        for_snakemake = 'running'
    elif raw.startswith('COMPLETED'):
        for_snakemake = 'success'
    else:
        for_snakemake = 'failed'

    return raw, for_snakemake


def _parse_submit_time(row):
    return _parse_datetime_col(row, 'Submit')


def _parse_start_time(row):
    return _parse_datetime_col(row, 'Start')


def _parse_end_time(row):
    return _parse_datetime_col(row, 'End')


def _parse_datetime_col(row, col):
    # yyyy-mm-ddThh:mm:ss
    raw = row.get(col)
    if not raw:
        return None

    return raw


def _parse_cpu_time(row):
    # 'mm:ss.millis'
    raw = row.get('TotalCPU')
    if not raw:
        return None

    return raw


def _parse_wallclock(row):
    # 'num_seconds'
    raw = row.get('ElapsedRaw')
    if not raw:
        return None

    return raw


def _parse_max_disk_read(row):
    return _parse_disk(row, 'MaxDiskRead')


def _parse_max_disk_write(row):
    return _parse_disk(row, 'MaxDiskWrite')


def _parse_disk(row, col):
    # '{float}M'
    raw = row.get(col)
    if not raw:
        return None

    return raw


def _parse_max_rss(row):
    return _parse_mem(row, 'MaxRss')


def _parse_max_vmem(row):
    return _parse_mem(row, 'MaxVMSize')


def _parse_mem(row, col):
    # '{int}K'
    raw = row.get(col)
    if not raw:
        return None

    return raw


def _parse_exit_code(row):
    # 'exitcode:signal_num'
    raw = row.get('ExitCode')
    if not raw:
        return None, None

    splits = raw.split(':')
    if len(splits) != 2:
        return raw, None

    return splits[0], splits[1]


def _overwrite_if_not_none_keys(source, dest, keys):
    for key in keys:
        value = source.get(key)
        _overwrite_if_not_none(key, value, dest)


def _overwrite_if_not_none(key, value, dest):
    if value is not None:
        dest[key] = value


def _add_if_not_none_keys(source, dest, keys):
    for key in keys:
        value = source.get(key)
        _add_if_not_none(key, value, dest)


def _add_if_not_none(key, value, dest):
    if value is not None and dest.get(key) is None:
        dest[key] = value
