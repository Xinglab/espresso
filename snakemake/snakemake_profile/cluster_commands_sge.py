import bs4


def submit_command(log_out, log_err, threads, time_hours, mem_mb,
                   mem_mb_per_thread, gpus, gpu_name, jobscript):
    command = ['qsub', '-o', log_out, '-e', log_err]
    if threads:
        command.extend(['-pe', 'smp', str(threads)])

    if mem_mb_per_thread:
        command.extend(['-l', 'h_vmem={}M'.format(mem_mb_per_thread)])

    command.append(jobscript)
    return command


def try_extract_job_id_from_submit_output(stdout):
    tokens = stdout.split()
    if len(tokens) < 3:
        return None, 'expected at least 3 tokens'

    if (tokens[0] != 'Your') or (tokens[1] != 'job'):
        return None, 'expected output to look like "Your job ..."'

    try:
        job_id = int(tokens[2])
    except ValueError as e:
        return None, 'could not parse {} as an int: {}'.format(tokens[2], e)

    return job_id, None


def status_command(job_id):
    return ['qstat', '-j', job_id, '-xml']


def try_extract_job_info_from_status_output(stdout, job_id, max_job_days):
    info = {'status': None, 'resource_usage': None}
    soup = bs4.BeautifulSoup(stdout, 'html.parser')

    unks = soup.find_all('unknown_jobs')
    if unks:
        names = unks[0].find_all('st_name')
        if names:
            if names[0].string == job_id:
                info['status'] = 'success'
                return info, None

    djobs = soup.find_all('djob_info')
    if djobs:
        numbers = djobs[0].find_all('jb_job_number')
        if numbers:
            if numbers[0].string == job_id:
                info['status'] = 'running'
                info['resource_usage'] = _extract_resource_usage(djobs[0])
                return info, None

    return None, 'unexpected output'


def _extract_resource_usage(soup):
    resources = _extract_resource_usage_components(soup)
    return ('wallclock: {wallclock}, cpu: {cpu}, io_wait: {io_wait},'
            ' max_vmem: {max_vmem}, max_rss: {max_rss}, vmem: {vmem},'
            ' rss: {rss}'.format(**resources))


def _extract_resource_usage_components(soup):
    resources = {
        'wallclock': None,
        'cpu': None,
        'io_wait': None,
        'vmem': None,
        'max_vmem': None,
        'rss': None,
        'max_rss': None,
    }
    usage_list = soup.find_all('jat_scaled_usage_list')
    if not usage_list:
        return resources

    wallclock_float = _extract_ua_by_name(usage_list[0], 'wallclock')
    if wallclock_float:
        resources['wallclock'] = '{:.2f}s'.format(wallclock_float)

    cpu_float = _extract_ua_by_name(usage_list[0], 'cpu')
    if cpu_float:
        resources['cpu'] = '{:.2f}s'.format(cpu_float)

    io_wait_float = _extract_ua_by_name(usage_list[0], 'iow')
    if io_wait_float:
        resources['io_wait'] = '{:.2f}s'.format(io_wait_float)

    bytes_per_gb = 1024**3
    vmem_float = _extract_ua_by_name(usage_list[0], 'vmem')
    if vmem_float:
        resources['vmem'] = '{:.2f}GB'.format(vmem_float / bytes_per_gb)

    max_vmem_float = _extract_ua_by_name(usage_list[0], 'maxvmem')
    if max_vmem_float:
        resources['max_vmem'] = '{:.2f}GB'.format(max_vmem_float /
                                                  bytes_per_gb)

    rss_float = _extract_ua_by_name(usage_list[0], 'rss')
    if rss_float:
        resources['rss'] = '{:.2f}GB'.format(rss_float / bytes_per_gb)

    max_rss_float = _extract_ua_by_name(usage_list[0], 'maxrss')
    if max_rss_float:
        resources['max_rss'] = '{:.2f}GB'.format(max_rss_float / bytes_per_gb)

    return resources


def _extract_ua_by_name(soup, ua_name):
    name_node = soup.find_all('ua_name', string=ua_name)
    if not name_node:
        return None

    value_node = name_node[0].find_next_siblings('ua_value')
    if not value_node:
        return None

    value_str = value_node[0].string
    return _try_parse_float(value_str)


def _try_parse_float(s):
    try:
        return float(s)
    except ValueError:
        return None
