def try_parse_float(string):
    try:
        return float(string), None
    except ValueError as e:
        return None, 'try_parse_float({}): {}'.format(string, e)


def try_parse_int(string):
    try:
        return int(string), None
    except ValueError as e:
        return None, 'try_parse_int({}): {}'.format(string, e)


def _expected_abundance_initial_headers():
    return ['transcript_ID', 'transcript_name', 'gene_ID']


def _parse_abundance_columns_from_line(line):
    expected_initial_headers = _expected_abundance_initial_headers()
    columns = line.strip().split('\t')
    initial_columns = columns[:len(expected_initial_headers)]
    sample_columns = columns[len(expected_initial_headers):]
    return initial_columns, sample_columns


def read_abundance_esp_header_line(line):
    expected_initial_headers = _expected_abundance_initial_headers()
    initial_headers, sample_headers = _parse_abundance_columns_from_line(line)
    if initial_headers != expected_initial_headers:
        raise Exception('expected abundance headers to start with: {}, but'
                        ' found: {}'.format(expected_initial_headers, line))

    return initial_headers, sample_headers


def read_abundance_esp_line(line, sample_headers):
    expected_initial_headers = _expected_abundance_initial_headers()
    values = dict()
    initial_values, sample_values = _parse_abundance_columns_from_line(line)
    for i, initial_value in enumerate(initial_values):
        values[expected_initial_headers[i]] = initial_value

    for i, sample_value_str in enumerate(sample_values):
        sample_value, error = try_parse_float(sample_value_str)
        if error:
            raise Exception(
                'could not parse sample value as float: line: {}, error: {}'.
                format(line, error))

        values[sample_headers[i]] = sample_value

    return values


def get_samples_from_abundance_esp(esp_path):
    with open(esp_path, 'rt') as in_f:
        line = in_f.readline()
        initial_columns, sample_columns = read_abundance_esp_header_line(line)
        return sample_columns


def read_chrom_sizes(chrom_sizes_path):
    chrom_names_to_sizes = dict()
    with open(chrom_sizes_path, 'rt') as in_f:
        for line in in_f:
            columns = line.strip().split('\t')
            chrom_name = columns[0]
            chrom_size, error = try_parse_int(columns[1])
            if error:
                raise Exception(
                    'could not parse chrom size as int: line: {}, error: {}'.
                    format(line, error))

            chrom_names_to_sizes[chrom_name] = chrom_size

    return chrom_names_to_sizes


def read_track_line(line):
    tokens = line.strip().split(' ')
    if tokens[0] != 'track':
        raise Exception(
            'track line should start with "track": {}'.format(line))

    options = dict()
    for token in tokens[1:]:
        split_token = token.split('=')
        if len(split_token) != 2:
            raise Exception(
                'track line attributes should be in the format key=value: {}'.
                format(line))

        key, value = split_token
        options[key] = value

    return options


def _parse_int_or_error(name, str_value, line):
    int_value, error = try_parse_int(str_value)
    if error:
        raise Exception('could not parse {} as int: line {}, error: {}'.format(
            name, line, error))

    return int_value


def _parse_list_of_ints_or_error(name, list_str, line):
    values = list()
    value_strs = list_str.split(',')
    for value_str in value_strs:
        value_int, error = try_parse_int(value_str)
        if error:
            raise Exception(
                'could not parse {} as int: line: {}, error: {}'.format(
                    name, line, error))

        values.append(value_int)

    return values


def read_bed_line(line):
    values = dict()
    columns = line.strip().split('\t')
    chrom = columns[0]
    chrom_start_str = columns[1]
    chrom_end_str = columns[2]
    name = columns[3]
    score_str = columns[4]
    strand = columns[5]
    thick_start_str = columns[6]
    thick_end_str = columns[7]
    item_rgb = columns[8]
    block_count_str = columns[9]
    block_sizes_str = columns[10]
    block_starts_str = columns[11]

    values['chrom'] = chrom
    chrom_start = _parse_int_or_error('chrom_start', chrom_start_str, line)
    values['chrom_start'] = chrom_start
    chrom_end = _parse_int_or_error('chrom_end', chrom_end_str, line)
    values['chrom_end'] = chrom_end
    values['name'] = name
    score = _parse_int_or_error('score', score_str, line)
    values['score'] = score
    values['strand'] = strand
    thick_start = _parse_int_or_error('thick_start', thick_start_str, line)
    values['thick_start'] = thick_start
    thick_end = _parse_int_or_error('thick_end', thick_end_str, line)
    values['thick_end'] = thick_end
    values['item_rgb'] = item_rgb
    block_count = _parse_int_or_error('block_count', block_count_str, line)
    values['block_count'] = block_count

    block_sizes = _parse_list_of_ints_or_error('block_size', block_sizes_str,
                                               line)
    block_starts = _parse_list_of_ints_or_error('block_start',
                                                block_starts_str, line)
    if len(block_sizes) != block_count or len(block_starts) != block_count:
        raise Exception(
            'block sizes: {} and starts: {} do not match count: {}'.format(
                len(block_sizes), len(block_starts), block_count))

    values['block_sizes'] = block_sizes
    values['block_starts'] = block_starts
    block_regions = list()
    for i, relative_block_start in enumerate(block_starts):
        block_size = block_sizes[i]
        block_start = values['chrom_start'] + relative_block_start
        block_end = block_start + block_size
        block_regions.append((block_start, block_end))

    values['block_regions'] = block_regions
    return values
