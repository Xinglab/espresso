import argparse
import os.path
import subprocess
import tempfile


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Compare ESPRESSO output files'))
    parser.add_argument(
        '--abun-a',
        required=True,
        help='the abundance.esp from the 1st set of output files')
    parser.add_argument(
        '--abun-b',
        required=True,
        help='the abundance.esp from the 2nd set of output files')
    parser.add_argument(
        '--gtf-a',
        required=True,
        help='the updated.gtf from the 1st set of output files')
    parser.add_argument(
        '--gtf-b',
        required=True,
        help='the updated.gtf from the 2nd set of output files')
    parser.add_argument(
        '--compat-a',
        required=False,
        help='the isoform.tsv from the 1st set of output files')
    parser.add_argument(
        '--compat-b',
        required=False,
        help='the isoform.tsv from the 2nd set of output files')
    parser.add_argument(
        '--sort-memory-buffer-size',
        default='2G',
        help='how much memory "sort" is allowed to use. Specified as {num_gb}G'
    )
    parser.add_argument(
        '--float-diff-threshold',
        type=float,
        default=0.1,
        help='how much a float value must differ to be reported')
    parser.add_argument('--tmp-dir',
                        required=True,
                        help='where to write temporary files')

    return parser.parse_args()


def sort_file(orig, output, sort_memory_buffer_size, temp_dir):
    sort_command = [
        'sort',
        orig,
        '--output',
        output,
        '--temporary-directory',
        temp_dir,
        '--buffer-size',
        sort_memory_buffer_size,
    ]
    env = {'LC_ALL': 'C'}
    subprocess.run(sort_command, env=env, check=True)


def compare_abundance(abun_a, abun_b, gtf_a, gtf_b, float_diff_threshold):
    with open(abun_a, 'rt') as abun_a_handle:
        with open(abun_b, 'rt') as abun_b_handle:
            with open(gtf_a, 'rt') as gtf_a_handle:
                with open(gtf_b, 'rt') as gtf_b_handle:
                    compare_abundance_with_handles(abun_a_handle,
                                                   abun_b_handle, gtf_a_handle,
                                                   gtf_b_handle,
                                                   float_diff_threshold)


def compare_matched_parsed_abundance_lines(parsed_a, parsed_b,
                                           float_diff_threshold):
    for i, a_val in enumerate(parsed_a['sample_values']):
        b_val = parsed_b['sample_values'][i]
        diff = a_val - b_val
        if abs(diff) > float_diff_threshold:
            print('transcript: {}, value_i: {}, a: {}, b: {}'.format(
                parsed_a['transcript_id'], i, a_val, b_val))

    if len(parsed_a['sample_values']) != len(parsed_b['sample_values']):
        print('transcript: {}, a_values: {}, b_values: {}'.format(
            parsed_a['transcript_id'], parsed_a['sample_values'],
            parsed_b['sample_values']))


def compare_unmatched_parsed_abundance_line(parsed, is_a, float_diff_threshold):
    for i, val in enumerate(parsed['sample_values']):
        if is_a:
            a_val = val
            b_val = 0
        else:
            a_val = 0
            b_val = val

        diff = a_val - b_val
        if abs(diff) > float_diff_threshold:
            print('transcript: {}, value_i: {}, a: {}, b: {}'.format(
                parsed['transcript_id'], i, a_val, b_val))


def read_novel_defs(gtf_handle, novel_chr):
    novel_defs = dict()
    offset = gtf_handle.tell()
    while True:
        line = gtf_handle.readline()
        columns = line.strip().split('\t')
        chr_name = columns[0]
        if chr_name != novel_chr:
            gtf_handle.seek(offset)
            break

        offset = gtf_handle.tell()
        strand = columns[1]
        exons_str = columns[2]
        is_novel = columns[5] == 'True'
        transcript_id = columns[6]

        if not is_novel:
            continue

        exon_strs = exons_str.split(';')
        exon_coords = list()
        for exon_str in exon_strs:
            coord_strs = exon_str.split(':')
            start = int(coord_strs[0])
            end = int(coord_strs[1])
            exon_coords.append((start, end))

        novel_defs[transcript_id] = {'exons': exon_coords, 'strand': strand}

    return novel_defs


def coords_from_definition(definition):
    coords = [definition['strand']]
    exons = definition['exons']
    if len(exons) == 1:
        # include the start and end for a single exon transcript
        coords.extend(exons[0])
        return tuple(coords)

    for i, exon in enumerate(exons):
        if i != 0:
            # do not include the start of the first exon
            coords.append(exon[0])
        if i != len(exons) - 1:
            # do not include the end of the last exon
            coords.append(exon[1])

    return tuple(coords)


def make_abuns_by_coords(abuns, defs):
    abuns_by_coords = dict()
    for abun in abuns:
        definition = defs.get(abun['transcript_id'])
        if not definition:
            raise Exception('No definition found for {}'.format(
                abun['transcript_id']))

        coords = coords_from_definition(definition)
        abuns_by_coords[coords] = abun

    return abuns_by_coords


def compare_novels(novel_defs_a, novel_defs_b, novel_abuns_a, novel_abuns_b,
                   float_diff_threshold):
    abuns_by_coords_a = make_abuns_by_coords(novel_abuns_a, novel_defs_a)
    abuns_by_coords_b = make_abuns_by_coords(novel_abuns_b, novel_defs_b)
    for coords_a, abun_a in abuns_by_coords_a.items():
        abun_b = abuns_by_coords_b.get(coords_a)
        if not abun_b:
            print('extra a transcript: {}'.format(abun_a['transcript_id']))
            compare_unmatched_parsed_abundance_line(abun_a, True, float_diff_threshold)
        else:
            del abuns_by_coords_b[coords_a]
            compare_matched_parsed_abundance_lines(abun_a, abun_b,
                                                   float_diff_threshold)

    for abun_b in abuns_by_coords_b.values():
        print('extra b transcript: {}'.format(abun_b['transcript_id']))
        compare_unmatched_parsed_abundance_line(abun_b, False, float_diff_threshold)


def print_extra_novels(name, novel_abuns, float_diff_threshold):
    is_a = name == 'a'
    for abun in novel_abuns:
        print('extra {} transcript: {}'.format(name, abun['transcript_id']))
        compare_unmatched_parsed_abundance_line(abun, is_a, float_diff_threshold)


def compare_abundance_with_handles(abun_a_handle, abun_b_handle, gtf_a_handle,
                                   gtf_b_handle, float_diff_threshold):
    line_a, line_b = compare_annotated_abundance(abun_a_handle, abun_b_handle,
                                                 float_diff_threshold)
    if not (line_a or line_b):
        # processed all lines
        return

    if not line_a:
        while line_b:
            parsed_b = parse_abundance_line(line_b)
            print('extra b transcript: {}'.format(parsed_b['transcript_id']))
            compare_unmatched_parsed_abundance_line(parsed_b, False, float_diff_threshold)
            line_b = abun_b_handle.readline()

        return

    if not line_b:
        while line_a:
            parsed_a = parse_abundance_line(line_a)
            print('extra a transcript: {}'.format(parsed_a['transcript_id']))
            compare_unmatched_parsed_abundance_line(parsed_a, True, float_diff_threshold)
            line_a = abun_a_handle.readline()

        return

    compare_novel_abundance(line_a, line_b, abun_a_handle, abun_b_handle,
                            gtf_a_handle, gtf_b_handle, float_diff_threshold)


def compare_annotated_abundance(abun_a_handle, abun_b_handle,
                                float_diff_threshold):
    line_a = abun_a_handle.readline()
    line_b = abun_b_handle.readline()
    while line_a and line_b:
        parsed_a = parse_abundance_line(line_a)
        parsed_b = parse_abundance_line(line_b)
        if parsed_a['is_novel'] or parsed_b['is_novel']:
            if parsed_a['is_novel'] and parsed_b['is_novel']:
                break

            if parsed_a['is_novel']:
                print('extra b transcript: {}'.format(
                    parsed_b['transcript_id']))
                compare_unmatched_parsed_abundance_line(parsed_b, False, float_diff_threshold)
                line_b = abun_b_handle.readline()

            if parsed_b['is_novel']:
                print('extra a transcript: {}'.format(
                    parsed_a['transcript_id']))
                compare_unmatched_parsed_abundance_line(parsed_a, True, float_diff_threshold)
                line_a = abun_a_handle.readline()
        elif ((parsed_a['transcript_id'] == parsed_b['transcript_id'])
              and (parsed_a['transcript_name'] == parsed_b['transcript_name'])
              and (parsed_a['gene_id'] == parsed_b['gene_id'])):
            compare_matched_parsed_abundance_lines(parsed_a, parsed_b,
                                                   float_diff_threshold)
            line_a = abun_a_handle.readline()
            line_b = abun_b_handle.readline()
        else:
            id_a = parsed_a['transcript_id']
            id_b = parsed_b['transcript_id']
            if id_a < id_b:
                print('extra a transcript: {}'.format(id_a))
                compare_unmatched_parsed_abundance_line(parsed_a, True, float_diff_threshold)
                line_a = abun_a_handle.readline()
            elif id_b < id_a:
                print('extra b transcript: {}'.format(id_b))
                compare_unmatched_parsed_abundance_line(parsed_b, False, float_diff_threshold)
                line_b = abun_b_handle.readline()
            else:
                print('differing transcript info for id: {}'.format(id_a))
                line_a = abun_a_handle.readline()
                line_b = abun_b_handle.readline()

    return line_a, line_b


def compare_novel_abundance(line_a, line_b, abun_a_handle, abun_b_handle,
                            gtf_a_handle, gtf_b_handle, float_diff_threshold):
    novel_chr_a = None
    novel_chr_b = None
    novel_abuns_a = None
    novel_abuns_b = None
    novel_defs_a = None
    novel_defs_b = None
    while line_a or line_b:
        if line_a:
            parsed_a = parse_abundance_line(line_a)
            if not parsed_a['is_novel']:
                raise Exception(
                    'unexpected annotated transcript when processing novel transcripts: {}'
                    .format(line_a))

            if novel_chr_a is None:
                novel_chr_a = parsed_a['chr']
                novel_abuns_a = [parsed_a]
                line_a = abun_a_handle.readline()
            elif parsed_a['chr'] == novel_chr_a:
                novel_abuns_a.append(parsed_a)
                line_a = abun_a_handle.readline()
            elif not novel_defs_a:
                novel_defs_a = read_novel_defs(gtf_a_handle, novel_chr_a)
        elif novel_chr_a and not novel_defs_a:
            # done reading all a lines
            novel_defs_a = read_novel_defs(gtf_a_handle, novel_chr_a)

        if line_b:
            parsed_b = parse_abundance_line(line_b)
            if not parsed_b['is_novel']:
                raise Exception(
                    'unexpected annotated transcript when processing novel transcripts: {}'
                    .format(line_b))

            if novel_chr_b is None:
                novel_chr_b = parsed_b['chr']
                novel_abuns_b = [parsed_b]
                line_b = abun_b_handle.readline()
            elif parsed_b['chr'] == novel_chr_b:
                novel_abuns_b.append(parsed_b)
                line_b = abun_b_handle.readline()
            elif not novel_defs_b:
                novel_defs_b = read_novel_defs(gtf_b_handle, novel_chr_b)
        elif novel_chr_b and not novel_defs_b:
            # done reading all b lines
            novel_defs_b = read_novel_defs(gtf_b_handle, novel_chr_b)

        if novel_defs_a and novel_defs_b:
            if novel_chr_a == novel_chr_b:
                compare_novels(novel_defs_a, novel_defs_b, novel_abuns_a,
                               novel_abuns_b, float_diff_threshold)
                novel_chr_a = None
                novel_abuns_a = None
                novel_defs_a = None
                novel_chr_b = None
                novel_abuns_b = None
                novel_defs_b = None
            elif novel_chr_a < novel_chr_b:
                print_extra_novels('a', novel_abuns_a, float_diff_threshold)
                novel_chr_a = None
                novel_abuns_a = None
                novel_defs_a = None
            else:
                print_extra_novels('b', novel_abuns_b, float_diff_threshold)
                novel_chr_b = None
                novel_abuns_b = None
                novel_defs_b = None

    if novel_abuns_a and not novel_defs_a:
        novel_defs_a = read_novel_defs(gtf_a_handle, novel_chr_a)

    if novel_abuns_b and not novel_defs_b:
        novel_defs_b = read_novel_defs(gtf_b_handle, novel_chr_b)

    if novel_defs_a and novel_defs_b:
        compare_novels(novel_defs_a, novel_defs_b, novel_abuns_a,
                       novel_abuns_b, float_diff_threshold)
    elif novel_defs_a:
        print_extra_novels('a', novel_abuns_a, float_diff_threshold)
    elif novel_defs_b:
        print_extra_novels('b', novel_abuns_b, float_diff_threshold)


def try_parse_float(string):
    try:
        return float(string)
    except ValueError:
        return None


def parse_abundance_line(line):
    columns = line.strip().split('\t')
    transcript_id = columns[0]
    transcript_name = columns[1]
    gene_id = columns[2]
    sample_float_strs = columns[3:]
    sample_floats = list()
    for string in sample_float_strs:
        val = try_parse_float(string)
        sample_floats.append(val)

    is_novel = False
    chr_name = None
    if transcript_id.startswith('ESPRESSO'):
        is_novel = True
        id_splits = transcript_id.split(':')
        chr_name = id_splits[1]

    return {
        'transcript_id': transcript_id,
        'transcript_name': transcript_name,
        'gene_id': gene_id,
        'sample_values': sample_floats,
        'is_novel': is_novel,
        'chr': chr_name,
    }


def get_temp_file_name(temp_dir):
    with tempfile.NamedTemporaryFile(dir=temp_dir,
                                     delete=False) as temp_file_handle:
        return temp_file_handle.name


def normalize_gtf(gtf, temp_dir, sort_memory_buffer_size):
    temp_file_name = get_temp_file_name(temp_dir)
    with open(gtf, 'rt') as in_handle:
        with open(temp_file_name, 'wt') as out_handle:
            normalize_gtf_with_handles(in_handle, out_handle)

    sorted_temp_file_name = get_temp_file_name(temp_dir)
    sort_file(temp_file_name, sorted_temp_file_name, sort_memory_buffer_size,
              temp_dir)

    os.remove(temp_file_name)
    return sorted_temp_file_name


def normalize_gtf_with_handles(in_handle, out_handle):
    current_transcript = {'id': None}
    for in_line in in_handle:
        if in_line.startswith('#'):
            continue

        columns = in_line.strip().split('\t')
        chr_name = columns[0]
        isoform_type = columns[1]
        feature_type = columns[2]
        start_str = columns[3]
        start = int(start_str)
        end_str = columns[4]
        end = int(end_str)
        # score = columns[5]
        strand = columns[6]
        # frame = columns[7]
        attributes_str = columns[8]
        attributes = parse_gtf_attributes_str(attributes_str)
        transcript_id = attributes.get('transcript_id')
        if not transcript_id:
            raise Exception('no transcript_id in gtf line: {}'.format(in_line))

        if current_transcript['id'] != transcript_id:
            if current_transcript['id']:
                write_transcript_gtf_line(current_transcript, out_handle)

            current_transcript = {'id': transcript_id}

        if feature_type == 'transcript':
            current_transcript['chr'] = chr_name
            current_transcript['is_novel'] = isoform_type == 'novel_isoform'
            current_transcript['start'] = start
            current_transcript['end'] = end
            current_transcript['strand'] = strand
        elif feature_type == 'exon':
            exon = (start, end)
            exons = current_transcript.get('exons')
            if not exons:
                exons = [exon]
                current_transcript['exons'] = exons
            else:
                exons.append(exon)
        else:
            raise Exception('unknown feature type in gtf: {}'.format(in_line))

    if current_transcript['id']:
        write_transcript_gtf_line(current_transcript, out_handle)


def remove_quotes(string):
    if len(string) < 2:
        return string

    if string[0] != string[-1]:
        return string

    if string[0] in ['"', "'"]:
        return string[1:-1]

    return string


def write_transcript_gtf_line(transcript, handle):
    sorted_exons = sorted(transcript['exons'])
    exon_strs = list()
    for exon in sorted_exons:
        exon_str = '{}:{}'.format(exon[0], exon[1])
        exon_strs.append(exon_str)

    exons_str = ';'.join(exon_strs)
    columns = [
        transcript['chr'],
        transcript['strand'],
        exons_str,
        str(transcript['start']),
        str(transcript['end']),
        str(transcript['is_novel']),
        remove_quotes(transcript['id']),
    ]
    handle.write('{}\n'.format('\t'.join(columns)))


def parse_gtf_attributes_str(attributes_str):
    attributes = dict()
    parts = attributes_str.strip().split(';')
    for part in parts:
        if part == '':
            continue

        key_and_value = part.strip().split(' ')
        if len(key_and_value) != 2:
            raise Exception(
                'unable to parse gtf attributes: {}'.format(attributes_str))

        key, value = key_and_value
        attributes[key] = value

    return attributes


def split_annotated_and_novel_abun(in_handle, annotated_handle, novel_handle):
    for i, line in enumerate(in_handle):
        columns = line.strip().split('\t')
        if i == 0:
            expected_headers = ['transcript_ID', 'transcript_name', 'gene_ID']
            if columns[:3] != expected_headers:
                raise Exception(
                    'unexpected headers in abundance file: {}'.format(line))

            continue

        transcript_id = columns[0]
        if transcript_id.startswith('ESPRESSO:'):
            colon_splits = transcript_id.split(':')
            espresso = colon_splits[0]
            chr_name = colon_splits[1]
            rest = ':'.join(colon_splits[2:])
            new_columns = [espresso, chr_name, rest]
            new_columns.extend(columns[1:])
            novel_handle.write('{}\n'.format('\t'.join(new_columns)))
        else:
            annotated_handle.write('{}\n'.format('\t'.join(columns)))


def reformat_and_append_novel_abun(append_handle, novel_handle):
    for novel_line in novel_handle:
        columns = novel_line.strip().split('\t')
        new_columns = list()
        combined_column = ':'.join(columns[:3])
        new_columns.append(combined_column)
        new_columns.extend(columns[3:])
        append_handle.write('{}\n'.format('\t'.join(new_columns)))


def remove_header_and_sort_abundance(orig, temp_dir, sort_memory_buffer_size):
    temp_file_name_1 = get_temp_file_name(temp_dir)
    temp_file_name_2 = get_temp_file_name(temp_dir)
    with open(orig, 'rt') as in_handle:
        with open(temp_file_name_1, 'wt') as annotated_handle:
            with open(temp_file_name_2, 'wt') as novel_handle:
                split_annotated_and_novel_abun(in_handle, annotated_handle,
                                               novel_handle)

    sorted_temp_file_name_1 = get_temp_file_name(temp_dir)
    sorted_temp_file_name_2 = get_temp_file_name(temp_dir)
    sort_file(temp_file_name_1, sorted_temp_file_name_1,
              sort_memory_buffer_size, temp_dir)
    sort_file(temp_file_name_2, sorted_temp_file_name_2,
              sort_memory_buffer_size, temp_dir)
    with open(sorted_temp_file_name_1, 'at') as combined_handle:
        with open(sorted_temp_file_name_2) as novel_handle:
            reformat_and_append_novel_abun(combined_handle, novel_handle)

    os.remove(temp_file_name_1)
    os.remove(temp_file_name_2)
    os.remove(sorted_temp_file_name_2)
    return sorted_temp_file_name_1


def split_compat_isoforms_with_handles(in_handle, annotated_handle,
                                       novel_handle):
    for line in in_handle:
        columns = line.strip().split('\t')
        parsed = parse_compat_line(line)
        if parsed['is_novel']:
            new_columns = [parsed['chr_name']] + columns
            novel_handle.write('{}\n'.format('\t'.join(new_columns)))
        else:
            annotated_handle.write(line)


def split_and_sort_compat_isoforms(in_path, temp_dir, sort_memory_buffer_size):
    temp_annotated_path = get_temp_file_name(temp_dir)
    temp_novel_path = get_temp_file_name(temp_dir)
    with open(in_path, 'rt') as in_handle:
        with open(temp_annotated_path, 'wt') as annotated_handle:
            with open(temp_novel_path, 'wt') as novel_handle:
                split_compat_isoforms_with_handles(in_handle, annotated_handle,
                                                   novel_handle)

    annotated_path = get_temp_file_name(temp_dir)
    novel_path = get_temp_file_name(temp_dir)
    sort_file(temp_annotated_path, annotated_path, sort_memory_buffer_size,
              temp_dir)
    sort_file(temp_novel_path, novel_path, sort_memory_buffer_size, temp_dir)
    os.remove(temp_annotated_path)
    os.remove(temp_novel_path)
    return annotated_path, novel_path


def parse_novel_compat_line(line):
    columns = line.strip().split('\t')
    # chr_name = columns[0]
    orig_line = '\t'.join(columns[1:])
    return parse_compat_line(orig_line)


def parse_compat_line(line):
    columns = line.strip().split('\t')
    read_id = columns[0]
    sample_name = columns[1]
    read_classification = columns[2]
    compat_isoforms_str = columns[3]
    compat_isoforms = compat_isoforms_str.strip(',').split(',')
    compat_isoforms.sort()
    is_novel = False
    chr_name = None
    for isoform in compat_isoforms:
        if isoform.startswith('ESPRESSO'):
            is_novel = True
            colon_splits = isoform.split(':')
            chr_name = colon_splits[1]
            break

    return {
        'read_id': read_id,
        'sample_name': sample_name,
        'read_classification': read_classification,
        'compat_isoforms': compat_isoforms,
        'is_novel': is_novel,
        'chr_name': chr_name,
    }


def compare_annotated_compat(compat_a, compat_b):
    with open(compat_a, 'rt') as a_handle:
        with open(compat_b, 'rt') as b_handle:
            line_a = a_handle.readline()
            line_b = b_handle.readline()
            while line_a and line_b:
                parsed_a = parse_compat_line(line_a)
                parsed_b = parse_compat_line(line_b)
                if parsed_a['read_id'] < parsed_b['read_id']:
                    print('extra a read_id: {}'.format(parsed_a['read_id']))
                    line_a = a_handle.readline()
                elif parsed_b['read_id'] < parsed_a['read_id']:
                    print('extra b read_id: {}'.format(parsed_b['read_id']))
                    line_b = b_handle.readline()
                else:
                    if parsed_a != parsed_b:
                        print('compat isoforms diff: {}, {}'.format(
                            parsed_a, parsed_b))

                    line_a = a_handle.readline()
                    line_b = b_handle.readline()

            while line_a:
                parsed_a = parse_compat_line(line_a)
                print('extra a read_id: {}'.format(parsed_a['read_id']))
                line_a = a_handle.readline()

            while line_b:
                parsed_b = parse_compat_line(line_b)
                print('extra b read_id: {}'.format(parsed_b['read_id']))
                line_b = b_handle.readline()


def normalize_isoforms_by_definitions(parsed, novel_defs):
    isoform_defs = list()
    for isoform in parsed['compat_isoforms']:
        if isoform.startswith('ESPRESSO'):
            definition = novel_defs.get(isoform)
            coords = coords_from_definition(definition)
            coord_str = ','.join([str(x) for x in coords])
            isoform_defs.append(coord_str)
        else:
            isoform_defs.append(isoform)

    isoform_defs.sort()
    return isoform_defs


def compare_parsed_novel_compat(parsed_a, parsed_b, novel_defs_a,
                                novel_defs_b):

    for key in ['read_id', 'sample_name', 'read_classification']:
        if parsed_a[key] != parsed_b[key]:
            print('compat isoforms diff {}: {}, {}'.format(
                key, parsed_a, parsed_b))
            return

    a_isoform_defs = normalize_isoforms_by_definitions(parsed_a, novel_defs_a)
    b_isoform_defs = normalize_isoforms_by_definitions(parsed_b, novel_defs_b)
    if a_isoform_defs != b_isoform_defs:
        print('compat_isoforms diff: {}, {}, {}'.format(
            parsed_a['read_id'], a_isoform_defs, b_isoform_defs))


def compare_novel_compat_with_handles(compat_a_handle, compat_b_handle,
                                      gtf_a_handle, gtf_b_handle):
    novel_chr_a = None
    novel_chr_b = None
    novel_defs_a = None
    novel_defs_b = None
    line_a = compat_a_handle.readline()
    line_b = compat_b_handle.readline()
    while line_a and line_b:
        parsed_a = parse_novel_compat_line(line_a)
        parsed_b = parse_novel_compat_line(line_b)
        if novel_chr_a != parsed_a['chr_name']:
            novel_chr_a = parsed_a['chr_name']
            novel_defs_a = read_novel_defs(gtf_a_handle, novel_chr_a)
        if novel_chr_b != parsed_b['chr_name']:
            novel_chr_b = parsed_b['chr_name']
            novel_defs_b = read_novel_defs(gtf_b_handle, novel_chr_b)

        if parsed_a['read_id'] < parsed_b['read_id']:
            print('extra a read_id: {}'.format(parsed_a['read_id']))
            line_a = compat_a_handle.readline()
        elif parsed_b['read_id'] < parsed_a['read_id']:
            print('extra b read_id: {}'.format(parsed_b['read_id']))
            line_b = compat_b_handle.readline()
        else:
            compare_parsed_novel_compat(parsed_a, parsed_b, novel_defs_a,
                                        novel_defs_b)
            line_a = compat_a_handle.readline()
            line_b = compat_b_handle.readline()

    while line_a:
        print('extra a read_id: {}'.format(parsed_a['read_id']))
        line_a = compat_a_handle.readline()

    while line_b:
        print('extra b read_id: {}'.format(parsed_b['read_id']))
        line_b = compat_b_handle.readline()


def compare_novel_compat(compat_a, compat_b, normed_gtf_a, normed_gtf_b):
    with open(compat_a, 'rt') as compat_a_handle:
        with open(compat_b, 'rt') as compat_b_handle:
            with open(normed_gtf_a, 'rt') as gtf_a_handle:
                with open(normed_gtf_b, 'rt') as gtf_b_handle:
                    compare_novel_compat_with_handles(compat_a_handle,
                                                      compat_b_handle,
                                                      gtf_a_handle,
                                                      gtf_b_handle)


def compare_compatible_isoforms(compat_a, compat_b, normed_gtf_a, normed_gtf_b,
                                temp_dir, sort_memory_buffer_size):
    annotated_compat_a, novel_compat_a = split_and_sort_compat_isoforms(
        compat_a, temp_dir, sort_memory_buffer_size)
    annotated_compat_b, novel_compat_b = split_and_sort_compat_isoforms(
        compat_b, temp_dir, sort_memory_buffer_size)

    compare_annotated_compat(annotated_compat_a, annotated_compat_b)
    compare_novel_compat(novel_compat_a, novel_compat_b, normed_gtf_a,
                         normed_gtf_b)


def compare_espresso_output(args, temp_dir):
    sorted_abun_a = remove_header_and_sort_abundance(
        args.abun_a, temp_dir, args.sort_memory_buffer_size)
    sorted_abun_b = remove_header_and_sort_abundance(
        args.abun_b, temp_dir, args.sort_memory_buffer_size)

    normed_gtf_a = normalize_gtf(args.gtf_a, temp_dir,
                                 args.sort_memory_buffer_size)
    normed_gtf_b = normalize_gtf(args.gtf_b, temp_dir,
                                 args.sort_memory_buffer_size)

    compare_abundance(sorted_abun_a, sorted_abun_b, normed_gtf_a, normed_gtf_b,
                      args.float_diff_threshold)

    if args.compat_a or args.compat_b:
        if not (args.compat_a and args.compat_b):
            raise Exception(
                'Only one of --compat-a or --compat-b was provided')

        compare_compatible_isoforms(args.compat_a, args.compat_b, normed_gtf_a,
                                    normed_gtf_b, temp_dir,
                                    args.sort_memory_buffer_size)


def main():
    args = parse_args()
    if os.path.exists(args.tmp_dir):
        if not os.path.isdir(args.tmp_dir):
            raise Exception('{} exists and is not a directory'.format(
                args.tmp_dir))
    else:
        os.makedirs(args.tmp_dir)

    with tempfile.TemporaryDirectory(dir=args.tmp_dir) as temp_dir:
        compare_espresso_output(args, temp_dir)

    print('finished')


if __name__ == '__main__':
    main()
