import argparse
import os.path
import subprocess
import sys
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
    parser.add_argument('--tmp-dir',
                        required=True,
                        help='where to write temporary files')
    parser.add_argument('--out-abun-tsv',
                        required=True,
                        help='where to write compared abundance')
    parser.add_argument('--out-compat-tsv',
                        required=False,
                        help='where to write compared compatible isoforms')

    return parser.parse_args()


def write_tsv_line(handle, columns):
    handle.write('{}\n'.format('\t'.join(columns)))


def read_tsv_line(line):
    return line.rstrip('\n').split('\t')


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


def compare_abundance(abun_a, abun_b, gtf_a, gtf_b, out_handle):
    with open(abun_a, 'rt') as abun_a_handle:
        with open(abun_b, 'rt') as abun_b_handle:
            with open(gtf_a, 'rt') as gtf_a_handle:
                with open(gtf_b, 'rt') as gtf_b_handle:
                    compare_abundance_with_handles(abun_a_handle,
                                                   abun_b_handle, gtf_a_handle,
                                                   gtf_b_handle, out_handle)


def output_abundance(parsed_a, parsed_b, sample_names, out_handle):
    num_samples = len(sample_names)
    id_a = None
    id_b = None
    values_a = [None] * num_samples
    values_b = [None] * num_samples
    if parsed_a:
        id_a = parsed_a['transcript_id']
        values_a = parsed_a['sample_values']
    if parsed_b:
        id_b = parsed_b['transcript_id']
        values_b = parsed_b['sample_values']

    columns = [id_a, id_b]
    for i, value_a in enumerate(values_a):
        value_b = values_b[i]
        columns.extend([value_a, value_b])

    write_tsv_line(out_handle, [str(x) for x in columns])


def read_novel_defs(gtf_handle, novel_chr):
    novel_defs = dict()
    offset = gtf_handle.tell()
    while True:
        line = gtf_handle.readline()
        if not line:
            break

        columns = read_tsv_line(line)
        chr_name = columns[0]
        if chr_name < novel_chr:
            continue
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
                   sample_names, out_handle):
    abuns_by_coords_a = make_abuns_by_coords(novel_abuns_a, novel_defs_a)
    abuns_by_coords_b = make_abuns_by_coords(novel_abuns_b, novel_defs_b)
    for coords_a, abun_a in abuns_by_coords_a.items():
        abun_b = abuns_by_coords_b.get(coords_a)
        if not abun_b:
            output_abundance(abun_a, None, sample_names, out_handle)
        else:
            del abuns_by_coords_b[coords_a]
            output_abundance(abun_a, abun_b, sample_names, out_handle)

    for abun_b in abuns_by_coords_b.values():
        output_abundance(None, abun_b, sample_names, out_handle)


def output_extra_novels(novel_abuns, is_a, sample_names, out_handle):
    for abun in novel_abuns:
        if is_a:
            output_abundance(abun, None, sample_names, out_handle)
        else:
            output_abundance(None, abun, sample_names, out_handle)


def read_samples_from_header(handle):
    line = handle.readline()
    headers = read_tsv_line(line)
    samples = headers[3:]
    return samples


def write_compat_headers(out_handle):
    columns = ['read_id', 'sample', 'compat_a', 'compat_b', 'both']
    write_tsv_line(out_handle, columns)


def output_annotated_compat_read(parsed_a, parsed_b, out_handle):
    output_novel_compat_read(parsed_a, parsed_b, None, None, out_handle)


def process_novel_compat_read_info(parsed, novel_defs):
    result = {
        'read_id': None,
        'sample': None,
        'isoforms': None,
        'coords_or_isoforms': None,
    }
    if not parsed:
        return result

    result['read_id'] = parsed['read_id']
    result['sample'] = parsed['sample_name']
    isoforms = parsed['compat_isoforms']
    if novel_defs:
        coords_or_isos = convert_novel_isoforms_to_coords(isoforms, novel_defs)
        combined_sorted = sorted(zip(isoforms, coords_or_isos),
                                 key=lambda pair: pair[1])
        result['isoforms'] = [pair[0] for pair in combined_sorted]
        result['coords_or_isoforms'] = [pair[1] for pair in combined_sorted]
    else:
        result['isoforms'] = isoforms
        result['coords_or_isoforms'] = isoforms

    return result


def get_shared_or_different_compat_isoforms(orig_a, converted_a, orig_b,
                                            converted_b, only_a, only_b, both):
    a_i = 0
    b_i = 0
    len_a = len(converted_a)
    len_b = len(converted_b)
    while a_i < len_a and b_i < len_b:
        iso_a = converted_a[a_i]
        iso_b = converted_b[b_i]
        if iso_a < iso_b:
            only_a.append([orig_a[a_i], iso_a])
            a_i += 1
        elif iso_b < iso_a:
            only_b.append([orig_b[b_i], iso_b])
            b_i += 1
        else:
            both.append(orig_a[a_i])
            a_i += 1
            b_i += 1

    while a_i < len_a:
        only_a.append([orig_a[a_i], converted_a[a_i]])
        a_i += 1

    while b_i < len_b:
        only_b.append([orig_b[b_i], converted_b[b_i]])
        b_i += 1


def output_novel_compat_read(parsed_a, parsed_b, novel_defs_a, novel_defs_b,
                             out_handle):
    details_a = process_novel_compat_read_info(parsed_a, novel_defs_a)
    details_b = process_novel_compat_read_info(parsed_b, novel_defs_b)

    # Each element of only_a (and only_b) is a list.
    # The 1st element of each inner list will be the isoform ID.
    # For novel isoforms the 2nd element will be the coordinates.
    only_a = list()
    only_b = list()
    both = list()
    orig_a = details_a['isoforms']
    converted_a = details_a['coords_or_isoforms']
    orig_b = details_b['isoforms']
    converted_b = details_b['coords_or_isoforms']
    if orig_a and orig_b:
        read_id = details_a['read_id']
        sample = details_a['sample']
        if details_a['sample'] != details_b['sample']:
            print('different samples for {}: {}, {}'.format(
                read_id, details_a['sample'], details_b['sample']),
                  file=sys.stderr)

        get_shared_or_different_compat_isoforms(orig_a, converted_a, orig_b,
                                                converted_b, only_a, only_b,
                                                both)
    elif orig_a:
        read_id = details_a['read_id']
        sample = details_a['sample']
        only_a = [[x] for x in orig_a]
    else:
        read_id = details_b['read_id']
        sample = details_b['sample']
        only_b = [[x] for x in orig_b]

    if not (only_a or only_b):
        return

    only_a_str = format_isoform_list_for_compat_output(only_a)
    only_b_str = format_isoform_list_for_compat_output(only_b)
    both_str = format_isoform_list_for_compat_output(both)
    columns = [read_id, sample, only_a_str, only_b_str, both_str]
    write_tsv_line(out_handle, columns)


def format_isoform_list_for_compat_output(isoforms):
    if not isoforms:
        return 'None'

    final_parts = list()
    for isoform_or_list in isoforms:
        if isinstance(isoform_or_list, list):
            if (((len(isoform_or_list) == 2)
                 and (isoform_or_list[0] == isoform_or_list[1]))):
                final_parts.append(isoform_or_list[0])
            else:
                final_parts.append('|'.join(isoform_or_list))
        else:
            final_parts.append(isoform_or_list)

    return ';'.join(final_parts)


def write_abundance_headers(out_handle, samples):
    columns = ['id_a', 'id_b']
    for sample in samples:
        columns.append('{}_a'.format(sample))
        columns.append('{}_b'.format(sample))

    write_tsv_line(out_handle, columns)


def compare_abundance_with_handles(abun_a_handle, abun_b_handle, gtf_a_handle,
                                   gtf_b_handle, out_handle):
    samples_a = read_samples_from_header(abun_a_handle)
    samples_b = read_samples_from_header(abun_b_handle)
    if samples_a != samples_b:
        print('sample names differ: {}, {}'.format(samples_a, samples_b),
              file=sys.stderr)

    write_abundance_headers(out_handle, samples_a)

    line_a, line_b = compare_annotated_abundance(abun_a_handle, abun_b_handle,
                                                 samples_a, out_handle)
    if not (line_a or line_b):
        # processed all lines
        return

    compare_novel_abundance(line_a, line_b, abun_a_handle, abun_b_handle,
                            gtf_a_handle, gtf_b_handle, samples_a, out_handle)


def compare_annotated_abundance(abun_a_handle, abun_b_handle, sample_names,
                                out_handle):
    line_a = abun_a_handle.readline()
    line_b = abun_b_handle.readline()
    while line_a and line_b:
        parsed_a = parse_abundance_line(line_a)
        parsed_b = parse_abundance_line(line_b)
        if parsed_a['is_novel'] or parsed_b['is_novel']:
            if parsed_a['is_novel'] and parsed_b['is_novel']:
                break

            if parsed_a['is_novel']:
                output_abundance(None, parsed_b, sample_names, out_handle)
                line_b = abun_b_handle.readline()

            if parsed_b['is_novel']:
                output_abundance(parsed_a, None, sample_names, out_handle)
                line_a = abun_a_handle.readline()
        elif parsed_a['transcript_id'] == parsed_b['transcript_id']:
            if (((parsed_a['transcript_name'] != parsed_b['transcript_name'])
                 or (parsed_a['gene_id'] != parsed_b['gene_id']))):
                print(
                    'differing transcript info for id: {} {} {} {} {}'.format(
                        parsed_a['transcript_id'], parsed_a['transcript_name'],
                        parsed_b['transcript_name'], parsed_a['gene_id'],
                        parsed_b['gene_id']),
                    file=sys.stderr)

            output_abundance(parsed_a, parsed_b, sample_names, out_handle)
            line_a = abun_a_handle.readline()
            line_b = abun_b_handle.readline()
        else:
            id_a = parsed_a['transcript_id']
            id_b = parsed_b['transcript_id']
            if id_a < id_b:
                output_abundance(parsed_a, None, sample_names, out_handle)
                line_a = abun_a_handle.readline()
            else:
                output_abundance(None, parsed_b, sample_names, out_handle)
                line_b = abun_b_handle.readline()

    both_novel = line_a and line_b
    both_finished = (not line_a) and (not line_b)
    if both_novel or both_finished:
        return line_a, line_b

    # Either a or b has no more lines.
    # Read the other file until a novel line or the end
    while line_a:
        parsed_a = parse_abundance_line(line_a)
        if parsed_a['is_novel']:
            break

        output_abundance(parsed_a, None, sample_names, out_handle)
        line_a = abun_a_handle.readline()

    while line_b:
        parsed_b = parse_abundance_line(line_b)
        if parsed_b['is_novel']:
            break

        output_abundance(None, parsed_b, sample_names, out_handle)
        line_b = abun_b_handle.readline()

    return line_a, line_b


def process_novel_abundance_line(details, abun_handle, gtf_handle):
    if details['line']:
        parsed = parse_abundance_line(details['line'])
        if not parsed['is_novel']:
            raise Exception(
                'unexpected annotated transcript when processing novel transcripts: {}'
                .format(details['line']))

        # novel_chr is set back to None after the chr has been processed
        if details['novel_chr'] is None:
            details['novel_chr'] = parsed['chr']
            details['novel_abuns'] = [parsed]
            details['line'] = abun_handle.readline()
        elif parsed['chr'] == details['novel_chr']:
            details['novel_abuns'].append(parsed)
            details['line'] = abun_handle.readline()
        elif details['novel_defs'] is None:
            # instead of moving to a new chr, load novel_defs_a to indicate
            # that novel_chr_a is ready to compare
            details['novel_defs'] = read_novel_defs(gtf_handle,
                                                    details['novel_chr'])
    elif details['novel_abuns'] and (not details['novel_defs']):
        # indicate ready to compare final novel_chr
        details['novel_defs'] = read_novel_defs(gtf_handle,
                                                details['novel_chr'])


def compare_novel_abundance(line_a, line_b, abun_a_handle, abun_b_handle,
                            gtf_a_handle, gtf_b_handle, sample_names,
                            out_handle):
    # novel_defs is set when all lines for the chr are loaded
    details_a = {
        'novel_chr': None,
        'novel_abuns': None,
        'novel_defs': None,
        'line': line_a,
    }
    details_b = {
        'novel_chr': None,
        'novel_abuns': None,
        'novel_defs': None,
        'line': line_b,
    }
    while details_a['line'] or details_b['line']:
        process_novel_abundance_line(details_a, abun_a_handle, gtf_a_handle)
        process_novel_abundance_line(details_b, abun_b_handle, gtf_b_handle)
        if details_a['novel_defs'] and details_b['novel_defs']:
            if details_a['novel_chr'] == details_b['novel_chr']:
                compare_novels(details_a['novel_defs'],
                               details_b['novel_defs'],
                               details_a['novel_abuns'],
                               details_b['novel_abuns'], sample_names,
                               out_handle)
                details_a['novel_chr'] = None
                details_a['novel_abuns'] = None
                details_a['novel_defs'] = None
                details_b['novel_chr'] = None
                details_b['novel_abuns'] = None
                details_b['novel_defs'] = None
            elif details_a['novel_chr'] < details_b['novel_chr']:
                output_extra_novels(details_a['novel_abuns'], True,
                                    sample_names, out_handle)
                details_a['novel_chr'] = None
                details_a['novel_abuns'] = None
                details_a['novel_defs'] = None
            else:
                output_extra_novels(details_b['novel_abuns'], False,
                                    sample_names, out_handle)
                details_b['novel_chr'] = None
                details_b['novel_abuns'] = None
                details_b['novel_defs'] = None
        elif details_a['novel_defs'] and (not details_b['line']):
            output_extra_novels(details_a['novel_abuns'], True, sample_names,
                                out_handle)
            details_a['novel_chr'] = None
            details_a['novel_abuns'] = None
            details_a['novel_defs'] = None
        elif details_b['novel_defs'] and (not details_a['line']):
            output_extra_novels(details_b['novel_abuns'], False, sample_names,
                                out_handle)
            details_b['novel_chr'] = None
            details_b['novel_abuns'] = None
            details_b['novel_defs'] = None

    # process each line one last time to potentially set novel_defs
    process_novel_abundance_line(details_a, abun_a_handle, gtf_a_handle)
    process_novel_abundance_line(details_b, abun_b_handle, gtf_b_handle)
    if details_a['novel_defs'] and details_b['novel_defs']:
        compare_novels(details_a['novel_defs'], details_b['novel_defs'],
                       details_a['novel_abuns'], details_b['novel_abuns'],
                       sample_names, out_handle)
    elif details_a['novel_defs']:
        output_extra_novels(details_a['novel_abuns'], True, sample_names,
                            out_handle)
    elif details_b['novel_defs']:
        output_extra_novels(details_b['novel_abuns'], False, sample_names,
                            out_handle)


def try_parse_float(string):
    try:
        return float(string)
    except ValueError:
        return None


def parse_abundance_line(line):
    columns = read_tsv_line(line)
    # The 1st column is added so that the abundance and
    # gtf files are sorted by chr.
    # For annotated transcripts chr_for_sort = None
    chr_for_sort = columns[0]
    transcript_id = columns[1]
    transcript_name = columns[2]
    gene_id = columns[3]
    sample_float_strs = columns[4:]
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

        columns = read_tsv_line(in_line)
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
    write_tsv_line(handle, columns)


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
    headers = None
    for i, line in enumerate(in_handle):
        columns = read_tsv_line(line)
        if i == 0:
            headers = columns
            expected_headers = ['transcript_ID', 'transcript_name', 'gene_ID']
            if headers[:3] != expected_headers:
                raise Exception(
                    'unexpected headers in abundance file: {}'.format(line))

            continue

        # chr_for_sort is used so that the sort order matches the sorted gtf
        transcript_id = columns[0]
        if transcript_id.startswith('ESPRESSO:'):
            id_splits = transcript_id.split(':')
            chr_for_sort = id_splits[1]
            write_tsv_line(novel_handle, [chr_for_sort] + columns)
        else:
            chr_for_sort = 'None'
            write_tsv_line(annotated_handle, [chr_for_sort] + columns)

    return headers


def combine_sorted_annotated_and_novel(annotated_handle, novel_handle,
                                       combined_handle, headers):
    write_tsv_line(combined_handle, headers)
    for annotated_line in annotated_handle:
        combined_handle.write(annotated_line)

    for novel_line in novel_handle:
        combined_handle.write(novel_line)


def sort_abundance(orig, temp_dir, sort_memory_buffer_size):
    temp_file_name_1 = get_temp_file_name(temp_dir)
    temp_file_name_2 = get_temp_file_name(temp_dir)
    with open(orig, 'rt') as in_handle:
        with open(temp_file_name_1, 'wt') as annotated_handle:
            with open(temp_file_name_2, 'wt') as novel_handle:
                headers = split_annotated_and_novel_abun(
                    in_handle, annotated_handle, novel_handle)

    sorted_temp_file_name_1 = get_temp_file_name(temp_dir)
    sorted_temp_file_name_2 = get_temp_file_name(temp_dir)
    sort_file(temp_file_name_1, sorted_temp_file_name_1,
              sort_memory_buffer_size, temp_dir)
    sort_file(temp_file_name_2, sorted_temp_file_name_2,
              sort_memory_buffer_size, temp_dir)

    final_temp_file_name = get_temp_file_name(temp_dir)
    with open(sorted_temp_file_name_1, 'rt') as annotated_handle:
        with open(sorted_temp_file_name_2, 'rt') as novel_handle:
            with open(final_temp_file_name, 'wt') as combined_handle:
                combine_sorted_annotated_and_novel(annotated_handle,
                                                   novel_handle,
                                                   combined_handle, headers)

    os.remove(temp_file_name_1)
    os.remove(temp_file_name_2)
    os.remove(sorted_temp_file_name_1)
    os.remove(sorted_temp_file_name_2)
    return final_temp_file_name


def split_compat_isoforms_with_handles(in_handle, annotated_handle,
                                       novel_handle):
    for line in in_handle:
        columns = read_tsv_line(line)
        parsed = parse_compat_line(line)
        if parsed['is_novel']:
            new_columns = [parsed['chr_name']] + columns
            write_tsv_line(novel_handle, new_columns)
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
    columns = read_tsv_line(line)
    # chr_name = columns[0]
    orig_line = '\t'.join(columns[1:])
    return parse_compat_line(orig_line)


def parse_compat_line(line):
    columns = read_tsv_line(line)
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


def compare_annotated_compat(compat_a, compat_b, out_handle):
    with open(compat_a, 'rt') as a_handle:
        with open(compat_b, 'rt') as b_handle:
            line_a = a_handle.readline()
            line_b = b_handle.readline()
            while line_a and line_b:
                parsed_a = parse_compat_line(line_a)
                parsed_b = parse_compat_line(line_b)
                if parsed_a['read_id'] < parsed_b['read_id']:
                    output_annotated_compat_read(parsed_a, None, out_handle)
                    line_a = a_handle.readline()
                elif parsed_b['read_id'] < parsed_a['read_id']:
                    output_annotated_compat_read(None, parsed_b, out_handle)
                    line_b = b_handle.readline()
                else:
                    if parsed_a != parsed_b:
                        output_annotated_compat_read(parsed_a, parsed_b,
                                                     out_handle)

                    line_a = a_handle.readline()
                    line_b = b_handle.readline()

            while line_a:
                parsed_a = parse_compat_line(line_a)
                output_annotated_compat_read(parsed_a, None, out_handle)
                line_a = a_handle.readline()

            while line_b:
                parsed_b = parse_compat_line(line_b)
                output_annotated_compat_read(None, parsed_b, out_handle)
                line_b = b_handle.readline()


def convert_novel_isoforms_to_coords(isoforms, novel_defs):
    isoform_defs = list()
    for isoform in isoforms:
        if isoform.startswith('ESPRESSO'):
            definition = novel_defs.get(isoform)
            coords = coords_from_definition(definition)
            coord_str = ','.join([str(x) for x in coords])
            isoform_defs.append(coord_str)
        else:
            isoform_defs.append(isoform)

    return isoform_defs


def compare_novel_compat_with_handles(compat_a_handle, compat_b_handle,
                                      gtf_a_handle, gtf_b_handle, out_handle):
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
            output_novel_compat_read(parsed_a, None, novel_defs_a, None,
                                     out_handle)
            line_a = compat_a_handle.readline()
        elif parsed_b['read_id'] < parsed_a['read_id']:
            output_novel_compat_read(None, parsed_b, None, novel_defs_b,
                                     out_handle)
            line_b = compat_b_handle.readline()
        else:
            output_novel_compat_read(parsed_a, parsed_b, novel_defs_a,
                                     novel_defs_b, out_handle)
            line_a = compat_a_handle.readline()
            line_b = compat_b_handle.readline()

    while line_a:
        parsed_a = parse_novel_compat_line(line_a)
        if novel_chr_a != parsed_a['chr_name']:
            novel_chr_a = parsed_a['chr_name']
            novel_defs_a = read_novel_defs(gtf_a_handle, novel_chr_a)

        output_novel_compat_read(parsed_a, None, novel_defs_a, None,
                                 out_handle)
        line_a = compat_a_handle.readline()

    while line_b:
        parsed_b = parse_novel_compat_line(line_b)
        if novel_chr_b != parsed_b['chr_name']:
            novel_chr_b = parsed_b['chr_name']
            novel_defs_b = read_novel_defs(gtf_b_handle, novel_chr_b)

        output_novel_compat_read(None, parsed_b, None, novel_defs_b,
                                 out_handle)
        line_b = compat_b_handle.readline()


def compare_novel_compat(compat_a, compat_b, normed_gtf_a, normed_gtf_b,
                         out_handle):
    with open(compat_a, 'rt') as compat_a_handle:
        with open(compat_b, 'rt') as compat_b_handle:
            with open(normed_gtf_a, 'rt') as gtf_a_handle:
                with open(normed_gtf_b, 'rt') as gtf_b_handle:
                    compare_novel_compat_with_handles(compat_a_handle,
                                                      compat_b_handle,
                                                      gtf_a_handle,
                                                      gtf_b_handle, out_handle)


def compare_compatible_isoforms(compat_a, compat_b, normed_gtf_a, normed_gtf_b,
                                temp_dir, sort_memory_buffer_size, out_handle):
    annotated_compat_a, novel_compat_a = split_and_sort_compat_isoforms(
        compat_a, temp_dir, sort_memory_buffer_size)
    annotated_compat_b, novel_compat_b = split_and_sort_compat_isoforms(
        compat_b, temp_dir, sort_memory_buffer_size)

    write_compat_headers(out_handle)
    compare_annotated_compat(annotated_compat_a, annotated_compat_b,
                             out_handle)
    compare_novel_compat(novel_compat_a, novel_compat_b, normed_gtf_a,
                         normed_gtf_b, out_handle)


def compare_espresso_output(args, temp_dir):
    sorted_abun_a = sort_abundance(args.abun_a, temp_dir,
                                   args.sort_memory_buffer_size)
    sorted_abun_b = sort_abundance(args.abun_b, temp_dir,
                                   args.sort_memory_buffer_size)

    normed_gtf_a = normalize_gtf(args.gtf_a, temp_dir,
                                 args.sort_memory_buffer_size)
    normed_gtf_b = normalize_gtf(args.gtf_b, temp_dir,
                                 args.sort_memory_buffer_size)

    with open(args.out_abun_tsv, 'wt') as abun_out_handle:
        compare_abundance(sorted_abun_a, sorted_abun_b, normed_gtf_a,
                          normed_gtf_b, abun_out_handle)

    if args.compat_a or args.compat_b:
        if not (args.compat_a and args.compat_b and args.out_compat_tsv):
            raise Exception('Either all or none of'
                            ' [--compat-a, --compat-b, --out-compat-tsv]'
                            ' is required')

        with open(args.out_compat_tsv, 'wt') as compat_out_handle:
            compare_compatible_isoforms(args.compat_a, args.compat_b,
                                        normed_gtf_a, normed_gtf_b, temp_dir,
                                        args.sort_memory_buffer_size,
                                        compat_out_handle)


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

    print('finished', file=sys.stderr)


if __name__ == '__main__':
    main()
