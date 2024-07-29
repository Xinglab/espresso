import argparse
import ctypes
import ctypes.util
import datetime
import json
import multiprocessing
import multiprocessing.sharedctypes
import os
import os.path
import queue
import shutil
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Create new .sam files using the original .sam files and'
                     ' the ESPRESSO corrections'))
    parser.add_argument('--samples-tsv',
                        required=True,
                        help='the samples.tsv used with ESPRESSO')
    parser.add_argument('--espresso-out-dir',
                        required=True,
                        help=('directory with ESPRESSO output'
                              ' (the *_read_final.txt files will be used)'))
    parser.add_argument(
        '--out-dir',
        required=True,
        help=('the path of the directory where output files (and temporary'
              ' files) will be written'))
    parser.add_argument('--fasta',
                        required=True,
                        help='the genome sequence used with ESPRESSO')
    parser.add_argument('--sort-memory-buffer-size',
                        default='2G',
                        help=('how much memory "sort" is allowed to use.'
                              ' Specified as {num_gb}G (default %(default)s)'))
    parser.add_argument('--alignment-batch-size',
                        type=int,
                        default=1000,
                        help=('how many alignments to process at a time'
                              ' (default %(default)s)'))
    parser.add_argument(
        '--libparasail-so-path',
        help=('The path to libparasail.so. If the path is not given,'
              ' the system paths and the directory of this script are'
              ' searched for the library'))
    parser.add_argument(
        '--progress-every-n',
        type=int,
        default=10_000,
        help=('print an update message after processing this number of'
              ' corrected alignments (default %(default)s)'))
    parser.add_argument('--num-threads',
                        type=int,
                        default=1,
                        help='How many threads to use (default %(default)s)')
    parser.add_argument(
        '--aligner-window-size',
        type=int,
        help=('A limit on the sequence length given to the aligner.'
              ' A smaller window can result in shorter running time,'
              ' but worse cigar strings'))

    return parser.parse_args()


def create_c_string(string):
    as_bytes = bytes(string, encoding='utf-8')
    return ctypes.create_string_buffer(as_bytes)


def print_with_timestamp(message, flush=True):
    now = datetime.datetime.now()
    # [Fri Jan 02 15:40:50 2060]
    now_string = now.strftime('[%a %b %d %H:%M:%S %Y]')

    print('{} {}'.format(now_string, message), flush=flush)


def parse_na_int_str(string):
    if string == 'NA':
        return None

    return int(string)


def parse_na_float_str(string):
    if string == 'NA':
        return None

    return float(string)


def parse_strand(string):
    if string == '0':
        return '+'
    if string == '1':
        return '-'

    # could be 'unknown'
    return string


def parse_sj_id(id_str):
    if id_str == 'NA':
        return {'chr': None, 'start': None, 'end': None, 'strand': None}

    parts = id_str.split(':')
    sj_chr = parts[0]
    sj_start = int(parts[1])
    sj_end = int(parts[2])
    sj_strand = parse_strand(parts[3])

    return {
        'chr': sj_chr,
        'start': sj_start,
        'end': sj_end,
        'strand': sj_strand
    }


def parse_read_final_read_id(in_path, handle, next_line):
    read_details = dict()
    if next_line:
        details_added = process_read_final_line(read_details, next_line,
                                                in_path)
        if not details_added:
            raise Exception('Unable to use line from {}: {}'.format(
                in_path, next_line))

    for line in handle:
        details_added = process_read_final_line(read_details, line, in_path)
        if not details_added:
            return {'read_details': read_details, 'next_line': line}

    return {'read_details': read_details, 'next_line': None}


# process_read_final_line returns True if it added to read_details.
# A False return indicates that the line is for a new read_id.
def process_read_final_line(read_details, line, in_path):
    columns = line.rstrip('\n').split('\t')
    line_read_id = columns[0]
    old_read_id = read_details.get('read_id')
    if old_read_id and line_read_id != old_read_id:
        return False

    read_details['read_id'] = line_read_id
    feature = columns[1]
    if feature == 'group_ID':
        read_details['group'] = int(columns[2])
        # ${n}_$read_ID_update
        read_details['group_updated'] = columns[3]
    elif feature == 'strand_isoform':
        read_details['strand_isoform'] = parse_strand(columns[2])
    elif feature == 'strand_read':
        read_details['strand_read'] = parse_strand(columns[2])
    elif feature == 'chr':
        read_details['chr'] = columns[2]
    elif feature == 'mapq':
        read_details['mapq'] = int(columns[2])
    elif feature == 'flag':
        read_details['flag'] = int(columns[2])
    elif feature == 'mapped_length_read':
        read_details['mapped_length'] = int(columns[2])
    elif feature == 'SJ':
        process_read_final_sj_feature(read_details, columns)
    elif feature in ['start', 'end']:
        endpoint_details = dict()
        read_details[feature] = endpoint_details
        endpoint_details['pos'] = int(columns[3]) - 1  # convert to 0-based
        endpoint_details['read_length_before'] = int(columns[4])
        endpoint_details['read_length_after'] = int(columns[5])

        endpoint_details['corrected_status'] = columns[6]
        # convert to 0-based
        endpoint_details['corrected_pos'] = int(columns[7]) - 1
        endpoint_details['corrected_read_length_before'] = parse_na_int_str(
            columns[8])
        endpoint_details['corrected_read_length_after'] = parse_na_int_str(
            columns[9])
    else:
        raise Exception('unrecognized feature {} in {}: {}'.format(
            feature, in_path, line))

    return True


def process_read_final_sj_feature(read_details, columns):
    sj_i = int(columns[2])
    sj_details = dict()
    if sj_i == 0:
        read_details['added_first'] = sj_details
    elif sj_i == -1:
        read_details['added_last'] = sj_details
    else:
        sj_i -= 1  # get 0-based index
        sjs = read_details.get('SJs')
        if not sjs:
            sjs = list()
            read_details['SJs'] = sjs

        while len(sjs) <= sj_i:
            sjs.append(None)

        sjs[sj_i] = sj_details

    parsed_sj_id = parse_sj_id(columns[3])
    sj_details['chr'] = parsed_sj_id['chr']
    # 0-based 1st position of intron
    sj_details['start'] = parsed_sj_id['start']
    # 1-based last position of intron
    sj_details['end'] = parsed_sj_id['end']
    sj_details['strand'] = parsed_sj_id['strand']
    sj_details['read_length_before_sj'] = parse_na_int_str(columns[4])
    sj_details['read_length_after_sj'] = parse_na_int_str(columns[5])

    sj_details['corrected_status'] = columns[6]
    parsed_corrected_sj_id = parse_sj_id(columns[7])
    sj_details['corrected_chr'] = parsed_corrected_sj_id['chr']
    sj_details['corrected_start'] = parsed_corrected_sj_id['start']
    sj_details['corrected_end'] = parsed_corrected_sj_id['end']
    sj_details['corrected_strand'] = parsed_corrected_sj_id['strand']
    sj_details['corrected_read_length_before_sj'] = parse_na_int_str(
        columns[8])
    sj_details['corrected_read_length_after_sj'] = parse_na_int_str(columns[9])
    sj_details['corrected_score'] = parse_na_float_str(columns[10])

    sj_details['is_annotated'] = columns[11] == 'yes'
    sj_details['corrected_is_annotated'] = columns[12] == 'yes'


def get_contig_order_from_fasta(fasta_path):
    print_with_timestamp('Reading {}'.format(fasta_path))
    contig_order = list()
    with open(fasta_path, 'rt') as handle:
        for line in handle:
            line = line.rstrip('\n')
            if (not contig_order) or line.startswith('>'):
                found_contig = parse_fasta_contig_line(line)
                contig_order.append(found_contig)

    return contig_order


def parse_samples_tsv(samples_tsv_path):
    orig_sam_or_bam_paths = list()
    with open(samples_tsv_path, 'rt') as handle:
        for line in handle:
            columns = line.rstrip('\n').split('\t')
            sam_or_bam_path = columns[0]
            abs_path = os.path.abspath(sam_or_bam_path)
            orig_sam_or_bam_paths.append(abs_path)

    return orig_sam_or_bam_paths


def index_sam_path(sam_path, index_path):
    index = dict()
    with open(sam_path, 'rt') as in_handle:
        offset = in_handle.tell()
        # Using a while loop and readline() to allow using tell()
        while True:
            line = in_handle.readline()
            if line == '':
                break

            previous_offset = offset
            offset = in_handle.tell()
            if line.startswith('@'):
                continue

            columns = line.rstrip('\n').split('\t')
            contig = columns[2]
            if contig in index:
                continue

            index[contig] = previous_offset

    with open(index_path, 'wt') as out_handle:
        index_json = json.dumps(index)
        out_handle.write('{}\n'.format(index_json))


def convert_bam_to_sam(bam_path, sam_path):
    command = ['samtools', 'view', '-h', '-o', sam_path, bam_path]
    subprocess.run(command, check=True)


# ESPRESSO may use a supplementary alignment and the full sequence for
# that read is only found in the primary alignment.  Ensure all
# entries for a read-id have the full sequence and quality.  That way
# the correct read can be found by coordinate and it will have the
# full sequence and quality.
def fix_sam_seq_and_qual(in_path, out_path):
    temp_path = '{}.tmp'.format(out_path)
    sort_sam_by_read_name(in_path, out_path)
    fix_seq_and_qual_for_name_sorted_sam(out_path, temp_path)
    sort_sam(temp_path, out_path)
    os.remove(temp_path)


def fix_seq_and_qual_for_name_sorted_sam(in_path, out_path):
    read_id = None
    rows_for_read = list()
    with open(in_path, 'rt') as in_handle:
        with open(out_path, 'wt') as out_handle:
            for line in in_handle:
                if line.startswith('@'):
                    out_handle.write(line)
                    continue

                columns = line.rstrip('\n').split('\t')
                line_read_id = columns[0]
                if (read_id is None) or (line_read_id != read_id):
                    if read_id is not None:
                        fix_seq_and_qual_for_read_rows(rows_for_read,
                                                       out_handle)

                    read_id = line_read_id
                    rows_for_read = list()

                rows_for_read.append(columns)

            if read_id is not None:
                fix_seq_and_qual_for_read_rows(rows_for_read, out_handle)


def fix_seq_and_qual_for_read_rows(rows, out_handle):
    seq = ''
    qual = ''
    for columns in rows:
        row_seq = columns[9]
        row_qual = columns[10]
        if len(row_seq) > len(seq):
            seq = row_seq
            qual = row_qual

    # The cigar must be modified to match the sequence length.
    # Otherwise samtools gives an error.
    fixed_cigar = '1M{}I'.format(len(seq) - 1)
    for columns in rows:
        columns[5] = fixed_cigar
        columns[9] = seq
        columns[10] = qual
        out_handle.write('{}\n'.format('\t'.join(columns)))


def index_sam_or_bam_path(new_basename, extension, orig_path, fixed_sam_path,
                          out_dir):
    print_with_timestamp('Indexing {}'.format(orig_path))
    if extension == '.sam':
        sam_path = orig_path
        fix_sam_seq_and_qual(sam_path, fixed_sam_path)
    elif extension == '.bam':
        sam_path = os.path.join(out_dir, 'orig_sams', new_basename)
        convert_bam_to_sam(orig_path, sam_path)
        fix_sam_seq_and_qual(sam_path, fixed_sam_path)
        os.remove(sam_path)
    else:
        raise Exception('Unrecognized file extension: {}'.format(orig_path))

    index_path = os.path.join(out_dir, 'indices',
                              '{}.index'.format(new_basename))
    index_sam_path(fixed_sam_path, index_path)
    index_details = {
        'orig_path': orig_path,
        'sam_path': fixed_sam_path,
        'index_path': index_path
    }
    index_details['current_contig'] = None
    index_details['current_offset'] = None
    index_details['loaded_alignments'] = None
    index_details['next_start_pos_i'] = None

    return {'new_basename': new_basename, 'index_details': index_details}


def try_get_from_queue_with_short_wait(in_queue):
    try:
        # Wait at most 1 second for value
        value = in_queue.get(True, 1)
    except queue.Empty:
        return None

    return value


def try_get_from_queue_without_wait(in_queue):
    try:
        value = in_queue.get(False)
    except queue.Empty:
        return None

    return value


def try_put_to_queue_with_short_wait(out_queue, value):
    try:
        # Wait at most 1 second
        out_queue.put(value, True, 1)
    except queue.Full:
        return False

    return True


def try_put_to_queue_without_wait(out_queue, value):
    try:
        out_queue.put(value, False)
    except queue.Full:
        return False

    return True


def index_sam_or_bam_path_thread(in_queue, out_queue):
    while True:
        arguments = try_get_from_queue_with_short_wait(in_queue)
        if arguments is None:
            return

        result = index_sam_or_bam_path(arguments['new_basename'],
                                       arguments['extension'],
                                       arguments['orig_path'],
                                       arguments['fixed_sam_path'],
                                       arguments['out_dir'])
        out_queue.put(result)


def run_index_sam_or_bam_paths_threads(thread_arguments, num_threads,
                                       indexed_paths):
    num_jobs = len(thread_arguments)
    threads = list()
    thread_outputs = multiprocessing.Queue(num_jobs)
    thread_inputs = multiprocessing.Queue(num_jobs)
    for arguments in thread_arguments:
        thread_inputs.put(arguments)

    for _ in range(num_threads):
        thread = multiprocessing.Process(target=index_sam_or_bam_path_thread,
                                         args=(thread_inputs, thread_outputs))
        threads.append(thread)
        thread.start()

    for _ in range(num_jobs):
        thread_result = thread_outputs.get()
        new_basename = thread_result['new_basename']
        index_details = thread_result['index_details']
        indexed_paths[new_basename] = index_details

    for thread in threads:
        thread.join()

    thread_inputs.close()
    thread_outputs.close()


def index_sam_or_bam_paths(orig_paths, out_dir, num_threads):
    indexed_paths = dict()
    out_sam_rename_path = os.path.join(out_dir, 'sam_file_name_mapping.txt')
    out_sam_renaming = dict()
    used_basenames = set()
    thread_arguments = list()
    for orig_path in orig_paths:
        basename = os.path.basename(orig_path)
        without_extension, extension = os.path.splitext(basename)
        extra_i = 0
        new_basename = '{}.sam'.format(without_extension)
        while new_basename in used_basenames:
            new_basename = '{}_{}.sam'.format(without_extension, extra_i)
            extra_i += 1

        used_basenames.add(new_basename)
        out_sam_renaming[orig_path] = new_basename
        fixed_sam_path = os.path.join(out_dir, 'orig_sams',
                                      'fixed_{}'.format(new_basename))
        thread_arguments.append({
            'new_basename': new_basename,
            'extension': extension,
            'orig_path': orig_path,
            'fixed_sam_path': fixed_sam_path,
            'out_dir': out_dir
        })

    run_index_sam_or_bam_paths_threads(thread_arguments, num_threads,
                                       indexed_paths)
    with open(out_sam_rename_path, 'wt') as rename_handle:
        for orig_path in sorted(out_sam_renaming):
            new_basename = out_sam_renaming[orig_path]
            rename_handle.write('{}\t{}\n'.format(orig_path, new_basename))

    return indexed_paths


def find_read_final_paths(espresso_dir):
    read_final_paths = list()
    path_suffix = '_read_final.txt'
    for dirpath, dirnames, filenames in os.walk(espresso_dir):
        for filename in filenames:
            if not filename.endswith(path_suffix):
                continue

            path = os.path.join(dirpath, filename)
            abs_path = os.path.abspath(path)
            read_final_paths.append(abs_path)

    return read_final_paths


def write_read_final_entries_for_sorting(in_path, in_handle, out_handle):
    next_line = None
    parsed_details = parse_read_final_read_id(in_path, in_handle, next_line)
    details = parsed_details['read_details']
    next_line = parsed_details['next_line']
    write_read_final_entry_for_sorting(details, out_handle)
    while next_line:
        parsed_details = parse_read_final_read_id(in_path, in_handle,
                                                  next_line)
        details = parsed_details['read_details']
        next_line = parsed_details['next_line']
        write_read_final_entry_for_sorting(details, out_handle)


def write_read_final_entry_for_sorting(read_details, out_handle):
    if not read_details:
        return

    json_string = json.dumps(read_details)
    start_details = read_details['start']
    start_pos = start_details['pos']
    out_handle.write('{}\t{}\n'.format(start_pos, json_string))


def parse_and_sort_read_final_paths(read_final_paths, out_dir,
                                    sort_memory_buffer_size):
    print_with_timestamp('Parsing read final files')
    read_finals_out_dir = os.path.join(out_dir, 'read_finals')
    unsorted_paths = dict()
    path_suffix = '_read_final.txt'
    for orig_path in read_final_paths:
        orig_basename = os.path.basename(orig_path)
        if not orig_basename.endswith(path_suffix):
            raise Exception('Expected {} to end with {}'.format(
                orig_path, path_suffix))

        contig = orig_basename[:-len(path_suffix)]
        unsorted_path = unsorted_paths.get(contig)
        if not unsorted_path:
            unsorted_path = os.path.join(
                read_finals_out_dir,
                '{}_unsorted{}'.format(contig, path_suffix))
            unsorted_paths[contig] = unsorted_path
            with open(unsorted_path, 'wt'):
                pass  # truncate file

        with open(orig_path, 'rt') as in_handle:
            with open(unsorted_path, 'at') as out_handle:
                write_read_final_entries_for_sorting(orig_path, in_handle,
                                                     out_handle)

    print_with_timestamp('Sorting read final files')
    sorted_paths = dict()
    sort_temp_dir = read_finals_out_dir
    for contig, unsorted_path in unsorted_paths.items():
        sorted_path = os.path.join(read_finals_out_dir,
                                   '{}_sorted{}'.format(contig, path_suffix))
        sorted_paths[contig] = sorted_path
        # sort numerically by the first column (start coordinate)
        sort_key = '1'
        sort_command = [
            'sort', '--numeric', '--key', sort_key, '--output', sorted_path,
            '--buffer-size', sort_memory_buffer_size, '--temporary-directory',
            sort_temp_dir, unsorted_path
        ]
        os.environ['LC_ALL'] = 'C'  # ensure standard sorting behavior
        subprocess.run(sort_command, check=True)
        os.remove(unsorted_path)

    return sorted_paths


def parse_contig_from_sq_header_line(line, orig_path):
    sq_parts = line.rstrip('\n').split('\t')
    if len(sq_parts) < 2:
        raise Exception('Unexpected @SQ format in {}: {}'.format(
            orig_path, line))

    contig_part = None
    sn_prefix = 'SN:'
    for sq_part in sq_parts:
        if sq_part.startswith(sn_prefix):
            contig_part = sq_part

    if not contig_part:
        raise Exception('Unexpected @SQ format in {}: {}'.format(
            orig_path, line))

    contig = contig_part[len(sn_prefix):]
    return contig


def sort_sq_headers(contig_order, sq_headers):
    sorted_sq_headers = list()
    for contig in contig_order:
        sq_header = sq_headers.get(contig)
        if not sq_header:
            continue

        sorted_sq_headers.append(sq_header)
        del sq_headers[contig]

    return sorted_sq_headers


def reorder_headers(orig_path, new_path, contig_order):
    non_sq_headers = list()
    sq_headers = dict()
    first_sq_index = None
    with open(orig_path, 'rt') as in_handle:
        for line in in_handle:
            if not line.startswith('@'):
                break

            if not line.startswith('@SQ'):
                non_sq_headers.append(line)
                continue

            if first_sq_index is None:
                first_sq_index = len(non_sq_headers)

            contig = parse_contig_from_sq_header_line(line, orig_path)
            if contig in sq_headers:
                raise Exception('duplicate @SQ contig in {}: {}'.format(
                    orig_path, line))

            sq_headers[contig] = line

    if first_sq_index is None:
        raise Exception('No @SQ headers found in {}'.format(orig_path))

    sorted_sq_headers = sort_sq_headers(contig_order, sq_headers)
    if sq_headers:
        raise Exception(
            'Extra @SQ contig(s) found in {} but not in fasta: {}'.format(
                orig_path, sq_headers))

    reordered_headers = non_sq_headers[:first_sq_index]
    reordered_headers.extend(sorted_sq_headers)
    reordered_headers.extend(non_sq_headers[first_sq_index:])
    with open(new_path, 'wt') as out_handle:
        for line in reordered_headers:
            out_handle.write(line)


def create_output_sams(indexed_orig_sam_paths, contig_order, out_dir):
    print_with_timestamp('Writing headers for output sams')
    output_sam_paths = dict()
    for new_basename, details in indexed_orig_sam_paths.items():
        sam_path = details['sam_path']
        out_path = os.path.join(out_dir, new_basename)
        reorder_headers(sam_path, out_path, contig_order)
        output_sam_paths[new_basename] = out_path

    return output_sam_paths


# read_contig reads the next contig sequence from the open
# fasta_handle.  The contig name from the handle is expected to
# match the contig argument.
def read_contig(fasta_details, contig):
    fasta_handle = fasta_details['handle']
    next_line = fasta_details['next_line']
    found_contig = None
    sequence = list()
    if next_line:
        found_contig = parse_fasta_contig_line(next_line)
        if found_contig != contig:
            raise Exception('Expected to find {} but found {}'.format(
                contig, found_contig))

    for line in fasta_handle:
        line = line.rstrip('\n')
        if found_contig is None:
            found_contig = parse_fasta_contig_line(line)
            if found_contig != contig:
                raise Exception('Expected to find {} but found {}'.format(
                    contig, found_contig))

            continue

        if line.startswith('>'):
            fasta_details['next_line'] = line
            break

        sequence.append(line)

    return ''.join(sequence)


def parse_fasta_contig_line(line):
    if not line.startswith('>'):
        raise Exception('Expected line to start with ">": {}'.format(line))

    line_parts = line.split()
    found_contig = line_parts[0][1:]
    return found_contig


def parse_sorted_read_final_line(line):
    first_tab_i = line.find('\t')
    if first_tab_i == -1:
        raise Exception('Expected to find a tab in {}'.format(line))

    read_final_json = line[(first_tab_i + 1):]
    read_final = json.loads(read_final_json)
    return read_final


def find_alignment_in_sam(sam_details, read_id, orig_start_pos, contig,
                          alignment_batch_size):
    if sam_details['current_contig'] != contig:
        with open(sam_details['index_path'], 'rt') as index_handle:
            index = json.load(index_handle)

        contig_offset = index.get(contig)
        if contig_offset is None:
            return None

        sam_details['current_contig'] = contig
        sam_details['current_offset'] = contig_offset
        load_alignments(sam_details, alignment_batch_size)

    while True:
        loaded_alignments = sam_details['loaded_alignments']
        if loaded_alignments:
            next_start_pos_i = sam_details['next_start_pos_i']
            for start_pos_i in range(next_start_pos_i, len(loaded_alignments)):
                start_pos_alignments = loaded_alignments[start_pos_i]
                found_pos = start_pos_alignments[0]['pos']
                if found_pos > orig_start_pos:
                    return None

                sam_details['next_start_pos_i'] = start_pos_i
                if found_pos < orig_start_pos:
                    continue

                for alignment in start_pos_alignments:
                    if alignment['read_id'] == read_id:
                        return alignment

            if found_pos == orig_start_pos:
                return None

        any_loaded = load_alignments(sam_details, alignment_batch_size)
        if not any_loaded:
            return None


def parse_alignment_line(line):
    alignment = dict()
    columns = line.rstrip('\n').split('\t')
    alignment['read_id'] = columns[0]
    alignment['flag'] = columns[1]
    alignment['contig'] = columns[2]
    alignment['pos'] = int(columns[3])
    alignment['mapq'] = columns[4]
    alignment['cigar'] = columns[5]
    alignment['contig_next'] = columns[6]
    alignment['pos_next'] = columns[7]
    alignment['template_len'] = columns[8]
    alignment['sequence'] = columns[9]
    alignment['quality'] = columns[10]
    alignment['tags'] = columns[11:]
    return alignment


def load_alignments(sam_details, alignment_batch_size):
    offset = sam_details['current_offset']
    sam_path = sam_details['sam_path']
    current_contig = sam_details['current_contig']
    alignments = list()
    current_start_pos = None
    total_loaded = 0
    previous_offset = offset
    with open(sam_path, 'rt') as handle:
        handle.seek(offset)
        # Using a while loop and readline() to allow using tell()
        while True:
            previous_offset = handle.tell()
            line = handle.readline()
            if line == '':
                break

            alignment = parse_alignment_line(line)
            if alignment['contig'] != current_contig:
                break

            if (((current_start_pos is None)
                 or (alignment['pos'] != current_start_pos))):
                if total_loaded >= alignment_batch_size:
                    break

                current_start_pos = alignment['pos']
                alignments.append(list())

            alignments[-1].append(alignment)
            total_loaded += 1

        sam_details['current_offset'] = previous_offset
        sam_details['loaded_alignments'] = alignments
        sam_details['next_start_pos_i'] = 0

    return len(alignments) > 0


def find_orig_alignment(read_final, indexed_orig_sam_paths,
                        alignment_batch_size):
    read_id = read_final['read_id']
    # convert to 1-based
    orig_start_pos = read_final['start']['pos'] + 1
    contig = read_final['chr']
    for new_basename, sam_details in indexed_orig_sam_paths.items():
        alignment = find_alignment_in_sam(sam_details, read_id, orig_start_pos,
                                          contig, alignment_batch_size)
        if alignment:
            alignment['new_basename'] = new_basename
            return alignment

    raise Exception(
        'Could not find original alignment for: {}'.format(read_final))


def update_intervals_for_junction(junction_details, ref_intervals,
                                  read_intervals):
    start = junction_details['corrected_start']
    if start is None:
        start = junction_details['start']
        end = junction_details['end']
    else:
        # corrected_start and corrected_end should both be non-None
        end = junction_details['corrected_end']

    # start is 0-based 1st position of intron
    ref_intervals[-1].append(start - 1)
    # end is 1-based last position of intron
    ref_intervals.append([end])

    read_before = junction_details['corrected_read_length_before_sj']
    if read_before is None:
        read_before = junction_details['read_length_before_sj']

    # The 0-based end position is the number of bases already used -1
    read_intervals[-1].append(read_before - 1)
    read_intervals.append([read_before])


def get_cigar_string_from_parts(parts):
    strings = list()
    for part in parts:
        operation, num = part
        strings.append('{}{}'.format(num, operation))

    return ''.join(strings)


def run_aligner(sequence, reference, parasail_api, aligner_window_size):
    # The alignment running time is roughly proportional to
    # len(sequence) * len(reference).
    # Also the alignment could have an error for long sequences due to
    # the score overflowing a 16 bit integer.
    # Limit the max sequence length used in a call to the aligner to:
    #  * safe_sequence_len (to avoid an error)
    #  * alignment_window_size (to speed things up)
    max_window_size = parasail_api.safe_sequence_len
    if (((aligner_window_size is not None)
         and (aligner_window_size < max_window_size))):
        max_window_size = aligner_window_size

    seq_len = len(sequence)
    ref_len = len(reference)
    longer_len = max([seq_len, ref_len])
    number_of_windows, extra = divmod(longer_len, max_window_size)
    if extra > 0:
        number_of_windows += 1

    if number_of_windows == 1:
        return run_aligner_on_window(sequence, reference, parasail_api)

    seq_window_i = 0
    ref_window_i = 0
    seq_increment, seq_extra = divmod(seq_len, number_of_windows)
    ref_increment, ref_extra = divmod(ref_len, number_of_windows)
    total_cigar_parts = list()
    for _ in range(number_of_windows):
        next_seq_window_i = seq_window_i + seq_increment
        if seq_extra > 0:
            seq_extra -= 1
            next_seq_window_i += 1

        next_ref_window_i = ref_window_i + ref_increment
        if ref_extra > 0:
            ref_extra -= 1
            next_ref_window_i += 1

        window_seq = sequence[seq_window_i:next_seq_window_i]
        window_ref = reference[ref_window_i:next_ref_window_i]
        seq_window_i = next_seq_window_i
        ref_window_i = next_ref_window_i
        window_cigar_parts = run_aligner_on_window(window_seq, window_ref,
                                                   parasail_api)
        total_cigar_parts.extend(window_cigar_parts)

    return total_cigar_parts


def run_aligner_on_window(sequence, reference, parasail_api):
    if len(sequence) == 0:
        if len(reference) == 0:
            raise Exception('Cannot align empty sequence and empty reference')

        deletion_len = len(reference)
        return [['D', deletion_len]]

    if len(reference) == 0:
        insertion_len = len(sequence)
        return [['I', insertion_len]]

    profile_p = parasail_api.profile_create_16(create_c_string(sequence),
                                               len(sequence),
                                               parasail_api.score_matrix_p)
    result_p = parasail_api.nw_trace_striped_profile_16(
        profile_p, create_c_string(reference), len(reference),
        parasail_api.gap_open, parasail_api.gap_extension)

    cigar_p = parasail_api.result_get_cigar(result_p,
                                            create_c_string(sequence),
                                            len(sequence),
                                            create_c_string(reference),
                                            len(reference),
                                            parasail_api.score_matrix_p)
    cigar_parts = list()
    for cigar_i in range(cigar_p[0].len):
        cigar_int = cigar_p[0].seq[cigar_i]
        cigar_op = parasail_api.cigar_decode_op(cigar_int).decode()
        cigar_num = parasail_api.cigar_decode_len(cigar_int)
        cigar_parts.append([cigar_op, cigar_num])

    parasail_api.cigar_free(cigar_p)

    query_end = parasail_api.result_get_end_query(result_p)
    ref_end = parasail_api.result_get_end_ref(result_p)
    parasail_api.result_free(result_p)
    parasail_api.profile_free(profile_p)
    if (query_end != (len(sequence) - 1)) or (ref_end != (len(reference) - 1)):
        raise Exception(
            'Alignment did not cover full sequence: {}/{}, {}/{}'.format(
                query_end + 1, len(sequence), ref_end + 1, len(reference)))

    return cigar_parts


def get_corrected_intervals(read_final):
    corrected_ref_intervals = list()
    corrected_read_intervals = list()
    corrected_start_pos = read_final['start']['corrected_pos']
    corrected_start_clipping = (
        read_final['start']['corrected_read_length_before'])
    if corrected_start_pos is None:
        corrected_start_pos = read_final['start']['pos']
        corrected_start_clipping = read_final['start']['read_length_before']

    corrected_ref_intervals.append([corrected_start_pos])
    corrected_read_intervals.append([corrected_start_clipping])

    added_first = read_final.get('added_first')
    if added_first:
        update_intervals_for_junction(added_first, corrected_ref_intervals,
                                      corrected_read_intervals)

    for junction in read_final.get('SJs', list()):
        update_intervals_for_junction(junction, corrected_ref_intervals,
                                      corrected_read_intervals)

    added_last = read_final.get('added_last')
    if added_last:
        update_intervals_for_junction(added_last, corrected_ref_intervals,
                                      corrected_read_intervals)

    corrected_end_pos = read_final['end']['corrected_pos']
    corrected_read_end = read_final['end']['corrected_read_length_before']
    if corrected_end_pos is None:
        corrected_end_pos = read_final['end']['pos']
        corrected_read_end = read_final['end']['read_length_before']

    corrected_ref_intervals[-1].append(corrected_end_pos)
    # The 0-based end position is the number of bases already used -1
    corrected_read_intervals[-1].append(corrected_read_end - 1)

    # TODO:
    # ESPRESSO can correct the 1st or last junction which could
    # result in parts of the read that were previously clipped
    # being used as the junction coordinate.
    # ESPRESSO adjusts the start coordinate and clipping for added
    # 1st and last SJs but not corrected 1st and last SJs.
    # For now, check for these cases and adjust the coordinates.
    previously_clipped_start_len = 0
    previously_clipped_end_len = 0
    if not added_first:
        first_start, first_end = corrected_read_intervals[0]
        if first_start > first_end:
            previously_clipped_start_len = first_start - first_end
            corrected_read_intervals[0][0] = first_end
            corrected_ref_intervals[0][0] = corrected_ref_intervals[0][1]

    if not added_last:
        last_start, last_end = corrected_read_intervals[-1]
        if last_start > last_end:
            previously_clipped_end_len = last_start - last_end
            corrected_read_intervals[-1][1] = last_start
            corrected_ref_intervals[-1][1] = corrected_ref_intervals[-1][0]

    return {
        'ref_intervals': corrected_ref_intervals,
        'read_intervals': corrected_read_intervals,
        'prev_clipped_start_len': previously_clipped_start_len,
        'prev_clipped_end_len': previously_clipped_end_len
    }


def collapse_cigar_parts_from_intervals(intervals):
    parts = list()
    for interval in intervals:
        for part in interval:
            add_cigar_part_and_maybe_collapse(parts, part)

    return parts


def add_cigar_part_and_maybe_collapse(parts, part, merge_equal_x=True):
    new_op, new_num = part
    if merge_equal_x and new_op in '=X':
        new_op = 'M'

    if len(parts) == 0:
        parts.append([new_op, new_num])
        return

    last_part_op, last_part_num = parts[-1]
    if last_part_op == new_op:
        combined_num = last_part_num + new_num
        parts[-1][1] = combined_num
        return

    parts.append([new_op, new_num])


def get_contig_seq_from_cache(contig_name, contig_cache):
    cached_contig_name = contig_cache['name']
    cached_contig_seq = contig_cache['sequence']
    if (((cached_contig_name is not None)
         and (cached_contig_name == contig_name))):
        return cached_contig_seq

    cache_queue = contig_cache['queue']
    while True:
        next_contig_info = try_get_from_queue_with_short_wait(cache_queue)
        if next_contig_info is None:
            raise Exception(
                'Could not find {} in contig cache'.format(contig_name))

        next_contig_name = next_contig_info['name']
        next_contig_seq = next_contig_info['sequence']
        if next_contig_name == contig_name:
            contig_cache['name'] = next_contig_name
            contig_cache['sequence'] = next_contig_seq
            return next_contig_seq


def check_if_interval_coords_increase(intervals):
    prev_coord = None
    for interval in intervals:
        for coord in interval:
            if (prev_coord is not None) and (coord < prev_coord):
                return False

            prev_coord = coord

    return True


def get_interval_cigars(corrected_ref_intervals, corrected_read_intervals,
                        contig_sequence, plus_strand_read_seq, read_final,
                        prev_clipped_start_len, prev_clipped_end_len,
                        parasail_api, aligner_window_size):
    # TODO:
    # ESPRESSO could correct an internal junction and end up
    # with the corrected junction coordinate being before/after
    # a previous/next SJ coordinate. This might only happen when
    # the previous/next SJ fails (is not passed or corrected).
    ref_coords_increase = check_if_interval_coords_increase(
        corrected_ref_intervals)
    read_coords_increase = check_if_interval_coords_increase(
        corrected_read_intervals)
    if not (ref_coords_increase and read_coords_increase):
        return None

    previous_ref_end = None
    junction_lengths = list()
    total_aligned_read_len = 0
    read_and_ref_sequences = list()
    for interval_i, ref_interval in enumerate(corrected_ref_intervals):
        ref_start, ref_end = ref_interval
        read_interval = corrected_read_intervals[interval_i]
        read_start, read_end = read_interval
        ref_sequence = contig_sequence[ref_start:(ref_end + 1)]
        read_sequence = plus_strand_read_seq[read_start:(read_end + 1)]
        total_aligned_read_len += len(read_sequence)
        read_and_ref_sequences.append([read_sequence, ref_sequence])
        if previous_ref_end:
            junction_lengths.append((ref_start - previous_ref_end) - 1)

        previous_ref_end = ref_end

    adjusted_mapped_length = read_final['mapped_length']
    adjusted_mapped_length += prev_clipped_start_len
    adjusted_mapped_length += prev_clipped_end_len
    if adjusted_mapped_length != total_aligned_read_len:
        raise Exception(
            'Length of sequence is inconsistent: {}, {}, {}'.format(
                adjusted_mapped_length, total_aligned_read_len,
                read_final['read_id']))

    interval_cigars = list()
    start_clipping = corrected_read_intervals[0][0]
    if start_clipping:
        interval_cigars.append([['S', start_clipping]])

    for read_and_ref_seq_i, read_and_ref_seq in enumerate(
            read_and_ref_sequences):
        read_seq, ref_seq = read_and_ref_seq
        if read_and_ref_seq_i > 0:
            junction_length = junction_lengths[read_and_ref_seq_i - 1]
            interval_cigars.append([['N', junction_length]])

        # TODO This should not happen? Raise an exception instead?
        if (len(read_seq) == 0) and (len(ref_seq) == 0):
            return None

        interval_cigar = run_aligner(read_seq, ref_seq, parasail_api,
                                     aligner_window_size)
        interval_cigars.append(interval_cigar)

    end_clipping = ((len(plus_strand_read_seq) - 1) -
                    corrected_read_intervals[-1][-1])
    if end_clipping:
        interval_cigars.append([['S', end_clipping]])

    return interval_cigars


def correct_alignment(read_final, orig_alignment, contig_name, contig_cache,
                      parasail_api, aligner_window_size):
    contig_sequence = get_contig_seq_from_cache(contig_name, contig_cache)

    alignment_is_minus_strand = read_final['strand_isoform'] == '-'
    # The sequence from the SAM file is always for the forward strand
    # since the read was mapped.
    plus_strand_read_seq = orig_alignment['sequence']

    corrected_interval_details = get_corrected_intervals(read_final)
    corrected_ref_intervals = corrected_interval_details['ref_intervals']
    corrected_read_intervals = corrected_interval_details['read_intervals']
    prev_clipped_start_len = (
        corrected_interval_details['prev_clipped_start_len'])
    prev_clipped_end_len = corrected_interval_details['prev_clipped_end_len']
    if len(corrected_ref_intervals) != len(corrected_read_intervals):
        raise Exception('Error finding exon intervals for read: {}'.format(
            read_final['read_id']))

    interval_cigars = get_interval_cigars(
        corrected_ref_intervals, corrected_read_intervals, contig_sequence,
        plus_strand_read_seq, read_final, prev_clipped_start_len,
        prev_clipped_end_len, parasail_api, aligner_window_size)
    if interval_cigars is None:
        final_cigar = '*'
    else:
        final_cigar_parts = collapse_cigar_parts_from_intervals(
            interval_cigars)
        final_cigar = get_cigar_string_from_parts(final_cigar_parts)

    corrected_start_pos = corrected_ref_intervals[0][0]
    corrected = dict()
    corrected['read_id'] = read_final['read_id']
    corrected['flag'] = 16 if alignment_is_minus_strand else 0
    corrected['contig'] = read_final['chr']
    corrected['pos'] = corrected_start_pos + 1  # convert to 1-based
    corrected['mapq'] = read_final['mapq']
    corrected['cigar'] = final_cigar
    corrected['contig_next'] = '*'
    corrected['pos_next'] = 0
    corrected['template_len'] = 0
    corrected['sequence'] = orig_alignment['sequence']
    corrected['quality'] = orig_alignment['quality']

    corrected['new_basename'] = orig_alignment['new_basename']

    return corrected


def flush_buffered_alignment_writes(buffered_writes_by_file_path):
    for sam_path, buffered_writes in buffered_writes_by_file_path.items():
        with open(sam_path, 'at') as handle:
            for line in buffered_writes:
                handle.write(line)

        buffered_writes_by_file_path[sam_path] = list()


def write_alignment(alignment, output_sam_paths, buffered_writes_by_file_path,
                    batch_size):
    new_basename = alignment['new_basename']
    sam_path = output_sam_paths[new_basename]
    buffered_writes = buffered_writes_by_file_path.get(sam_path)
    if buffered_writes is None:
        buffered_writes = list()
        buffered_writes_by_file_path[sam_path] = buffered_writes

    columns = list()
    columns.append(alignment['read_id'])
    columns.append(alignment['flag'])
    columns.append(alignment['contig'])
    columns.append(alignment['pos'])
    columns.append(alignment['mapq'])
    columns.append(alignment['cigar'])
    columns.append(alignment['contig_next'])
    columns.append(alignment['pos_next'])
    columns.append(alignment['template_len'])
    columns.append(alignment['sequence'])
    columns.append(alignment['quality'])
    line = '{}\n'.format('\t'.join([str(x) for x in columns]))
    buffered_writes.append(line)

    if len(buffered_writes) >= batch_size:
        with open(sam_path, 'at') as handle:
            for line in buffered_writes:
                handle.write(line)

        buffered_writes_by_file_path[sam_path] = list()


def sort_output_sams_thread(in_queue):
    while True:
        sam_path = try_get_from_queue_with_short_wait(in_queue)
        if sam_path is None:
            return

        print_with_timestamp('Sorting: {}'.format(sam_path))
        temp_path = '{}.tmp'.format(sam_path)
        sort_sam(sam_path, temp_path)
        shutil.move(temp_path, sam_path)


def sort_output_sams(output_sam_paths, num_threads):
    num_jobs = len(output_sam_paths)
    threads = list()
    thread_inputs = multiprocessing.Queue(num_jobs)
    for sam_path in output_sam_paths.values():
        thread_inputs.put(sam_path)

    for _ in range(num_threads):
        thread = multiprocessing.Process(target=sort_output_sams_thread,
                                         args=(thread_inputs, ))
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    thread_inputs.close()


def sort_sam(in_path, out_path):

    command = [
        'samtools', 'sort', '-o', out_path, '--output-fmt', 'SAM', in_path
    ]
    subprocess.run(command, check=True)


def sort_sam_by_read_name(in_path, out_path):
    command = [
        'samtools', 'sort', '-n', '-o', out_path, '--output-fmt', 'SAM',
        in_path
    ]
    subprocess.run(command, check=True)


def create_out_directories(out_dir):
    if os.path.exists(out_dir):
        if not (os.path.isdir(out_dir) and len(os.listdir(out_dir)) == 0):
            raise Exception('{} already exists'.format(out_dir))
    else:
        print_with_timestamp('mkdir -p {}'.format(out_dir))
        os.makedirs(out_dir)

    for dir_name in ['orig_sams', 'indices', 'read_finals']:
        os.mkdir(os.path.join(out_dir, dir_name))


def cleanup_temp_files(out_dir):
    print_with_timestamp('Cleaning up temporary files')
    for dir_name in ['orig_sams', 'indices', 'read_finals']:
        path = os.path.join(out_dir, dir_name)
        if not os.path.exists(path):
            continue

        print_with_timestamp('rm -r {}'.format(path))
        shutil.rmtree(path)


def correct_alignments_from_main_thread(result_to_queue, existing_arguments,
                                        inputs, outputs, contig_cache,
                                        parasail_api, aligner_window_size):
    if result_to_queue is not None:
        if not try_put_to_queue_with_short_wait(outputs, result_to_queue):
            return result_to_queue, existing_arguments

        result_to_queue = None

    while True:
        if existing_arguments is not None:
            arguments = existing_arguments
            existing_arguments = None
        else:
            arguments = try_get_from_queue_with_short_wait(inputs)
            if arguments is None:
                return result_to_queue, existing_arguments

        corrected_alignment = correct_alignment(arguments['read_final'],
                                                arguments['orig_alignment'],
                                                arguments['contig_name'],
                                                contig_cache, parasail_api,
                                                aligner_window_size)
        if not try_put_to_queue_with_short_wait(outputs, corrected_alignment):
            return corrected_alignment, existing_arguments


def create_contig_cache(contig_queue):
    return {'sequence': None, 'name': None, 'queue': contig_queue}


def correct_alignment_thread(inputs, outputs, signals, contig_queue,
                             parasail_api, aligner_window_size):
    contig_cache = create_contig_cache(contig_queue)
    while True:
        signal = try_get_from_queue_without_wait(signals)
        if signal is not None:
            return

        arguments = try_get_from_queue_with_short_wait(inputs)
        if arguments is None:
            continue

        corrected_alignment = correct_alignment(arguments['read_final'],
                                                arguments['orig_alignment'],
                                                arguments['contig_name'],
                                                contig_cache, parasail_api,
                                                aligner_window_size)
        outputs.put(corrected_alignment)


def write_alignments_from_main_thread(result_to_queue, inputs, progress,
                                      output_sam_paths,
                                      buffered_writes_by_file_path, batch_size,
                                      progress_every_n):
    while True:
        if result_to_queue is not None:
            corrected_alignment = result_to_queue
            result_to_queue = None
        else:
            corrected_alignment = try_get_from_queue_with_short_wait(inputs)
            if corrected_alignment is None:
                break

        write_alignment(corrected_alignment, output_sam_paths,
                        buffered_writes_by_file_path, batch_size)
        progress.value += 1
        if (progress.value % progress_every_n) == 0:
            print_with_timestamp('Corrected {} alignments'.format(
                progress.value))

    flush_buffered_alignment_writes(buffered_writes_by_file_path)
    return result_to_queue


def write_alignment_thread(inputs, progress, signals, output_sam_paths,
                           buffered_writes_by_file_path, batch_size,
                           progress_every_n):
    while True:
        signal = try_get_from_queue_without_wait(signals)
        if signal is not None:
            flush_buffered_alignment_writes(buffered_writes_by_file_path)
            return

        corrected_alignment = try_get_from_queue_with_short_wait(inputs)
        if corrected_alignment is None:
            continue

        write_alignment(corrected_alignment, output_sam_paths,
                        buffered_writes_by_file_path, batch_size)
        progress.value += 1
        if (progress.value % progress_every_n) == 0:
            print_with_timestamp('Corrected {} alignments'.format(
                progress.value))


def process_corrections_from_main_thread(
        correct_arguments, writer_progress, num_corrections_started,
        has_writer_thread, writer_thread, aligner_threads, aligner_inputs,
        aligner_outputs, output_sam_paths, buffered_writes_by_file_path,
        alignment_batch_size, progress_every_n, contig_cache, parasail_api,
        aligner_window_size):
    # In order to avoid a deadlock, the main thread may need to
    # hold onto a corrected alignment while it waits for the writer thread
    # (or works as the writer thread itself)
    result_to_queue = None
    while ((writer_progress.value != num_corrections_started)
           or (result_to_queue is not None)
           or (correct_arguments is not None)):
        if not has_writer_thread:
            result_to_queue = write_alignments_from_main_thread(
                result_to_queue, aligner_outputs, writer_progress,
                output_sam_paths, buffered_writes_by_file_path,
                alignment_batch_size, progress_every_n)
        else:
            raise_exception_if_thread_error(writer_thread)

        result_to_queue, correct_arguments = (
            correct_alignments_from_main_thread(result_to_queue,
                                                correct_arguments,
                                                aligner_inputs,
                                                aligner_outputs, contig_cache,
                                                parasail_api,
                                                aligner_window_size))
        for thread in aligner_threads:
            raise_exception_if_thread_error(thread)


# The main thread is the reader thread.
# There can also be a writer thread and many aligner threads.
# The main thread will also function as writer and aligner if necessary.
def decide_threads_per_task(num_threads):
    if num_threads <= 1:
        num_aligner_threads = 0
        has_writer_thread = False
    elif num_threads == 2:
        num_aligner_threads = 1
        has_writer_thread = False
    else:
        num_aligner_threads = num_threads - 2
        has_writer_thread = True

    return num_aligner_threads, has_writer_thread


def create_contig_queues(num_queues, contig_order):
    queues = list()
    num_contigs = len(contig_order)
    for _ in range(num_queues):
        contig_queue = multiprocessing.Queue(num_contigs)
        queues.append(contig_queue)

    return queues


def start_aligner_threads(num_aligner_threads, aligner_inputs, aligner_outputs,
                          aligner_signals, contig_queues, parasail_api,
                          aligner_window_size):
    aligner_threads = list()
    for thread_i in range(num_aligner_threads):
        contig_queue = contig_queues[thread_i]
        aligner_thread = multiprocessing.Process(
            target=correct_alignment_thread,
            args=(aligner_inputs, aligner_outputs, aligner_signals,
                  contig_queue, parasail_api, aligner_window_size))
        aligner_threads.append(aligner_thread)
        aligner_thread.start()

    return aligner_threads


def maybe_start_writer_thread(has_writer_thread, aligner_outputs,
                              writer_progress, writer_signal, output_sam_paths,
                              buffered_writes_by_file_path,
                              alignment_batch_size, progress_every_n):
    if not has_writer_thread:
        return None

    writer_thread = multiprocessing.Process(
        target=write_alignment_thread,
        args=(aligner_outputs, writer_progress, writer_signal,
              output_sam_paths, buffered_writes_by_file_path,
              alignment_batch_size, progress_every_n))
    writer_thread.start()
    return writer_thread


def correct_alignments_cleanup(aligner_threads, aligner_inputs,
                               aligner_outputs, aligner_signals,
                               has_writer_thread, writer_thread, writer_signal,
                               contig_queues):
    writer_signal.put(True)
    for _ in aligner_threads:
        aligner_signals.put(True)

    # Threads will not be joinable until they flush all data to the
    # queue. Throw away any data on the queues to make room for data
    # that needs to be flushed.
    open_queues = [aligner_inputs, aligner_outputs]
    open_queues.extend(contig_queues)
    for open_queue in open_queues:
        while True:
            value = try_get_from_queue_with_short_wait(open_queue)
            if value is None:
                break

    for aligner_thread in aligner_threads:
        aligner_thread.join()

    if has_writer_thread:
        writer_thread.join()


def correct_alignments(fasta_path, contig_order, sorted_read_final_paths,
                       indexed_orig_sam_paths, output_sam_paths,
                       alignment_batch_size, parasail_api, aligner_window_size,
                       progress_every_n, num_threads):
    num_aligner_threads, has_writer_thread = decide_threads_per_task(
        num_threads)

    aligner_inputs = multiprocessing.Queue(alignment_batch_size)
    aligner_outputs = multiprocessing.Queue(alignment_batch_size)
    aligner_signals = multiprocessing.Queue(num_aligner_threads)
    contig_queues = create_contig_queues(num_aligner_threads + 1, contig_order)
    writer_progress = multiprocessing.sharedctypes.Value(ctypes.c_int64, 0)
    writer_signal = multiprocessing.Queue(1)
    buffered_writes_by_file_path = dict()
    aligner_threads = start_aligner_threads(num_aligner_threads,
                                            aligner_inputs, aligner_outputs,
                                            aligner_signals, contig_queues,
                                            parasail_api, aligner_window_size)
    writer_thread = maybe_start_writer_thread(has_writer_thread,
                                              aligner_outputs, writer_progress,
                                              writer_signal, output_sam_paths,
                                              buffered_writes_by_file_path,
                                              alignment_batch_size,
                                              progress_every_n)
    try:
        correct_alignments_main(
            fasta_path, contig_order, contig_queues, sorted_read_final_paths,
            indexed_orig_sam_paths, alignment_batch_size, aligner_inputs,
            writer_progress, has_writer_thread, writer_thread, aligner_threads,
            aligner_outputs, output_sam_paths, buffered_writes_by_file_path,
            progress_every_n, parasail_api, aligner_window_size)
    finally:
        correct_alignments_cleanup(aligner_threads, aligner_inputs,
                                   aligner_outputs, aligner_signals,
                                   has_writer_thread, writer_thread,
                                   writer_signal, contig_queues)


def correct_alignments_main(fasta_path, contig_order, contig_queues,
                            sorted_read_final_paths, indexed_orig_sam_paths,
                            alignment_batch_size, aligner_inputs,
                            writer_progress, has_writer_thread, writer_thread,
                            aligner_threads, aligner_outputs, output_sam_paths,
                            buffered_writes_by_file_path, progress_every_n,
                            parasail_api, aligner_window_size):
    main_thread_contig_queue = contig_queues[-1]
    num_corrections_started = 0
    contig_cache = create_contig_cache(main_thread_contig_queue)
    with open(fasta_path, 'rt') as fasta_handle:
        fasta_details = {'handle': fasta_handle, 'next_line': None}
        for contig_name in contig_order:
            contig_sequence = read_contig(fasta_details, contig_name)
            contig_cache_entry = {
                'name': contig_name,
                'sequence': contig_sequence
            }
            for contig_queue in contig_queues:
                contig_queue.put(contig_cache_entry)

            # There might not be reads for every contig in the fasta
            contig_read_final_path = sorted_read_final_paths.get(contig_name)
            if contig_read_final_path is None:
                continue

            with open(contig_read_final_path, 'rt') as read_final_handle:
                previous_line = None
                for read_final_line in read_final_handle:
                    # There could be duplicate read final lines if the
                    # output is from the snakemake which copies files
                    # from c_work_dir to q_work_dir.
                    if (((previous_line is not None)
                         and (read_final_line == previous_line))):
                        continue

                    previous_line = read_final_line
                    read_final = parse_sorted_read_final_line(read_final_line)
                    orig_alignment = find_orig_alignment(
                        read_final, indexed_orig_sam_paths,
                        alignment_batch_size)
                    correct_arguments = {
                        'read_final': read_final,
                        'orig_alignment': orig_alignment,
                        'contig_name': contig_name
                    }
                    num_corrections_started += 1
                    if try_put_to_queue_with_short_wait(
                            aligner_inputs, correct_arguments):
                        correct_arguments = None

                    # If the aligner_inputs queue is full, or if there is
                    # a full batch of alignments waiting to be written,
                    # then the main thread can do other work before
                    # adding any new alignments to the queue.
                    num_unfinished = num_corrections_started - writer_progress.value
                    if (((correct_arguments is not None)
                         or (num_unfinished >= alignment_batch_size))):
                        process_corrections_from_main_thread(
                            correct_arguments, writer_progress,
                            num_corrections_started, has_writer_thread,
                            writer_thread, aligner_threads, aligner_inputs,
                            aligner_outputs, output_sam_paths,
                            buffered_writes_by_file_path, alignment_batch_size,
                            progress_every_n, contig_cache, parasail_api,
                            aligner_window_size)

    process_corrections_from_main_thread(
        correct_arguments, writer_progress, num_corrections_started,
        has_writer_thread, writer_thread, aligner_threads, aligner_inputs,
        aligner_outputs, output_sam_paths, buffered_writes_by_file_path,
        alignment_batch_size, progress_every_n, contig_cache, parasail_api,
        aligner_window_size)

    print_with_timestamp('Corrected {} alignments'.format(
        writer_progress.value))


def raise_exception_if_thread_error(thread):
    if not thread.is_alive():
        raise Exception('{} exited early'.format(thread.name))


class ParasailResultStruct(ctypes.Structure):
    _fields_ = [('score', ctypes.c_int32), ('end_query', ctypes.c_int32),
                ('end_ref', ctypes.c_int32), ('flag', ctypes.c_int32),
                ('extra', ctypes.c_void_p)]


class ParasailCigarStruct(ctypes.Structure):
    _fields_ = [('seq', ctypes.POINTER(ctypes.c_uint32)),
                ('len', ctypes.c_int32), ('beg_query', ctypes.c_int32),
                ('beg_ref', ctypes.c_int32)]


class ParasailMatrixStruct(ctypes.Structure):
    _fields_ = [('name', ctypes.POINTER(ctypes.c_char)),
                ('matrix', ctypes.POINTER(ctypes.c_int32)),
                ('mapper', ctypes.POINTER(ctypes.c_int32)),
                ('size', ctypes.c_int32), ('max', ctypes.c_int32),
                ('min', ctypes.c_int32),
                ('user_matrix', ctypes.POINTER(ctypes.c_int32)),
                ('type', ctypes.c_int32), ('length', ctypes.c_int32),
                ('alphabet', ctypes.POINTER(ctypes.c_char)),
                ('query', ctypes.POINTER(ctypes.c_char))]


class ParasailProfileStruct(ctypes.Structure):
    _fields_ = [('s1', ctypes.POINTER(ctypes.c_char)),
                ('s1Len', ctypes.c_int32),
                ('matrix', ctypes.POINTER(ParasailMatrixStruct)),
                ('profile8', ctypes.c_void_p), ('profile16', ctypes.c_void_p),
                ('profile32', ctypes.c_void_p), ('profile64', ctypes.c_void_p),
                ('free', ctypes.c_void_p), ('stop', ctypes.c_int32)]


def find_parasail_library_path():
    so_path = ctypes.util.find_library('parasail')
    if so_path is not None:
        return so_path

    try:
        script_path = __file__
    except NameError:
        return None

    abs_script_path = os.path.abspath(script_path)
    abs_script_dir = os.path.dirname(abs_script_path)
    return os.path.join(abs_script_dir, 'libparasail.so')


class ParasailApi:
    def __init__(self, so_path):
        if so_path:
            so_path = os.path.abspath(so_path)
        else:
            so_path = find_parasail_library_path()

        self._parasail = ctypes.cdll.LoadLibrary(so_path)
        self._expose_c_api()
        self._set_parameters()

    def _expose_c_api(self):
        self.cigar_decode_len = self._parasail.parasail_cigar_decode_len
        self.cigar_decode_len.argtypes = [ctypes.c_uint32]
        self.cigar_decode_len.restype = ctypes.c_uint32

        self.cigar_decode_op = self._parasail.parasail_cigar_decode_op
        self.cigar_decode_op.argtypes = [ctypes.c_uint32]
        self.cigar_decode_op.restype = ctypes.c_char

        self.cigar_free = self._parasail.parasail_cigar_free
        self.cigar_free.argtypes = [ctypes.POINTER(ParasailCigarStruct)]
        self.cigar_free.restype = None

        self.matrix_create = self._parasail.parasail_matrix_create
        self.matrix_create.argtypes = [
            ctypes.POINTER(ctypes.c_char), ctypes.c_int32, ctypes.c_int32
        ]
        self.matrix_create.restype = ctypes.POINTER(ParasailMatrixStruct)

        self.matrix_free = self._parasail.parasail_matrix_free
        self.matrix_free.argtypes = [ctypes.POINTER(ParasailMatrixStruct)]
        self.matrix_free.restype = None

        self.matrix_set_value = self._parasail.parasail_matrix_set_value
        self.matrix_set_value.argtypes = [
            ctypes.POINTER(ParasailMatrixStruct), ctypes.c_int32,
            ctypes.c_int32, ctypes.c_int32
        ]
        self.matrix_set_value.restype = None

        self.nw_trace_striped_profile_16 = (
            self._parasail.parasail_nw_trace_striped_profile_16)
        self.nw_trace_striped_profile_16.argtypes = [
            ctypes.POINTER(ParasailProfileStruct),
            ctypes.POINTER(ctypes.c_char), ctypes.c_int32, ctypes.c_int32,
            ctypes.c_int32
        ]
        self.nw_trace_striped_profile_16.restype = ctypes.POINTER(
            ParasailResultStruct)

        self.profile_create_16 = self._parasail.parasail_profile_create_16
        self.profile_create_16.argtypes = [
            ctypes.POINTER(ctypes.c_char), ctypes.c_int32,
            ctypes.POINTER(ParasailMatrixStruct)
        ]
        self.profile_create_16.restype = ctypes.POINTER(ParasailProfileStruct)

        self.profile_free = self._parasail.parasail_profile_free
        self.profile_free.argtypes = [ctypes.POINTER(ParasailProfileStruct)]
        self.profile_free.restype = None

        self.result_free = self._parasail.parasail_result_free
        self.result_free.argtypes = [ctypes.POINTER(ParasailResultStruct)]
        self.result_free.restype = None

        self.result_get_cigar = self._parasail.parasail_result_get_cigar
        self.result_get_cigar.argtypes = [
            ctypes.POINTER(ParasailResultStruct),
            ctypes.POINTER(ctypes.c_char), ctypes.c_int32,
            ctypes.POINTER(ctypes.c_char), ctypes.c_int32,
            ctypes.POINTER(ParasailMatrixStruct)
        ]
        self.result_get_cigar.restype = ctypes.POINTER(ParasailCigarStruct)

        self.result_get_end_query = self._parasail.parasail_result_get_end_query
        self.result_get_end_query.argtypes = [
            ctypes.POINTER(ParasailResultStruct)
        ]
        self.result_get_end_query.restype = ctypes.c_int32

        self.result_get_end_ref = self._parasail.parasail_result_get_end_ref
        self.result_get_end_ref.argtypes = [
            ctypes.POINTER(ParasailResultStruct)
        ]
        self.result_get_end_ref.restype = ctypes.c_int32

        self.result_get_score = self._parasail.parasail_result_get_score
        self.result_get_score.argtypes = [ctypes.POINTER(ParasailResultStruct)]
        self.result_get_score.restype = ctypes.c_int32

    def _compute_safe_sequence_length(self):
        # * The 16 bit alignment function seems to have a score limit of:
        #   [-(2**15), 2**15] = [-32768, 32768].
        # * For 2 sequences being aligned (x and y) the algorithm conceptually
        #   has a table with dimensions len(x) by len(y)
        # * A bound on the score is (len(x) + len(y)) * the_biggest_parameter
        biggest_param = max([
            abs(self.match),
            abs(self.mismatch),
            abs(self.gap_open),
            abs(self.gap_extension)
        ])
        score_cap = 2**15
        safe_length = score_cap / (2 * biggest_param)
        safe_length = int(safe_length)
        return safe_length

    def _set_parameters(self):
        self.match = 3  # espresso: uses 5
        self.mismatch = -2  # espresso: -4
        self.gap_open = 3  # espresso: 8
        self.gap_extension = 1  # espresso: 6
        self.safe_sequence_len = self._compute_safe_sequence_length()

        self.alphabet = 'ACGTN'
        self.ambiguous_base = 'N'
        self.score_matrix_p = self.matrix_create(
            create_c_string(self.alphabet), self.match, self.mismatch)
        for alpha_i, alpha in enumerate(self.alphabet):
            if alpha == self.ambiguous_base:
                for alpha_j in range(len(self.alphabet)):
                    self.matrix_set_value(self.score_matrix_p, alpha_i,
                                          alpha_j, 0)
                    self.matrix_set_value(self.score_matrix_p, alpha_j,
                                          alpha_i, 0)

    def cleanup_parameters(self):
        self.matrix_free(self.score_matrix_p)


def main():
    args = parse_args()

    abs_out_dir = os.path.abspath(args.out_dir)
    abs_fasta = os.path.abspath(args.fasta)
    abs_samples_tsv = os.path.abspath(args.samples_tsv)
    abs_espresso_out_dir = os.path.abspath(args.espresso_out_dir)

    parasail_api = ParasailApi(args.libparasail_so_path)

    create_out_directories(abs_out_dir)
    contig_order = get_contig_order_from_fasta(abs_fasta)
    orig_sam_or_bam_paths = parse_samples_tsv(abs_samples_tsv)
    indexed_orig_sam_paths = index_sam_or_bam_paths(orig_sam_or_bam_paths,
                                                    abs_out_dir,
                                                    args.num_threads)
    orig_read_final_paths = find_read_final_paths(abs_espresso_out_dir)
    sorted_read_final_paths = parse_and_sort_read_final_paths(
        orig_read_final_paths, abs_out_dir, args.sort_memory_buffer_size)
    output_sam_paths = create_output_sams(indexed_orig_sam_paths, contig_order,
                                          abs_out_dir)
    correct_alignments(abs_fasta, contig_order, sorted_read_final_paths,
                       indexed_orig_sam_paths, output_sam_paths,
                       args.alignment_batch_size, parasail_api,
                       args.aligner_window_size, args.progress_every_n,
                       args.num_threads)
    sort_output_sams(output_sam_paths, args.num_threads)
    parasail_api.cleanup_parameters()
    cleanup_temp_files(abs_out_dir)


if __name__ == '__main__':
    main()
