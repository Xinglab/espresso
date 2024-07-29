import argparse
import os
import os.path
import shutil

OTHER_CHR_NAME = 'other'


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Split the ESPRESSO_S output to efficiently distribute'
                     ' reads among ESPRESSO_C jobs'))
    parser.add_argument('--orig-work-dir',
                        required=True,
                        help='the output directory of ESPRESSO_S')
    parser.add_argument(
        '--new-base-dir',
        required=True,
        help='a directory to write new output for the new split of ESPRESSO_C')
    parser.add_argument(
        '--target-reads-per-c',
        type=int,
        required=True,
        help='how many reads should be allocated to a single ESPRESSO_C job')
    parser.add_argument('--num-threads-per-c',
                        type=int,
                        default=1,
                        help='how many threads will each ESPRESSO_C job have')
    parser.add_argument('--genome-fasta',
                        required=True,
                        help='the .fa file to use as input to ESPRESSO')

    return parser.parse_args()


def chr_name_from_chr_groups(chr_name, chr_groups):
    if chr_name in chr_groups:
        return chr_name

    return OTHER_CHR_NAME


def sam_list_path_from_c_dir_path(path):
    return os.path.join(path, 'sam.list3')


def sj_group_all_from_work_dir(path):
    return os.path.join(path, 'SJ_group_all.fa')


def sj_simplified_from_work_dir_and_chr(path, chr_name, chr_groups):
    adjusted_chr_name = chr_name_from_chr_groups(chr_name, chr_groups)
    sj_simplified_name = '{}_SJ_simplified.list'.format(adjusted_chr_name)
    return os.path.join(path, sj_simplified_name)


# ESPRESSO_Q just needs the samples.tsv to have at least
# 1 entry for each c dir and 1 entry for each sample.
def add_new_samples_entries(new_c_dir_i, samples, new_samples_handle):
    for sample in samples:
        fake_sam_path = '{}_{}.sam'.format(new_c_dir_i, sample)
        new_samples_handle.write('{}\t{}\t{}\n'.format(fake_sam_path, sample,
                                                       new_c_dir_i))


def split_files_to_new_c_dirs(sample_by_c_dir, new_base_dir, new_samples_path,
                              target_reads_per_c, num_threads_per_c,
                              chr_groups, fasta_dir, fasta_index, genome_fasta,
                              sorted_copy_dir, indices_by_path):
    new_c_dir_details = {'all_dirs': list(), 'partial_dirs': list()}
    target_per_thread = target_reads_per_c // num_threads_per_c
    orig_sj_group_path = sj_group_all_from_work_dir(sorted_copy_dir)
    sj_group_index = indices_by_path[orig_sj_group_path]
    read_count_details = get_read_counts_by_chr_by_group(
        sample_by_c_dir, sorted_copy_dir, indices_by_path)
    read_counts_by_chr_by_group = (
        read_count_details['read_counts_by_chr_by_group'])
    c_dirs_by_group = read_count_details['c_dirs_by_group']
    sorted_chrs, sorted_groups_by_chr = sort_chrs_and_groups(
        read_counts_by_chr_by_group)
    for chr_name in sorted_chrs:
        for group in sorted_groups_by_chr[chr_name]:
            c_dirs = c_dirs_by_group[group]
            split_files_to_new_c_dirs_for_group(
                chr_name, group, c_dirs, new_base_dir, target_reads_per_c,
                target_per_thread, chr_groups, fasta_dir, fasta_index,
                genome_fasta, sorted_copy_dir, indices_by_path,
                new_c_dir_details, orig_sj_group_path, sj_group_index,
                sample_by_c_dir)

        cleanup_partials_after_group(new_c_dir_details, target_reads_per_c)

    with open(new_samples_path, 'wt') as new_samples_handle:
        for details in new_c_dir_details['all_dirs']:
            add_new_samples_entries(details['dir_i'], details['samples'],
                                    new_samples_handle)

    new_c_sub_dir_paths = list()
    for details in new_c_dir_details['all_dirs']:
        new_c_sub_dir_paths.append(details['sub_dir_path'])

    return new_c_sub_dir_paths


def cleanup_partials_after_group(new_c_dir_details, target_reads_per_c):
    new_partial_dirs = list()
    for partial in new_c_dir_details['partial_dirs']:
        if partial['read_count'] < target_reads_per_c:
            new_partial_dirs.append(partial)

    new_c_dir_details['partial_dirs'] = new_partial_dirs


def get_available_new_c_dir(new_c_dir_details, num_partials_out_for_group,
                            new_base_dir):
    partial_dirs = new_c_dir_details['partial_dirs']
    if not partial_dirs or (len(partial_dirs) <= num_partials_out_for_group):
        new_i = len(new_c_dir_details['all_dirs'])
        new_dir_path = os.path.join(new_base_dir, str(new_i))
        new_sub_dir_path = os.path.join(new_dir_path, '0')
        new_dir = {
            'dir_i': new_i,
            'chrs': list(),
            'read_count': 0,
            'read_count_by_group': dict(),
            'dir_path': new_dir_path,
            'sub_dir_path': new_sub_dir_path,
            'samples': set(),
        }
        new_c_dir_details['all_dirs'].append(new_dir)
        partial_dirs.append(new_dir)
        os.makedirs(new_sub_dir_path)
        return new_dir

    return partial_dirs[num_partials_out_for_group]


def split_files_to_new_c_dirs_for_group(
        chr_name, group, c_dirs, new_base_dir, target_reads_per_c,
        target_per_thread, chr_groups, fasta_dir, fasta_index, genome_fasta,
        sorted_copy_dir, indices_by_path, new_c_dir_details,
        orig_sj_group_path, sj_group_index, sample_by_c_dir):
    num_partials_out_for_group = 0
    current_new_c_dir_details = None
    for c_dir in c_dirs:
        sample = sample_by_c_dir[c_dir]
        c_dir_path = os.path.join(sorted_copy_dir, c_dir)
        orig_sam_list_path = sam_list_path_from_c_dir_path(c_dir_path)
        orig_sj_simplified_path = sj_simplified_from_work_dir_and_chr(
            sorted_copy_dir, chr_name, chr_groups)
        sam_index = indices_by_path[orig_sam_list_path]
        sj_simplified_index = indices_by_path[orig_sj_simplified_path]
        remaining_reads_in_c = sam_index[group]['line_count']
        with open(orig_sam_list_path, 'rt') as in_sam_handle:
            in_sam_handle.seek(sam_index[group]['start'])
            while remaining_reads_in_c:
                if current_new_c_dir_details is None:
                    current_new_c_dir_details = get_available_new_c_dir(
                        new_c_dir_details, num_partials_out_for_group,
                        new_base_dir)
                    append_sj_group_lines_for_group(group, orig_sj_group_path,
                                                    sj_group_index,
                                                    current_new_c_dir_details)
                    append_sj_simplified_lines_for_group(
                        group, chr_name, chr_groups, orig_sj_simplified_path,
                        sj_simplified_index, current_new_c_dir_details)
                    if chr_name not in current_new_c_dir_details['chrs']:
                        current_new_c_dir_details['chrs'].append(chr_name)
                        append_fasta_lines_for_chr(chr_name, genome_fasta,
                                                   fasta_index,
                                                   current_new_c_dir_details,
                                                   fasta_dir)

                remaining_to_target = (target_reads_per_c -
                                       current_new_c_dir_details['read_count'])
                new_count_by_group = (
                    current_new_c_dir_details['read_count_by_group'])
                count_by_group_value = new_count_by_group.get(group, 0)
                remaining_per_thread = target_per_thread - count_by_group_value
                remaining_for_new = min(remaining_to_target,
                                        remaining_per_thread)
                count_to_write = min(remaining_reads_in_c, remaining_for_new)
                append_sam_lines_for_group(count_to_write, in_sam_handle,
                                           current_new_c_dir_details, sample)
                remaining_reads_in_c -= count_to_write
                current_new_c_dir_details['read_count'] += count_to_write
                new_count_by_group[group] = (count_by_group_value +
                                             count_to_write)
                if remaining_for_new <= remaining_reads_in_c:
                    current_new_c_dir_details = None
                    num_partials_out_for_group += 1


def append_sj_group_lines_for_group(group, orig_sj_group_path, sj_group_index,
                                    current_new_c_dir_details):
    new_dir = current_new_c_dir_details['dir_path']
    new_sj_group_path = sj_group_all_from_work_dir(new_dir)
    # A read group may not have any lines in SJ_group_all.
    # Still open the file to create an empty file.
    group_index = sj_group_index.get(group)
    with open(new_sj_group_path, 'at') as out_handle:
        if group_index is None:
            return

        with open(orig_sj_group_path, 'rt') as in_handle:
            in_handle.seek(group_index['start'])
            for line_i in range(group_index['line_count']):
                line = in_handle.readline()
                out_handle.write(line)


def append_sj_simplified_lines_for_group(group, chr_name, chr_groups,
                                         orig_sj_simplified_path,
                                         sj_simplified_index,
                                         current_new_c_dir_details):
    new_dir = current_new_c_dir_details['dir_path']
    new_sj_simplified_path = sj_simplified_from_work_dir_and_chr(
        new_dir, chr_name, chr_groups)
    # A read group may not have any lines in its SJ_simplified.list.
    # Still open the file to create an empty file.
    group_index = sj_simplified_index.get(group)
    with open(new_sj_simplified_path, 'at') as out_handle:
        if group_index is None:
            return

        with open(orig_sj_simplified_path, 'rt') as in_handle:
            in_handle.seek(group_index['start'])
            for line_i in range(group_index['line_count']):
                line = in_handle.readline()
                out_handle.write(line)


def append_fasta_lines_for_chr(chr_name, genome_fasta, fasta_index,
                               current_new_c_dir_details, fasta_dir):
    dir_i = current_new_c_dir_details['dir_i']
    new_fasta_path = os.path.join(fasta_dir, '{}.fa'.format(dir_i))
    chr_index = fasta_index[chr_name]
    with open(new_fasta_path, 'at') as out_handle:
        with open(genome_fasta, 'rt') as in_handle:
            in_handle.seek(chr_index['start'])
            for line_i in range(chr_index['line_count']):
                line = in_handle.readline()
                out_handle.write(line)


def append_sam_lines_for_group(count_to_write, in_sam_handle,
                               current_new_c_dir_details, sample):
    current_new_c_dir_details['samples'].add(sample)
    sub_dir_path = current_new_c_dir_details['sub_dir_path']
    new_sam_path = sam_list_path_from_c_dir_path(sub_dir_path)
    with open(new_sam_path, 'at') as out_sam_handle:
        for line_i in range(count_to_write):
            line = in_sam_handle.readline()
            out_sam_handle.write(line)


def sort_chrs_and_groups(read_counts_by_chr_by_group):
    unsorted_chrs_and_counts = list()
    sorted_groups_by_chr = dict()
    for chr_name, counts_by_group in read_counts_by_chr_by_group.items():
        unsorted_groups_and_counts = list()
        for group, count in counts_by_group.items():
            unsorted_groups_and_counts.append((group, count))

        sorted_groups_and_counts = sorted(unsorted_groups_and_counts,
                                          key=lambda pair: pair[1],
                                          reverse=True)
        largest_count = sorted_groups_and_counts[0][1]
        unsorted_chrs_and_counts.append((chr_name, largest_count))
        sorted_groups_by_chr[chr_name] = [
            pair[0] for pair in sorted_groups_and_counts
        ]

    sorted_chrs_and_counts = sorted(unsorted_chrs_and_counts,
                                    key=lambda pair: pair[1],
                                    reverse=True)
    sorted_chrs = [pair[0] for pair in sorted_chrs_and_counts]
    return sorted_chrs, sorted_groups_by_chr


def get_read_counts_by_chr_by_group(sample_by_c_dir, sorted_copy_dir,
                                    indices_by_path):
    read_counts_by_chr_by_group = dict()
    c_dirs_by_group = dict()
    for c_dir_name in sample_by_c_dir.keys():
        c_dir_path = os.path.join(sorted_copy_dir, c_dir_name)
        sam_list_path = sam_list_path_from_c_dir_path(c_dir_path)
        index = indices_by_path[sam_list_path]
        for group, values in index.items():
            c_dirs = c_dirs_by_group.get(group)
            if not c_dirs:
                c_dirs = list()
                c_dirs_by_group[group] = c_dirs

            c_dirs.append(c_dir_name)
            chr_name = values['chr']
            line_count = values['line_count']
            read_counts_by_group = read_counts_by_chr_by_group.get(chr_name)
            if not read_counts_by_group:
                read_counts_by_group = dict()
                read_counts_by_chr_by_group[chr_name] = read_counts_by_group

            old_count = read_counts_by_group.get(group, 0)
            read_counts_by_group[group] = old_count + line_count

    for c_dirs in c_dirs_by_group.values():
        c_dirs.sort(key=int)

    return {
        'read_counts_by_chr_by_group': read_counts_by_chr_by_group,
        'c_dirs_by_group': c_dirs_by_group
    }


def copy_sj_list_lines(new_c_sub_dirs, orig_work_dir, sample_by_c_dir):
    # Write the lines from the original sj.list files to the new_c_sub_dirs.
    # The distribution among dirs doesn't matter because ESPRESSO_Q will combine them.
    # Each new dir needs at least 1 line to avoid a warning.
    # Put 1 line in each dir (starting from last) and the rest in the first dir.
    new_sj_list_handle = None
    new_c_sub_dir_i = len(new_c_sub_dirs) - 1
    new_c_sub_dir = new_c_sub_dirs[new_c_sub_dir_i]
    try:
        new_sj_list_handle = open(os.path.join(new_c_sub_dir, 'sj.list'), 'wt')
        for c_dir in sample_by_c_dir.keys():
            sj_list_path = os.path.join(orig_work_dir, c_dir, 'sj.list')
            with open(sj_list_path, 'rt') as in_handle:
                for line in in_handle:
                    new_sj_list_handle.write(line)
                    if new_c_sub_dir_i != 0:
                        new_c_sub_dir_i -= 1
                        new_c_sub_dir = new_c_sub_dirs[new_c_sub_dir_i]
                        new_sj_list_handle.close()
                        new_sj_list_handle = open(
                            os.path.join(new_c_sub_dir, 'sj.list'), 'wt')
    finally:
        if new_sj_list_handle:
            new_sj_list_handle.close()


def parse_orig_samples(samples_path):
    sample_by_c_dir = dict()
    with open(samples_path, 'rt') as handle:
        for line in handle:
            columns = line.strip().split('\t')
            sam, sample, c_dir_num = columns
            old_sample_name = sample_by_c_dir.get(c_dir_num)
            if old_sample_name is not None and old_sample_name != sample:
                raise Exception(
                    'c dir has multiple sample names: {} {} {} {}'.format(
                        c_dir_num, sam, sample, old_sample_name))

            sample_by_c_dir[c_dir_num] = sample

    return sample_by_c_dir


def sort_by_chr_read_orig(orig_path, line_to_chr_and_group, length_by_group,
                          temp_sorting_dir, temp_files_by_group_number):
    prev_chr = None
    temp_chr_handle = None
    group_start_offset = None
    prev_group = None
    group_line_count = None
    try:
        with open(orig_path, 'rt') as orig_handle:
            while True:
                offset = orig_handle.tell()
                line = orig_handle.readline()
                if not line:
                    if prev_group is not None:
                        group_length = offset - group_start_offset
                        length_by_group[prev_group] = [
                            group_length, group_line_count
                        ]

                    break

                chr_name, group = line_to_chr_and_group(line)
                if chr_name is None:
                    # line goes with previous
                    pass
                elif group != prev_group:
                    if prev_group is not None:
                        group_length = offset - group_start_offset
                        length_by_group[prev_group] = [
                            group_length, group_line_count
                        ]

                    group_start_offset = offset
                    group_line_count = 0
                    prev_group = group
                    if chr_name != prev_chr:
                        prev_chr = chr_name
                        if temp_chr_handle:
                            temp_chr_handle.close()

                        temp_chr_path = os.path.join(temp_sorting_dir,
                                                     chr_name)
                        temp_files_by_group_number[group] = temp_chr_path
                        temp_chr_handle = open(temp_chr_path, 'wt')

                temp_chr_handle.write(line)
                group_line_count += 1
    finally:
        if temp_chr_handle:
            temp_chr_handle.close()


def sort_by_chr(orig_path, new_path, line_to_chr_and_group, temp_sorting_dir):
    index = dict()
    length_by_group = dict()
    temp_files_by_group_number = dict()
    sort_by_chr_read_orig(orig_path, line_to_chr_and_group, length_by_group,
                          temp_sorting_dir, temp_files_by_group_number)

    sorted_group_numbers = sorted(list(temp_files_by_group_number.keys()))
    with open(new_path, 'wt') as out_handle:
        for group_number in sorted_group_numbers:
            temp_file_path = temp_files_by_group_number[group_number]
            with open(temp_file_path, 'rt') as in_handle:
                while True:
                    line = in_handle.readline()
                    if not line:
                        break

                    chr_name, group = line_to_chr_and_group(line)
                    length, line_count = length_by_group[group]
                    start_offset = out_handle.tell()
                    group_end_offset = start_offset + length
                    index[group] = {
                        'start': start_offset,
                        'end': group_end_offset,
                        'line_count': line_count,
                        'chr': chr_name,
                    }
                    out_handle.write(line)
                    for line_i in range(line_count - 1):
                        line = in_handle.readline()
                        out_handle.write(line)

            os.remove(temp_file_path)

    return index


def sj_group_all_line_to_chr_and_group(line):
    if not line.startswith('>'):
        return None, None

    colon_split = line.split(':')
    group_str = colon_split[-2]
    group = int(group_str)
    chr_name = colon_split[0][1:]
    return chr_name, group


def sam_list_line_to_chr_and_group(line):
    columns = line.strip().split('\t')
    group = int(columns[0])
    chr_name = columns[6]
    return chr_name, group


def sj_simplified_list_line_to_chr_and_group(line):
    if not line.startswith('SJ_cluster'):
        return None, None

    columns = line.strip().split('\t')
    group = int(columns[1])
    chr_name = columns[4]
    return chr_name, group


def copy_sort_and_index_some_orig_files(orig_work_dir, new_base_dir):
    result = dict()
    sorted_copy_dir = os.path.join(new_base_dir, 'sorted_copies')
    os.makedirs(sorted_copy_dir)
    result['directory'] = sorted_copy_dir

    indices_by_path = dict()
    result['indices_by_path'] = indices_by_path

    temp_sorting_dir = os.path.join(sorted_copy_dir, 'temp_sorting_dir')
    os.makedirs(temp_sorting_dir)
    orig_sj_group_path = sj_group_all_from_work_dir(orig_work_dir)
    sorted_sj_group_path = sj_group_all_from_work_dir(sorted_copy_dir)
    sj_group_index = sort_by_chr(orig_sj_group_path, sorted_sj_group_path,
                                 sj_group_all_line_to_chr_and_group,
                                 temp_sorting_dir)
    indices_by_path[sorted_sj_group_path] = sj_group_index

    orig_c_dirs = list()
    orig_sj_simplified_lists = list()
    orig_file_names = os.listdir(orig_work_dir)
    for file_name in orig_file_names:
        path = os.path.join(orig_work_dir, file_name)
        if os.path.isdir(path) and file_name.isdigit():
            orig_c_dirs.append(file_name)
        elif file_name.endswith('SJ_simplified.list'):
            orig_sj_simplified_lists.append(file_name)

    for orig_c_dir_name in orig_c_dirs:
        orig_c_dir_path = os.path.join(orig_work_dir, orig_c_dir_name)
        new_c_dir_path = os.path.join(sorted_copy_dir, orig_c_dir_name)
        os.makedirs(new_c_dir_path)
        orig_sam_list_path = sam_list_path_from_c_dir_path(orig_c_dir_path)
        new_sam_list_path = sam_list_path_from_c_dir_path(new_c_dir_path)
        sam_list_index = sort_by_chr(orig_sam_list_path, new_sam_list_path,
                                     sam_list_line_to_chr_and_group,
                                     temp_sorting_dir)
        indices_by_path[new_sam_list_path] = sam_list_index

    for sj_list_name in orig_sj_simplified_lists:
        orig_path = os.path.join(orig_work_dir, sj_list_name)
        new_path = os.path.join(sorted_copy_dir, sj_list_name)
        sj_list_index = sort_by_chr(orig_path, new_path,
                                    sj_simplified_list_line_to_chr_and_group,
                                    temp_sorting_dir)
        indices_by_path[new_path] = sj_list_index

    return result


def index_fasta(fasta):
    index = dict()
    chr_name = None
    line_count = 0
    with open(fasta, 'rt') as handle:
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line:
                if chr_name is not None:
                    index[chr_name]['line_count'] = line_count

                break

            if line.startswith('>'):
                if chr_name is not None:
                    index[chr_name]['line_count'] = line_count

                line_count = 0
                chr_name = line[1:].split()[0]
                if chr_name in index:
                    raise Exception('found {} multiple times in {}'.format(
                        chr_name, fasta))

                index[chr_name] = {
                    'start': offset,
                    'line_count': 0,
                }

            line_count += 1

    return index


def get_chr_groups_from_s_output(orig_work_dir):
    chr_groups = list()
    file_names = os.listdir(orig_work_dir)
    expected_suffix = '_SJ_simplified.list'
    for file_name in file_names:
        if not file_name.endswith(expected_suffix):
            continue

        chr_name = file_name[:-len(expected_suffix)]
        chr_groups.append(chr_name)

    return chr_groups


def create_output_dir(dir_path):
    if os.path.exists(dir_path):
        if not os.path.isdir(dir_path):
            raise Exception('new output dir already exists and is a file')

        if os.listdir(dir_path):
            raise Exception('new output dir already exists and is not empty')

        return

    os.makedirs(dir_path)


def split_espresso_s_output_for_c(orig_work_dir, new_base_dir,
                                  target_reads_per_c, num_threads_per_c,
                                  genome_fasta):
    create_output_dir(new_base_dir)
    chr_groups = get_chr_groups_from_s_output(orig_work_dir)
    fasta_dir = os.path.join(new_base_dir, 'fastas')
    os.makedirs(fasta_dir)
    fasta_index = index_fasta(genome_fasta)

    copied_details = copy_sort_and_index_some_orig_files(
        orig_work_dir, new_base_dir)
    sorted_copy_dir = copied_details['directory']
    indices_by_path = copied_details['indices_by_path']

    orig_samples_path = os.path.join(orig_work_dir, 'samples.tsv.updated')
    new_samples_path = os.path.join(new_base_dir, 'samples.tsv.updated')
    sample_by_c_dir = parse_orig_samples(orig_samples_path)
    temp_bam_to_sam_dir = os.path.join(sorted_copy_dir, 'temp_bam_to_sam')
    os.makedirs(temp_bam_to_sam_dir)

    new_c_sub_dirs = split_files_to_new_c_dirs(
        sample_by_c_dir, new_base_dir, new_samples_path, target_reads_per_c,
        num_threads_per_c, chr_groups, fasta_dir, fasta_index, genome_fasta,
        sorted_copy_dir, indices_by_path)
    copy_sj_list_lines(new_c_sub_dirs, orig_work_dir, sample_by_c_dir)
    shutil.rmtree(sorted_copy_dir)


def main():
    args = parse_args()
    split_espresso_s_output_for_c(args.orig_work_dir, args.new_base_dir,
                                  args.target_reads_per_c,
                                  args.num_threads_per_c, args.genome_fasta)


if __name__ == '__main__':
    main()
