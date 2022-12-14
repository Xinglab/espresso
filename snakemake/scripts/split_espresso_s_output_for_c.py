import argparse
import os
import os.path
import shutil
import subprocess

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
    parser.add_argument('--genome-fasta',
                        required=True,
                        help='the .fa file to use as input to ESPRESSO')
    parser.add_argument(
        '--sort-memory-buffer-size',
        default='2G',
        help='how much memory "sort" is allowed to use. Specified as {num_gb}G'
    )

    return parser.parse_args()


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


def process_fasta_sequence_with_handles(in_handle, out_handle):
    in_line = in_handle.readline()
    while in_line and (not in_line.startswith('>')):
        out_handle.write(in_line)
        in_line = in_handle.readline()

    return in_line


def chr_name_from_chr_groups(chr_name, chr_groups):
    if chr_name in chr_groups:
        return chr_name

    return OTHER_CHR_NAME


def split_fasta_path_from_chr(chr_name, fasta_dir):
    return os.path.join(fasta_dir, 'split_{}.fa'.format(chr_name))


def split_fasta_by_chr_groups(genome_fasta, chr_groups, fasta_dir):
    with open(genome_fasta, 'rt') as in_handle:
        in_line = in_handle.readline()
        while in_line:
            if in_line.startswith('>'):
                chr_name = in_line[1:].split()[0]
                chr_name = chr_name_from_chr_groups(chr_name, chr_groups)
                split_fasta_path = split_fasta_path_from_chr(
                    chr_name, fasta_dir)
                with open(split_fasta_path, 'at') as out_handle:
                    out_handle.write(in_line)
                    in_line = process_fasta_sequence_with_handles(
                        in_handle, out_handle)
            else:
                in_line = in_handle.readline()


def delete_split_fastas(fasta_dir):
    file_names = os.listdir(fasta_dir)
    for file_name in file_names:
        if file_name.startswith('split_') and file_name.endswith('.fa'):
            path = os.path.join(fasta_dir, file_name)
            os.remove(path)


def parse_orig_samples(samples_path):
    info_by_c_dir = dict()
    with open(samples_path, 'rt') as handle:
        for line in handle:
            columns = line.strip().split('\t')
            sam, sample, c_dir_num = columns
            c_dir_info = info_by_c_dir.get(c_dir_num)
            if not c_dir_info:
                c_dir_info = dict()
                info_by_c_dir[c_dir_num] = c_dir_info

            old_sample_name = c_dir_info.get(sam)
            if not old_sample_name:
                c_dir_info[sam] = sample
            elif old_sample_name != sample:
                raise Exception(
                    'c dir has conflicting info for sam: {} {} {}'.format(
                        sam, old_sample_name, sample))

    return info_by_c_dir


def get_read_id_to_sample(info_by_c_dir):
    read_id_to_sample = dict()
    for sample_by_sam in info_by_c_dir.values():
        for sam, sample in sample_by_sam.items():
            with open(sam, 'rt') as handle:
                for line in handle:
                    if line.startswith('@'):
                        continue

                    columns = line.strip().split('\t')
                    read_id = columns[0]
                    old_sample = read_id_to_sample.get(read_id)
                    if old_sample and old_sample != sample:
                        raise Exception(
                            'conflicting sample names for read id: {} {} {}'.
                            format(read_id, old_sample, sample))

                    read_id_to_sample[read_id] = sample

    return read_id_to_sample


def sam_list_path_from_c_dir_path(path):
    return os.path.join(path, 'sam.list3')


def sj_group_all_from_work_dir(path):
    return os.path.join(path, 'SJ_group_all.fa')


def initialize_state_by_c_dir(orig_info_by_c_dir, sorted_copy_dir,
                              orig_work_dir):
    state_by_c_dir = dict()
    for c_dir in orig_info_by_c_dir:
        state = dict()
        state_by_c_dir[c_dir] = state
        sam_list_path = sam_list_path_from_c_dir_path(
            os.path.join(sorted_copy_dir, c_dir))
        state['sam_list_path'] = sam_list_path
        state['sj_list_path'] = os.path.join(orig_work_dir, c_dir, 'sj.list')
        state['sam_offset'] = 0
        with open(sam_list_path, 'rt') as handle:
            line = handle.readline()
            columns = line.strip().split('\t')
            group_num = int(columns[0])
            state['group'] = group_num

    return state_by_c_dir


def get_next_group_and_c_dirs(state_by_orig_c_dir):
    group = None
    c_dirs = list()
    for c_dir, state in state_by_orig_c_dir.items():
        dir_group = state['group']
        if dir_group is None:
            # no more groups in this c_dir
            continue

        if group is None:
            group = dir_group
            c_dirs.append(c_dir)
            continue

        if dir_group < group:
            group = dir_group
            c_dirs = [c_dir]
        elif dir_group == group:
            c_dirs.append(c_dir)

    # Ensure that the reads are processed in a consistent order
    # when multiple c_dirs have reads for the same group.
    c_dirs.sort()
    return group, c_dirs


def copy_sj_simplified_group_info(group, in_handle, out_handle):
    offset_of_next_group = in_handle.tell()
    writing = False
    # need to use .readline() instead of "for line in ..."
    # to allow using .tell()
    while True:
        line = in_handle.readline()
        if not line:
            break

        columns = line.strip().split('\t')
        is_sj_cluster_line = columns[0] == 'SJ_cluster'
        if not is_sj_cluster_line:
            # line goes with previous
            if writing:
                out_handle.write(line)

            offset_of_next_group = in_handle.tell()
            continue

        found_group = int(columns[1])
        if found_group > group:
            break

        offset_of_next_group = in_handle.tell()
        if found_group < group:
            continue

        out_handle.write(line)
        writing = True

    in_handle.seek(offset_of_next_group)


def copy_sj_group_all_info(group, in_handle, out_handle):
    offset_of_next_group = in_handle.tell()
    writing = False
    # need to use .readline() instead of "for line in ..."
    # to allow using .tell()
    while True:
        line = in_handle.readline()
        if not line:
            break

        chr_name, found_group = sj_group_all_line_to_chr_and_group(line)
        if chr_name is None:
            # line goes with previous
            if writing:
                out_handle.write(line)

            offset_of_next_group = in_handle.tell()
            continue

        if found_group > group:
            break

        offset_of_next_group = in_handle.tell()
        if found_group < group:
            continue

        out_handle.write(line)
        writing = True

    in_handle.seek(offset_of_next_group)


def split_read_ids_from_temp_file(temp_path, sam_dir, new_c_dir_i,
                                  new_samples_handle, sort_memory_buffer_size):
    temp_sorted_path = os.path.join(sam_dir,
                                    '{}_read_id.sorted'.format(new_c_dir_i))
    sort_command = [
        'sort',
        temp_path,
        '--output',
        temp_sorted_path,
        '--temporary-directory',
        sam_dir,
        '--buffer-size',
        sort_memory_buffer_size,
    ]
    subprocess.run(sort_command, check=True)
    os.remove(temp_path)

    current_sample = None
    out_handle = None
    try:
        with open(temp_sorted_path, 'rt') as in_handle:
            for line in in_handle:
                columns = line.strip().split('\t')
                sample = columns[0]
                read_id = columns[1]
                if current_sample != sample:
                    current_sample = sample
                    if out_handle:
                        out_handle.close()

                    out_path = os.path.join(
                        sam_dir, '{}_{}.sam'.format(new_c_dir_i, sample))
                    out_handle = open(out_path, 'wt')
                    new_samples_handle.write('{}\t{}\t{}\n'.format(
                        out_path, sample, new_c_dir_i))

                # ESPRESSO_Q expects out_handle to be a sam file where
                # read_id is the first of many tab separated columns.
                # Write a tab to allow ESPRESSO_Q to parse the read_id
                out_handle.write('{}\t\n'.format(read_id))
    finally:
        if out_handle:
            out_handle.close()

    os.remove(temp_sorted_path)


def process_orig_c_dirs(orig_work_dir, new_base_dir, target_reads_per_c,
                        chr_groups, fasta_dir, sam_dir, sorted_copy_dir,
                        sort_memory_buffer_size):
    orig_samples_path = os.path.join(orig_work_dir, 'samples.tsv.updated')
    new_samples_path = os.path.join(new_base_dir, 'samples.tsv.updated')
    orig_info_by_c_dir = parse_orig_samples(orig_samples_path)
    read_id_to_sample = get_read_id_to_sample(orig_info_by_c_dir)

    new_c_sub_dirs = list()
    new_c_dir_i = None
    new_c_dir = None
    new_c_sub_dir = None
    current_read_total = 0
    current_chr = None
    state_by_orig_c_dir = initialize_state_by_c_dir(orig_info_by_c_dir,
                                                    sorted_copy_dir,
                                                    orig_work_dir)
    new_samples_handle = None
    new_sam_list_handle = None
    new_sj_group_handle = None
    old_sj_group_handle = None
    new_sj_simplified_handle = None
    old_sj_simplified_handle = None
    old_sj_simplified_chr = None
    old_sj_simplified_other_chr_offset = None
    new_fasta_handle = None
    have_copied_other_chr_fasta = False
    temp_read_id_handle = None
    temp_read_id_path = None
    try:
        new_samples_handle = open(new_samples_path, 'wt')
        old_sj_group_handle = open(sj_group_all_from_work_dir(sorted_copy_dir),
                                   'rt')
        while True:
            if (((new_c_dir_i is None)
                 or (current_read_total >= target_reads_per_c))):
                if new_c_dir_i is None:
                    new_c_dir_i = 0
                else:
                    new_sam_list_handle.close()
                    new_sam_list_handle = None
                    new_sj_group_handle.close()
                    new_sj_group_handle = None
                    new_fasta_handle.close()
                    new_fasta_handle = None
                    temp_read_id_handle.close()
                    temp_read_id_handle = None
                    split_read_ids_from_temp_file(temp_read_id_path, sam_dir,
                                                  new_c_dir_i,
                                                  new_samples_handle,
                                                  sort_memory_buffer_size)
                    new_c_dir_i += 1

                current_read_total = 0
                current_chr = None
                new_c_dir = os.path.join(new_base_dir, str(new_c_dir_i))
                os.makedirs(new_c_dir)
                new_c_sub_dir = os.path.join(new_c_dir, '0')
                os.makedirs(new_c_sub_dir)
                new_c_sub_dirs.append(new_c_sub_dir)
                new_sam_list_handle = open(
                    sam_list_path_from_c_dir_path(new_c_sub_dir), 'wt')
                new_sj_group_handle = open(
                    sj_group_all_from_work_dir(new_c_dir), 'wt')
                new_fasta_path = os.path.join(fasta_dir,
                                              '{}.fa'.format(new_c_dir_i))
                new_fasta_handle = open(new_fasta_path, 'wt')
                temp_read_id_path = os.path.join(
                    sam_dir, '{}_read_id.tmp'.format(new_c_dir_i))
                temp_read_id_handle = open(temp_read_id_path, 'wt')

            current_group, c_dirs = get_next_group_and_c_dirs(
                state_by_orig_c_dir)
            did_once_per_group = False
            if current_group is None:
                # All groups have been read.
                # Do final processing of the current new_c_dir.
                new_sam_list_handle.close()
                new_sam_list_handle = None
                new_sj_group_handle.close()
                new_sj_group_handle = None
                new_fasta_handle.close()
                new_fasta_handle = None
                temp_read_id_handle.close()
                temp_read_id_handle = None
                if current_read_total == 0:
                    # remove extra empty new_c_dir
                    new_c_sub_dirs.pop()
                    os.remove(new_fasta_path)
                    os.remove(temp_read_id_path)
                    shutil.rmtree(new_c_dir)
                else:
                    split_read_ids_from_temp_file(temp_read_id_path, sam_dir,
                                                  new_c_dir_i,
                                                  new_samples_handle,
                                                  sort_memory_buffer_size)
                break

            for c_dir in c_dirs:
                dir_state = state_by_orig_c_dir[c_dir]
                sam_list_path = dir_state['sam_list_path']
                with open(sam_list_path, 'rt') as in_sam_handle:
                    offset = dir_state['sam_offset']
                    in_sam_handle.seek(offset)
                    # need to use .readline() instead of "for line in ..."
                    # to allow using .tell()
                    while True:
                        line = in_sam_handle.readline()
                        if not line:
                            # finished reading the file
                            dir_state['group'] = None
                            dir_state['sam_offset'] = None
                            break

                        columns = line.strip().split('\t')
                        group_num = int(columns[0])
                        if group_num != current_group:
                            dir_state['group'] = group_num
                            dir_state['sam_offset'] = offset
                            break

                        offset = in_sam_handle.tell()
                        read_id = columns[2]
                        sample = read_id_to_sample.get(read_id)
                        if not sample:
                            raise Exception(
                                'read id with unknown sample: {}'.format(
                                    read_id))

                        temp_read_id_handle.write('{}\t{}\n'.format(
                            sample, read_id))
                        new_sam_list_handle.write(line)
                        current_read_total += 1

                        if not did_once_per_group:
                            # handle info that is shared for all c_dirs in a group
                            did_once_per_group = True
                            chr_name = columns[5]
                            chr_name = chr_name_from_chr_groups(
                                chr_name, chr_groups)
                            sj_simplified_name = '{}_SJ_simplified.list'.format(
                                chr_name)

                            if chr_name != current_chr:
                                current_chr = chr_name
                                split_fasta_path = split_fasta_path_from_chr(
                                    chr_name, fasta_dir)
                                # Many different contigs are under OTHER_CHR_NAME.
                                # Make sure to only copy the "other" fasta once.
                                if current_chr == OTHER_CHR_NAME:
                                    if not have_copied_other_chr_fasta:
                                        have_copied_other_chr_fasta = True
                                    else:
                                        split_fasta_path = None

                                if split_fasta_path:
                                    with open(split_fasta_path,
                                              'rt') as split_fasta_handle:
                                        for split_fasta_line in split_fasta_handle:
                                            new_fasta_handle.write(
                                                split_fasta_line)

                                if new_sj_simplified_handle:
                                    new_sj_simplified_handle.close()
                                    new_sj_simplified_handle = None

                                new_sj_open_mode = 'wt'
                                # Many different contigs are under OTHER_CHR_NAME.
                                # Resume writing where previous "other" group left off.
                                if ((current_chr == OTHER_CHR_NAME
                                     and old_sj_simplified_other_chr_offset)):
                                    new_sj_open_mode = 'at'

                                new_sj_simplified_handle = open(
                                    os.path.join(new_c_dir,
                                                 sj_simplified_name),
                                    new_sj_open_mode)

                            if old_sj_simplified_chr != chr_name:
                                if old_sj_simplified_handle:
                                    # Record offset for OTHER_CHR_NAME to allow
                                    # resuming later
                                    if old_sj_simplified_chr == OTHER_CHR_NAME:
                                        old_sj_simplified_other_chr_offset = (
                                            old_sj_simplified_handle.tell())

                                    old_sj_simplified_handle.close()
                                    old_sj_simplified_handle = None

                                old_sj_simplified_chr = chr_name
                                old_sj_simplified_handle = open(
                                    os.path.join(sorted_copy_dir,
                                                 sj_simplified_name), 'rt')
                                # Many different contigs are under OTHER_CHR_NAME.
                                # Resume reading where previous "other" group left off.
                                if ((current_chr == OTHER_CHR_NAME
                                     and old_sj_simplified_other_chr_offset)):
                                    old_sj_simplified_handle.seek(
                                        old_sj_simplified_other_chr_offset)

                            copy_sj_simplified_group_info(
                                current_group, old_sj_simplified_handle,
                                new_sj_simplified_handle)
                            copy_sj_group_all_info(current_group,
                                                   old_sj_group_handle,
                                                   new_sj_group_handle)
    finally:
        for handle in [
                new_samples_handle, new_sam_list_handle, new_sj_group_handle,
                old_sj_group_handle, new_sj_simplified_handle,
                old_sj_simplified_handle, new_fasta_handle, temp_read_id_handle
        ]:
            if handle:
                handle.close()

    copy_sj_list_lines(new_c_sub_dirs, state_by_orig_c_dir)


def copy_sj_list_lines(new_c_sub_dirs, state_by_orig_c_dir):
    # Write the lines from the original sj.list files to the new_c_sub_dirs.
    # The distribution among dirs doesn't matter because ESPRESSO_Q will combine them.
    # Each new dir needs at least 1 line to avoid a warning.
    # Put 1 line in each dir (starting from last) and the rest in the first dir.
    new_sj_list_handle = None
    new_c_sub_dir_i = len(new_c_sub_dirs) - 1
    new_c_sub_dir = new_c_sub_dirs[new_c_sub_dir_i]
    try:
        new_sj_list_handle = open(os.path.join(new_c_sub_dir, 'sj.list'), 'wt')
        for state in state_by_orig_c_dir.values():
            with open(state['sj_list_path'], 'rt') as in_handle:
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


def sort_by_chr(orig_path, new_path, line_to_chr_and_group, temp_sorting_dir):
    temp_files_by_group_number = dict()
    prev_chr = None
    temp_chr_handle = None
    try:
        with open(orig_path, 'rt') as orig_handle:
            for line in orig_handle:
                chr_name, group = line_to_chr_and_group(line)
                if chr_name is None:
                    # line goes with previous
                    pass
                elif chr_name != prev_chr:
                    prev_chr = chr_name
                    if temp_chr_handle:
                        temp_chr_handle.close()

                    temp_chr_path = os.path.join(temp_sorting_dir, chr_name)
                    temp_files_by_group_number[group] = temp_chr_path
                    temp_chr_handle = open(temp_chr_path, 'wt')

                temp_chr_handle.write(line)
    finally:
        if temp_chr_handle:
            temp_chr_handle.close()

    sorted_group_numbers = sorted(list(temp_files_by_group_number.keys()))
    with open(new_path, 'wt') as out_handle:
        for group_number in sorted_group_numbers:
            temp_file_path = temp_files_by_group_number[group_number]
            with open(temp_file_path, 'rt') as in_handle:
                for line in in_handle:
                    out_handle.write(line)

            os.remove(temp_file_path)


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
    chr_name = columns[5]
    return chr_name, group


def sj_simplified_list_line_to_chr_and_group(line):
    if not line.startswith('SJ_cluster'):
        return None, None

    columns = line.strip().split('\t')
    group = int(columns[1])
    chr_name = columns[4]
    return chr_name, group


def copy_some_orig_files_and_sort(orig_work_dir, new_base_dir):
    sorted_copy_dir = os.path.join(new_base_dir, 'sorted_copies')
    os.makedirs(sorted_copy_dir)
    temp_sorting_dir = os.path.join(sorted_copy_dir, 'temp_sorting_dir')
    os.makedirs(temp_sorting_dir)
    orig_sj_group_path = sj_group_all_from_work_dir(orig_work_dir)
    sorted_sj_group_path = sj_group_all_from_work_dir(sorted_copy_dir)
    sort_by_chr(orig_sj_group_path, sorted_sj_group_path,
                sj_group_all_line_to_chr_and_group, temp_sorting_dir)
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
        sort_by_chr(orig_sam_list_path, new_sam_list_path,
                    sam_list_line_to_chr_and_group, temp_sorting_dir)

    for sj_list_name in orig_sj_simplified_lists:
        orig_path = os.path.join(orig_work_dir, sj_list_name)
        new_path = os.path.join(sorted_copy_dir, sj_list_name)
        name_split = sj_list_name.split('_')
        chr_name = name_split[0]
        # Groups are already sorted within a chr.
        # Only OTHER_CHR_NAME has multiple chrs so other files can be
        # copied directly without sorting.
        if chr_name == OTHER_CHR_NAME:
            sort_by_chr(orig_path, new_path,
                        sj_simplified_list_line_to_chr_and_group,
                        temp_sorting_dir)
        else:
            shutil.copy(orig_path, new_path)

    return sorted_copy_dir


def split_espresso_s_output_for_c(orig_work_dir, new_base_dir,
                                  target_reads_per_c, genome_fasta,
                                  sort_memory_buffer_size):
    create_output_dir(new_base_dir)
    chr_groups = get_chr_groups_from_s_output(orig_work_dir)
    fasta_dir = os.path.join(new_base_dir, 'fastas')
    os.makedirs(fasta_dir)
    split_fasta_by_chr_groups(genome_fasta, chr_groups, fasta_dir)
    sam_dir = os.path.join(new_base_dir, 'sams')
    os.makedirs(sam_dir)
    sorted_copy_dir = copy_some_orig_files_and_sort(orig_work_dir,
                                                    new_base_dir)
    process_orig_c_dirs(orig_work_dir, new_base_dir, target_reads_per_c,
                        chr_groups, fasta_dir, sam_dir, sorted_copy_dir,
                        sort_memory_buffer_size)
    shutil.rmtree(sorted_copy_dir)
    delete_split_fastas(fasta_dir)


def main():
    args = parse_args()
    split_espresso_s_output_for_c(args.orig_work_dir, args.new_base_dir,
                                  args.target_reads_per_c, args.genome_fasta,
                                  args.sort_memory_buffer_size)


if __name__ == '__main__':
    main()
