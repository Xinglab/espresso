import argparse
import os
import os.path
import shutil


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Combine the ESPRESSO_C results orgainized the way that'
                     ' ESPRESSO_Q expects'))
    parser.add_argument(
        '--c-base-work-dir',
        required=True,
        help=
        'the directory which contains a sub-directory for each ESPRESSO_C result'
    )
    parser.add_argument('--new-base-dir',
                        required=True,
                        help='a directory to write new output for ESPRESSO_Q')

    return parser.parse_args()


def create_output_dir(dir_path):
    if os.path.exists(dir_path):
        if not os.path.isdir(dir_path):
            raise Exception('new output dir already exists and is a file')

        if os.listdir(dir_path):
            raise Exception('new output dir already exists and is not empty')

        return

    os.makedirs(dir_path)


def parse_c_dir_nums_from_samples_file(samples_path):
    c_dir_nums = set()
    with open(samples_path, 'rt') as handle:
        for line in handle:
            columns = line.strip().split('\t')
            c_dir_num_str = columns[2]
            c_dir_num = int(c_dir_num_str)
            c_dir_nums.add(c_dir_num)

    return sorted(c_dir_nums)


def find_read_final_names(c_dir):
    read_final_names = list()
    file_names = os.listdir(c_dir)
    for file_name in file_names:
        if file_name.endswith('_read_final.txt'):
            read_final_names.append(file_name)

    return read_final_names


def combine_espresso_c_output_for_q(c_base_work_dir, new_base_dir):
    create_output_dir(new_base_dir)
    samples_file_name = 'samples.tsv.updated'
    orig_samples_path = os.path.join(c_base_work_dir, samples_file_name)
    new_samples_path = os.path.join(new_base_dir, samples_file_name)
    shutil.copy(orig_samples_path, new_samples_path)
    c_dir_nums = parse_c_dir_nums_from_samples_file(new_samples_path)
    for c_dir_num in c_dir_nums:
        orig_c_dir_path = os.path.join(c_base_work_dir, str(c_dir_num), '0')
        new_c_dir_path = os.path.join(new_base_dir, str(c_dir_num))
        os.makedirs(new_c_dir_path)
        sj_list_file_name = 'sj.list'
        orig_sj_list_path = os.path.join(orig_c_dir_path, sj_list_file_name)
        new_sj_list_path = os.path.join(new_c_dir_path, sj_list_file_name)
        shutil.copy(orig_sj_list_path, new_sj_list_path)
        read_final_names = find_read_final_names(orig_c_dir_path)
        for read_final_name in read_final_names:
            orig_read_final_path = os.path.join(orig_c_dir_path,
                                                read_final_name)
            new_read_final_path = os.path.join(new_c_dir_path, read_final_name)
            shutil.copy(orig_read_final_path, new_read_final_path)


def main():
    args = parse_args()
    combine_espresso_c_output_for_q(args.c_base_work_dir, args.new_base_dir)


if __name__ == '__main__':
    main()
