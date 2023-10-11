import argparse
import os
import os.path


def parse_args():
    parser = argparse.ArgumentParser(
        description='Combine the summary files from different C steps')
    parser.add_argument(
        '--c-base-work-dir',
        required=True,
        help='the working directory which contains the summary files')
    parser.add_argument('--out-path',
                        required=True,
                        help='where to write the combined summary file')

    return parser.parse_args()


def find_summary_files(c_dir):
    summary_paths = list()
    for dir_path, dir_names, file_names in os.walk(c_dir):
        for file_name in file_names:
            path = os.path.join(dir_path, file_name)
            if file_name == 'espresso_c_summary.txt':
                summary_paths.append(path)

    return summary_paths


def parse_value(value_str):
    if value_str.isdigit():
        return int(value_str)

    return float(value_str)


def combine_summary_files(summary_paths, out_path):
    combined_results = list()
    with open(out_path, 'wt') as out_handle:
        for summary_i, summary_path in enumerate(summary_paths):
            with open(summary_path, 'rt') as in_handle:
                for line_i, line in enumerate(in_handle):
                    if line_i == 0:
                        # The 1st line is the command used to run ESPRESSO C
                        out_handle.write(line)
                        continue

                    line = line.rstrip('\n')
                    key, value_str = line.split(': ')
                    value = parse_value(value_str)
                    if summary_i == 0:
                        combined_results.append([key, value])
                    else:
                        result_entry = combined_results[line_i - 1]
                        result_entry[1] += value

        for key, value in combined_results:
            out_handle.write('{}: {}\n'.format(key, value))


def main():
    args = parse_args()
    summary_paths = find_summary_files(args.c_base_work_dir)
    combine_summary_files(summary_paths, args.out_path)


if __name__ == '__main__':
    main()
