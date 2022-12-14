import argparse
import os
import os.path


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Find fastq files and combine them into a single file'))
    parser.add_argument('--input-dir',
                        required=True,
                        help='path to a directory with .fastq files')
    parser.add_argument('--out-path',
                        required=True,
                        help='path to write the final output .fastq file')

    return parser.parse_args()


def combine_fastq_files(input_dir, out_handle):
    file_names = os.listdir(input_dir)
    file_names.sort()
    for file_name in file_names:
        file_path = os.path.join(input_dir, file_name)
        # some guppy versions output .fastq files under pass/ and fail/ directories
        if os.path.isdir(file_path) and file_name in ['pass', 'fail']:
            combine_fastq_files(file_path, out_handle)
            continue

        # other guppy versions output .fastq files at the top level directory
        if os.path.isfile(file_path) and file_name.endswith('.fastq'):
            with open(file_path, 'rb') as in_handle:
                for chunk in in_handle:
                    out_handle.write(chunk)


def main():
    args = parse_args()
    with open(args.out_path, 'wb') as out_handle:
        combine_fastq_files(args.input_dir, out_handle)


if __name__ == '__main__':
    main()
