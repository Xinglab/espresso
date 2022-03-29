import argparse
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Wrapper to write the tsv needed to run ESPRESSO_S'))
    parser.add_argument(
        '--sams', help='comma separated paths of sam files for each input')
    parser.add_argument('--sample-names',
                        help='comma separated sample names for each input')
    parser.add_argument('--out-tsv',
                        help='the path of the sample tsv to write')
    parser.add_argument('--command',
                        nargs=argparse.REMAINDER,
                        help='the espresso_s command to add the sample tsv to')

    return parser.parse_args()


def write_sample_tsv(sams, sample_names, file_name):
    with open(file_name, 'wt') as f_handle:
        for i, sam in enumerate(sams):
            sample_name = sample_names[i]
            columns = '\t'.join([sam, sample_name])
            f_handle.write('{}\n'.format(columns))


def run_espresso_s(args):
    write_sample_tsv(args.sams.split(','), args.sample_names.split(','),
                     args.out_tsv)
    espresso_command = args.command + ['-L', args.out_tsv]
    subprocess.run(espresso_command, check=True)


def main():
    args = parse_args()
    run_espresso_s(args)


if __name__ == '__main__':
    main()
