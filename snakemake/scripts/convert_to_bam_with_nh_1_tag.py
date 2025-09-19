import argparse
import os
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add NH:i:1 and convert SAM to BAM")
    parser.add_argument('--in-sam', required=True, help='input .sam file')
    parser.add_argument('--out-bam', required=True, help='output .bam file')
    parser.add_argument(
        '--tmp-suffix',
        default='.tmp',
        help='suffix added to input sam file when adding NH tag')

    return parser.parse_args()


def convert_to_bam_with_nh_1_tag(in_sam, out_bam, tmp_suffix):
    tmp_path = '{}{}'.format(in_sam, tmp_suffix)
    nh_1 = 'NH:i:1'
    with open(in_sam, 'rt') as in_handle:
        with open(tmp_path, 'wt') as out_handle:
            for line in in_handle:
                if line.startswith('@'):
                    out_handle.write(line)
                    continue

                line = line.rstrip('\n')
                out_handle.write(line)
                out_handle.write('\t{}\n'.format(nh_1))

    command = ['samtools', 'view', '-h', '-o', out_bam, tmp_path]
    subprocess.run(command, check=True)
    os.remove(tmp_path)


def main():
    args = parse_args()
    convert_to_bam_with_nh_1_tag(args.in_sam, args.out_bam, args.tmp_suffix)


if __name__ == '__main__':
    main()
