import argparse
import itertools
import os.path
import subprocess
import tempfile

import numpy as np

import file_reader


def parse_args():
    parser = argparse.ArgumentParser(description=(
        'Create a bigwig for each sample based on the isoform counts'))
    parser.add_argument('--abundance-esp',
                        help='the *_abundance.esp file output by ESPRESSO',
                        required=True)
    parser.add_argument('--chrom-sizes',
                        help='the .chrom.sizes file based on the genome fasta',
                        required=True)
    parser.add_argument('--annotation-bed',
                        help='the .bed file based on the gtf',
                        required=True)
    parser.add_argument('--output-dir',
                        help='where to write the output files',
                        required=True)

    return parser.parse_args()


def run_bedgraph_to_bigwig(bedgraph, chrom_sizes, bigwig_path):
    bedgraph_to_bigwig_cmd = [
        'bedGraphToBigWig', bedgraph, chrom_sizes, bigwig_path
    ]
    print(bedgraph_to_bigwig_cmd)
    subprocess.run(bedgraph_to_bigwig_cmd, check=True)


def get_chr_coverage(abundance_by_transcript, annotation_bed, chrom_name,
                     chrom_size):
    chr_coverage = np.zeros(chrom_size, dtype=float)
    with open(annotation_bed, 'rt') as annotation_handle:
        for i, line in enumerate(annotation_handle):
            if i == 0:
                file_reader.read_track_line(line)
                continue

            values = file_reader.read_bed_line(line)
            if values['chrom'] != chrom_name:
                continue

            abundance = abundance_by_transcript.get(values['name'])
            if abundance is None:
                continue

            for block_start, block_end in values['block_regions']:
                chr_coverage[block_start:block_end] += abundance

    return chr_coverage


def get_abundance_by_transcript(abundance_esp, sample):
    abundance_by_transcript = dict()
    with open(abundance_esp, 'rt') as esp_handle:
        for i, line in enumerate(esp_handle):
            if i == 0:
                initial_headers, sample_headers = (
                    file_reader.read_abundance_esp_header_line(line))
                continue

            values = file_reader.read_abundance_esp_line(line, sample_headers)
            transcript_id = values['transcript_ID']
            abundance_by_transcript[transcript_id] = values[sample]

    return abundance_by_transcript


def make_bedgraph_for_chrom(abundance_esp, annotation_bed, sample, chrom_name,
                            chrom_size, output_path):
    abundance_by_transcript = get_abundance_by_transcript(
        abundance_esp, sample)
    chr_coverage = get_chr_coverage(abundance_by_transcript, annotation_bed,
                                    chrom_name, chrom_size)

    values_with_repeat_counts = list()
    for value, values in itertools.groupby(chr_coverage):
        repeat_count = len(list(values))
        values_with_repeat_counts.append((value, repeat_count))

    with open(output_path, 'wt') as out_handle:
        chr_i = 0
        for value, repeat_count in values_with_repeat_counts:
            if value > 0:
                block_start = chr_i
                block_end = block_start + repeat_count
                columns = [chrom_name, block_start, block_end, value]
                str_columns = [str(x) for x in columns]
                out_handle.write('{}\n'.format('\t'.join(str_columns)))

            chr_i += repeat_count


def make_chrom_bedgraphs(chrom_names_to_sizes, sample, abundance_esp,
                         annotation_bed, output_dir):
    chrom_bedgraphs = list()
    for chrom_name, chrom_size in chrom_names_to_sizes.items():
        output_name = '{}_{}.bed'.format(sample, chrom_name)
        output_path = os.path.join(output_dir, output_name)
        make_bedgraph_for_chrom(abundance_esp, annotation_bed, sample,
                                chrom_name, chrom_size, output_path)
        chrom_bedgraphs.append(output_path)

    return chrom_bedgraphs


def merge_chrom_bedgraphs(chrom_bedgraphs, output_path, temp_dir_path, sample):
    combine_bedgraphs_cmd = ['cat'] + chrom_bedgraphs
    unsorted_bedgraph_path = os.path.join(temp_dir_path, 'unsorted.bed')
    with open(unsorted_bedgraph_path, 'wb') as unsorted_handle:
        print(combine_bedgraphs_cmd)
        subprocess.run(combine_bedgraphs_cmd,
                       check=True,
                       stdout=unsorted_handle)

    green_rgb = '33,189,52'
    header_values = [
        'type=bedGraph', 'color={}'.format(green_rgb),
        'name=\"{}\"'.format(sample), 'visibility=pack'
    ]
    header_str = 'track {}'.format(' '.join(header_values))
    with open(output_path, 'wt') as output_handle:
        output_handle.write('{}\n'.format(header_str))

    # bedGraphToBigWig requires sorted input
    # 1st sort alphabetically on column 1 (chrom): -k1,1
    # then sort numerically on column 2 (start): -k2,2n
    sort_bedgraph_cmd = ['sort', '-k1,1', '-k2,2n', unsorted_bedgraph_path]
    with open(output_path, 'ab') as output_handle:
        print(sort_bedgraph_cmd)
        subprocess.run(sort_bedgraph_cmd, check=True, stdout=output_handle)


def main():
    args = parse_args()
    chrom_names_to_sizes = file_reader.read_chrom_sizes(args.chrom_sizes)
    samples = file_reader.get_samples_from_abundance_esp(args.abundance_esp)
    for sample in samples:
        merged_bedgraph_name = '{}.bed'.format(sample)
        merged_bedgraph_path = os.path.join(args.output_dir,
                                            merged_bedgraph_name)
        bigwig_name = '{}.bw'.format(sample)
        bigwig_path = os.path.join(args.output_dir, bigwig_name)
        with tempfile.TemporaryDirectory(dir=args.output_dir) as temp_dir_path:
            chrom_bedgraphs = make_chrom_bedgraphs(chrom_names_to_sizes,
                                                   sample, args.abundance_esp,
                                                   args.annotation_bed,
                                                   temp_dir_path)
            merge_chrom_bedgraphs(chrom_bedgraphs, merged_bedgraph_path,
                                  temp_dir_path, sample)
            run_bedgraph_to_bigwig(merged_bedgraph_path, args.chrom_sizes,
                                   bigwig_path)


if __name__ == '__main__':
    main()
