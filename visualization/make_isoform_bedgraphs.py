import argparse
import os
import os.path

import file_reader


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Create bedgraphs for isoforms in the target gene that'
                     ' have sufficient counts'))
    parser.add_argument('--abundance-esp',
                        help='the *_abundance.esp file output by ESPRESSO',
                        required=True)
    parser.add_argument('--annotation-bed',
                        help='the .bed file based on the gtf',
                        required=True)
    parser.add_argument(
        '--target-gene',
        help=('the name of the gene to visualize. transcripts with name like'
              ' {target-gene}-{number} or gene_id like {target-gene}.* will'
              ' have output generated. Use the gene_id to match novel'
              ' isoforms output by ESPRESSO'),
        required=True)
    parser.add_argument(
        '--minimum-count',
        help=('only isoforms where the (normalized) count for a sample meets'
              ' the minimum count are considered'),
        type=float,
        required=True)
    parser.add_argument('--output-dir',
                        help='where to write the output files',
                        required=True)

    return parser.parse_args()


def matches_target_gene(target_gene, transcript_name, gene_id):
    name_splits = transcript_name.split('-')
    if ((len(name_splits) == 2 and name_splits[0] == target_gene
         and name_splits[1].isdigit())):
        return True

    gene_id_splits = gene_id.split('.')
    return gene_id_splits[0] == target_gene


def get_abundance_by_transcript(abundance_esp, target_sample, target_gene,
                                minimum_count):
    abundance_by_transcript = dict()
    with open(abundance_esp, 'rt') as esp_handle:
        for i, line in enumerate(esp_handle):
            if i == 0:
                initial_headers, sample_headers = (
                    file_reader.read_abundance_esp_header_line(line))
                continue

            values = file_reader.read_abundance_esp_line(line, sample_headers)

            found_minimum = False
            for sample in sample_headers:
                if values[sample] >= minimum_count:
                    found_minimum = True
                    break

            if not found_minimum:
                continue

            transcript_id = values['transcript_ID']
            transcript_name = values['transcript_name']
            gene_id = values['gene_ID']
            if matches_target_gene(target_gene, transcript_name, gene_id):
                abundance_by_transcript[transcript_id] = values[target_sample]

    return abundance_by_transcript


def make_isoform_bedgraphs_for_sample(sample, gene, abundance_by_transcript,
                                      annotation_bed, output_dir):
    with open(annotation_bed, 'rt') as annotation_handle:
        for i, line in enumerate(annotation_handle):
            if i == 0:
                file_reader.read_track_line(line)
                continue

            values = file_reader.read_bed_line(line)
            transcript_id = values['name']
            abundance = abundance_by_transcript.get(transcript_id)
            if abundance is None:
                continue

            rounded_abundance_str = '{:.2f}'.format(abundance)

            output_name = '{}_{}_{}.bed'.format(sample, gene, transcript_id)
            output_path = os.path.join(output_dir, output_name)
            is_novel = transcript_id.startswith('ESPRESSO:')
            red_rgb = '203,32,55'
            blue_rgb = '0,0,178'
            color = red_rgb if is_novel else blue_rgb
            with open(output_path, 'wt') as isoform_file_handle:
                header_values = [
                    'type=bedGraph',
                    'name="{}_{}"'.format(sample, transcript_id),
                    'visibility=dense', 'color={}'.format(color)
                ]
                header_str = 'track {}'.format(' '.join(header_values))
                isoform_file_handle.write('{}\n'.format(header_str))
                for region_start, region_end in values['block_regions']:
                    columns = [
                        values['chrom'], region_start, region_end,
                        rounded_abundance_str
                    ]
                    column_strs = [str(x) for x in columns]
                    isoform_file_handle.write('{}\n'.format(
                        '\t'.join(column_strs)))


def main():
    args = parse_args()
    isoform_output_dir = os.path.join(args.output_dir, 'target_genes')
    os.makedirs(isoform_output_dir, exist_ok=True)

    samples = file_reader.get_samples_from_abundance_esp(args.abundance_esp)
    for sample in samples:
        abundance_by_transcript = get_abundance_by_transcript(
            args.abundance_esp, sample, args.target_gene, args.minimum_count)
        make_isoform_bedgraphs_for_sample(sample, args.target_gene,
                                          abundance_by_transcript,
                                          args.annotation_bed,
                                          isoform_output_dir)


if __name__ == '__main__':
    main()
