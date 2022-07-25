import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Create a bed file for the isoforms in the gtf'))
    parser.add_argument('--updated-gtf',
                        help='the *_updated.gtf file output by ESPRESSO',
                        required=True)
    parser.add_argument('--descriptive-name',
                        help='used as a label in the visualization',
                        required=True)
    parser.add_argument('--output-bed',
                        help='where to write the output file',
                        required=True)

    return parser.parse_args()


def write_bed_line(out_bed, chrom, strand, trans_start, trans_end, trans_id,
                   exons):
    # bed start coords are 0-based
    trans_start -= 1
    thick_start = trans_start
    thick_end = trans_start
    item_rgb = '0'
    score = '0'
    block_count = len(exons)
    block_sizes = list()
    block_starts = list()
    for exon in exons:
        # bed start coords are 0-based
        exon_start = exon['start'] - 1
        exon_end = exon['end']
        block_sizes.append(exon_end - exon_start)
        block_starts.append(exon_start - trans_start)

    block_sizes = ','.join([str(x) for x in block_sizes])
    block_starts = ','.join([str(x) for x in block_starts])
    columns = [
        chrom, trans_start, trans_end, trans_id, score, strand, thick_start,
        thick_end, item_rgb, block_count, block_sizes, block_starts
    ]
    str_columns = [str(x) for x in columns]
    out_bed.write('{}\n'.format('\t'.join(str_columns)))


def remove_quotes(string):
    if len(string) < 2:
        return string

    if string[0] != '"':
        return string

    if string[-1] != '"':
        raise Exception('remove_quotes({}): no trailing quote'.format(string))

    return string[1:-1]


def get_transcript_id_from_attribute_string(string):
    attribute_pairs = string.split(';')
    for attribute_pair in attribute_pairs:
        attribute_pair = attribute_pair.strip()
        first_space = attribute_pair.find(' ')
        if first_space <= 0:
            continue

        key = attribute_pair[:first_space]
        value = attribute_pair[first_space + 1:]
        if key != 'transcript_id':
            continue

        transcript_id = remove_quotes(value)
        return transcript_id

    return None


def gtf_to_bed_with_file_handles(in_gtf, out_bed, descriptive_name):
    # header for visualizing with a genome browser
    header = 'track name="{}" visibility=dense'.format(descriptive_name)
    out_bed.write('{}\n'.format(header))

    transcript_exons = list()
    transcript_chrom = ''
    transcript_strand = ''
    transcript_start = ''
    transcript_end = ''
    transcript_id = ''
    for line in in_gtf:
        in_columns = line.strip().split('\t')
        if in_columns[0].startswith('#'):
            # skip comment lines
            continue

        in_seq = in_columns[0]
        # in_source = in_columns[1]
        in_feature = in_columns[2]
        in_start = in_columns[3]
        in_end = in_columns[4]
        # in_score = in_columns[5]
        in_strand = in_columns[6]
        # in_frame = in_columns[7]
        in_attributes = in_columns[8]

        if in_feature == 'transcript':
            if transcript_exons:
                write_bed_line(out_bed, transcript_chrom, transcript_strand,
                               transcript_start, transcript_end, transcript_id,
                               transcript_exons)

            transcript_exons = list()
            transcript_chrom = in_seq
            transcript_strand = in_strand
            transcript_start = int(in_start)
            transcript_end = int(in_end)
            transcript_id = get_transcript_id_from_attribute_string(
                in_attributes)
        elif in_feature == 'exon':
            transcript_exons.append({
                'start': int(in_start),
                'end': int(in_end)
            })

    # write the last transcript in the file
    if transcript_exons:
        write_bed_line(out_bed, transcript_chrom, transcript_strand,
                       transcript_start, transcript_end, transcript_id,
                       transcript_exons)


def main():
    args = parse_args()
    with open(args.updated_gtf, 'rt') as in_gtf:
        with open(args.output_bed, 'wt') as out_bed:
            gtf_to_bed_with_file_handles(in_gtf, out_bed,
                                         args.descriptive_name)


if __name__ == '__main__':
    main()
