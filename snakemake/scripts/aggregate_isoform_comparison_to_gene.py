import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Aggregate isoform abundance comparison to gene level'))
    parser.add_argument('--isoform-tsv',
                        required=True,
                        help='isoform output from compare_espresso_output.py')
    parser.add_argument('--gene-tsv',
                        required=True,
                        help='where to write gene level output')
    parser.add_argument(
        '--gtf-a',
        required=True,
        help='the updated.gtf from the 1st set of output files')
    parser.add_argument(
        '--gtf-b',
        required=True,
        help='the updated.gtf from the 2nd set of output files')

    return parser.parse_args()


def write_tsv_line(handle, columns):
    handle.write('{}\n'.format('\t'.join(columns)))


def read_tsv_line(line):
    return line.rstrip('\n').split('\t')


def parse_gtf_attributes_str(attributes_str):
    attributes = dict()
    parts = attributes_str.strip().split(';')
    for part in parts:
        if part == '':
            continue

        key_and_value = part.strip().split(' ')
        if len(key_and_value) != 2:
            raise Exception(
                'unable to parse gtf attributes: {}'.format(attributes_str))

        key, value = key_and_value
        attributes[key] = value

    return attributes


def remove_quotes(string):
    if len(string) < 2:
        return string

    if string[0] != string[-1]:
        return string

    if string[0] in ['"', "'"]:
        return string[1:-1]

    return string


def read_gtf_path(path):
    isoform_to_gene = dict()
    with open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue

            columns = read_tsv_line(line)
            feature = columns[2]
            if feature != 'transcript':
                continue

            attributes_str = columns[8]
            attributes = parse_gtf_attributes_str(attributes_str)
            transcript = attributes['transcript_id']
            transcript = remove_quotes(transcript)
            gene = attributes['gene_id']
            gene = remove_quotes(gene)
            isoform_to_gene[transcript] = gene

    return isoform_to_gene


def update_gene_abundance(gene, samples, row, abun_by_gene_by_sample):
    abun_by_sample = abun_by_gene_by_sample.get(gene)
    if not abun_by_sample:
        abun_by_sample = dict()
        abun_by_gene_by_sample[gene] = abun_by_sample

    for sample in samples:
        old_count = abun_by_sample.get(sample, 0)
        new_count = float(row[sample])
        abun_by_sample[sample] = old_count + new_count


def aggregate_isoforms_to_genes(in_path, out_path, isoform_to_gene_a,
                                isoform_to_gene_b):
    abun_by_gene_by_sample = dict()
    with open(out_path, 'wt') as out_handle:
        with open(in_path, 'rt') as in_handle:
            for line_i, line in enumerate(in_handle):
                columns = read_tsv_line(line)
                if line_i == 0:
                    headers = columns
                    if headers[:2] != ['id_a', 'id_b']:
                        raise Exception(
                            'unexpected headers: {}'.format(headers))

                    # headers are alternating sample_a, sample_b
                    samples_a = list()
                    samples_b = list()
                    is_a = True
                    for header in headers[2:]:
                        if is_a:
                            samples_a.append(header)
                            is_a = False
                        else:
                            samples_b.append(header)
                            is_a = True

                    out_handle.write(line)
                    continue

                row = dict(zip(headers, columns))
                id_a = row['id_a']
                id_b = row['id_b']
                if id_a != 'None':
                    gene_a_str = isoform_to_gene_a[id_a]
                    genes_a = gene_a_str.split(',')
                    for gene_a in genes_a:
                        update_gene_abundance(gene_a, samples_a, row,
                                              abun_by_gene_by_sample)

                if id_b != 'None':
                    gene_b_str = isoform_to_gene_b[id_b]
                    genes_b = gene_b_str.split(',')
                    for gene_b in genes_b:
                        update_gene_abundance(gene_b, samples_b, row,
                                              abun_by_gene_by_sample)

        genes = sorted(abun_by_gene_by_sample.keys())
        for gene in genes:
            abun_by_sample = abun_by_gene_by_sample[gene]
            out_columns = list()
            for header in headers:
                if header in ['id_a', 'id_b']:
                    out_columns.append(gene)
                    continue

                out_columns.append(abun_by_sample.get(header, 'None'))

            write_tsv_line(out_handle, [str(x) for x in out_columns])


def main():
    args = parse_args()
    isoform_to_gene_a = read_gtf_path(args.gtf_a)
    isoform_to_gene_b = read_gtf_path(args.gtf_b)
    aggregate_isoforms_to_genes(args.isoform_tsv, args.gene_tsv,
                                isoform_to_gene_a, isoform_to_gene_b)
    print('finished')


if __name__ == '__main__':
    main()
