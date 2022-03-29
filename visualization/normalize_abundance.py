import argparse

import file_reader


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Normalize the raw counts in the abundance file'))
    parser.add_argument('--abundance-esp',
                        help='the *_abundance.esp file output by ESPRESSO',
                        required=True)
    parser.add_argument('--output-path',
                        help='where to write the output file',
                        required=True)

    return parser.parse_args()


def get_totals_by_sample(in_esp_path):
    with open(in_esp_path, 'rt') as in_esp:
        for line_i, line in enumerate(in_esp):
            if line_i == 0:
                initial_headers, sample_headers = (
                    file_reader.read_abundance_esp_header_line(line))
                sample_totals = [0] * len(sample_headers)
                continue

            values = file_reader.read_abundance_esp_line(line, sample_headers)
            for i, sample in enumerate(sample_headers):
                sample_totals[i] += values[sample]

    return sample_totals


def write_columns(out_f, columns):
    str_columns = list()
    for col in columns:
        if isinstance(col, float):
            str_val = '{:.2f}'.format(col)
        else:
            str_val = str(col)

        str_columns.append(str_val)

    out_f.write('{}\n'.format('\t'.join(str_columns)))


def write_normalized_esp(totals_by_sample, in_esp_path, out_path):
    with open(in_esp_path, 'rt') as in_esp:
        with open(out_path, 'wt') as out_esp:
            for line_i, line in enumerate(in_esp):
                if line_i == 0:
                    initial_headers, sample_headers = (
                        file_reader.read_abundance_esp_header_line(line))
                    write_columns(out_esp, initial_headers + sample_headers)
                    continue

                values = file_reader.read_abundance_esp_line(
                    line, sample_headers)
                cpms = list()
                for sample_i, sample in enumerate(sample_headers):
                    sample_val = values[sample]
                    total = totals_by_sample[sample_i]
                    cpms.append((sample_val * 1e6) / total)

                initial_columns = [values[h] for h in initial_headers]
                write_columns(out_esp, initial_columns + cpms)


def main():
    args = parse_args()
    totals_by_sample = get_totals_by_sample(args.abundance_esp)
    write_normalized_esp(totals_by_sample, args.abundance_esp,
                         args.output_path)


if __name__ == '__main__':
    main()
