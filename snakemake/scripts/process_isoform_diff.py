import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description=('prepare isoform diff data for plotting'))
    parser.add_argument('--in-diff', required=True,
                        help='path to a file with the output from compare_espresso_output.py')
    parser.add_argument('--out-tsv', required=True,
                        help='path to write output')

    return parser.parse_args()


def read_diffs(in_diff_path):
    diffs = dict()
    extras = {'a': list(), 'b': list()}
    with open(in_diff_path, 'rt') as handle:
        for line in handle:
            transcript_prefix = 'transcript: '
            extra_a_prefix = 'extra a transcript: '
            extra_b_prefix = 'extra b transcript: '
            extra_a_read_prefix = 'extra a read_id: '
            extra_b_read_prefix = 'extra b read_id: '
            compat_diff_prefix = 'compat isoforms diff: '
            compat_diff_underscore_prefix = 'compat_isoforms diff: '
            compat_diff_read_class_prefix = 'compat isoforms diff read_classification: '
            if line.startswith(transcript_prefix):
                without_prefix = line.strip()[len(transcript_prefix):]
                comma_split = without_prefix.split(',')
                transcript_id = comma_split[0]
                kv_pairs = comma_split[1:]
                value_i = a_diff = b_diff = None
                for kv_pair in kv_pairs:
                    k, v = [x.strip() for x in kv_pair.split(': ')]
                    if k == 'value_i':
                        value_i = v
                    elif k == 'a':
                        a_diff = v
                    elif k == 'b':
                        b_diff = v
                    else:
                        raise Exception('unknown value in transcript line: {}'.format(line))

                if value_i is None:
                    raise Exception('missing value_i in {}'.format(line))
                if a_diff is None:
                    raise Exception('missing a_diff in {}'.format(line))
                if b_diff is None:
                    raise Exception('missing b_diff in {}'.format(line))

                trans_diff = diffs.get(transcript_id)
                if trans_diff is None:
                    trans_diff = dict()
                    diffs[transcript_id] = trans_diff

                trans_diff[value_i] = {'a': float(a_diff), 'b': float(b_diff)}
            elif line.startswith(extra_a_prefix):
                transcript_id = line.strip().split()[-1]
                extras['a'].append(transcript_id)
            elif line.startswith(extra_b_prefix):
                transcript_id = line.strip().split()[-1]
                extras['b'].append(transcript_id)
            elif (line.startswith(extra_a_read_prefix)
                  or line.startswith(extra_b_read_prefix)
                  or line.startswith(compat_diff_prefix)
                  or line.startswith(compat_diff_underscore_prefix)
                  or line.startswith(compat_diff_read_class_prefix)
                  or line == 'finished\n'):
                pass  # ignoring these lines
            else:
                raise Exception('unknown line format: {}'.format(line))

    return diffs, extras


def write_tsv_line(handle, columns):
    handle.write('{}\n'.format('\t'.join([str(x) for x in columns])))


def write_diffs(diffs, out_tsv_path):
    with open(out_tsv_path, 'wt') as handle:
        headers = ['type', 'id', 'abs_diff', 'percent_diff', 'count', 'smaller', 'larger']
        write_tsv_line(handle, headers)
        for transcript, details in diffs.items():
            max_abs = -1
            max_percent = -1
            count = 0
            for value_i, a_and_b in details.items():
                count += 1
                a_diff = a_and_b['a']
                b_diff = a_and_b['b']
                if a_diff > b_diff:
                    higher = a_diff
                    lower = b_diff
                else:
                    higher = b_diff
                    lower = a_diff

                abs_diff = round(higher - lower, ndigits=2)
                if lower == 0:
                    percent_diff = -1
                else:
                    percent_diff = round(higher / lower, ndigits=2)

                if abs_diff > max_abs:
                    max_abs = abs_diff
                if percent_diff > max_percent:
                    max_percent = percent_diff

                id_str = '{}:{}'.format(transcript, value_i)
                write_tsv_line(handle, ['sample', id_str, abs_diff, percent_diff, 1, lower, higher])

            id_str = transcript
            write_tsv_line(handle, ['transcript', id_str, max_abs, max_percent, count, 0, 0])


def main():
    args = parse_args()
    diffs, extras = read_diffs(args.in_diff)
    print('transcripts only in a: {}'.format(len(extras['a'])))
    print('transcripts only in b: {}'.format(len(extras['b'])))
    write_diffs(diffs, args.out_tsv)


if __name__ == '__main__':
    main()
