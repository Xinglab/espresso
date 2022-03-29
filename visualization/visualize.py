import argparse
import os
import os.path
import subprocess
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Generate files for visualizing ESPRESSO output'))
    parser.add_argument('--genome-fasta',
                        help='the fasta input to ESPRESSO',
                        required=True)
    parser.add_argument('--updated-gtf',
                        help='the *_updated.gtf file output by ESPRESSO',
                        required=True)
    parser.add_argument('--abundance-esp',
                        help='the *_abundance.esp file output by ESPRESSO',
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
    parser.add_argument('--normalize-counts-to-cpm',
                        help='convert raw counts to counts per million',
                        action='store_true')
    parser.add_argument(
        '--descriptive-name',
        help='used as a label in the visualization and for filenames',
        required=True)
    parser.add_argument('--output-dir',
                        help='where to write visualization files',
                        required=True)

    return parser.parse_args()


def get_python_executable():
    # Try to get the absolute path of the executable for the running
    # Python interpreter.
    python_executable = sys.executable
    if not python_executable:
        # Fallback
        print(
            'Absolute path for current Python interpreter not found.'
            ' Using "python" without a full path to run scripts',
            file=sys.stderr)
        python_executable = 'python'

    return python_executable


def get_script_dir():
    script_path = os.path.realpath(__file__)
    return os.path.dirname(script_path)


def basename_without_extension(path):
    basename = os.path.basename(path)
    without_ext, ext = os.path.splitext(basename)
    return without_ext


def run_gtf_to_bed(python_executable, script_dir, updated_gtf, output_dir,
                   descriptive_name):
    gtf_to_bed_script_path = os.path.join(script_dir, 'gtf_to_bed.py')
    gtf_without_ext = basename_without_extension(updated_gtf)
    bed_basename = '{}.bed'.format(gtf_without_ext)
    bed_path = os.path.join(output_dir, bed_basename)

    gtf_to_bed_cmd = [
        python_executable, gtf_to_bed_script_path, '--updated-gtf',
        updated_gtf, '--output-bed', bed_path, '--descriptive-name',
        descriptive_name
    ]
    print(gtf_to_bed_cmd)
    subprocess.run(gtf_to_bed_cmd, check=True)
    return bed_path


def run_normalize_abundance(python_executable, script_dir, abundance_esp,
                            output_dir, descriptive_name):
    normalize_abundance_script_path = os.path.join(script_dir,
                                                   'normalize_abundance.py')
    output_path = os.path.join(output_dir,
                               '{}_cpm.esp'.format(descriptive_name))
    normalize_abundance_cmd = [
        python_executable, normalize_abundance_script_path, '--abundance-esp',
        abundance_esp, '--output-path', output_path
    ]
    print(normalize_abundance_cmd)
    subprocess.run(normalize_abundance_cmd, check=True)
    return output_path


def index_fasta(fasta):
    index_path = '{}.fai'.format(fasta)
    if os.path.exists(index_path):
        return

    index_cmd = ['samtools', 'faidx', fasta]
    print(index_cmd)
    subprocess.run(index_cmd, check=True)


def run_fasta_to_chrom_sizes(fasta, output_dir):
    fasta_without_ext = basename_without_extension(fasta)
    two_bit_basename = '{}.2bit'.format(fasta_without_ext)
    two_bit_path = os.path.join(output_dir, two_bit_basename)
    fa_to_two_bit_cmd = ['faToTwoBit', '-long', fasta, two_bit_path]
    print(fa_to_two_bit_cmd)
    subprocess.run(fa_to_two_bit_cmd, check=True)

    chrom_sizes_basename = '{}.chrom.sizes'.format(fasta_without_ext)
    chrom_sizes_path = os.path.join(output_dir, chrom_sizes_basename)
    two_bit_info_cmd = ['twoBitInfo', two_bit_path, chrom_sizes_path]
    print(two_bit_info_cmd)
    subprocess.run(two_bit_info_cmd, check=True)
    return chrom_sizes_path


def run_make_sample_bigwigs(python_executable, script_dir, abundance_esp,
                            chrom_sizes, annotation_bed, output_dir):
    make_sample_bigwigs_script_path = os.path.join(script_dir,
                                                   'make_sample_bigwigs.py')
    make_sample_bigwigs_cmd = [
        python_executable, make_sample_bigwigs_script_path, '--abundance-esp',
        abundance_esp, '--chrom-sizes', chrom_sizes, '--annotation-bed',
        annotation_bed, '--output-dir', output_dir
    ]
    print(make_sample_bigwigs_cmd)
    subprocess.run(make_sample_bigwigs_cmd, check=True)


def run_make_isoform_bedgraphs(python_executable, script_dir, abundance_esp,
                               annotation_bed, target_gene, minimum_count,
                               output_dir):
    make_isoform_bedgraphs_script_path = os.path.join(
        script_dir, 'make_isoform_bedgraphs.py')
    make_isoform_bedgraphs_cmd = [
        python_executable, make_isoform_bedgraphs_script_path,
        '--abundance-esp', abundance_esp, '--annotation-bed', annotation_bed,
        '--target-gene', target_gene, '--minimum-count',
        str(minimum_count), '--output-dir', output_dir
    ]
    print(make_isoform_bedgraphs_cmd)
    subprocess.run(make_isoform_bedgraphs_cmd, check=True)


def main():
    args = parse_args()
    python_executable = get_python_executable()
    script_dir = get_script_dir()
    os.makedirs(args.output_dir, exist_ok=True)

    index_fasta(args.genome_fasta)
    annotation_bed = run_gtf_to_bed(python_executable, script_dir,
                                    args.updated_gtf, args.output_dir,
                                    args.descriptive_name)

    normalized_esp = args.abundance_esp
    if args.normalize_counts_to_cpm:
        normalized_esp = run_normalize_abundance(python_executable, script_dir,
                                                 args.abundance_esp,
                                                 args.output_dir,
                                                 args.descriptive_name)

    chrom_sizes = run_fasta_to_chrom_sizes(args.genome_fasta, args.output_dir)
    run_make_sample_bigwigs(python_executable, script_dir, normalized_esp,
                            chrom_sizes, annotation_bed, args.output_dir)
    run_make_isoform_bedgraphs(python_executable, script_dir, normalized_esp,
                               annotation_bed, args.target_gene,
                               args.minimum_count, args.output_dir)


if __name__ == '__main__':
    main()
