import copy
import math
import os
import os.path
import random
import shutil
import subprocess
import sys
import unittest

INTRON_START = 'GT'
INTRON_END = 'AG'
DEFAULT_EXON_LENGTH = (100, 200)
DEFAULT_INTRON_LENGTH = (100, 1_000)
DEFAULT_INTERGENIC_LENGTH = (1_000, 10_000)


def run_command_with_log(command, log_path, check=True):
    with open(log_path, 'wb') as log_handle:
        process = subprocess.run(command,
                                 check=check,
                                 stdout=log_handle,
                                 stderr=subprocess.STDOUT)

    return process


def run_command_with_output_and_log(command, out_path, log_path):
    with open(out_path, 'wb') as out_handle:
        with open(log_path, 'wb') as log_handle:
            process = subprocess.run(command,
                                     check=True,
                                     stdout=out_handle,
                                     stderr=log_handle)

    return process


def sort_sam_with_temp_file(unsorted_sam, temp_file, log_path):
    command = ['samtools', 'sort', '-o', temp_file, unsorted_sam]
    run_command_with_log(command, log_path)
    shutil.move(temp_file, unsorted_sam)


def remove_then_create_directories(directories):
    for directory in directories:
        if os.path.exists(directory):
            shutil.rmtree(directory)

        os.makedirs(directory)


def _get_python_executable():
    exe = sys.executable
    if not exe:
        return 'python'

    return exe


def write_tsv_line(handle, columns):
    handle.write('{}\n'.format('\t'.join([str(x) for x in columns])))


def out_file_prefix_based_on_params(num_cutoff='2', ratio_cutoff='0'):
    return 'samples_N{}_R{}'.format(num_cutoff, ratio_cutoff)


def get_espresso_novel_id(chrom, group, num):
    return 'ESPRESSO:{}:{}:{}'.format(chrom, group, num)


def write_fasta_from_chroms(fasta_path, chroms):
    with open(fasta_path, 'wt') as handle:
        for chrom in chroms:
            handle.write('>{}\n'.format(chrom.name))
            handle.write('{}\n'.format(chrom.get_sequence_and_set_coords()))


def write_gtf_from_chroms(gtf_path, chroms):
    with open(gtf_path, 'wt') as handle:
        for chrom in chroms:
            for gene in chrom.genes:
                write_gtf_from_gene(chrom.name, gene, handle)


def write_gtf_from_gene(chrom_name, gene, handle):
    source = "test"
    gene_start = gene.exons[0].start + 1  # 1-based
    gene_end = gene.exons[-1].end + 1  # 1-based
    score = '.'  # no score
    strand = gene.strand
    frame = '.'  # no frame
    gene_attributes = [
        'gene_id "{}"'.format(gene.id), 'gene_name "{}"'.format(gene.name)
    ]
    gene_attributes_str = '; '.join(gene_attributes)
    gene_columns = [
        chrom_name, source, 'gene', gene_start, gene_end, score, strand, frame,
        gene_attributes_str
    ]
    write_tsv_line(handle, gene_columns)
    for isoform in gene.isoforms:
        isoform_start = isoform.exons[0].start + 1
        isoform_end = isoform.exons[-1].end + 1
        isoform_attributes = gene_attributes + [
            'transcript_id "{}"'.format(isoform.id),
            'transcript_name "{}"'.format(isoform.name)
        ]
        isoform_attributes_str = '; '.join(isoform_attributes)
        isoform_columns = [
            chrom_name, source, 'transcript', isoform_start, isoform_end,
            score, strand, frame, isoform_attributes_str
        ]
        write_tsv_line(handle, isoform_columns)
        for exon_i, exon in enumerate(isoform.exons):
            exon_start = exon.start + 1
            exon_end = exon.end + 1
            exon_attributes = isoform_attributes + [
                'exon_number {}'.format(exon_i + 1)
            ]
            exon_attributes_str = '; '.join(exon_attributes)
            exon_columns = [
                chrom_name, source, 'exon', exon_start, exon_end, score,
                strand, frame, exon_attributes_str
            ]
            write_tsv_line(handle, exon_columns)


def write_bed_from_chroms(bed_path, chroms):
    with open(bed_path, 'wt') as handle:
        for chrom in chroms:
            for gene in chrom.genes:
                write_bed_from_gene(chrom.name, gene, handle)


def write_bed_from_gene(chrom_name, gene, handle):
    feature_name = 'sj'
    score = 0
    for isoform in gene.isoforms:
        for exon_i, exon in enumerate(isoform.exons):
            if exon_i == 0:
                continue

            upstream_exon = isoform.exons[exon_i - 1]
            sj_start = upstream_exon.end + 1  # 0-based intron start
            sj_end = exon.start  # 1-based intron end
            sj_columns = [
                chrom_name, sj_start, sj_end, feature_name, score, gene.strand
            ]
            write_tsv_line(handle, sj_columns)


def write_sam_from_alignments(sam_path,
                              chroms,
                              alignments,
                              read_id_offset,
                              use_m_cigar_op=True):
    chrom_names_sequences_and_lengths = list()
    for chrom in chroms:
        chrom_name = chrom.name
        chrom_seq = chrom.get_sequence_and_set_coords()
        chrom_length = len(chrom_seq)
        chrom_names_sequences_and_lengths.append({
            'name': chrom_name,
            'seq': chrom_seq,
            'length': chrom_length
        })

    chrom_names_sequences_and_lengths.sort(key=lambda details: details['name'])
    chrom_name_to_i = dict()
    for details_i, details in enumerate(chrom_names_sequences_and_lengths):
        chrom_name_to_i[details['name']] = details_i

    # reads are not paired
    ref_next = '*'
    pos_next = '0'
    # not computing template length
    template_length = '0'
    with open(sam_path, 'wt') as handle:
        for details in chrom_names_sequences_and_lengths:
            header_columns = [
                '@SQ', 'SN:{}'.format(details['name']),
                'LN:{}'.format(details['length'])
            ]
            write_tsv_line(handle, header_columns)

        for alignment in alignments:
            if alignment.read_id is not None:
                read_id = alignment.read_id
            else:
                read_id = 'read-{}'.format(read_id_offset)
                read_id_offset += 1

            chrom_i = chrom_name_to_i[alignment.chrom_name]
            chrom_details = chrom_names_sequences_and_lengths[chrom_i]
            if alignment.is_secondary:
                alignment.sequence = '*'
                alignment.quality = '*'
            else:
                alignment.set_sequence_from_cigar_and_chr_seq(
                    chrom_details['seq'])
                alignment.set_default_quality_from_sequence()
            if use_m_cigar_op:
                cigar_string = (
                    alignment.cigar.str_with_match_and_mismatch_as_m())
            else:
                cigar_string = (
                    alignment.cigar.str_with_unmodified_match_and_mismatch())

            columns = [
                read_id,
                alignment.flag(), alignment.chrom_name, alignment.start,
                alignment.mapq, cigar_string, ref_next, pos_next,
                template_length, alignment.sequence, alignment.quality
            ]
            write_tsv_line(handle, columns)

    return read_id_offset


def sam_flag_int(has_multiple_segments=False,
                 all_segments_aligned=True,
                 segment_unmapped=False,
                 next_segment_unmapped=False,
                 is_reversed=False,
                 next_segment_is_reversed=False,
                 is_first_segment=True,
                 is_last_segment=False,
                 is_secondary=False,
                 failed_filter=False,
                 is_duplicate=False,
                 is_supplementary=False):
    flags_int = 0
    if has_multiple_segments:
        flags_int += 1
    if all_segments_aligned:
        flags_int += 2
    if segment_unmapped:
        flags_int += 4
    if next_segment_unmapped:
        flags_int += 8
    if is_reversed:
        flags_int += 16
    if next_segment_is_reversed:
        flags_int += 32
    if is_first_segment:
        flags_int += 64
    if is_last_segment:
        flags_int += 128
    if is_secondary:
        flags_int += 256
    if failed_filter:
        flags_int += 512
    if is_duplicate:
        flags_int += 1024
    if is_supplementary:
        flags_int += 2048

    return flags_int


def perfect_alignment_for_gene_and_exon_numbers(chr_name, gene, exon_numbers):
    alignment = Alignment()
    alignment.chrom_name = chr_name
    alignment.strand = gene.strand
    start = None
    length_to_skip = 0
    for exon_i, exon in enumerate(gene.exons):
        exon_length = len(exon.sequence)
        if exon_i in exon_numbers:
            if start is None:
                start = exon.start + 1
                alignment.start = start

            if length_to_skip:
                alignment.cigar.add_skip(length_to_skip)
                length_to_skip = 0

            alignment.cigar.add_match(exon_length)
        elif start is not None:
            length_to_skip += exon_length

        if (start is not None) and (exon_i != (len(gene.exons) - 1)):
            intron = gene.introns[exon_i]
            length_to_skip += len(intron.sequence)

    return alignment


def _get_junction_i_from_cigar(cigar, junction_i):
    seen_junctions = 0
    junction_op_i = None
    for op_i, op in enumerate(cigar.operations):
        if op.char == 'N':
            if seen_junctions == junction_i:
                junction_op_i = op_i
                break

            seen_junctions += 1

    if junction_op_i is None:
        raise Exception('did not find junction_i={} in {}'.format(
            junction_i, cigar))

    return junction_op_i


def _get_exon_i_from_cigar(cigar, exon_i):
    seen_exons = 0
    exon_op_i = None
    expecting_exon = True
    for op_i, op in enumerate(cigar.operations):
        if op.char in ['M', 'X', '=']:
            if expecting_exon:
                if seen_exons == exon_i:
                    exon_op_i = op_i
                    break

                seen_exons += 1
                expecting_exon = False

        if op.char == 'N':
            expecting_exon = True

    if exon_op_i is None:
        raise Exception('did not find exon_i={} in {}'.format(exon_i, cigar))

    return exon_op_i


def add_insertion_at_junction_number(alignment, junction_i=0, insert_len=1):
    junction_op_i = _get_junction_i_from_cigar(alignment.cigar, junction_i)
    new_op = CigarOp(insert_len, 'I')
    alignment.cigar.operations.insert(junction_op_i, new_op)


def add_mismatch_into_intron_number(alignment, junction_i=0, match_len=1):
    junction_op_i = _get_junction_i_from_cigar(alignment.cigar, junction_i)
    new_op = CigarOp(match_len, 'X')
    # reduce the 'N' operation so that only the start of the intron is changed
    alignment.cigar.operations[junction_op_i].num -= match_len
    alignment.cigar.operations.insert(junction_op_i, new_op)


def add_mismatch_in_exon(alignment,
                         exon_i=0,
                         offset_from_start=1,
                         mismatch_len=1):
    exon_op_i = _get_exon_i_from_cigar(alignment.cigar, exon_i)
    existing_op = alignment.cigar.operations[exon_op_i]
    if existing_op.num <= (mismatch_len + offset_from_start):
        raise Exception('Cigar operation not long enough to add mismatch of'
                        ' length {} offset by {} in {}'.format(
                            mismatch_len, offset_from_start, existing_op))

    if offset_from_start > 0:
        existing_op.num -= offset_from_start
        before_op = CigarOp(offset_from_start, existing_op.char)
        alignment.cigar.operations.insert(exon_op_i, before_op)
        exon_op_i += 1

    existing_op.num -= mismatch_len
    mismatch_op = CigarOp(mismatch_len, 'X')
    alignment.cigar.operations.insert(exon_op_i, mismatch_op)


def add_deletion_in_exon(alignment,
                         exon_i=0,
                         offset_from_start=1,
                         deletion_length=1):
    if exon_i == 0:
        junction_op_i = -1
    else:
        junction_op_i = _get_junction_i_from_cigar(alignment.cigar, exon_i - 1)

    exon_match_op_i = junction_op_i + 1
    while exon_match_op_i < len(alignment.cigar.operations):
        if alignment.cigar.operations[exon_match_op_i].char == 'M':
            break

        exon_match_op_i += 1

    exon_match_op = alignment.cigar.operations[exon_match_op_i]
    orig_match_len = exon_match_op.num
    exon_match_op.num = offset_from_start
    new_del_op = CigarOp(deletion_length, 'D')
    new_match_length = orig_match_len - (offset_from_start + deletion_length)
    new_match_op = CigarOp(new_match_length, 'M')
    alignment.cigar.operations.insert(exon_match_op_i + 1, new_del_op)
    alignment.cigar.operations.insert(exon_match_op_i + 2, new_match_op)


def adjust_start_of_junction(alignment, junction_i=0, offset=1):
    junction_op_i = _get_junction_i_from_cigar(alignment.cigar, junction_i)
    junction_op = alignment.cigar.operations[junction_op_i]
    junction_op.num -= offset
    # assume cigar is MNM(NM...)
    pre_match_op_i = junction_op_i - 1
    pre_match_op = alignment.cigar.operations[pre_match_op_i]
    pre_match_op.num += offset


def trim_alignment_start(alignment, num_nt):
    alignment.cigar.operations[0].num -= num_nt
    alignment.start += num_nt


def trim_alignment_end(alignment, num_nt):
    alignment.cigar.operations[-1].num -= num_nt


def append_copies(values, to_copy, num_copies):
    for _ in range(num_copies):
        values.append(copy.deepcopy(to_copy))


def parse_float_with_NA(float_str):
    if float_str == 'NA':
        return math.nan

    return float(float_str)


class Cigar(object):
    def __init__(self):
        self.operations = list()

    def copy(self):
        copied = Cigar()
        copied.operations = [op.copy() for op in self.operations]
        return copied

    def str_with_match_and_mismatch_as_m(self):
        # 'X' and '=' operations are written as 'M' to be
        # compatible with minimap output.
        # Adjacent 'M' operations are combined together.
        parts = list()
        combined_m = CigarOp(0, 'M')
        for op in self.operations:
            if op.char in ['M', 'X', '=']:
                combined_m.num += op.num
                continue

            if combined_m.num:
                parts.append(str(combined_m))
                combined_m.num = 0

            parts.append(str(op))

        if combined_m.num:
            parts.append(str(combined_m))

        return ''.join(parts)

    def str_with_unmodified_match_and_mismatch(self):
        parts = list()
        for op in self.operations:
            parts.append(str(op))

        return ''.join(parts)

    def convert_alignment_match_to_sequence_match(self):
        for op in self.operations:
            if op.char == 'M':
                op.char = '='

    def __str__(self):
        return self.str_with_match_and_mismatch_as_m()

    def add_match(self, length):
        match = CigarOp(length, 'M')
        self.operations.append(match)

    def add_skip(self, length):
        skip = CigarOp(length, 'N')
        self.operations.append(skip)


class CigarOp(object):
    def __init__(self, num, char):
        self.num = num
        self.char = char

    def copy(self):
        copied = CigarOp(self.num, self.char)
        return copied

    def __str__(self):
        return '{}{}'.format(self.num, self.char)


class Alignment(object):
    def __init__(self):
        self.chrom_name = None
        self.start = None  # 1-based
        self.cigar = Cigar()
        self.mapq = 30  # default to 30 which is 0.999 accuracy
        # even for - strand the sequence written to the SAM is the + strand
        self.sequence = None
        self.quality = None
        self.strand = '+'
        self.read_id = None
        self.is_secondary = False

    def copy(self):
        copied = Alignment()
        copied.chrom_name = self.chrom_name
        copied.start = self.start
        copied.cigar = self.cigar.copy()
        copied.mapq = self.mapq
        if self.sequence is not None:
            copied.sequence = self.sequence[:]

        if self.quality is not None:
            copied.quality = self.quality[:]

        copied.strand = self.strand
        if self.read_id is not None:
            copied.read_id = self.read_id[:]

        copied.is_secondary = self.is_secondary

    def set_sequence_from_cigar_and_chr_seq(self, chr_seq):
        sequence_parts = list()
        position = self.start - 1
        for op in self.cigar.operations:
            if op.char == 'M':
                next_pos = position + op.num
                sequence_parts.append(chr_seq[position:next_pos])
                position = next_pos
            if op.char == 'I':
                # insert sequence is all 'N' for Nucleic acid
                sequence_parts.append('N' * op.num)
            if op.char == 'D':
                position = position + op.num
            if op.char == 'N':
                position = position + op.num
            if op.char == 'S':
                # soft clip sequence is all 'N'
                sequence_parts.append('N' * op.num)
            if op.char == 'H':
                pass  # no change for hard clipping
            if op.char == 'P':
                pass  # no change for padding
            if op.char == '=':
                next_pos = position + op.num
                sequence_parts.append(chr_seq[position:next_pos])
                position = next_pos
            if op.char == 'X':
                next_pos = position + op.num
                mismatch_seq = list(chr_seq[position:next_pos])
                for mismatch_i in range(op.num):
                    mismatch_seq[mismatch_i] = self._mismatch_base(
                        mismatch_seq[mismatch_i])

                sequence_parts.append(''.join(mismatch_seq))
                position = next_pos

        self.sequence = ''.join(sequence_parts)

    # arbitrary mismatch values
    def _mismatch_base(self, base):
        if base in ['A', 'a']:
            return 'C'
        if base in ['C', 'c']:
            return 'G'
        if base in ['G', 'g']:
            return 'T'
        if base in ['T', 't']:
            return 'A'

        raise Exception('could not mismatch: {}'.format(base))

    def set_default_quality_from_sequence(self):
        # 'E' is about 0.9998 accuracy
        self.quality = 'E' * len(self.sequence)

    def flag(self):
        is_minus_strand = self.strand == '-'
        return sam_flag_int(is_reversed=is_minus_strand,
                            is_secondary=self.is_secondary)


class Region(object):
    def __init__(self, sequence):
        self.start = None  # 0-based
        self.end = None  # 0-based
        self.sequence = sequence
        self.set_coords_from_zero()

    def set_coords_from_position(self, pos):
        self.start = pos
        pos += len(self.sequence)
        self.end = pos - 1
        return pos

    def set_coords_from_zero(self):
        return self.set_coords_from_position(0)

    def copy(self):
        copied = Region(self.sequence)
        copied.start = self.start
        copied.end = self.end
        return copied


class Isoform(object):
    def __init__(self):
        self.id = None
        self.name = None
        self.exons = list()

    def copy(self):
        copied = Isoform()
        copied.id = self.id
        copied.name = self.name
        copied.exons = [exon.copy() for exon in self.exons]
        return copied


class Gene(object):
    def __init__(self):
        self.id = None
        self.name = None
        self.strand = '+'
        self.exons = list()
        self.introns = list()
        self.isoforms = list()

    def copy(self):
        copied = Gene()
        copied.id = self.id
        copied.name = self.name
        copied.strand = self.strand
        copied.exons = [exon.copy() for exon in self.exons]
        copied.introns = [intron.copy() for intron in self.introns]
        copied.isoforms = [isoform.copy() for isoform in self.isoforms]
        return copied

    def get_sequence_and_set_coords(self, start):
        parts = list()
        pos = start
        # handle the last exon separately
        for exon_i, exon in enumerate(self.exons[:-1]):
            pos = exon.set_coords_from_position(pos)
            parts.append(exon.sequence)
            intron = self.introns[exon_i]
            pos = intron.set_coords_from_position(pos)
            parts.append(intron.sequence)

        last_exon = self.exons[-1]
        last_exon.set_coords_from_position(pos)
        parts.append(last_exon.sequence)
        return ''.join(parts)

    def add_isoform(self, isoform):
        self.isoforms.append(isoform)

    def add_isoform_for_exon_numbers(self, exon_numbers, isoform_suffix):
        isoform = Isoform()
        isoform.exons = list()
        for exon_i in exon_numbers:
            isoform.exons.append(self.exons[exon_i])

        isoform.id = '{}T{}'.format(self.id, isoform_suffix)
        isoform.name = '{}T{}'.format(self.name, isoform_suffix)
        self.add_isoform(isoform)

    def set_intron_start_end_bases(self,
                                   exon_i=0,
                                   offset_from_start=1,
                                   length=4):
        exon = self.exons[exon_i]
        intron_end_i = offset_from_start + (length - 1)
        before_seq = exon.sequence[:offset_from_start]
        intron_seq = exon.sequence[offset_from_start:intron_end_i + 1]
        after_seq = exon.sequence[intron_end_i + 1:]
        intron_seq_mid = intron_seq[len(INTRON_START):-len(INTRON_END)]
        new_intron_seq = '{}{}{}'.format(INTRON_START, intron_seq_mid,
                                         INTRON_END)
        new_exon_seq = '{}{}{}'.format(before_seq, new_intron_seq, after_seq)
        exon.sequence = new_exon_seq

    def set_intron_start_bases(self, intron_i=0, offset_from_start=2):
        intron = self.introns[intron_i]
        before_seq = intron.sequence[:offset_from_start]
        after_seq = intron.sequence[offset_from_start + len(INTRON_START):]
        new_intron_seq = '{}{}{}'.format(before_seq, INTRON_START, after_seq)
        intron.sequence = new_intron_seq


class Chromosome(object):
    def __init__(self):
        self.name = None
        self.genes = list()
        self.intergenic_regions = list()

    def get_sequence_and_set_coords(self):
        pos = 0
        parts = list()
        # TODO ESPRESSO does not correctly handle isoforms that
        # start within the first 2 coordinates of the chromosome.
        # For now, just add 'NN' at the start of the chromosome.
        parts.append('NN')
        pos += 2
        # handle the last gene separately
        for gene_i, gene in enumerate(self.genes[:-1]):
            gene_sequence = gene.get_sequence_and_set_coords(pos)
            parts.append(gene_sequence)
            pos = gene.exons[-1].end + 1
            intergenic = self.intergenic_regions[gene_i]
            pos = intergenic.set_coords_from_position(pos)
            parts.append(intergenic.sequence)

        last_gene = self.genes[-1]
        last_sequence = last_gene.get_sequence_and_set_coords(pos)
        parts.append(last_sequence)
        return ''.join(parts)

    def copy(self):
        copied = Chromosome()
        copied.name = self.name
        copied.genes = [gene.copy() for gene in self.genes]
        copied.intergenic_regions = [
            region.copy() for region in self.intergenic_regions
        ]
        return copied


class BaseTest(unittest.TestCase):
    def setUp(self):
        random.seed(1)  # get consistent "random" values
        self._base_test_dir = os.path.dirname(__file__)
        self._espresso_dir = os.path.dirname(self._base_test_dir)
        self._espresso_src_dir = os.path.join(self._espresso_dir, 'src')
        self._test_data_dir = os.path.join(self._espresso_dir, 'test_data')
        self._visualization_src_dir = os.path.join(self._espresso_dir,
                                                   'visualization')
        self._snakemake_script_dir = os.path.join(self._espresso_dir,
                                                  'snakemake', 'scripts')
        self._espresso_s = os.path.join(self._espresso_src_dir,
                                        'ESPRESSO_S.pl')
        self._espresso_c = os.path.join(self._espresso_src_dir,
                                        'ESPRESSO_C.pl')
        self._espresso_q = os.path.join(self._espresso_src_dir,
                                        'ESPRESSO_Q.pl')
        self._visualization_py = os.path.join(self._visualization_src_dir,
                                              'visualize.py')
        self._split_s_py = os.path.join(self._snakemake_script_dir,
                                        'split_espresso_s_output_for_c.py')
        self._combine_c_py = os.path.join(
            self._snakemake_script_dir, 'combine_espresso_c_output_for_q.py')
        self._py_executable = _get_python_executable()

    def _assert_within_x_percent_or_y(self, actual, expected, x_percent, y):
        offset = x_percent * expected * 0.01
        offset = max(offset, y)
        lower = expected - offset
        upper = expected + offset
        if lower <= actual <= upper:
            return

        self.fail('{} not within {}% or {} of {} [{}, {}]'.format(
            actual, x_percent, y, expected, lower, upper))

    def _assert_does_not_exist_or_is_empty(self, paths):
        for path in paths:
            if not os.path.exists(path):
                continue

            stat_result = os.stat(path)
            self.assertEqual(stat_result.st_size,
                             0,
                             msg='{} is not empty'.format(path))

    def _assert_exists_and_non_empty(self, paths):
        for path in paths:
            if not os.path.exists(path):
                self.fail('{} does not exist'.format(path))

            stat_result = os.stat(path)
            self.assertGreater(stat_result.st_size,
                               0,
                               msg='{} is empty'.format(path))

    def _assert_does_not_exist(self, path):
        if os.path.exists(path):
            self.fail('{} exists'.format(path))

    def _get_random_range(self, low, high):
        return random.randint(low, high)

    def _get_random_nucleotides(self, count):
        choices = ['A', 'C', 'G', 'T']
        return ''.join(random.choices(choices, k=count))

    def _get_random_gene(self, num_exons=4):
        gene = Gene()
        if num_exons < 1:
            raise Exception('num_exons < 1: {}'.format(num_exons))

        first_exon = self._get_random_exon()
        gene.exons.append(first_exon)
        for _ in range(num_exons - 1):
            intron = self._get_random_intron()
            gene.introns.append(intron)
            exon = self._get_random_exon()
            gene.exons.append(exon)

        return gene

    def _get_random_exon(self,
                         min_length=DEFAULT_EXON_LENGTH[0],
                         max_length=DEFAULT_EXON_LENGTH[1]):
        length = self._get_random_range(min_length, max_length)
        seq = self._get_random_nucleotides(length)
        exon = Region(seq)
        return exon

    def _get_random_intron(self,
                           min_length=DEFAULT_INTRON_LENGTH[0],
                           max_length=DEFAULT_INTRON_LENGTH[1]):
        length = self._get_random_range(min_length, max_length)
        middle_length = length - (len(INTRON_START) + len(INTRON_END))
        middle = self._get_random_nucleotides(middle_length)
        intron = Region(INTRON_START + middle + INTRON_END)
        return intron

    def _get_random_intergenic_region(self,
                                      min_length=DEFAULT_INTERGENIC_LENGTH[0],
                                      max_length=DEFAULT_INTERGENIC_LENGTH[1]):
        length = self._get_random_range(min_length, max_length)
        seq = self._get_random_nucleotides(length)
        intergenic = Region(seq)
        return intergenic

    def _make_chromosome_from_genes(self, genes):
        chrom = Chromosome()
        chrom.genes = genes
        for _ in genes[1:]:
            intergenic = self._get_random_intergenic_region()
            chrom.intergenic_regions.append(intergenic)

        chrom.get_sequence_and_set_coords()
        return chrom

    def _parse_abundance(self, abundance_path, sample_names):
        expected_headers = (['transcript_ID', 'transcript_name', 'gene_ID'] +
                            sample_names)
        rows = dict()
        with open(abundance_path, 'rt') as abundance_handle:
            for i, line in enumerate(abundance_handle):
                columns = line.strip().split('\t')
                if i == 0:
                    self.assertEqual(columns, expected_headers)
                    continue

                self.assertEqual(len(columns), len(expected_headers))
                row = dict()
                transcript_id = columns[0]
                row['transcript_name'] = columns[1]
                row['gene_id'] = columns[2]
                for sample_i, sample_name in enumerate(sample_names):
                    row[sample_name] = float(columns[3 + sample_i])

                rows[transcript_id] = row

        return rows

    def _parse_attributes(self, attributes_str):
        attributes = dict()
        attribute_pairs = attributes_str.split(';')
        for attribute_pair in attribute_pairs:
            attribute_pair = attribute_pair.strip()
            first_space = attribute_pair.find(' ')
            if first_space <= 0:
                continue

            key = attribute_pair[:first_space]
            value = attribute_pair[first_space + 1:]
            if value[0] == '"' and value[-1] == '"':
                # remove quotes
                value = value[1:-1]

            attributes[key] = value

        return attributes

    def _parse_updated_gtf(self, gtf_path):
        isoforms = dict()
        with open(gtf_path, 'rt') as gtf_handle:
            for line in gtf_handle:
                columns = line.strip().split('\t')
                if columns[0].startswith('#'):
                    continue

                seq = columns[0]
                # source = columns[1]
                feature = columns[2]
                start_str = columns[3]
                end_str = columns[4]
                # score = columns[5]
                strand = columns[6]
                # frame = columns[7]
                attributes_str = columns[8]
                attributes = self._parse_attributes(attributes_str)

                if feature != 'exon':
                    continue

                transcript_id = attributes['transcript_id']
                exon_number = int(attributes['exon_number'])
                exon_i = exon_number - 1
                isoform_details = isoforms.get(transcript_id)
                if not isoform_details:
                    isoform_details = {'transcript_id': transcript_id}
                    isoforms[transcript_id] = isoform_details
                    isoform_details['isoform'] = Isoform()
                    isoform_details['isoform'].id = transcript_id
                    isoform_details['chr'] = seq
                    isoform_details['strand'] = strand

                while len(isoform_details['isoform'].exons) < exon_number:
                    isoform_details['isoform'].exons.append(None)

                # 1-based
                start = int(start_str)
                end = int(end_str)
                length = (end - start) + 1
                exon = Region(length * 'N')
                exon.set_coords_from_position(start - 1)
                isoform_details['isoform'].exons[exon_i] = exon

        index = dict()
        for isoform in isoforms.values():
            key = self._get_key_from_isoform_details(isoform)
            index[key] = isoform

        return {'isoforms': isoforms, 'index': index}

    def _parse_sj_id_string(self, id_str):
        parts = id_str.split(':')
        sj_chr = parts[0]
        sj_start = int(parts[1])
        sj_end = int(parts[2])
        sj_strand = self._get_strand_from_int_str(parts[3])

        return {
            'chr': sj_chr,
            'start': sj_start,
            'end': sj_end,
            'strand': sj_strand
        }

    def _parse_sj_simplified(self, sj_simplified_path):
        sjs_from_simplified = dict()
        sj_clusters = dict()
        with open(sj_simplified_path, 'rt') as handle:
            for line in handle:
                columns = line.strip().split('\t')
                if columns[0] == 'SJ_cluster':
                    read_group = int(columns[1])
                    cluster_index = int(columns[2])
                    # cluster_sort_index = int(columns[3])
                    cluster_chr = columns[4]
                    cluster_start = int(columns[5])
                    cluster_end = int(columns[6])
                    cluster_key = (cluster_chr, cluster_start, cluster_end)
                    sj_clusters[cluster_key] = {
                        'chr': cluster_chr,
                        'start': cluster_start,
                        'end': cluster_end,
                        'read_group': read_group,
                        'cluster_index': cluster_index
                    }
                    continue

                read_group = int(columns[0])
                sj_id_details = self._parse_sj_id_string(columns[1])
                # columns[2] is empty
                # SJ start and end are (0-based intron_start, 1-based intron_end)
                sj_start = int(columns[3])
                sj_end = int(columns[4])
                strand_str = columns[5]
                strand = self._get_strand_from_int_str(strand_str)

                perfect_count = int(columns[6])
                all_count = int(columns[7])
                upstream_2nt = columns[8]
                downstream_2nt = columns[9]
                # tag = columns[10]
                is_putative = columns[11] == 'yes'
                is_annotated = columns[12] == 'yes'
                is_high_confidence = columns[13] == '1'
                cluster_index = int(columns[14])
                sj_key = (sj_id_details['chr'], sj_start, sj_end)
                sjs_from_simplified[sj_key] = {
                    'chr': sj_id_details['chr'],
                    'start': sj_start,
                    'end': sj_end,
                    'read_group': read_group,
                    'cluster_index': cluster_index,
                    'strand': strand,
                    'perfect_count': perfect_count,
                    'all_count': all_count,
                    'upstream_nt': upstream_2nt,
                    'downstream_nt': downstream_2nt,
                    'is_annotated': is_annotated,
                    'is_high_confidence': is_high_confidence
                }

        return (sjs_from_simplified, sj_clusters)

    def _parse_read_final(self, read_final_path):
        parsed = dict()
        with open(read_final_path, 'rt') as handle:
            for line in handle:
                columns = line.strip().split('\t')
                read_id = columns[0]
                read_details = parsed.get(read_id)
                if not read_details:
                    read_details = dict()
                    parsed[read_id] = read_details

                feature = columns[1]
                if feature == 'group_ID':
                    read_details['read_group'] = int(columns[2])
                    read_details['read_id_update'] = columns[3]
                elif feature == 'strand_isoform':
                    read_details[
                        'strand_isoform'] = self._get_strand_from_int_str(
                            columns[2])
                elif feature == 'strand_read':
                    read_details[
                        'strand_read'] = self._get_strand_from_int_str(
                            columns[2])
                elif feature == 'chr':
                    read_details['chr'] = columns[2]
                elif feature == 'mapq':
                    read_details['mapq'] = int(columns[2])
                elif feature == 'flag':
                    read_details['flag'] = int(columns[2])
                elif feature == 'SJ':
                    read_sjs = read_details.get('SJs')
                    if not read_sjs:
                        read_sjs = list()
                        read_details['SJs'] = read_sjs

                    junction_i = int(columns[2]) - 1  # 0-based
                    while len(read_sjs) <= junction_i:
                        read_sjs.append(None)

                    sj_id_str = columns[3]
                    sj_details = self._parse_sj_id_string(sj_id_str)
                    read_length_used_up_to_sj = int(columns[4])
                    read_length_remaining_after_sj = int(columns[5])
                    corrected_status = columns[6]
                    corrected_id_str = columns[7]
                    if corrected_id_str == 'NA':
                        # if not corrected than use the uncorrected value
                        corrected_sj_details = self._parse_sj_id_string(
                            sj_id_str)
                    else:
                        corrected_sj_details = self._parse_sj_id_string(
                            corrected_id_str)

                    corrected_length_up_to = parse_float_with_NA(columns[8])
                    corrected_length_remaining = parse_float_with_NA(
                        columns[9])
                    corrected_score = parse_float_with_NA(columns[10])
                    is_annotated = columns[11] == 'yes'
                    corrected_is_annotated = columns[12] == 'yes'
                    read_sjs[junction_i] = {
                        'id_str': sj_id_str,
                        'chr': sj_details['chr'],
                        'start': sj_details['start'],
                        'end': sj_details['end'],
                        'strand': sj_details['strand'],
                        'read_length_used_up_to_sj': read_length_used_up_to_sj,
                        'read_length_remaining_after_sj':
                        read_length_remaining_after_sj,
                        'corrected_id_str': corrected_id_str,
                        'corrected_chr': corrected_sj_details['chr'],
                        'corrected_start': corrected_sj_details['start'],
                        'corrected_end': corrected_sj_details['end'],
                        'corrected_strand': corrected_sj_details['strand'],
                        'corrected_length_up_to': corrected_length_up_to,
                        'corrected_length_remaining':
                        corrected_length_remaining,
                        'corrected_status': corrected_status,
                        'corrected_score': corrected_score,
                        'is_annotated': is_annotated,
                        'corrected_is_annotated': corrected_is_annotated,
                    }
                elif feature in ['start', 'end']:
                    # columns[2] == 'NA'
                    position = int(columns[3]) - 1  # 0-based
                    length_up_to = int(columns[4])
                    length_remaining = int(columns[5])
                    corrected_status = columns[6]
                    corrected_position = int(columns[7])
                    corrected_length_up_to = int(columns[8])
                    corrected_length_remaining = int(columns[9])
                    # columns[10,11,12] == 'NA'
                    read_details[feature] = {
                        'position': position,
                        'length_up_to': length_up_to,
                        'length_remaining': length_remaining,
                        'corrected_status': corrected_status,
                        'corrected_position': corrected_position,
                        'corrected_length_up_to': corrected_length_up_to,
                        'corrected_length_remaining':
                        corrected_length_remaining,
                    }
                elif feature == 'mapped_length_read':
                    read_details['mapped_length_read'] = int(columns[2])
                else:
                    raise Exception('unknown feature: {} in {}'.format(
                        feature, read_final_path))

        return parsed

    def _get_strand_from_int_str(self, int_str):
        if int_str == '0':
            return '+'
        if int_str == '1':
            return '-'

        return int_str

    def _get_key_from_isoform_details(self, isoform):
        key = list()
        key.append(isoform['chr'])
        key.append(isoform['strand'])
        for exon in isoform['isoform'].exons:
            key.append(exon.start)
            key.append(exon.end)

        return tuple(key)

    def _get_key_from_chr_strand_and_exons(self, chrom, strand, exons):
        key = list()
        key.append(chrom)
        key.append(strand)
        for exon in exons:
            key.append(exon.start)
            key.append(exon.end)

        return tuple(key)

    def _get_transcript_ids_from_detected(self, expected_isoforms,
                                          detected_isoforms):
        expected_by_id = dict()
        for expected_isoform in expected_isoforms:
            expected_key = self._get_key_from_chr_strand_and_exons(
                expected_isoform['chr'], expected_isoform['gene'].strand,
                expected_isoform['exons'])
            found = detected_isoforms['index'].get(expected_key)
            if not found:
                self.fail('could not find: {}'.format(expected_key))

            existing_id = expected_isoform.get('transcript_id')
            if existing_id:
                self.assertEqual(existing_id, found['transcript_id'])
            else:
                expected_isoform['transcript_id'] = found['transcript_id']

            expected_by_id[found['transcript_id']] = expected_isoform

        return expected_by_id

    def _check_abundance_row_values(self, expected_by_id, abundance_rows,
                                    sample_names):
        for transcript_id, abundance_details in abundance_rows.items():
            expected = expected_by_id[transcript_id]
            if expected.get('transcript_name'):
                self.assertEqual(abundance_details['transcript_name'],
                                 expected['transcript_name'])

            if not expected['is_novel']:
                self.assertEqual(abundance_details['gene_id'],
                                 expected['gene'].id)

            for name in sample_names:
                self.assertAlmostEqual(abundance_details[name],
                                       expected['abundance'][name],
                                       msg='transcript: {}, sample: {}'.format(
                                           transcript_id, name))
