import os
import os.path
import unittest

import tests.base_test


class AlignmentsBaseTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_dir = os.path.dirname(__file__)
        # set in subclasses
        self._test_name = None
        self._data_dir = None
        self._out_dir = None
        self._log_dir = None
        self._work_dir = None

        self._out_file_prefix = (
            tests.base_test.out_file_prefix_based_on_params())
        self._chromosomes = None
        self._alignment_chrom = None
        self._fasta = None
        self._gtf = None
        self._sam = None
        self._sample_names = None
        self._samples_tsv = None
        self._samples_updated = None
        self._abundance = None
        self._updated_gtf = None
        self._s_log = None

    def _run_test(self):
        self._initialize_dirs()
        self._create_test_data()
        self._create_samples_tsv()
        self._run_espresso_s()
        self._check_s_output()
        if self._should_expect_s_to_fail():
            return

        self._run_espresso_c()
        self._check_c_output()
        self._run_espresso_q()
        self._check_q_output()

    def _initialize_dirs(self):
        directories = [self._data_dir, self._out_dir, self._log_dir]
        tests.base_test.remove_then_create_directories(directories)

    def _create_test_data(self):
        self._fasta = os.path.join(self._data_dir, 'test.fasta')
        self._create_test_fasta(self._fasta)
        self._gtf = os.path.join(self._data_dir, 'test.gtf')
        self._create_test_gtf(self._gtf)
        self._sam = os.path.join(self._data_dir, 'test.sam')
        self._create_test_sam(self._sam)

    def _create_test_fasta(self, fasta_path):
        self._chromosomes = list()
        chr_id = 'chr1'
        gene_id = 'ENSG1'
        gene_name = 'GENE1'
        # 1 gene with 4 exons
        gene = self._get_random_gene(num_exons=4)
        gene.id = gene_id
        gene.name = gene_name
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id
        self._chromosomes.append(chromosome)

        if self._rename_chr_in_sam():
            # Use 1 instead of chr1 as name in .sam file
            self._alignment_chrom = chromosome.copy()
            self._alignment_chrom.name = '1'
        else:
            self._alignment_chrom = chromosome

        tests.base_test.write_fasta_from_chroms(fasta_path, self._chromosomes)

    def _create_test_gtf(self, gtf_path):
        gene = self._chromosomes[0].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2], '0')
        gene.add_isoform_for_exon_numbers([0, 2], '1')

        tests.base_test.write_gtf_from_chroms(gtf_path, self._chromosomes)

    def _get_alignments(self):
        alignments = list()
        chrom = self._alignment_chrom
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.append_copies(alignments, alignment, 2)
        return alignments

    def _get_expected_isoforms(self):
        expected_isoforms = list()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[0].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[0].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        return expected_isoforms

    def _create_test_sam(self, sam_path):
        alignments = self._get_alignments()

        use_m_cigar_op = self._use_m_cigar_op()
        read_id_offset = 0
        read_id_offset = tests.base_test.write_sam_from_alignments(
            sam_path, [self._alignment_chrom],
            alignments,
            read_id_offset,
            use_m_cigar_op=use_m_cigar_op)
        temp_file = '{}.tmp.sam'.format(sam_path)
        sort_log = os.path.join(self._log_dir, 'sort.log')
        tests.base_test.sort_sam_with_temp_file(sam_path, temp_file, sort_log)

    def _create_samples_tsv(self):
        self._samples_tsv = os.path.join(self._data_dir, 'samples.tsv')
        self._sample_names = [self._test_name]
        with open(self._samples_tsv, 'wt') as handle:
            columns = [self._sam, self._test_name]
            tests.base_test.write_tsv_line(handle, columns)

    def _should_expect_s_to_fail(self):
        return False

    def _use_m_cigar_op(self):
        return True

    def _rename_chr_in_sam(self):
        return False

    def _run_espresso_s(self):
        command = [
            'perl', self._espresso_s, '-L', self._samples_tsv, '-F',
            self._fasta, '-O', self._work_dir, '-A', self._gtf
        ]

        self._s_log = os.path.join(self._log_dir, 'espresso_s.log')
        check = not self._should_expect_s_to_fail()
        process = tests.base_test.run_command_with_log(command,
                                                       self._s_log,
                                                       check=check)
        if self._should_expect_s_to_fail():
            self.assertNotEqual(process.returncode, 0)

    def _check_s_output(self):
        self._samples_updated = os.path.join(self._work_dir,
                                             'samples.tsv.updated')
        sj_simplified_list = os.path.join(
            self._work_dir,
            '{}_SJ_simplified.list'.format(self._chromosomes[0].name))
        all_sjs = os.path.join(self._work_dir, 'SJ_group_all.fa')
        sam_list = os.path.join(self._work_dir, '0', 'sam.list3')
        sj_list = os.path.join(self._work_dir, '0', 'sj.list')

        if self._should_expect_s_to_fail():
            expected_error = 'chr name (1) not recognized'
            with open(self._s_log, 'rt') as log_handle:
                for line in log_handle:
                    if expected_error in line:
                        break
                else:
                    self.fail('did not find {} in {}'.format(
                        expected_error, self._s_log))
        else:
            self._assert_exists_and_non_empty([
                self._samples_updated, sj_simplified_list, all_sjs, sam_list,
                sj_list
            ])

    def _run_espresso_c(self):
        threads = '2'
        sample_i = '0'  # only 1 sample
        command = [
            'perl', self._espresso_c, '-I', self._work_dir, '-F', self._fasta,
            '-X', sample_i, '-T', threads
        ]
        c_log = os.path.join(self._log_dir, 'espresso_c.log')
        tests.base_test.run_command_with_log(command, c_log)

    def _check_c_output(self):
        sample_i = '0'
        chrom = self._chromosomes[0]
        read_final = os.path.join(self._work_dir, sample_i,
                                  '{}_read_final.txt'.format(chrom.name))
        self._assert_exists_and_non_empty([read_final])
        parsed_read_final = self._parse_read_final(read_final)
        self.assertGreater(len(parsed_read_final), 0)

    def _run_espresso_q(self):
        command = [
            'perl', self._espresso_q, '-L', self._samples_updated, '-A',
            self._gtf
        ]

        q_log = os.path.join(self._log_dir, 'espresso_q.log')
        tests.base_test.run_command_with_log(command, q_log)

    def _check_q_output(self):
        gtf_basename = '{}_updated.gtf'.format(self._out_file_prefix)
        self._updated_gtf = os.path.join(self._work_dir, gtf_basename)
        abundance_basename = '{}_abundance.esp'.format(self._out_file_prefix)
        self._abundance = os.path.join(self._work_dir, abundance_basename)
        self._assert_exists_and_non_empty([self._updated_gtf, self._abundance])
        detected_isoforms = self._parse_updated_gtf(self._updated_gtf)
        abundance_rows = self._parse_abundance(self._abundance,
                                               self._sample_names)

        expected_isoforms = self._get_expected_isoforms()
        expected_by_id = self._get_transcript_ids_from_detected(
            expected_isoforms, detected_isoforms)
        self.assertEqual(len(abundance_rows), len(expected_isoforms))
        self._check_abundance_row_values(expected_by_id, abundance_rows,
                                         self._sample_names)


class ChrNameMismatchTest(AlignmentsBaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'chr_name_mismatch'
        self._data_dir = os.path.join(self._test_dir, 'chr_name_mismatch_data')
        self._out_dir = os.path.join(self._test_dir, 'chr_name_mismatch_out')
        self._log_dir = os.path.join(self._test_dir, 'chr_name_mismatch_logs')
        self._work_dir = os.path.join(self._out_dir, 'work_dir')

    def test(self):
        self._run_test()

    def _should_expect_s_to_fail(self):
        return True

    def _rename_chr_in_sam(self):
        return True


class CigarFormatTest(AlignmentsBaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'cigar_format'
        self._data_dir = os.path.join(self._test_dir, 'cigar_format_data')
        self._out_dir = os.path.join(self._test_dir, 'cigar_format_out')
        self._log_dir = os.path.join(self._test_dir, 'cigar_format_logs')
        self._work_dir = os.path.join(self._out_dir, 'work_dir')

    def test(self):
        self._run_test()

    def _use_m_cigar_op(self):
        return False

    def _get_alignments(self):
        alignments = super()._get_alignments()
        chrom = self._alignment_chrom
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2])
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=0,
                                             offset_from_start=10,
                                             mismatch_len=3)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2])
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=0,
                                             mismatch_len=5)
        tests.base_test.append_copies(alignments, alignment, 2)

        for alignment in alignments:
            alignment.cigar.convert_alignment_match_to_sequence_match()

        return alignments

    def _get_expected_isoforms(self):
        expected_isoforms = super()._get_expected_isoforms()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['abundance'] = {self._test_name: 4}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[0].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[0].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)
        return expected_isoforms


class SecondaryAlignmentTest(AlignmentsBaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'secondary_alignment'
        self._data_dir = os.path.join(self._test_dir,
                                      'secondary_alignment_data')
        self._out_dir = os.path.join(self._test_dir, 'secondary_alignment_out')
        self._log_dir = os.path.join(self._test_dir,
                                     'secondary_alignment_logs')
        self._work_dir = os.path.join(self._out_dir, 'work_dir')

    def test(self):
        self._run_test()

    def _get_alignments(self):
        alignments = super()._get_alignments()
        chrom = self._alignment_chrom
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2])
        alignment.read_id = 'first_read_id'
        tests.base_test.append_copies(alignments, alignment, 1)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        alignment.read_id = 'first_read_id'
        alignment.is_secondary = True
        tests.base_test.append_copies(alignments, alignment, 1)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        alignment.read_id = 'second_read_id'
        alignment.is_secondary = True
        tests.base_test.append_copies(alignments, alignment, 1)

        return alignments

    def _get_expected_isoforms(self):
        expected_isoforms = super()._get_expected_isoforms()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['abundance'] = {self._test_name: 1}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[0].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[0].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)
        return expected_isoforms


class MissingSequenceTest(AlignmentsBaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'missing_sequence'
        self._data_dir = os.path.join(self._test_dir,
                                      '{}_data'.format(self._test_name))
        self._out_dir = os.path.join(self._test_dir,
                                     '{}_out'.format(self._test_name))
        self._log_dir = os.path.join(self._test_dir,
                                     '{}_logs'.format(self._test_name))
        self._work_dir = os.path.join(self._out_dir, 'work_dir')

    def test(self):
        self._run_test()

    def _get_alignments(self):
        alignments = super()._get_alignments()
        chrom = self._alignment_chrom
        chrom_seq = chrom.get_sequence_and_set_coords()
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2])
        tests.base_test.append_copies(alignments, alignment, 2)
        # Trim off 1 base from sequence and quality
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        alignment.set_default_quality_from_sequence()
        alignment.sequence = alignment.sequence[:-1]
        alignment.quality = alignment.quality[:-1]
        # Replace last matched base with a hard clip
        alignment.cigar.operations[-1].num -= 1
        alignment.cigar.add_hard_clip(1)
        tests.base_test.append_copies(alignments, alignment, 3)
        # No sequence
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2])
        alignment.sequence = '*'
        alignment.quality = '*'
        tests.base_test.append_copies(alignments, alignment, 4)

        # like above but for [0, 2, 3]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        alignment.set_default_quality_from_sequence()
        alignment.sequence = alignment.sequence[:-1]
        alignment.quality = alignment.quality[:-1]
        alignment.cigar.operations[-1].num -= 1
        alignment.cigar.add_hard_clip(1)
        tests.base_test.append_copies(alignments, alignment, 3)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2, 3])
        alignment.sequence = '*'
        alignment.quality = '*'
        tests.base_test.append_copies(alignments, alignment, 4)

        return alignments

    def _get_expected_isoforms(self):
        expected_isoforms = super()._get_expected_isoforms()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[0].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[0].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)

        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[0].name, 0,
                                                  0))
        expected_isoforms.append(expected_isoform)

        return expected_isoforms


if __name__ == '__main__':
    unittest.main(verbosity=2)
