import os
import os.path
import random
import unittest

import tests.base_test


class SortedInputBaseTest(tests.base_test.BaseTest):
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
        chr_id_template = 'chr{}'
        gene_id_template = 'ENSG{}'
        gene_name_template = 'GENE{}'
        # 2 chrs each with 2 genes each with 4 exons
        gene_template_i = 0
        for chr_i in [1, 2]:
            genes = list()
            for gene_i in [1, 2]:
                gene_template_i += 1
                gene = self._get_random_gene(num_exons=4)
                gene.id = gene_id_template.format(gene_template_i)
                gene.name = gene_name_template.format(gene_template_i)
                genes.append(gene)

            chromosome = self._make_chromosome_from_genes(genes)
            chromosome.name = chr_id_template.format(chr_i)
            self._chromosomes.append(chromosome)

        tests.base_test.write_fasta_from_chroms(fasta_path, self._chromosomes)

    def _create_test_gtf(self, gtf_path):
        gene = self._chromosomes[0].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3], '0')
        gene = self._chromosomes[0].genes[1]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3], '0')

        gene = self._chromosomes[1].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3], '0')
        gene = self._chromosomes[1].genes[1]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3], '0')

        tests.base_test.write_gtf_from_chroms(gtf_path, self._chromosomes)

    def _get_alignments(self):
        alignments = list()
        for chrom in self._chromosomes:
            for gene in chrom.genes:
                alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
                    chrom.name, gene, [0, 1, 2, 3])
                tests.base_test.append_copies(alignments, alignment, 2)
                alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
                    chrom.name, gene, [0, 1, 2])
                tests.base_test.append_copies(alignments, alignment, 2)
                alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
                    chrom.name, gene, [1, 2, 3])
                tests.base_test.append_copies(alignments, alignment, 2)
                alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
                    chrom.name, gene, [0, 2, 3])
                tests.base_test.append_copies(alignments, alignment, 2)

        return alignments

    def _get_expected_isoforms(self):
        expected_isoforms = list()
        novel_i = 0
        for chrom in self._chromosomes:
            for gene in chrom.genes:
                expected_isoform = dict()
                expected_isoform['chr'] = chrom.name
                expected_isoform['gene'] = gene
                expected_isoform['exons'] = [gene.exons[0]]
                expected_isoform['exons'].append(gene.exons[1])
                expected_isoform['exons'].append(gene.exons[2])
                expected_isoform['exons'].append(gene.exons[3])
                expected_isoform['abundance'] = {self._test_name: 6}
                expected_isoform['is_novel'] = False
                expected_isoform['transcript_id'] = gene.isoforms[0].id
                expected_isoform['transcript_name'] = gene.isoforms[0].name
                expected_isoforms.append(expected_isoform)

                expected_isoform = dict()
                expected_isoform['chr'] = chrom.name
                expected_isoform['gene'] = gene
                expected_isoform['exons'] = [gene.exons[0]]
                expected_isoform['exons'].append(gene.exons[2])
                expected_isoform['exons'].append(gene.exons[3])
                expected_isoform['abundance'] = {self._test_name: 2}
                expected_isoform['is_novel'] = True
                expected_isoform['transcript_id'] = (
                    tests.base_test.get_espresso_novel_id(
                        chrom.name, novel_i, 0))
                novel_i += 1
                expected_isoforms.append(expected_isoform)

        return expected_isoforms

    def _create_test_sam(self, sam_path):
        alignments = self._get_alignments()
        random.shuffle(alignments)
        read_id_offset = 0
        read_id_offset += tests.base_test.write_sam_from_alignments(
            sam_path, self._chromosomes, alignments, read_id_offset)
        if self._should_sort_alignments():
            temp_file = '{}.tmp.sam'.format(sam_path)
            sort_log = os.path.join(self._log_dir, 'sort.log')
            tests.base_test.sort_sam_with_temp_file(sam_path, temp_file,
                                                    sort_log)

    def _create_samples_tsv(self):
        self._samples_tsv = os.path.join(self._data_dir, 'samples.tsv')
        self._sample_names = [self._test_name]
        with open(self._samples_tsv, 'wt') as handle:
            columns = [self._sam, self._test_name]
            tests.base_test.write_tsv_line(handle, columns)

    def _should_expect_s_to_fail(self):
        return False

    def _should_sort_alignments(self):
        return True

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
        sj_simplified_lists = list()
        for chrom in self._chromosomes:
            sj_simplified_lists.append(
                os.path.join(self._work_dir,
                             '{}_SJ_simplified.list'.format(chrom.name)))

        all_sjs = os.path.join(self._work_dir, 'SJ_group_all.fa')
        sam_list = os.path.join(self._work_dir, '0', 'sam.list3')
        sj_list = os.path.join(self._work_dir, '0', 'sj.list')

        if self._should_expect_s_to_fail():
            expected_error = '{} is not sorted'.format(self._sam)
            with open(self._s_log, 'rt') as log_handle:
                for line in log_handle:
                    if expected_error in line:
                        break
                else:
                    self.fail('did not find {} in {}'.format(
                        expected_error, self._s_log))
        else:
            self._assert_exists_and_non_empty(
                [self._samples_updated, all_sjs, sam_list, sj_list] +
                sj_simplified_lists)

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
        for chrom in self._chromosomes:
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


class SortedInputTest(SortedInputBaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'sorted'
        self._data_dir = os.path.join(self._test_dir, 'sorted_data')
        self._out_dir = os.path.join(self._test_dir, 'sorted_out')
        self._log_dir = os.path.join(self._test_dir, 'sorted_logs')
        self._work_dir = os.path.join(self._out_dir, 'work_dir')

    def test(self):
        self._run_test()


class UnsortedInputTest(SortedInputBaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'unsorted'
        self._data_dir = os.path.join(self._test_dir, 'unsorted_data')
        self._out_dir = os.path.join(self._test_dir, 'unsorted_out')
        self._log_dir = os.path.join(self._test_dir, 'unsorted_logs')
        self._work_dir = os.path.join(self._out_dir, 'work_dir')

    def test(self):
        self._run_test()

    def _should_expect_s_to_fail(self):
        return True

    def _should_sort_alignments(self):
        return False


if __name__ == '__main__':
    unittest.main(verbosity=2)
