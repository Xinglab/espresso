import os
import os.path
import unittest

import tests.base_test


class ReadFiltersTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'read_filters'
        self._test_dir = os.path.dirname(__file__)
        self._data_dir = os.path.join(self._test_dir, 'data')
        self._out_dir = os.path.join(self._test_dir, 'out')
        self._log_dir = os.path.join(self._test_dir, 'logs')
        self._work_dir = os.path.join(self._out_dir, 'work_dir')
        self._out_file_prefix = (
            tests.base_test.out_file_prefix_based_on_params())
        self._chr_name = 'chr1'
        self._chromosome = None
        self._fasta = None
        self._gtf = None
        self._sam = None
        self._sample_names = None
        self._samples_tsv = None
        self._samples_updated = None
        self._abundance = None
        self._updated_gtf = None

    def test(self):
        self._initialize_dirs()
        self._create_test_data()
        configs = list()
        configs.append({
            'name': 'default',
            'chrM': None,
            'mapq': None,
        })
        configs.append({
            'name': 'no_chr1',
            'chrM': self._chr_name,
            'mapq': None,
        })
        configs.append({
            'name': 'mapq_3',
            'chrM': None,
            'mapq': 3,
        })
        for config in configs:
            print('(config={})'.format(config['name']))
            # keep the generated files in data_dir
            self._clean_out_and_log_dirs()
            self._create_samples_tsv()
            self._run_espresso_s(config['chrM'], config['mapq'])
            self._check_s_output(config['name'])
            if config['name'] == 'no_chr1':
                continue

            self._run_espresso_c()
            self._check_c_output(config['name'])
            self._run_espresso_q()
            self._check_q_output(config['name'])

    def _initialize_dirs(self):
        directories = [self._data_dir, self._out_dir, self._log_dir]
        tests.base_test.remove_then_create_directories(directories)

    def _clean_out_and_log_dirs(self):
        directories = [self._out_dir, self._log_dir]
        tests.base_test.remove_then_create_directories(directories)

    def _create_test_data(self):
        self._fasta = os.path.join(self._data_dir, 'test.fasta')
        self._create_test_fasta(self._fasta)
        self._gtf = os.path.join(self._data_dir, 'test.gtf')
        self._create_test_gtf(self._gtf)
        self._sam = os.path.join(self._data_dir, 'test.sam')
        self._create_test_sam(self._sam)

    def _create_test_fasta(self, fasta_path):
        genes = list()
        # 1 gene with 3 exons and an SE event
        gene = self._get_random_gene(num_exons=3)
        gene.id = 'ENSG01'
        gene.name = 'GENE01'
        genes.append(gene)

        self._chromosome = self._make_chromosome_from_genes(genes)
        self._chromosome.name = self._chr_name
        chroms = [self._chromosome]
        tests.base_test.write_fasta_from_chroms(fasta_path, chroms)

    def _create_test_gtf(self, gtf_path):
        gene = self._chromosome.genes[0]
        # All exons
        gene.add_isoform_for_exon_numbers([0, 1, 2], '0')

        # Skip middle exon
        gene.add_isoform_for_exon_numbers([0, 2], '1')
        tests.base_test.write_gtf_from_chroms(gtf_path, [self._chromosome])

    def _create_test_sam(self, sam_path):
        alignments = list()
        chroms = [self._chromosome]
        read_id_offset = 0
        gene = self._chromosome.genes[0]
        # 1 read for each isoform at each mapping quality level
        for mapping_quality in range(1, 10):
            alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
                self._chr_name, gene, [0, 1, 2])
            alignment.mapq = mapping_quality
            alignments.append(alignment)
            alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
                self._chr_name, gene, [0, 2])
            alignment.mapq = mapping_quality
            alignments.append(alignment)

        read_id_offset += tests.base_test.write_sam_from_alignments(
            sam_path, chroms, alignments, read_id_offset)
        temp_file = '{}.tmp.sam'.format(sam_path)
        sort_log = os.path.join(self._log_dir, 'sort_1.log')
        tests.base_test.sort_sam_with_temp_file(sam_path, temp_file, sort_log)

    def _create_samples_tsv(self):
        self._samples_tsv = os.path.join(self._data_dir, 'samples.tsv')
        self._sample_names = [self._test_name]
        with open(self._samples_tsv, 'wt') as handle:
            columns = [self._sam, self._test_name]
            tests.base_test.write_tsv_line(handle, columns)

    def _run_espresso_s(self, chrm, mapq):
        command = [
            'perl', self._espresso_s, '-L', self._samples_tsv, '-F',
            self._fasta, '-O', self._work_dir, '-A', self._gtf
        ]

        if chrm:
            command.extend(['--chrM', chrm])

        if mapq:
            command.extend(['--mapq_cutoff', str(mapq)])

        s_log = os.path.join(self._log_dir, 'espresso_s.log')
        tests.base_test.run_command_with_log(command, s_log)

    def _check_s_output(self, config_name):
        self._samples_updated = os.path.join(self._work_dir,
                                             'samples.tsv.updated')
        sj_simplified_list = os.path.join(
            self._work_dir, '{}_SJ_simplified.list'.format(self._chr_name))
        all_sjs = os.path.join(self._work_dir, 'SJ_group_all.fa')
        sam_list = os.path.join(self._work_dir, '0', 'sam.list3')
        sj_list = os.path.join(self._work_dir, '0', 'sj.list')

        self._assert_exists_and_non_empty([self._samples_updated])

        if config_name == 'no_chr1':
            self._assert_does_not_exist_or_is_empty(
                [sj_simplified_list, all_sjs, sam_list, sj_list])
        else:
            self._assert_exists_and_non_empty(
                [sj_simplified_list, all_sjs, sam_list, sj_list])

    def _run_espresso_c(self):
        threads = '2'
        sample_i = '0'  # only 1 sample
        command = [
            'perl', self._espresso_c, '-I', self._work_dir, '-F', self._fasta,
            '-X', sample_i, '-T', threads
        ]
        c_log = os.path.join(self._log_dir, 'espresso_c.log')
        tests.base_test.run_command_with_log(command, c_log)

    def _check_c_output(self, config_name):
        sample_i = '0'
        read_final = os.path.join(self._work_dir, sample_i,
                                  '{}_read_final.txt'.format(self._chr_name))
        self._assert_exists_and_non_empty([read_final])
        parsed_read_final = self._parse_read_final(read_final)
        if config_name == 'default':
            self.assertEqual(len(parsed_read_final), 18)
        if config_name == 'no_chr1':
            self.assertEqual(len(parsed_read_final), 0)
        if config_name == 'mapq_3':
            self.assertEqual(len(parsed_read_final), 14)

    def _run_espresso_q(self):
        command = [
            'perl', self._espresso_q, '-L', self._samples_updated, '-A',
            self._gtf
        ]

        q_log = os.path.join(self._log_dir, 'espresso_q.log')
        tests.base_test.run_command_with_log(command, q_log)

    def _check_q_output(self, config_name):
        gtf_basename = '{}_updated.gtf'.format(self._out_file_prefix)
        self._updated_gtf = os.path.join(self._work_dir, gtf_basename)
        abundance_basename = '{}_abundance.esp'.format(self._out_file_prefix)
        self._abundance = os.path.join(self._work_dir, abundance_basename)
        self._assert_exists_and_non_empty([self._updated_gtf, self._abundance])
        abundance_rows = self._parse_abundance(self._abundance,
                                               self._sample_names)
        detected_isoforms = self._parse_updated_gtf(self._updated_gtf)

        expected_isoforms = list()
        # All exons
        expected_isoform = dict()
        expected_isoform['chr'] = self._chr_name
        expected_isoform['gene'] = self._chromosome.genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        if config_name == 'default':
            expected_isoform['abundance'] = {self._test_name: 9}
        if config_name == 'no_chr1':
            expected_isoform['abundance'] = {self._test_name: 0}
        if config_name == 'mapq_3':
            expected_isoform['abundance'] = {self._test_name: 7}

        expected_isoform['is_novel'] = False
        expected_isoform['transcript_name'] = 'GENE01T0'
        expected_isoforms.append(expected_isoform)

        # Skip middle
        expected_isoform = dict()
        expected_isoform['chr'] = self._chr_name
        expected_isoform['gene'] = self._chromosome.genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        if config_name == 'default':
            expected_isoform['abundance'] = {self._test_name: 9}
        if config_name == 'no_chr1':
            expected_isoform['abundance'] = {self._test_name: 0}
        if config_name == 'mapq_3':
            expected_isoform['abundance'] = {self._test_name: 7}

        expected_isoform['is_novel'] = False
        expected_isoform['transcript_name'] = 'GENE01T1'
        expected_isoforms.append(expected_isoform)

        expected_by_id = self._get_transcript_ids_from_detected(
            expected_isoforms, detected_isoforms)
        self.assertEqual(len(abundance_rows), len(expected_isoforms))
        self._check_abundance_row_values(expected_by_id, abundance_rows,
                                         self._sample_names)


if __name__ == '__main__':
    unittest.main(verbosity=2)
