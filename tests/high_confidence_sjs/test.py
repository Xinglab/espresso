import os
import os.path
import unittest

import tests.base_test


class HighConfidenceSjsTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'high_confidence_sjs'
        self._test_dir = os.path.dirname(__file__)
        self._data_dir = os.path.join(self._test_dir, 'data')
        self._out_dir = os.path.join(self._test_dir, 'out')
        self._log_dir = os.path.join(self._test_dir, 'logs')
        self._work_dir = os.path.join(self._out_dir, 'work_dir')
        self._out_file_prefix = None
        self._chr_name = 'chr1'
        self._chromosome = None
        self._fasta = None
        self._gtf = None
        self._bed = None
        self._sams = None
        self._sample_names = None
        self._samples_tsv = None
        self._samples_updated = None
        self._isoform_tsv = None
        self._updated_gtf = None
        self._abundance = None

    def test(self):
        self._initialize_dirs()
        self._create_test_data()
        configs = list()
        # Do 1 restrictive run.
        # Then do additional runs that make certain parameters more permissive
        configs.append({
            'name': 'strict',
            'gtf': None,
            'bed': None,
            'num_cutoff': '4',
            'ratio_cutoff': '0.6'
        })
        configs.append({
            'name': 'num',
            'gtf': None,
            'bed': None,
            'num_cutoff': '2',
            'ratio_cutoff': '0.6'
        })
        configs.append({
            'name': 'ratio',
            'gtf': None,
            'bed': None,
            'num_cutoff': '4',
            'ratio_cutoff': '0.1'
        })
        configs.append({
            'name': 'num_and_ratio',
            'gtf': None,
            'bed': None,
            'num_cutoff': '2',
            'ratio_cutoff': '0.1'
        })
        configs.append({
            'name': 'gtf',
            'gtf': self._gtf,
            'bed': None,
            'num_cutoff': '4',
            'ratio_cutoff': '0.6'
        })
        configs.append({
            'name': 'bed',
            'gtf': None,
            'bed': self._bed,
            'num_cutoff': '4',
            'ratio_cutoff': '0.6'
        })
        configs.append({
            'name': 'bed_and_num_and_ratio',
            'gtf': None,
            'bed': self._bed,
            'num_cutoff': '1',
            'ratio_cutoff': '0.1'
        })
        for config in configs:
            print('(config={})'.format(config['name']))
            # keep the generated files in data_dir
            self._clean_out_and_log_dirs()
            gtf = config['gtf']
            bed = config['bed']
            num_cutoff = config['num_cutoff']
            ratio_cutoff = config['ratio_cutoff']
            self._out_file_prefix = (
                tests.base_test.out_file_prefix_based_on_params(
                    num_cutoff=num_cutoff, ratio_cutoff=ratio_cutoff))
            self._create_samples_tsv()
            self._run_espresso_s(gtf, bed, num_cutoff, ratio_cutoff)
            self._check_s_output(config['name'])
            self._run_espresso_c()
            self._check_c_output(config['name'])
            self._run_espresso_q(gtf, num_cutoff, ratio_cutoff)
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
        self._bed = os.path.join(self._data_dir, 'test.bed')
        self._create_test_bed(self._bed)
        self._sams = list()
        self._sams.append(os.path.join(self._data_dir, 'test_1.sam'))
        self._sams.append(os.path.join(self._data_dir, 'test_2.sam'))
        self._create_test_sams(self._sams)

    def _create_test_fasta(self, fasta_path):
        genes = list()
        # 3 genes each with 4 exons.
        # Each gene is testing a different param:
        #   GENE00. num_cutoff and ratio_cutoff
        #   GENE01. gtf
        #   GENE02. bed
        # Each gene has 2 isoforms: {[0,1,2,3], [0,1,3]}
        # The skipping junction will require the relevant parameter(s) to be set
        for gene_i in range(3):
            gene = self._get_random_gene(num_exons=4)
            gene.id = 'ENSG0{}'.format(gene_i)
            gene.name = 'GENE0{}'.format(gene_i)
            genes.append(gene)

        self._chromosome = self._make_chromosome_from_genes(genes)
        self._chromosome.name = self._chr_name
        chroms = [self._chromosome]
        tests.base_test.write_fasta_from_chroms(fasta_path, chroms)

    def _create_test_gtf(self, gtf_path):
        gene = self._chromosome.genes[1]
        gene.add_isoform_for_exon_numbers([0, 1, 3], '0')
        tests.base_test.write_gtf_from_chroms(gtf_path, [self._chromosome])

    def _create_test_bed(self, bed_path):
        isoform = tests.base_test.Isoform()
        gene = self._chromosome.genes[2]
        isoform.exons = [
            gene.exons[0],
            gene.exons[1],
            gene.exons[3],
        ]
        isoform.id = '{}T0'.format(gene.id)
        isoform.name = '{}T0'.format(gene.name)
        gene.add_isoform(isoform)
        with open(bed_path, 'wt') as handle:
            tests.base_test.write_bed_from_gene(self._chromosome.name, gene,
                                                handle)

    def _create_test_sams(self, sam_paths):
        first_sam = sam_paths[0]
        # For each gene do 5 reads of isoform [0,1,2,3] with no deletions or mismatch.
        # Then add reads for isoform [0,1,3] for each gene so that these reads are
        # treated differently based on the parameter that is being tested with that gene.
        alignments = list()
        chroms = [self._chromosome]
        read_id_offset = 0
        # 5 perfect reads of isoform [0,1,2,3] for each gene
        for gene in self._chromosome.genes:
            alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
                self._chr_name, gene, [0, 1, 2, 3])
            tests.base_test.append_copies(alignments, alignment, 5)

        # 2 perfect reads of [0,1,3] for gene_0 to recognize the junction at --read_num_cutoff=2
        gene_0 = self._chromosome.genes[0]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            self._chr_name, gene_0, [0, 1, 3])
        tests.base_test.append_copies(alignments, alignment, 2)

        # 3 reads with correct SJ chain for gene_0, but not perfect.
        # There are (2 + 3 = 5) total reads for the SJ chain.
        # The ratio of perfect reads in (2/5 = 0.4).
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            self._chr_name, gene_0, [0, 1, 3])
        tests.base_test.add_insertion_at_junction_number(alignment,
                                                         junction_i=1,
                                                         insert_len=3)
        tests.base_test.append_copies(alignments, alignment, 3)

        # 2 reads with the SJ moved slightly for gene_0.
        # When the SJ becomes high-confidence the read should get corrected.
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            self._chr_name, gene_0, [0, 1, 3])
        tests.base_test.add_mismatch_into_intron_number(alignment,
                                                        junction_i=1,
                                                        match_len=3)
        tests.base_test.append_copies(alignments, alignment, 2)

        # Genes 1 and 2 have imperfect reads for the skipping isoform.
        # Gene 1 can rely on the annotated isoform from the gtf
        for gene_i in [1, 2]:
            gene = self._chromosome.genes[gene_i]
            # 2 reads with insertions
            alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
                self._chr_name, gene, [0, 1, 3])
            tests.base_test.add_insertion_at_junction_number(alignment,
                                                             junction_i=1,
                                                             insert_len=3)
            tests.base_test.append_copies(alignments, alignment, 2)

            # 2 reads with SJ moved
            alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
                self._chr_name, gene, [0, 1, 3])
            tests.base_test.add_mismatch_into_intron_number(alignment,
                                                            junction_i=1,
                                                            match_len=3)
            tests.base_test.append_copies(alignments, alignment, 2)

        # Gene 2 needs at least 1 perfect read since the bed file
        # annotates SJs but not the full isoform
        gene_2 = self._chromosome.genes[2]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            self._chr_name, gene_2, [0, 1, 3])
        alignments.append(alignment)

        read_id_offset += tests.base_test.write_sam_from_alignments(
            first_sam, chroms, alignments, read_id_offset)
        temp_file = '{}.tmp.sam'.format(first_sam)
        sort_log = os.path.join(self._log_dir, 'sort_1.log')
        tests.base_test.sort_sam_with_temp_file(first_sam, temp_file, sort_log)

        # This second_sam has no perfect reads for gene_0 to test that
        # perfect reads can influence other samples.
        second_sam = sam_paths[1]
        alignments = list()
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            self._chr_name, gene_0, [0, 1, 3])
        tests.base_test.add_mismatch_into_intron_number(alignment,
                                                        junction_i=1,
                                                        match_len=3)
        tests.base_test.append_copies(alignments, alignment, 2)
        read_id_offset += tests.base_test.write_sam_from_alignments(
            second_sam, chroms, alignments, read_id_offset)
        temp_file = '{}.tmp.sam'.format(second_sam)
        sort_log = os.path.join(self._log_dir, 'sort_2.log')
        tests.base_test.sort_sam_with_temp_file(second_sam, temp_file,
                                                sort_log)

    def _create_samples_tsv(self):
        self._samples_tsv = os.path.join(self._data_dir, 'samples.tsv')
        self._sample_names = list()
        with open(self._samples_tsv, 'wt') as handle:
            for sam_i, sam in enumerate(self._sams):
                sample_name = '{}_{}'.format(self._test_name, sam_i)
                self._sample_names.append(sample_name)
                columns = [sam, sample_name]
                tests.base_test.write_tsv_line(handle, columns)

    def _run_espresso_s(self, gtf, bed, num_cutoff, ratio_cutoff):
        command = [
            'perl', self._espresso_s, '-L', self._samples_tsv, '-F',
            self._fasta, '-O', self._work_dir, '--read_num_cutoff', num_cutoff,
            '--read_ratio_cutoff', ratio_cutoff
        ]
        if gtf:
            command.extend(['-A', gtf])

        if bed:
            command.extend(['--SJ_bed', bed])

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

        self._assert_exists_and_non_empty([
            self._samples_updated, sj_simplified_list, all_sjs, sam_list,
            sj_list
        ])

        sjs_from_simplified, sj_clusters = self._parse_sj_simplified(
            sj_simplified_list)
        expected_sjs = dict()
        unformatted_expected_values = [
            (self._chromosome.genes[0], ((0, 0), (1, 0)), 14, 14, False, True),
            (self._chromosome.genes[0], ((1, 0), (2, 0)), 5, 5, False, True),
            (self._chromosome.genes[0], ((2, 0), (3, 0)), 5, 5, False, True),
            (self._chromosome.genes[0], ((1, 3), (3, 0)), 0, 4, False, False),
            (self._chromosome.genes[1], ((1, 0), (2, 0)), 5, 5, False, True),
            (self._chromosome.genes[1], ((2, 0), (3, 0)), 5, 5, False, True),
            (self._chromosome.genes[1], ((1, 3), (3, 0)), 0, 2, False, False),
            (self._chromosome.genes[2], ((1, 0), (2, 0)), 5, 5, False, True),
            (self._chromosome.genes[2], ((2, 0), (3, 0)), 5, 5, False, True),
            (self._chromosome.genes[2], ((1, 3), (3, 0)), 0, 2, False, False),
        ]

        if config_name in ['num_and_ratio', 'bed_and_num_and_ratio']:
            unformatted_expected_values.append(
                (self._chromosome.genes[0], ((1, 0), (3, 0)), 2, 5, False,
                 True))
        else:
            unformatted_expected_values.append(
                (self._chromosome.genes[0], ((1, 0), (3, 0)), 2, 5, False,
                 False))

        if config_name == 'gtf':
            unformatted_expected_values.extend([
                (self._chromosome.genes[1], ((0, 0), (1, 0)), 9, 9, True,
                 True),
                (self._chromosome.genes[1], ((1, 0), (3, 0)), 0, 2, True,
                 True),
            ])
        else:
            unformatted_expected_values.extend([
                (self._chromosome.genes[1], ((0, 0), (1, 0)), 9, 9, False,
                 True),
                (self._chromosome.genes[1], ((1, 0), (3, 0)), 0, 2, False,
                 False),
            ])

        if config_name in ['bed', 'bed_and_num_and_ratio']:
            unformatted_expected_values.extend([
                (self._chromosome.genes[2], ((0, 0), (1, 0)), 10, 10, True,
                 True),
                (self._chromosome.genes[2], ((1, 0), (3, 0)), 1, 3, True,
                 True),
            ])
        else:
            unformatted_expected_values.extend([
                (self._chromosome.genes[2], ((0, 0), (1, 0)), 10, 10, False,
                 True),
                (self._chromosome.genes[2], ((1, 0), (3, 0)), 1, 3, False,
                 False),
            ])

        for unformatted in unformatted_expected_values:
            gene = unformatted[0]
            exon_indices_and_offsets = unformatted[1]
            perfect_count = unformatted[2]
            all_count = unformatted[3]
            is_annotated = unformatted[4]
            is_high_confidence = unformatted[5]
            start_index, start_offset = exon_indices_and_offsets[0]
            end_index, end_offset = exon_indices_and_offsets[1]
            first_exon = gene.exons[start_index]
            second_exon = gene.exons[end_index]
            sj_key = (self._chr_name, first_exon.end + 1 + start_offset,
                      second_exon.start + end_offset)
            expected_values = {
                'chr': self._chr_name,
                'start': sj_key[1],
                'end': sj_key[2],
                'perfect_count': perfect_count,
                'all_count': all_count,
                'is_annotated': is_annotated,
                'is_high_confidence': is_high_confidence,
            }
            # if the junction is offset then the intron motif may not be GTAG
            if not (start_offset or end_offset):
                expected_values['strand'] = gene.strand

            expected_sjs[sj_key] = expected_values

        self.assertEqual(len(sjs_from_simplified), len(expected_sjs))
        for expected_sj_key, expected_sj_details in expected_sjs.items():
            found = sjs_from_simplified[expected_sj_key]
            for detail, value in expected_sj_details.items():
                self.assertEqual(found[detail],
                                 value,
                                 msg='{}[{}] was {} but expected {}'.format(
                                     expected_sj_key, detail, found[detail],
                                     value))

    def _run_espresso_c(self):
        for sample_i in range(len(self._sample_names)):
            threads = '2'
            command = [
                'perl', self._espresso_c, '-I', self._work_dir, '-F',
                self._fasta, '-X',
                str(sample_i), '-T', threads
            ]
            c_log = os.path.join(self._log_dir,
                                 'espresso_c_{}.log'.format(sample_i))
            tests.base_test.run_command_with_log(command, c_log)

    def _check_c_output(self, config_name):
        gene_skip_coords = list()
        gene_moved_skip_coords = list()
        for gene_i in [0, 1, 2]:
            gene = self._chromosome.genes[gene_i]
            # (0-based first nt in intron, 1-based last nt in intron)
            gene_skip_coords.append([
                (gene.exons[0].end + 1, gene.exons[1].start),
                (gene.exons[1].end + 1, gene.exons[3].start),
            ])
            # The skipping junction is moved 3 nt into the intron
            gene_moved_skip_coords.append([
                (gene.exons[0].end + 1, gene.exons[1].start),
                (gene.exons[1].end + 4, gene.exons[3].start),
            ])

        for sample_i in range(len(self._sample_names)):
            read_final = os.path.join(
                self._work_dir, str(sample_i),
                '{}_read_final.txt'.format(self._chr_name))
            self._assert_exists_and_non_empty([read_final])
            parsed_read_final = self._parse_read_final(read_final)
            if sample_i == 0:
                self.assertEqual(len(parsed_read_final), 31)
            else:
                self.assertEqual(len(parsed_read_final), 2)

            for read_id, details in parsed_read_final.items():
                # Check the reads for the skipping isoform.
                # With the correct config the skipping junction should pass.
                # With the correct config the moved junction should be corrected.
                sj_coords = [(sj['start'], sj['end']) for sj in details['SJs']]
                handled = False
                for gene_i, correct_configs in [
                    (0, ['num_and_ratio', 'bed_and_num_and_ratio']),
                    (1, ['gtf']), (2, ['bed', 'bed_and_num_and_ratio'])
                ]:
                    if sj_coords not in [
                            gene_skip_coords[gene_i],
                            gene_moved_skip_coords[gene_i]
                    ]:
                        continue

                    self.assertEqual(details['SJs'][0]['corrected_status'],
                                     'pass')
                    if config_name in correct_configs:
                        if sj_coords == gene_moved_skip_coords[gene_i]:
                            self.assertEqual(
                                details['SJs'][1]['corrected_status'],
                                'corrected')
                        else:
                            self.assertEqual(
                                details['SJs'][1]['corrected_status'], 'pass')
                    else:
                        self.assertEqual(details['SJs'][1]['corrected_status'],
                                         'fail')

                    handled = True
                    break

                # All other junctions should 'pass'
                if not handled:
                    for sj_details in details['SJs']:
                        self.assertEqual(sj_details['corrected_status'],
                                         'pass')

                self.assertEqual(details['start']['corrected_status'], 'pass')
                self.assertEqual(details['end']['corrected_status'], 'pass')

    def _run_espresso_q(self, gtf, num_cutoff, ratio_cutoff):
        isoform_basename = '{}_isoform.tsv'.format(self._out_file_prefix)
        self._isoform_tsv = os.path.join(self._work_dir, isoform_basename)
        command = [
            'perl', self._espresso_q, '-L', self._samples_updated,
            '--read_num_cutoff', num_cutoff, '--read_ratio_cutoff',
            ratio_cutoff
        ]
        if gtf:
            command.extend(['-A', gtf])

        q_log = os.path.join(self._log_dir, 'espresso_q.log')
        tests.base_test.run_command_with_log(command, q_log)

    def _check_q_output(self, config_name):
        gtf_basename = '{}_updated.gtf'.format(self._out_file_prefix)
        self._updated_gtf = os.path.join(self._work_dir, gtf_basename)
        abundance_basename = '{}_abundance.esp'.format(self._out_file_prefix)
        self._abundance = os.path.join(self._work_dir, abundance_basename)
        self._assert_does_not_exist(self._isoform_tsv)
        self._assert_exists_and_non_empty([self._updated_gtf, self._abundance])
        abundance_rows = self._parse_abundance(self._abundance,
                                               self._sample_names)
        detected_isoforms = self._parse_updated_gtf(self._updated_gtf)

        expected_isoforms = list()
        for gene_i in [0, 1, 2]:
            expected_isoform = dict()
            expected_isoform['chr'] = self._chr_name
            expected_isoform['gene'] = self._chromosome.genes[gene_i]
            expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
            expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
            expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
            expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
            expected_isoform['abundance'] = {
                self._sample_names[0]: 5,
                self._sample_names[1]: 0
            }
            expected_isoform['is_novel'] = True
            expected_isoforms.append(expected_isoform)

            # skipping isoform [0,1,3]
            expected_isoform = dict()
            expected_isoform['chr'] = self._chr_name
            expected_isoform['gene'] = self._chromosome.genes[gene_i]
            expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
            expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
            expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
            expected_isoform['abundance'] = None
            expected_isoform['is_novel'] = True

            if ((config_name in ['num_and_ratio', 'bed_and_num_and_ratio'])
                    and (gene_i == 0)):
                expected_isoform['abundance'] = {
                    self._sample_names[0]: 7,
                    self._sample_names[1]: 2
                }
            if (config_name == 'gtf') and (gene_i == 1):
                expected_isoform['abundance'] = {
                    self._sample_names[0]: 4,
                    self._sample_names[1]: 0
                }
                expected_isoform['is_novel'] = False
                expected_isoform['transcript_name'] = 'GENE01T0'
            # The bed annotation is only for the junctions not the isoform.
            # Each junction in the novel isoform still needs to pass read num and ratio.
            if (config_name == 'bed_and_num_and_ratio') and (gene_i == 2):
                expected_isoform['abundance'] = {
                    self._sample_names[0]: 5,
                    self._sample_names[1]: 0
                }

            if expected_isoform['abundance'] is not None:
                expected_isoforms.append(expected_isoform)

        expected_by_id = self._get_transcript_ids_from_detected(
            expected_isoforms, detected_isoforms)
        self.assertEqual(len(abundance_rows), len(expected_isoforms))
        self._check_abundance_row_values(expected_by_id, abundance_rows,
                                         self._sample_names)


if __name__ == '__main__':
    unittest.main(verbosity=2)
