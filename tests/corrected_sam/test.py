import os
import os.path
import unittest

import tests.base_test


class CorrectedSamTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_dir = os.path.dirname(__file__)
        self._test_name = 'corrected_sam'
        self._data_dir = os.path.join(self._test_dir, 'data')
        self._out_dir = os.path.join(self._test_dir, 'out')
        self._log_dir = os.path.join(self._test_dir, 'logs')
        self._work_dir = os.path.join(self._out_dir, 'work_dir')
        self._corrected_out_dir = os.path.join(self._out_dir, 'corrected')

        self._out_file_prefix = (
            tests.base_test.out_file_prefix_based_on_params())
        self._chromosomes = None
        self._fasta = None
        self._gtf = None
        self._sams = None
        self._sample_names = None
        self._samples_tsv = None
        self._samples_updated = None
        self._abundance = None
        self._updated_gtf = None

    def test(self):
        self._initialize_dirs()
        self._create_test_data()
        self._create_samples_tsv()
        self._run_espresso_s()
        self._check_s_output()
        self._run_espresso_c()
        self._check_c_output()
        self._run_espresso_q()
        self._check_q_output()
        self._create_corrected_sams()
        self._check_corrected_sams()

    def _initialize_dirs(self):
        directories = [self._data_dir, self._out_dir, self._log_dir]
        tests.base_test.remove_then_create_directories(directories)

    def _create_test_data(self):
        self._fasta = os.path.join(self._data_dir, 'test.fasta')
        self._create_test_fasta(self._fasta)
        self._gtf = os.path.join(self._data_dir, 'test.gtf')
        self._create_test_gtf(self._gtf)
        self._sams = list()
        # Using the order 'sub_dir', '1', '2' to match sorted order of
        # samples in abundance output file.
        self._sams.append(os.path.join(self._data_dir, 'sub_dir', 'test.sam'))
        self._sams.append(os.path.join(self._data_dir, 'test.sam'))
        self._sams.append(os.path.join(self._data_dir, 'test_2.sam'))
        self._create_test_sams(self._sams)

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

        chr_id = 'chr2'
        gene_id = 'ENSG2'
        gene_name = 'GENE2'
        # 1 gene with 4 exons
        gene = self._get_random_gene(num_exons=4)
        gene.id = gene_id
        gene.name = gene_name
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id
        self._chromosomes.append(chromosome)

        # Test having a contig with no reads
        chr_id = 'chr3'
        gene_id = 'ENSG3'
        gene_name = 'GENE3'
        # 1 gene with 4 exons
        gene = self._get_random_gene(num_exons=4)
        gene.id = gene_id
        gene.name = gene_name
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id
        self._chromosomes.append(chromosome)

        tests.base_test.write_fasta_from_chroms(fasta_path, self._chromosomes)

    def _create_test_gtf(self, gtf_path):
        gene = self._chromosomes[0].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3], '0')

        gene = self._chromosomes[1].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3], '0')

        tests.base_test.write_gtf_from_chroms(gtf_path, self._chromosomes)

    def _get_alignments_1(self):
        alignments = list()
        chrom = self._chromosomes[0]
        chrom_seq = chrom.get_sequence_and_set_coords()
        exon_0_length = len(chrom.genes[0].exons[0].sequence)
        exon_3_length = len(chrom.genes[0].exons[3].sequence)
        # added 1st
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.trim_alignment_start(alignment, exon_0_length - 25)
        junc_op = tests.base_test.remove_junction_operation(alignment,
                                                            junction_i=0)
        alignment.start += junc_op.num
        alignment.sequence = alignment.sequence[(exon_0_length - 25):]
        tests.base_test.append_copies(alignments, alignment, 1)

        # added last
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.trim_alignment_end(alignment, exon_3_length - 25)
        tests.base_test.remove_junction_operation(alignment, junction_i=2)
        alignment.sequence = alignment.sequence[:-(exon_3_length - 25)]
        tests.base_test.append_copies(alignments, alignment, 1)

        # corrected internal SJ
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=5)
        tests.base_test.add_deletion_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=20,
                                             deletion_length=5)
        tests.base_test.append_copies(alignments, alignment, 1)

        # clipping
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.trim_alignment_start(alignment, exon_0_length - 25)
        tests.base_test.add_soft_clipping_at_alignment_start(
            alignment, exon_0_length - 25)
        tests.base_test.trim_alignment_end(alignment, exon_3_length - 25)
        tests.base_test.add_soft_clipping_at_alignment_end(
            alignment, exon_3_length - 25)
        tests.base_test.append_copies(alignments, alignment, 1)

        # added 1st, last, corrected internal, clipping
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.trim_alignment_start(alignment, exon_0_length - 25)
        tests.base_test.add_soft_clipping_at_alignment_start(
            alignment, exon_0_length - 25)
        tests.base_test.trim_alignment_end(alignment, exon_3_length - 25)
        tests.base_test.add_soft_clipping_at_alignment_end(
            alignment, exon_3_length - 25)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=5)
        tests.base_test.add_deletion_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=20,
                                             deletion_length=5)
        tests.base_test.remove_junction_operation(alignment, junction_i=2)
        junc_op = tests.base_test.remove_junction_operation(alignment,
                                                            junction_i=0)
        alignment.start += junc_op.num
        tests.base_test.append_copies(alignments, alignment, 1)

        chrom = self._chromosomes[1]
        chrom_seq = chrom.get_sequence_and_set_coords()
        exon_0_length = len(chrom.genes[0].exons[0].sequence)
        exon_3_length = len(chrom.genes[0].exons[3].sequence)
        # added 1st
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.trim_alignment_start(alignment, exon_0_length - 25)
        junc_op = tests.base_test.remove_junction_operation(alignment,
                                                            junction_i=0)
        alignment.start += junc_op.num
        alignment.sequence = alignment.sequence[(exon_0_length - 25):]
        tests.base_test.append_copies(alignments, alignment, 1)

        # added last
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.trim_alignment_end(alignment, exon_3_length - 25)
        tests.base_test.remove_junction_operation(alignment, junction_i=2)
        alignment.sequence = alignment.sequence[:-(exon_3_length - 25)]
        tests.base_test.append_copies(alignments, alignment, 1)

        # corrected internal SJ
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=5)
        tests.base_test.add_deletion_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=20,
                                             deletion_length=5)
        tests.base_test.append_copies(alignments, alignment, 1)

        # clipping
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.trim_alignment_start(alignment, exon_0_length - 25)
        tests.base_test.add_soft_clipping_at_alignment_start(
            alignment, exon_0_length - 25)
        tests.base_test.trim_alignment_end(alignment, exon_3_length - 25)
        tests.base_test.add_soft_clipping_at_alignment_end(
            alignment, exon_3_length - 25)
        tests.base_test.append_copies(alignments, alignment, 1)

        # added 1st, last, corrected internal, clipping
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.trim_alignment_start(alignment, exon_0_length - 25)
        tests.base_test.add_soft_clipping_at_alignment_start(
            alignment, exon_0_length - 25)
        tests.base_test.trim_alignment_end(alignment, exon_3_length - 25)
        tests.base_test.add_soft_clipping_at_alignment_end(
            alignment, exon_3_length - 25)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=5)
        tests.base_test.add_deletion_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=20,
                                             deletion_length=5)
        tests.base_test.remove_junction_operation(alignment, junction_i=2)
        junc_op = tests.base_test.remove_junction_operation(alignment,
                                                            junction_i=0)
        alignment.start += junc_op.num
        tests.base_test.append_copies(alignments, alignment, 1)

        return alignments

    def _get_alignments_2(self):
        alignments = list()
        chrom = self._chromosomes[0]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)

        chrom = self._chromosomes[1]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)

        return alignments

    def _get_alignments_sub_dir(self):
        alignments = list()
        chrom = self._chromosomes[0]
        chrom_seq = chrom.get_sequence_and_set_coords()
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        insert_len = 10
        tests.base_test.add_insertion_in_exon(alignment,
                                              exon_i=1,
                                              offset_from_start=20,
                                              insertion_length=insert_len)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.replace_insertion_sequence(alignment, insert_len)
        tests.base_test.append_copies(alignments, alignment, 2)

        chrom = self._chromosomes[1]
        chrom_seq = chrom.get_sequence_and_set_coords()
        exon_0_length = len(chrom.genes[0].exons[0].sequence)
        exon_1_length = len(chrom.genes[0].exons[1].sequence)
        exon_2_length = len(chrom.genes[0].exons[2].sequence)
        exon_3_length = len(chrom.genes[0].exons[3].sequence)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)

        # corrected 1st SJ that uses the originally clipped start
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.trim_alignment_start(alignment, exon_0_length - 30)
        tests.base_test.trim_alignment_end(alignment, 30)
        # Add insertion in sequence
        tests.base_test.add_insertion_in_exon(alignment,
                                              exon_i=1,
                                              offset_from_start=10,
                                              insertion_length=10)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        # Use insertion length as match in first exon
        alignment.cigar.operations.pop(3)
        alignment.cigar.operations[0].num += 10
        alignment.start -= 10
        # replace most of the (increased) inital match with clipping
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 35)
        alignment.cigar.operations[1].num -= 35
        # move the start to account for clipping
        alignment.start += 35
        # move the start back to offset the 1st SJ start
        alignment.start -= 10
        # increase the length of the 1st junction to correct the 1st SJ end
        alignment.cigar.operations[2].num += 10
        tests.base_test.append_copies(alignments, alignment, 1)

        # corrected last SJ that uses the originally clipped end
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.trim_alignment_start(alignment, 30)
        tests.base_test.trim_alignment_end(alignment, exon_3_length - 30)
        # Add insertion in sequence
        tests.base_test.add_insertion_in_exon(
            alignment,
            exon_i=2,
            offset_from_start=(exon_2_length - 10),
            insertion_length=10)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        # Use insertion length as match in last exon
        alignment.cigar.operations.pop(5)
        alignment.cigar.operations[-1].num += 10
        # replace most of the (increased) final match with clipping
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 35)
        alignment.cigar.operations[-2].num -= 35
        # Increase the length of the last junction to offset the SJ end
        alignment.cigar.operations[-3].num += 10
        tests.base_test.append_copies(alignments, alignment, 1)

        # corrected internal SJ with a read position before an earlier SJ
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.trim_alignment_end(alignment, exon_3_length - 30)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        # Move the 2nd SJ forward slightly
        alignment.cigar.operations[2].num += 5
        alignment.cigar.operations[4].num -= 5
        # Insert an incorrect SJ right before the actual 2nd SJ
        fake_junc_op = tests.base_test.CigarOp(2, 'N')
        after_junc_match_op = tests.base_test.CigarOp(2, 'M')
        alignment.cigar.operations.insert(3, fake_junc_op)
        alignment.cigar.operations.insert(4, after_junc_match_op)
        # Remove matches from the end to account for the newly added matches
        alignment.cigar.operations[-1].num -= 2
        # Adjust the 2nd SJ len to avoid offsetting later SJs
        alignment.cigar.operations[5].num -= 4
        tests.base_test.append_copies(alignments, alignment, 1)

        return alignments

    def _get_expected_isoforms(self):
        expected_isoforms = list()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {
            self._sample_names[0]: 2,
            self._sample_names[1]: 5,
            self._sample_names[2]: 2
        }
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[0].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[0].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)

        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[1].name
        expected_isoform['gene'] = self._chromosomes[1].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {
            self._sample_names[0]: 4,
            self._sample_names[1]: 5,
            self._sample_names[2]: 2
        }
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[1].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[1].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)

        return expected_isoforms

    def _create_test_sams(self, sam_paths):
        sam_path_sub_dir = sam_paths[0]
        sam_path_1 = sam_paths[1]
        sam_path_2 = sam_paths[2]
        alignments_1 = self._get_alignments_1()
        read_id_offset = 0
        read_id_offset = tests.base_test.write_sam_from_alignments(
            sam_path_1, self._chromosomes, alignments_1, read_id_offset)

        alignments_2 = self._get_alignments_2()
        read_id_offset = tests.base_test.write_sam_from_alignments(
            sam_path_2, self._chromosomes, alignments_2, read_id_offset)

        os.mkdir(os.path.join(self._data_dir, 'sub_dir'))
        alignments_sub_dir = self._get_alignments_sub_dir()
        read_id_offset = tests.base_test.write_sam_from_alignments(
            sam_path_sub_dir, self._chromosomes, alignments_sub_dir,
            read_id_offset)

        temp_file_1 = '{}.tmp.sam'.format(sam_path_1)
        sort_log_1 = os.path.join(self._log_dir, 'sort_1.log')
        tests.base_test.sort_sam_with_temp_file(sam_path_1, temp_file_1,
                                                sort_log_1)

        temp_file_2 = '{}.tmp.sam'.format(sam_path_2)
        sort_log_2 = os.path.join(self._log_dir, 'sort_2.log')
        tests.base_test.sort_sam_with_temp_file(sam_path_2, temp_file_2,
                                                sort_log_2)

        temp_file_sub_dir = '{}.tmp.sam'.format(sam_path_sub_dir)
        sort_log_sub_dir = os.path.join(self._log_dir, 'sort_sub_dir.log')
        tests.base_test.sort_sam_with_temp_file(sam_path_sub_dir,
                                                temp_file_sub_dir,
                                                sort_log_sub_dir)

    def _create_samples_tsv(self):
        self._samples_tsv = os.path.join(self._data_dir, 'samples.tsv')
        self._sample_names = ['sub_dir', '1', '2']
        with open(self._samples_tsv, 'wt') as handle:
            for sam_i, sam in enumerate(self._sams):
                columns = [sam, self._sample_names[sam_i]]
                tests.base_test.write_tsv_line(handle, columns)

    def _run_espresso_s(self):
        command = [
            'perl', self._espresso_s, '-L', self._samples_tsv, '-F',
            self._fasta, '-O', self._work_dir, '-A', self._gtf
        ]

        s_log = os.path.join(self._log_dir, 'espresso_s.log')
        tests.base_test.run_command_with_log(command, s_log)

    def _check_s_output(self):
        self._samples_updated = os.path.join(self._work_dir,
                                             'samples.tsv.updated')
        sj_simplified_lists = list()
        for chrom in self._chromosomes[:-1]:
            sj_simplified_lists.append(
                os.path.join(self._work_dir,
                             '{}_SJ_simplified.list'.format(chrom.name)))
        all_sjs = os.path.join(self._work_dir, 'SJ_group_all.fa')
        sam_lists = list()
        sj_lists = list()
        for sample_i in range(len(self._sample_names)):
            sam_list = os.path.join(self._work_dir, str(sample_i), 'sam.list3')
            sj_list = os.path.join(self._work_dir, str(sample_i), 'sj.list')
            sam_lists.append(sam_list)
            sj_lists.append(sj_list)

        expected_files = ([self._samples_updated, all_sjs] +
                          sj_simplified_lists + sam_lists + sj_lists)
        self._assert_exists_and_non_empty(expected_files)

    def _run_espresso_c(self):
        threads = '2'
        for sample_i in range(len(self._sample_names)):
            command = [
                'perl', self._espresso_c, '-I', self._work_dir, '-F',
                self._fasta, '-X',
                str(sample_i), '-T', threads
            ]
            c_log = os.path.join(self._log_dir,
                                 'espresso_c_{}.log'.format(sample_i))
            tests.base_test.run_command_with_log(command, c_log)

    def _check_c_output(self):
        for sample_i in range(len(self._sample_names)):
            for chrom in self._chromosomes[:-1]:
                read_final = os.path.join(
                    self._work_dir, str(sample_i),
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

    def _create_corrected_sams(self):
        num_threads = 4
        command = [
            self._py_executable, self._create_corrected_sam_py,
            '--samples-tsv', self._samples_tsv, '--espresso-out-dir',
            self._work_dir, '--out-dir', self._corrected_out_dir, '--fasta',
            self._fasta, '--libparasail-so-path', self._libparasail_so_path,
            '--num-threads',
            str(num_threads)
        ]
        corrected_sams_log = os.path.join(self._log_dir,
                                          'create_corrected_sams.log')
        tests.base_test.run_command_with_log(command, corrected_sams_log)

    def _check_corrected_sams(self):
        sam_mapping_path = os.path.join(self._corrected_out_dir,
                                        'sam_file_name_mapping.txt')
        sam_mapping = self._parse_corrected_sam_mapping(sam_mapping_path)
        corrected_sams = list()
        for old_sam in self._sams:
            abs_old_sam = os.path.abspath(old_sam)
            self.assertIn(abs_old_sam, sam_mapping)
            corrected_name = sam_mapping[abs_old_sam]
            corrected_sam = os.path.join(self._corrected_out_dir,
                                         corrected_name)
            corrected_sams.append(corrected_sam)

        self._assert_exists_and_non_empty(corrected_sams)
        expected_alignments = self._get_expected_corrected_alignments()
        for sam_i, corrected_sam in enumerate(corrected_sams):
            with open(corrected_sam, 'rt') as handle:
                for line in handle:
                    columns = line.rstrip('\n').split('\t')
                    alignment = self._parse_alignment_columns(columns)
                    if not alignment:
                        continue

                    self.assertIn(alignment['read_id'], expected_alignments)
                    expected_align = expected_alignments[alignment['read_id']]
                    del expected_alignments[alignment['read_id']]
                    self.assertEqual(sam_i, expected_align['sam_i'])
                    for key, value in expected_align.items():
                        if key == 'sam_i':
                            continue

                        self.assertIn(key, alignment)
                        self.assertEqual(alignment[key],
                                         value,
                                         msg='key: {}, read_id: {}'.format(
                                             key, alignment['read_id']))

        self.assertEqual(expected_alignments, dict())

    def _get_expected_corrected_alignments(self):
        expected = dict()
        chrom = self._chromosomes[0]
        chrom_seq = chrom.get_sequence_and_set_coords()
        exon_0_length = len(chrom.genes[0].exons[0].sequence)
        exon_3_length = len(chrom.genes[0].exons[3].sequence)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        read_0_cigar = alignment.cigar.copy()
        read_0_cigar.operations[0].num -= (exon_0_length - 25)
        expected['read-0'] = {
            'sam_i': 1,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1 + (exon_0_length - 25),
            'cigar': str(read_0_cigar),
            'sequence': alignment.sequence[(exon_0_length - 25):]
        }
        read_1_cigar = alignment.cigar.copy()
        read_1_cigar.operations[-1].num -= (exon_3_length - 25)
        expected['read-1'] = {
            'sam_i': 1,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(read_1_cigar),
            'sequence': alignment.sequence[:-(exon_3_length - 25)]
        }
        expected['read-2'] = {
            'sam_i': 1,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(alignment.cigar),
            'sequence': alignment.sequence
        }
        read_3_cigar = alignment.cigar.copy()
        read_3_cigar.operations[0].num -= (exon_0_length - 25)
        read_3_cigar.operations[-1].num -= (exon_3_length - 25)
        read_3_start_clip_op = tests.base_test.CigarOp(exon_0_length - 25, 'S')
        read_3_end_clip_op = tests.base_test.CigarOp(exon_3_length - 25, 'S')
        read_3_cigar.operations.insert(0, read_3_start_clip_op)
        read_3_cigar.operations.append(read_3_end_clip_op)
        expected['read-3'] = {
            'sam_i': 1,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1 + (exon_0_length - 25),
            'cigar': str(read_3_cigar),
            'sequence': alignment.sequence
        }
        expected['read-4'] = {
            'sam_i': 1,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(alignment.cigar),
            'sequence': alignment.sequence
        }
        expected['read-10'] = {
            'sam_i': 2,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(alignment.cigar),
            'sequence': alignment.sequence
        }
        expected['read-11'] = {
            'sam_i': 2,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(alignment.cigar),
            'sequence': alignment.sequence
        }
        read_14_15_cigar = alignment.cigar.copy()
        old_match_len = read_14_15_cigar.operations[2].num
        read_14_15_cigar.operations[2].num = 20
        read_14_15_cigar.operations.insert(3, tests.base_test.CigarOp(10, 'I'))
        read_14_15_cigar.operations.insert(
            4, tests.base_test.CigarOp(old_match_len - 20, 'M'))
        read_14_15_alignment = alignment.copy()
        tests.base_test.add_insertion_in_exon(read_14_15_alignment,
                                              exon_i=1,
                                              offset_from_start=20,
                                              insertion_length=10)
        read_14_15_alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.replace_insertion_sequence(read_14_15_alignment, 10)
        expected['read-14'] = {
            'sam_i': 0,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(read_14_15_cigar),
            'sequence': read_14_15_alignment.sequence
        }
        expected['read-15'] = {
            'sam_i': 0,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(read_14_15_cigar),
            'sequence': read_14_15_alignment.sequence
        }

        chrom = self._chromosomes[1]
        chrom_seq = chrom.get_sequence_and_set_coords()
        exon_0_length = len(chrom.genes[0].exons[0].sequence)
        exon_1_length = len(chrom.genes[0].exons[1].sequence)
        exon_2_length = len(chrom.genes[0].exons[2].sequence)
        exon_3_length = len(chrom.genes[0].exons[3].sequence)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        read_5_cigar = alignment.cigar.copy()
        read_5_cigar.operations[0].num -= (exon_0_length - 25)
        expected['read-5'] = {
            'sam_i': 1,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1 + (exon_0_length - 25),
            'cigar': str(read_5_cigar),
            'sequence': alignment.sequence[(exon_0_length - 25):]
        }
        read_6_cigar = alignment.cigar.copy()
        read_6_cigar.operations[-1].num -= (exon_3_length - 25)
        expected['read-6'] = {
            'sam_i': 1,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(read_6_cigar),
            'sequence': alignment.sequence[:-(exon_3_length - 25)]
        }
        expected['read-7'] = {
            'sam_i': 1,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(alignment.cigar),
            'sequence': alignment.sequence
        }
        read_8_cigar = alignment.cigar.copy()
        read_8_cigar.operations[0].num -= (exon_0_length - 25)
        read_8_cigar.operations[-1].num -= (exon_3_length - 25)
        read_8_start_clip_op = tests.base_test.CigarOp(exon_0_length - 25, 'S')
        read_8_end_clip_op = tests.base_test.CigarOp(exon_3_length - 25, 'S')
        read_8_cigar.operations.insert(0, read_8_start_clip_op)
        read_8_cigar.operations.append(read_8_end_clip_op)
        expected['read-8'] = {
            'sam_i': 1,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1 + (exon_0_length - 25),
            'cigar': str(read_8_cigar),
            'sequence': alignment.sequence
        }
        expected['read-9'] = {
            'sam_i': 1,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(alignment.cigar),
            'sequence': alignment.sequence
        }
        expected['read-12'] = {
            'sam_i': 2,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(alignment.cigar),
            'sequence': alignment.sequence
        }
        expected['read-13'] = {
            'sam_i': 2,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(alignment.cigar),
            'sequence': alignment.sequence
        }
        expected['read-16'] = {
            'sam_i': 0,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(alignment.cigar),
            'sequence': alignment.sequence
        }
        expected['read-17'] = {
            'sam_i': 0,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': str(alignment.cigar),
            'sequence': alignment.sequence
        }
        read_18_cigar = alignment.cigar.copy()
        read_18_cigar.operations[0].num = 1
        read_18_start_clip_op = tests.base_test.CigarOp(29, 'S')
        read_18_cigar.operations.insert(0, read_18_start_clip_op)
        read_18_cigar.operations[3].num -= 10
        read_18_pre_insertion_match_op = tests.base_test.CigarOp(10, 'M')
        read_18_insertion_op = tests.base_test.CigarOp(10, 'I')
        read_18_cigar.operations.insert(3, read_18_insertion_op)
        read_18_cigar.operations.insert(3, read_18_pre_insertion_match_op)
        read_18_cigar.operations[-1].num -= 30
        read_18_seq = alignment.sequence[(exon_0_length - 30):(exon_0_length +
                                                               10)]
        read_18_seq += 'N' * 10
        read_18_seq += alignment.sequence[(exon_0_length + 10):-30]
        expected['read-18'] = {
            'sam_i': 0,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].end + 1,
            'cigar': str(read_18_cigar),
            'sequence': read_18_seq
        }
        read_19_cigar = alignment.cigar.copy()
        read_19_cigar.operations[-1].num = 1
        read_19_end_clip_op = tests.base_test.CigarOp(29, 'S')
        read_19_cigar.operations.append(read_19_end_clip_op)
        # Final CigarOps will be MNMNMIMNMS
        read_19_cigar.operations[4].num -= 10
        read_19_post_insertion_match_op = tests.base_test.CigarOp(10, 'M')
        read_19_insertion_op = tests.base_test.CigarOp(10, 'I')
        read_19_cigar.operations.insert(5, read_19_post_insertion_match_op)
        read_19_cigar.operations.insert(5, read_19_insertion_op)
        read_19_cigar.operations[0].num -= 30
        read_19_seq_up_to_insert = (exon_0_length + exon_1_length +
                                    (exon_2_length - 10))
        read_19_seq = alignment.sequence[30:read_19_seq_up_to_insert]
        read_19_seq += 'N' * 10
        read_19_seq += alignment.sequence[read_19_seq_up_to_insert:-(
            exon_3_length - 30)]
        expected['read-19'] = {
            'sam_i': 0,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 31,
            'cigar': str(read_19_cigar),
            'sequence': read_19_seq
        }
        read_20_seq = alignment.sequence[:-(exon_3_length - 30)]
        expected['read-20'] = {
            'sam_i': 0,
            'contig': chrom.name,
            'pos': chrom.genes[0].exons[0].start + 1,
            'cigar': '*',
            'sequence': read_20_seq
        }

        return expected

    def _parse_alignment_columns(self, columns):
        if columns[0].startswith('@'):
            return None

        alignment = dict()
        alignment['read_id'] = columns[0]
        alignment['flag'] = columns[1]
        alignment['contig'] = columns[2]
        alignment['pos'] = int(columns[3])
        alignment['mapq'] = columns[4]
        alignment['cigar'] = columns[5]
        alignment['contig_next'] = columns[6]
        alignment['pos_next'] = columns[7]
        alignment['template_len'] = columns[8]
        alignment['sequence'] = columns[9]
        alignment['quality'] = columns[10]
        alignment['tags'] = columns[11:]
        return alignment

    def _parse_corrected_sam_mapping(self, path):
        mapping = dict()
        with open(path, 'rt') as handle:
            for line in handle:
                columns = line.rstrip('\n').split('\t')
                orig_path, new_name = columns
                mapping[orig_path] = new_name

        return mapping


if __name__ == '__main__':
    unittest.main(verbosity=2)
