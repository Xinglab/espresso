import os
import os.path
import unittest

import tests.base_test


class IsoformAssignmentTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'isoform_assignment'
        self._test_dir = os.path.dirname(__file__)
        self._data_dir = os.path.join(self._test_dir, 'data')
        self._out_dir = os.path.join(self._test_dir, 'out')
        self._log_dir = os.path.join(self._test_dir, 'logs')
        self._work_dir = os.path.join(self._out_dir, 'work_dir')
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
        chr_i = 1
        gene_i = 1
        # {:02d} pads with 0: 1 -> 01
        chr_id_format = 'chr{:02d}'
        gene_id_format = 'ENSG{:02d}'
        gene_name_format = 'GENE{:02d}'
        # chr1: check --cont_del_max 30
        # gene has annotated exons: [0, 1]
        # reads with length 30 deletion in exon 0
        # reads with length 31 deletion in exon 1
        # longer deletion should lead to detection of new isoform
        gene = self._get_random_gene(num_exons=2)
        gene.set_intron_start_end_bases(exon_i=0,
                                        offset_from_start=20,
                                        length=30)
        gene.set_intron_start_end_bases(exon_i=1,
                                        offset_from_start=20,
                                        length=31)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr2: check FSM ISM NIC NNC classification
        # gene has annotated exons: [0, 1, 2, 3], [0, 1, 3]
        # ISM: [0, 1], [2, 3]
        # NIC: [0, 2, 3]
        # NNC: [0, 1, 2, 3], but exon 2 is 10 nt longer
        gene = self._get_random_gene(num_exons=4)
        gene.set_intron_start_bases(intron_i=2, offset_from_start=10)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr3: Like chr2, but with reads for each annotated isoform
        # FSM: [0, 1, 2, 3], [0, 1, 3]
        # ISM: [0, 1], [2, 3]
        # NIC: [0, 2, 3]
        # NNC: [0, 1, 2, 3], but exon 2 is 10 nt longer
        # EM init: ([0,1,2,3], 3.66), ([0,1,3], 2.66), ([0,2,3], 3), ([0,1,2,3], 2.66)
        # EM 1st pass: ([0,1,2,3], 3.92), ([0,1,3], 2.59), ([0,2,3], 2.9), ([0,1,2,3], 2.59)
        # EM 2nd pass: ([0,1,2,3], 4.01), ([0,1,3], 2.57), ([0,2,3], 2.85), ([0,1,2,3], 2.57)
        gene = self._get_random_gene(num_exons=4)
        gene.set_intron_start_bases(intron_i=2, offset_from_start=10)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr4: check EM assignment without perfect reads
        # gene has annotated exons: [0, 1, 2],
        #                           [0, 1, 2] (but shortened by 50nt at start and end)
        # no perfect reads for either annotated isoform (endpoints off by 25nt)
        # reads can get assigned to either
        gene = self._get_random_gene(num_exons=3)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr5: check EM assignment with perfect reads
        # same annotated isoforms as chr4
        # perfect reads for both annotated isoforms
        # reads get assigned proportional to perfect reads
        gene = self._get_random_gene(num_exons=3)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr6: reads get assigned to longest compatible novel SJ chain
        # gene has annotated exons: [0, 1, 2, 3]
        # perfect reads for all 2 exon pairs: [0,1], [0,2], [0,3], [1,2], [1,3], [2,3]
        # perfect reads just skipping 1 or 2: [0,2,3], [0,1,3]
        # Compatible chains:
        # [0,1,2,3] -> [0,1,2,3]
        # [0,1] -> [0,1,3], [0,1,2,3]
        # [0,2] -> [0,2,3]
        # [0,3] -> no longer chain
        # [1,2] -> [0,1,2,3]
        # [1,3] -> [0,1,3]
        # [2,3] -> [0,2,3], [0,1,2,3]
        # [0,2,3] -> no longer chain
        # [0,1,3] -> no longer chain
        # EM init: ([0,1,2,3], 6), ([0,3], 2), ([0,2,3], 5), ([0,1,3], 5)
        # EM 1st pass: ([0,1,2,3], 6.18), ([0,3], 2), ([0,2,3], 4.91), ([0,1,3], 4.91)
        # EM 2nd pass: ([0,1,2,3], 6.22), ([0,3], 2), ([0,2,3], 4.89), ([0,1,3], 4.89)
        gene = self._get_random_gene(num_exons=4)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr7: single exon gene
        gene = self._get_random_gene(num_exons=1)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr8: reads for an SJ chain should not be compatible with
        # isoforms if the read goes --SJ_dist into an intron.
        # gene has annotated exons [0, 1, 2, 3]
        # novel isoform extends exon 1 by 10nt and skips exon 2 [0, 1', 3]
        # [0, 1] could be assigned to either isoform
        # [0, 1] and 35 into intron at end of read: no isoforms
        # [0, 1] and 25 into intron at end of read: novel isoform
        # [0, 1] and 15 into intron at end of read: either isoform
        # EM init: ([0,1,2,3], 4), ([0,1',3], 2), ([0,1], (6 either, 2 novel))
        # EM 1st pass: ([0,1,2,3], 8), ([0,1',3], 6)
        # EM 2nd pass: ([0,1,2,3], 7.42), ([0,1',3], 6.58)
        # EM 3rd pass: ([0,1,2,3], 7.07), ([0,1',3], 6.93)
        # EM 4th pass: ([0,1,2,3], 7.03), ([0,1',3], 6.97)
        gene = self._get_random_gene(num_exons=4)
        gene.set_intron_start_bases(intron_i=1, offset_from_start=10)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr9: similar to chr8, but where the novel isoform shortens exon 1 by 10nt
        # [0, 1] could be assigned to either isoform
        # [0, 1] and 25 into intron at end of read: no isoforms
        # [0, 1] and 15 into intron at end of read: annotated isoform
        # [0, 1] and 5 into intron at end of read: either isoform
        gene = self._get_random_gene(num_exons=4)
        gene.set_intron_start_bases(intron_i=1, offset_from_start=-10)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        tests.base_test.write_fasta_from_chroms(fasta_path, self._chromosomes)

    def _create_test_gtf(self, gtf_path):
        gene = self._chromosomes[0].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1], '0')

        gene = self._chromosomes[1].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3], '0')
        gene.add_isoform_for_exon_numbers([0, 1, 3], '1')

        gene = self._chromosomes[2].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3], '0')
        gene.add_isoform_for_exon_numbers([0, 1, 3], '1')

        gene = self._chromosomes[3].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2], '0')
        gene.add_isoform_for_exon_numbers([0, 1, 2], '1')
        # Remove 50nt from start of exon 0
        exon = gene.isoforms[-1].exons[0].copy()
        gene.isoforms[-1].exons[0] = exon
        exon.start += 50
        exon.sequence = exon.sequence[50:]
        # Remove 50nt from end of exon 2
        exon = gene.isoforms[-1].exons[2].copy()
        gene.isoforms[-1].exons[2] = exon
        exon.end -= 50
        exon.sequence = exon.sequence[:-50]

        gene = self._chromosomes[4].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2], '0')
        gene.add_isoform_for_exon_numbers([0, 1, 2], '1')
        # Remove 50nt from start of exon 0
        exon = gene.isoforms[-1].exons[0].copy()
        gene.isoforms[-1].exons[0] = exon
        exon.start += 50
        exon.sequence = exon.sequence[50:]
        # Remove 50nt from end of exon 2
        exon = gene.isoforms[-1].exons[2].copy()
        gene.isoforms[-1].exons[2] = exon
        exon.end -= 50
        exon.sequence = exon.sequence[:-50]

        gene = self._chromosomes[5].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3], '0')

        gene = self._chromosomes[6].genes[0]
        gene.add_isoform_for_exon_numbers([0], '0')

        gene = self._chromosomes[7].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3], '0')

        gene = self._chromosomes[8].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3], '0')

        tests.base_test.write_gtf_from_chroms(gtf_path, self._chromosomes)

    def _create_test_sam(self, sam_path):
        alignments = list()
        read_id_offset = 0
        chrom = self._chromosomes[0]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_deletion_in_exon(alignment,
                                             exon_i=0,
                                             offset_from_start=20,
                                             deletion_length=30)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_deletion_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=20,
                                             deletion_length=31)
        tests.base_test.append_copies(alignments, alignment, 2)

        chrom = self._chromosomes[1]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=2,
                                                 offset=10)
        tests.base_test.append_copies(alignments, alignment, 2)

        chrom = self._chromosomes[2]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=2,
                                                 offset=10)
        tests.base_test.append_copies(alignments, alignment, 2)

        chrom = self._chromosomes[3]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.trim_alignment_start(alignment, 25)
        tests.base_test.trim_alignment_end(alignment, 25)
        tests.base_test.append_copies(alignments, alignment, 10)

        chrom = self._chromosomes[4]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.trim_alignment_start(alignment, 25)
        tests.base_test.trim_alignment_end(alignment, 25)
        tests.base_test.append_copies(alignments, alignment, 10)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.trim_alignment_start(alignment, 50)
        tests.base_test.trim_alignment_end(alignment, 50)
        tests.base_test.append_copies(alignments, alignment, 8)

        chrom = self._chromosomes[5]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 3])
        tests.base_test.append_copies(alignments, alignment, 2)

        chrom = self._chromosomes[6]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0])
        tests.base_test.append_copies(alignments, alignment, 2)

        chrom = self._chromosomes[7]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 4)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 3])
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=10)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.append_copies(alignments, alignment, 4)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_match_at_end_of_alignment(alignment, 35)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_match_at_end_of_alignment(alignment, 25)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_match_at_end_of_alignment(alignment, 15)
        tests.base_test.append_copies(alignments, alignment, 2)

        chrom = self._chromosomes[8]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 3])
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=-10)
        tests.base_test.append_copies(alignments, alignment, 4)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.append_copies(alignments, alignment, 4)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_match_at_end_of_alignment(alignment, 25)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_match_at_end_of_alignment(alignment, 15)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_match_at_end_of_alignment(alignment, 5)
        tests.base_test.append_copies(alignments, alignment, 2)

        read_id_offset += tests.base_test.write_sam_from_alignments(
            sam_path, self._chromosomes, alignments, read_id_offset)
        temp_file = '{}.tmp.sam'.format(sam_path)
        sort_log = os.path.join(self._log_dir, 'sort_1.log')
        tests.base_test.sort_sam_with_temp_file(sam_path, temp_file, sort_log)

    def _create_samples_tsv(self):
        self._samples_tsv = os.path.join(self._data_dir, 'samples.tsv')
        self._sample_names = [self._test_name]
        with open(self._samples_tsv, 'wt') as handle:
            columns = [self._sam, self._test_name]
            tests.base_test.write_tsv_line(handle, columns)

    def _run_espresso_s(self):
        command = [
            'perl', self._espresso_s, '-L', self._samples_tsv, '-F',
            self._fasta, '-O', self._work_dir, '-A', self._gtf,
            '--cont_del_max', '30'
        ]

        s_log = os.path.join(self._log_dir, 'espresso_s.log')
        tests.base_test.run_command_with_log(command, s_log)

    def _check_s_output(self):
        self._samples_updated = os.path.join(self._work_dir,
                                             'samples.tsv.updated')
        sj_simplified_lists = list()
        for chrom in self._chromosomes:
            sj_simplified_lists.append(
                os.path.join(self._work_dir,
                             '{}_SJ_simplified.list'.format(chrom.name)))

        # chr7 is a single exon isoform
        sj_lists_without_7 = sj_simplified_lists[:]
        chr_7_sj_simplified = sj_lists_without_7.pop(6)

        all_sjs = os.path.join(self._work_dir, 'SJ_group_all.fa')
        sam_list = os.path.join(self._work_dir, '0', 'sam.list3')
        sj_list = os.path.join(self._work_dir, '0', 'sj.list')

        self._assert_exists_and_non_empty(
            [self._samples_updated, all_sjs, sam_list, sj_list] +
            sj_lists_without_7)
        self._assert_does_not_exist_or_is_empty([chr_7_sj_simplified])

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
            self._gtf, '--SJ_dist', '20'
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

        expected_isoforms = list()
        # chr1: [0, 1]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[0].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[0].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # deletion splits exon 1: [0, 1_a, 1_b]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        # 1_a is first 20 nt
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[1].copy())
        orig_start = expected_isoform['exons'][-1].start
        expected_isoform['exons'][-1].end = orig_start + 19
        expected_isoform['exons'][-1].sequence = (
            expected_isoform['exons'][-1].sequence[:20])
        # 1_b is rest of exon after length 31 deletion
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[1].copy())
        expected_isoform['exons'][-1].start = orig_start + 51
        expected_isoform['exons'][-1].sequence = (
            expected_isoform['exons'][-1].sequence[51:])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[0].name, 0,
                                                  0))
        expected_isoforms.append(expected_isoform)

        # chr2: [0, 2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[1].name
        expected_isoform['gene'] = self._chromosomes[1].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 4}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[1].name, 1,
                                                  1))
        expected_isoforms.append(expected_isoform)
        # [0, 1, 2, 3] with exon 2 10 nt longer
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[1].name
        expected_isoform['gene'] = self._chromosomes[1].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[2].copy())
        orig_end = expected_isoform['exons'][-1].end
        expected_isoform['exons'][-1].end += 10
        chrom_sequence = self._chromosomes[1].get_sequence_and_set_coords()
        extra_sequence = chrom_sequence[orig_end + 1:orig_end + 11]
        expected_isoform['exons'][-1].sequence += extra_sequence
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 4}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[1].name, 1,
                                                  0))
        expected_isoforms.append(expected_isoform)

        # chr3: [0, 1, 2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[2].name
        expected_isoform['gene'] = self._chromosomes[2].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 4.01}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[2].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[2].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 1, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[2].name
        expected_isoform['gene'] = self._chromosomes[2].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 2.57}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[2].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[2].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)
        # [0, 2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[2].name
        expected_isoform['gene'] = self._chromosomes[2].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 2.85}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[2].name, 2,
                                                  1))
        expected_isoforms.append(expected_isoform)
        # [0, 1, 2, 3] with exon 2 10 nt longer
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[2].name
        expected_isoform['gene'] = self._chromosomes[2].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[2].copy())
        orig_end = expected_isoform['exons'][-1].end
        expected_isoform['exons'][-1].end += 10
        chrom_sequence = self._chromosomes[2].get_sequence_and_set_coords()
        extra_sequence = chrom_sequence[orig_end + 1:orig_end + 11]
        expected_isoform['exons'][-1].sequence += extra_sequence
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 2.57}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[2].name, 2,
                                                  0))
        expected_isoforms.append(expected_isoform)

        # chr4: [0, 1, 2]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[3].name
        expected_isoform['gene'] = self._chromosomes[3].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['abundance'] = {self._test_name: 5}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[3].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[3].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 1, 2] (shortened by 50nt at start and end)
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[3].name
        expected_isoform['gene'] = self._chromosomes[3].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0].copy()]
        expected_isoform['exons'][0].start += 50
        expected_isoform['exons'][0].sequence = (
            expected_isoform['exons'][0].sequence[50:])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[2].copy())
        expected_isoform['exons'][-1].end -= 50
        expected_isoform['exons'][-1].sequence = (
            expected_isoform['exons'][0].sequence[:-50])
        expected_isoform['abundance'] = {self._test_name: 5}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[3].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[3].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)

        # chr5: [0, 1, 2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[4].name
        expected_isoform['gene'] = self._chromosomes[4].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['abundance'] = {self._test_name: 4}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[4].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[4].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 1, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[4].name
        expected_isoform['gene'] = self._chromosomes[4].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0].copy()]
        expected_isoform['exons'][0].start += 50
        expected_isoform['exons'][0].sequence = (
            expected_isoform['exons'][0].sequence[50:])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[2].copy())
        expected_isoform['exons'][-1].end -= 50
        expected_isoform['exons'][-1].sequence = (
            expected_isoform['exons'][0].sequence[:-50])
        expected_isoform['abundance'] = {self._test_name: 16}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[4].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[4].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)

        # chr6: [0, 1, 2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[5].name
        expected_isoform['gene'] = self._chromosomes[5].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 6.22}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[5].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[5].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[5].name
        expected_isoform['gene'] = self._chromosomes[5].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[5].name, 5,
                                                  0))
        expected_isoforms.append(expected_isoform)
        # [0, 2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[5].name
        expected_isoform['gene'] = self._chromosomes[5].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 4.89}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[5].name, 5,
                                                  2))
        expected_isoforms.append(expected_isoform)
        # [0, 1, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[5].name
        expected_isoform['gene'] = self._chromosomes[5].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 4.89}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[5].name, 5,
                                                  1))
        expected_isoforms.append(expected_isoform)

        # chr7: [0]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[6].name
        expected_isoform['gene'] = self._chromosomes[6].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[6].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[6].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)

        # chr8: [0, 1, 2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[7].name
        expected_isoform['gene'] = self._chromosomes[7].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 7.03}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[7].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[7].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 1, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[7].name
        expected_isoform['gene'] = self._chromosomes[7].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[1].copy())
        orig_end = expected_isoform['exons'][-1].end
        expected_isoform['exons'][-1].end += 10
        chrom_sequence = self._chromosomes[7].get_sequence_and_set_coords()
        extra_sequence = chrom_sequence[orig_end + 1:orig_end + 11]
        expected_isoform['exons'][-1].sequence += extra_sequence
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 6.97}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[7].name, 7,
                                                  0))
        expected_isoforms.append(expected_isoform)

        # chr9: [0, 1, 2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[8].name
        expected_isoform['gene'] = self._chromosomes[8].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 6.97}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[8].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[8].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 1, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[8].name
        expected_isoform['gene'] = self._chromosomes[8].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[1].copy())
        expected_isoform['exons'][-1].end -= 10
        expected_isoform['exons'][-1].sequence = (
            expected_isoform['exons'][-1].sequence[:-10])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 7.03}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[8].name, 8,
                                                  0))
        expected_isoforms.append(expected_isoform)

        expected_by_id = self._get_transcript_ids_from_detected(
            expected_isoforms, detected_isoforms)
        self.assertEqual(len(abundance_rows), len(expected_isoforms))
        self._check_abundance_row_values(expected_by_id, abundance_rows,
                                         self._sample_names)


if __name__ == '__main__':
    unittest.main(verbosity=2)
