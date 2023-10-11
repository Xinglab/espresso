import os
import os.path
import unittest

import tests.base_test


class IsoformAssignmentTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'isoform_assignment'
        self._test_dir = os.path.dirname(__file__)
        self._data_dir = os.path.join(self._test_dir,
                                      '{}_data'.format(self._test_name))
        self._out_dir = os.path.join(self._test_dir,
                                     '{}_out'.format(self._test_name))
        self._log_dir = os.path.join(self._test_dir,
                                     '{}_logs'.format(self._test_name))
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
        self._compat_isoforms = None
        self._expected_compat_isoforms = None

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
        # isoforms if the read goes --internal_boundary_limit into an intron.
        # gene has annotated exons [0, 1, 2, 3]
        # novel isoform extends exon 1 by 10nt and skips exon 2 [0, 1', 3]
        # [0, 1] could be assigned to either isoform
        # [0, 1] and 20 into intron at end of read: no isoforms
        # [0, 1] and 14 into intron at end of read: novel isoform
        # [0, 1] and 4 into intron at end of read: either isoform
        # Median end for [0, 1] chain needs to be compatible with both isoforms so
        # that the individual reads can be checked against both isoforms
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
        # [0, 1] only assigned to annotated isoform
        # [0, 1] and 10 into intron at end of read: no isoforms
        # [0, 1] and 4 into intron at end of read: annotated isoform
        # [0, 1] with exon 1 shortened by 6: either isoform
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
        self._expected_compat_isoforms = list()
        read_id_offset = 0
        chrom = self._chromosomes[0]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_deletion_in_exon(alignment,
                                             exon_i=0,
                                             offset_from_start=20,
                                             deletion_length=30)
        compat = {
            'class': 'FSM',
            'isoforms': [chrom.genes[0].isoforms[0].id],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_deletion_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=20,
                                             deletion_length=31)
        compat = {
            'class': 'NNC',
            'isoforms':
            [tests.base_test.get_espresso_novel_id(chrom.name, 0, 0)],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)

        chrom = self._chromosomes[1]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        compat = {
            'class': 'ISM',
            'isoforms':
            [tests.base_test.get_espresso_novel_id(chrom.name, 1, 0)],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        compat = {
            'class': 'ISM',
            'isoforms':
            [tests.base_test.get_espresso_novel_id(chrom.name, 1, 1)],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2, 3])
        compat = {
            'class': 'NIC',
            'isoforms':
            [tests.base_test.get_espresso_novel_id(chrom.name, 1, 1)],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=2,
                                                 offset=10)
        compat = {
            'class': 'NNC',
            'isoforms':
            [tests.base_test.get_espresso_novel_id(chrom.name, 1, 0)],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)

        chrom = self._chromosomes[2]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        compat = {
            'class': 'FSM',
            'isoforms': [chrom.genes[0].isoforms[0].id],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 3])
        compat = {
            'class': 'FSM',
            'isoforms': [chrom.genes[0].isoforms[1].id],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        compat = {
            'class':
            'ISM',
            'isoforms': [
                chrom.genes[0].isoforms[0].id, chrom.genes[0].isoforms[1].id,
                tests.base_test.get_espresso_novel_id(chrom.name, 2, 0)
            ],
            'sample':
            self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        compat = {
            'class':
            'ISM',
            'isoforms': [
                chrom.genes[0].isoforms[0].id,
                tests.base_test.get_espresso_novel_id(chrom.name, 2, 1)
            ],
            'sample':
            self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2, 3])
        compat = {
            'class': 'NIC',
            'isoforms':
            [tests.base_test.get_espresso_novel_id(chrom.name, 2, 1)],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=2,
                                                 offset=10)
        compat = {
            'class': 'NNC',
            'isoforms':
            [tests.base_test.get_espresso_novel_id(chrom.name, 2, 0)],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
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

        # Add minimal compat isoform info for alignments on chrs 4->6
        compat = {'sample': self._test_name}
        num_without_compat_info = (len(alignments) -
                                   len(self._expected_compat_isoforms))
        self._expected_compat_isoforms.extend([compat] *
                                              num_without_compat_info)
        chrom = self._chromosomes[6]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0])
        compat = {
            'class': 'single-exon',
            'isoforms': [chrom.genes[0].isoforms[0].id],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
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
        tests.base_test.add_match_at_end_of_alignment(alignment, 20)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_match_at_end_of_alignment(alignment, 14)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_match_at_end_of_alignment(alignment, 4)
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
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_match_at_end_of_alignment(alignment, 10)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.add_match_at_end_of_alignment(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1])
        tests.base_test.trim_alignment_end(alignment, 6)
        tests.base_test.append_copies(alignments, alignment, 7)

        read_id_offset += tests.base_test.write_sam_from_alignments(
            sam_path, self._chromosomes, alignments, read_id_offset)
        temp_file = '{}.tmp.sam'.format(sam_path)
        sort_log = os.path.join(self._log_dir, 'sort.log')
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
        compat_basename = '{}_compatible_isoform.tsv'.format(
            self._out_file_prefix)
        self._compat_isoforms = os.path.join(self._work_dir, compat_basename)
        command = [
            'perl', self._espresso_q, '-L', self._samples_updated, '-A',
            self._gtf, '--SJ_dist', '20', '-V', self._compat_isoforms
        ]

        q_log = os.path.join(self._log_dir, 'espresso_q.log')
        tests.base_test.run_command_with_log(command, q_log)

    def _check_q_output(self):
        gtf_basename = '{}_updated.gtf'.format(self._out_file_prefix)
        self._updated_gtf = os.path.join(self._work_dir, gtf_basename)
        abundance_basename = '{}_abundance.esp'.format(self._out_file_prefix)
        self._abundance = os.path.join(self._work_dir, abundance_basename)
        self._assert_exists_and_non_empty(
            [self._updated_gtf, self._abundance, self._compat_isoforms])
        detected_isoforms = self._parse_updated_gtf(self._updated_gtf)
        abundance_rows = self._parse_abundance(self._abundance,
                                               self._sample_names)
        compat_isoforms = self._parse_compat_isoforms(self._compat_isoforms)
        self._check_compat_isoforms(self._expected_compat_isoforms,
                                    compat_isoforms)

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
        expected_isoform['abundance'] = {self._test_name: 10.15}
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
        expected_isoform['abundance'] = {self._test_name: 6.85}
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


class ReadEndpointsTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'read_endpoints'
        self._test_dir = os.path.dirname(__file__)
        self._data_dir = os.path.join(self._test_dir,
                                      '{}_data'.format(self._test_name))
        self._out_dir = os.path.join(self._test_dir,
                                     '{}_out'.format(self._test_name))
        self._log_dir = os.path.join(self._test_dir,
                                     '{}_logs'.format(self._test_name))
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
        self._compat_isoforms = None
        self._expected_compat_isoforms = None

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
        # chr1: check added 1st junction
        # annotated [0, 1, 2], [0, 1, 2, 3]
        # read with only X nt in exon 0
        # alignment only has [1,2] junction with match past start of exon 1
        # the [0,1] junction can be added depending on:
        #   number of bases aligned to exon 0
        #   number of inserted or deleted bases in exon 1
        # Uncorrected alignments will be ISM but
        # be incompatible if more than --internal_boundary_limit into intron.
        gene = self._get_random_gene(num_exons=4)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr2: like chr 1 but for last junction
        gene = self._get_random_gene(num_exons=4)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr3: NNC reads
        # Chain endpoints are based on median start,end of reads for that chain.
        # A chain is incompatible if the read endpoint goes through a
        # terminal exon boundary by 35 or an internal exon boundary by 6.
        # exon 2 is not in gtf.
        # [1,2] and [2,3] pass all endpoint checks
        #
        # Use a 5 exon gene but use only the middle 3 exons so that
        # there is room to adjust start and end positions past the
        # middle exons.
        gene = self._get_random_gene(num_exons=5)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr4: like chr3 but [1,2] and [2,3] fail start checks
        gene = self._get_random_gene(num_exons=5)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr5: like chr3 but [1,2] and [2,3] fail end checks
        gene = self._get_random_gene(num_exons=5)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr6: NIC can be substring of NNC
        # annotated: [0, 1, 4], [0, 2, 4]
        # NNC: [1, 2, 3]
        # NIC: [1, 2]
        gene = self._get_random_gene(num_exons=5)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr7: like chr6 but NIC fails endpoint checks
        gene = self._get_random_gene(num_exons=5)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr8: NIC substring of NIC
        # annotated: [0, 1, 2, 3, 4, 5]
        # NIC: [1, 3, 4], [1, 3]
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr9: like chr8 but fails endpoint check
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr10: ISM can be substring of NIC or NNC
        # annotated: [0, 1, 2, 3, 4, 5]
        # NIC: [1, 3, 4]
        # ISM: [3, 4]
        # Individual reads of ISM chain will have their endpoints checked
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr11: like chr10 but ISM fails endpoint checks
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr12: ISM substring of a longest novel ISM
        # annotated: [0, 1, 2, 3, 4, 5] (no read)
        #            [0, 2, 3, 4] with 4 extended by 8nt (no read)
        # NIC: [2, 3, 5]
        # ISM: [1, 2, 3, 4] but with exon 4 extended by 4nt.
        #                   Longest novel ISM
        # ISM: [2, 3, 4] but with exon 4 extended by 8nt.
        #                Substring of longest ISM.
        # ISM: [2, 3] substring of NIC and both ISM.
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr13: like chr12 but longest ISM fails start endpoint check
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr14: like chr12 but there is a read for [0, 2, 3, 4]
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr15: like chr12 but there is a read for [0, 1, 2, 3, 4, 5]
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr16: like chr12 but there is a read for both annotated isoforms
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr17: Filtered chain compatible with annotated and novel
        # annotated: [0, 1, 2, 3, 4, 5]
        # NNC: [1, 2, 3, 4] with exon 1 start shortened by 34nt
        #                   and exon 2 end extended by 10nt
        #                   and exon 4 end extended by 10nt
        # filtered: [1, 2, 3, 4] where 2->3 junction is filtered
        # filtered: [2, 3] where 2->3 junction is filtered (no unfiltered SJ)
        gene = self._get_random_gene(num_exons=6)
        gene.set_intron_start_bases(intron_i=2, offset_from_start=10)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr18: like chr17 but filtered chain fails endpoint check for annotated
        gene = self._get_random_gene(num_exons=6)
        gene.set_intron_start_bases(intron_i=2, offset_from_start=10)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr19: like chr17 but filtered chain fails endpoint check for novel
        gene = self._get_random_gene(num_exons=6)
        gene.set_intron_start_bases(intron_i=2, offset_from_start=10)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr20: Chain with 1st or last SJ=x has endpoints checked
        # annotated: [0, 1, 2, 3, 4, 5]
        # NIC: [0, 1, 2, 3, 5]
        # filtered: [1, 2, 3] where 1->2 junction is filtered
        #           [1, 2, 3] where 2->3 junction is filtered
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr21: like chr20 but where endpoint checks fail
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr22: Single exon isoform
        # annotated: [1], [1] with endpoints shifted by 20nt
        gene = self._get_random_gene(num_exons=3)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr23: novel chain individual read endpoints are not checked
        # annotated: [0, 2, 4, 6] [1, 3, 5]
        # NIC: [1, 2, 3, 4], [1, 3, 4], [3, 4]
        gene = self._get_random_gene(num_exons=7)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr24: FSM reads can exceed endpoints.
        #        Filtered chains match against FSM chain endpoints.
        #        ISM chains match against FSM isoform endpoints.
        # annotated: [1, 2, 3, 4, 5], and with exon 1 trimmed by 40nt
        # ISM: [1, 2, 3, 4], [2, 3, 4, 5]
        # filtered: [1, 2, 3, 4, 5] 3->4 fail
        # filtered: [1, 2, 3, 4, 5] 2->3 fail
        gene = self._get_random_gene(num_exons=7)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr25: ISM chains have individual read endpoints checked
        #        even if there is only 1 compatible chain
        # annotated: [0, 1, 2, 3, 4, 5] (no read)
        #            [0, 1, 3, 4, 5]
        # ISM: [1, 2, 3, 4]
        #      [1, 2, 3]
        #      [1, 3]
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        tests.base_test.write_fasta_from_chroms(fasta_path, self._chromosomes)

    def _create_test_gtf(self, gtf_path):
        gene = self._chromosomes[0].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2], '0')
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3], '1')

        gene = self._chromosomes[1].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2], '0')
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3], '1')

        gene = self._chromosomes[2].genes[0]
        gene.add_isoform_for_exon_numbers([1, 3], '0')

        gene = self._chromosomes[3].genes[0]
        gene.add_isoform_for_exon_numbers([1, 3], '0')

        gene = self._chromosomes[4].genes[0]
        gene.add_isoform_for_exon_numbers([1, 3], '0')

        gene = self._chromosomes[5].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 4], '0')
        gene.add_isoform_for_exon_numbers([0, 2, 4], '1')

        gene = self._chromosomes[6].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 4], '0')
        gene.add_isoform_for_exon_numbers([0, 2, 4], '1')

        gene = self._chromosomes[7].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')

        gene = self._chromosomes[8].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')

        gene = self._chromosomes[9].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')

        gene = self._chromosomes[10].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')

        gene = self._chromosomes[11].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')
        gene.add_isoform_for_exon_numbers([0, 2, 3, 4], '1')
        # Add 8nt to end of exon 4
        chr_seq = self._chromosomes[11].get_sequence_and_set_coords()
        exon = gene.isoforms[-1].exons[3].copy()
        gene.isoforms[-1].exons[3] = exon
        tests.base_test.extend_exon_end(exon, 8, chr_seq)

        gene = self._chromosomes[12].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')
        gene.add_isoform_for_exon_numbers([0, 2, 3, 4], '1')
        # Add 8nt to end of exon 4
        chr_seq = self._chromosomes[12].get_sequence_and_set_coords()
        exon = gene.isoforms[-1].exons[3].copy()
        gene.isoforms[-1].exons[3] = exon
        tests.base_test.extend_exon_end(exon, 8, chr_seq)

        gene = self._chromosomes[13].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')
        gene.add_isoform_for_exon_numbers([0, 2, 3, 4], '1')
        # Add 8nt to end of exon 4
        chr_seq = self._chromosomes[13].get_sequence_and_set_coords()
        exon = gene.isoforms[-1].exons[3].copy()
        gene.isoforms[-1].exons[3] = exon
        tests.base_test.extend_exon_end(exon, 8, chr_seq)

        gene = self._chromosomes[14].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')
        gene.add_isoform_for_exon_numbers([0, 2, 3, 4], '1')
        # Add 8nt to end of exon 4
        chr_seq = self._chromosomes[14].get_sequence_and_set_coords()
        exon = gene.isoforms[-1].exons[3].copy()
        gene.isoforms[-1].exons[3] = exon
        tests.base_test.extend_exon_end(exon, 8, chr_seq)

        gene = self._chromosomes[15].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')
        gene.add_isoform_for_exon_numbers([0, 2, 3, 4], '1')
        # Add 8nt to end of exon 4
        chr_seq = self._chromosomes[15].get_sequence_and_set_coords()
        exon = gene.isoforms[-1].exons[3].copy()
        gene.isoforms[-1].exons[3] = exon
        tests.base_test.extend_exon_end(exon, 8, chr_seq)

        gene = self._chromosomes[16].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')

        gene = self._chromosomes[17].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')

        gene = self._chromosomes[18].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')

        gene = self._chromosomes[19].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')

        gene = self._chromosomes[20].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')

        gene = self._chromosomes[21].genes[0]
        gene.add_isoform_for_exon_numbers([1], '0')
        gene.add_isoform_for_exon_numbers([1], '1')
        # Remove 20nt at start and add 20nt at end
        chr_seq = self._chromosomes[21].get_sequence_and_set_coords()
        exon = gene.isoforms[-1].exons[0].copy()
        gene.isoforms[-1].exons[0] = exon
        tests.base_test.trim_exon_start(exon, 20)
        tests.base_test.extend_exon_end(exon, 20, chr_seq)

        gene = self._chromosomes[22].genes[0]
        gene.add_isoform_for_exon_numbers([0, 2, 4, 6], '0')
        gene.add_isoform_for_exon_numbers([1, 3, 5], '1')

        gene = self._chromosomes[23].genes[0]
        gene.add_isoform_for_exon_numbers([1, 2, 3, 4, 5], '0')
        gene.add_isoform_for_exon_numbers([1, 2, 3, 4, 5], '1')
        exon = gene.isoforms[-1].exons[0].copy()
        gene.isoforms[-1].exons[0] = exon
        tests.base_test.trim_exon_start(exon, 40)

        gene = self._chromosomes[24].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')
        gene.add_isoform_for_exon_numbers([0, 1, 3, 4, 5], '1')

        tests.base_test.write_gtf_from_chroms(gtf_path, self._chromosomes)

    def _create_test_sam(self, sam_path):
        alignments = list()
        self._expected_compat_isoforms = list()
        read_id_offset = 0
        self._add_chr1_alignments(alignments)
        self._add_chr2_alignments(alignments)
        self._add_chr3_alignments(alignments)
        self._add_chr4_alignments(alignments)
        self._add_chr5_alignments(alignments)
        self._add_chr6_alignments(alignments)
        self._add_chr7_alignments(alignments)
        self._add_chr8_alignments(alignments)
        self._add_chr9_alignments(alignments)
        self._add_chr10_alignments(alignments)
        self._add_chr11_alignments(alignments)
        self._add_chr12_alignments(alignments)
        self._add_chr13_alignments(alignments)
        self._add_chr14_alignments(alignments)
        self._add_chr15_alignments(alignments)
        self._add_chr16_alignments(alignments)
        # Add minimal compat isoform info for alignments on chrs 1->16
        compat = {'sample': self._test_name}
        self._expected_compat_isoforms.extend([compat] * len(alignments))
        self._add_chr17_alignments(alignments)
        self._add_chr18_alignments(alignments)
        self._add_chr19_alignments(alignments)
        self._add_chr20_alignments(alignments)
        self._add_chr21_alignments(alignments)
        self._add_chr22_alignments(alignments)
        self._add_chr23_alignments(alignments)
        self._add_chr24_alignments(alignments)
        # Add minimal compat isoform info for alignments on chrs 18->24
        compat = {'sample': self._test_name}
        num_without_compat_info = (len(alignments) -
                                   len(self._expected_compat_isoforms))
        self._expected_compat_isoforms.extend([compat] *
                                              num_without_compat_info)
        self._add_chr25_alignments(alignments)

        read_id_offset += tests.base_test.write_sam_from_alignments(
            sam_path, self._chromosomes, alignments, read_id_offset)
        temp_file = '{}.tmp.sam'.format(sam_path)
        sort_log = os.path.join(self._log_dir, 'sort.log')
        tests.base_test.sort_sam_with_temp_file(sam_path, temp_file, sort_log)

    def _add_chr1_alignments(self, alignments):
        chrom = self._chromosomes[0]
        chrom_seq = chrom.get_sequence_and_set_coords()
        exon_0_length = len(chrom.genes[0].exons[0].sequence)
        exon_1_length = len(chrom.genes[0].exons[1].sequence)
        exon_1_15_percent = int(exon_1_length * 0.15)
        exon_1_30_percent = int(exon_1_length * 0.3)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        # 25nt in exon 0. no insertion or deletion
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        # only 25nt into exon 0
        tests.base_test.trim_alignment_start(alignment, exon_0_length - 25)
        # clipping is ignored
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        # set sequence for "real" alignment, then misalign
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        # drop 1st N operation
        junc_op = tests.base_test.remove_junction_operation(alignment,
                                                            junction_i=0)
        alignment.start += junc_op.num
        tests.base_test.append_copies(alignments, alignment, 2)

        # 35nt in exon 0. no insertion or deletion
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.trim_alignment_start(alignment, exon_0_length - 35)
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        junc_op = tests.base_test.remove_junction_operation(alignment,
                                                            junction_i=0)
        alignment.start += junc_op.num
        tests.base_test.append_copies(alignments, alignment, 2)

        # 25nt in exon 0. 15% insertion in exon 1
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.trim_alignment_start(alignment, exon_0_length - 25)
        tests.base_test.add_insertion_in_exon(
            alignment,
            exon_i=1,
            offset_from_start=25,
            insertion_length=exon_1_15_percent)
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        junc_op = tests.base_test.remove_junction_operation(alignment,
                                                            junction_i=0)
        alignment.start += junc_op.num
        tests.base_test.append_copies(alignments, alignment, 2)

        # 25nt in exon 0. 30% insertion in exon 1
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.trim_alignment_start(alignment, exon_0_length - 25)
        tests.base_test.add_insertion_in_exon(
            alignment,
            exon_i=1,
            offset_from_start=25,
            insertion_length=exon_1_30_percent)
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        junc_op = tests.base_test.remove_junction_operation(alignment,
                                                            junction_i=0)
        alignment.start += junc_op.num
        tests.base_test.append_copies(alignments, alignment, 2)

        # 25nt in exon 0. 15% deletion in exon 1
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.trim_alignment_start(alignment, exon_0_length - 25)
        tests.base_test.add_deletion_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=25,
                                             deletion_length=exon_1_15_percent)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        junc_op = tests.base_test.remove_junction_operation(alignment,
                                                            junction_i=0)
        alignment.start += junc_op.num
        tests.base_test.append_copies(alignments, alignment, 2)

        # 25nt in exon 0. 30% deletion in exon 1
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.trim_alignment_start(alignment, exon_0_length - 25)
        tests.base_test.add_deletion_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=25,
                                             deletion_length=exon_1_30_percent)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        junc_op = tests.base_test.remove_junction_operation(alignment,
                                                            junction_i=0)
        alignment.start += junc_op.num
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr2_alignments(self, alignments):
        chrom = self._chromosomes[1]
        chrom_seq = chrom.get_sequence_and_set_coords()
        exon_1_length = len(chrom.genes[0].exons[1].sequence)
        exon_2_length = len(chrom.genes[0].exons[2].sequence)
        exon_1_15_percent = int(exon_1_length * 0.15)
        exon_1_30_percent = int(exon_1_length * 0.3)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3])
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        # 25nt in exon 2. no insertion or deletion
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        # only 25nt into exon 2
        tests.base_test.trim_alignment_end(alignment, exon_2_length - 25)
        # clipping is ignored
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        # set sequence for "real" alignment, then misalign
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        # drop last N operation
        tests.base_test.remove_junction_operation(alignment, junction_i=1)
        tests.base_test.append_copies(alignments, alignment, 2)

        # 35nt in exon 2. no insertion or deletion
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.trim_alignment_end(alignment, exon_2_length - 35)
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.remove_junction_operation(alignment, junction_i=1)
        tests.base_test.append_copies(alignments, alignment, 2)

        # 25nt in exon 2. 15% insertion in exon 1
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.trim_alignment_end(alignment, exon_2_length - 25)
        tests.base_test.add_insertion_in_exon(
            alignment,
            exon_i=1,
            offset_from_start=25,
            insertion_length=exon_1_15_percent)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.remove_junction_operation(alignment, junction_i=1)
        tests.base_test.append_copies(alignments, alignment, 2)

        # 25nt in exon 2. 30% insertion in exon 1
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.trim_alignment_end(alignment, exon_2_length - 25)
        tests.base_test.add_insertion_in_exon(
            alignment,
            exon_i=1,
            offset_from_start=25,
            insertion_length=exon_1_30_percent)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.remove_junction_operation(alignment, junction_i=1)
        tests.base_test.append_copies(alignments, alignment, 2)

        # 25nt in exon 2. 15% deletion in exon 1
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.trim_alignment_end(alignment, exon_2_length - 25)
        tests.base_test.add_deletion_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=25,
                                             deletion_length=exon_1_15_percent)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.remove_junction_operation(alignment, junction_i=1)
        tests.base_test.append_copies(alignments, alignment, 2)

        # 25nt in exon 2. 30% deletion in exon 1
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2])
        tests.base_test.trim_alignment_end(alignment, exon_2_length - 25)
        tests.base_test.add_deletion_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=25,
                                             deletion_length=exon_1_30_percent)
        alignment.set_sequence_from_cigar_and_chr_seq(chrom_seq)
        tests.base_test.remove_junction_operation(alignment, junction_i=1)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr3_alignments(self, alignments):
        chrom = self._chromosomes[2]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.trim_alignment_start(alignment, 5)
        tests.base_test.extend_alignment_end(alignment, 5)
        tests.base_test.append_copies(alignments, alignment, 1)
        tests.base_test.trim_alignment_start(alignment, 25)
        tests.base_test.extend_alignment_end(alignment, 25)
        tests.base_test.append_copies(alignments, alignment, 3)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.trim_alignment_start(alignment, 40)
        tests.base_test.extend_alignment_end(alignment, 5)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.trim_alignment_start(alignment, 40)
        tests.base_test.extend_alignment_end(alignment, 5)
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 5)
        tests.base_test.trim_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 5)
        tests.base_test.trim_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr4_alignments(self, alignments):
        chrom = self._chromosomes[3]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 40)
        tests.base_test.trim_alignment_end(alignment, 5)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 5)
        tests.base_test.trim_alignment_end(alignment, 5)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 10)
        tests.base_test.trim_alignment_end(alignment, 5)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 5)
        tests.base_test.trim_alignment_end(alignment, 5)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr5_alignments(self, alignments):
        chrom = self._chromosomes[4]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.trim_alignment_start(alignment, 5)
        tests.base_test.extend_alignment_end(alignment, 10)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.trim_alignment_start(alignment, 5)
        tests.base_test.extend_alignment_end(alignment, 5)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.trim_alignment_start(alignment, 5)
        tests.base_test.extend_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.trim_alignment_start(alignment, 5)
        tests.base_test.extend_alignment_end(alignment, 5)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr6_alignments(self, alignments):
        chrom = self._chromosomes[5]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 20)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 20)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr7_alignments(self, alignments):
        chrom = self._chromosomes[6]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 40)
        tests.base_test.extend_alignment_end(alignment, 8)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 5)
        tests.base_test.extend_alignment_end(alignment, 5)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr8_alignments(self, alignments):
        chrom = self._chromosomes[7]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 3, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 20)
        tests.base_test.trim_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 20)
        tests.base_test.trim_alignment_end(alignment, 5)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr9_alignments(self, alignments):
        chrom = self._chromosomes[8]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 3, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 40)
        tests.base_test.trim_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 5)
        tests.base_test.trim_alignment_end(alignment, 5)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr10_alignments(self, alignments):
        chrom = self._chromosomes[9]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 3, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [3, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [3, 4])
        tests.base_test.extend_alignment_start(alignment, 40)
        tests.base_test.extend_alignment_end(alignment, 40)
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [3, 4])
        tests.base_test.trim_alignment_start(alignment, 40)
        tests.base_test.trim_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr11_alignments(self, alignments):
        chrom = self._chromosomes[10]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 3, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [3, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [3, 4])
        tests.base_test.extend_alignment_start(alignment, 40)
        tests.base_test.extend_alignment_end(alignment, 40)
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [3, 4])
        tests.base_test.extend_alignment_start(alignment, 40)
        tests.base_test.trim_alignment_end(alignment, 40)
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr12_alignments(self, alignments):
        chrom = self._chromosomes[11]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3, 5])
        tests.base_test.append_copies(alignments, alignment, 6)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3, 4])
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 36)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr13_alignments(self, alignments):
        chrom = self._chromosomes[12]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3, 5])
        tests.base_test.append_copies(alignments, alignment, 6)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        tests.base_test.extend_alignment_start(alignment, 10)
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3, 4])
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 36)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr14_alignments(self, alignments):
        chrom = self._chromosomes[13]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2, 3, 4])
        tests.base_test.extend_alignment_end(alignment, 8)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3, 5])
        tests.base_test.append_copies(alignments, alignment, 6)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3, 4])
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 36)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr15_alignments(self, alignments):
        chrom = self._chromosomes[14]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3, 5])
        tests.base_test.append_copies(alignments, alignment, 6)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3, 4])
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 36)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr16_alignments(self, alignments):
        chrom = self._chromosomes[15]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 2, 3, 4])
        tests.base_test.extend_alignment_end(alignment, 8)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3, 5])
        tests.base_test.append_copies(alignments, alignment, 6)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3, 4])
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 36)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr17_alignments(self, alignments):
        chrom = self._chromosomes[16]
        exon_2_length = len(chrom.genes[0].exons[2].sequence)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        compat = {
            'class': 'FSM',
            'isoforms': [chrom.genes[0].isoforms[0].id],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        tests.base_test.trim_alignment_start(alignment, 34)
        tests.base_test.extend_alignment_end(alignment, 10)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=10)
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        compat = {
            'class': 'NNC',
            'isoforms':
            [tests.base_test.get_espresso_novel_id(chrom.name, 16, 0)],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        mismatch_len = 20
        offset = exon_2_length - mismatch_len
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=offset,
                                             mismatch_len=mismatch_len)
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=2,
                                             offset_from_start=0,
                                             mismatch_len=mismatch_len)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=1)
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        compat = {
            'class':
            'NCD',
            'isoforms': [
                chrom.genes[0].isoforms[0].id,
                tests.base_test.get_espresso_novel_id(chrom.name, 16, 0)
            ],
            'sample':
            self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 6)
        tests.base_test.append_copies(alignments, alignment, 6)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3])
        mismatch_len = 20
        offset = exon_2_length - mismatch_len
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=0,
                                             offset_from_start=offset,
                                             mismatch_len=mismatch_len)
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=0,
                                             mismatch_len=mismatch_len)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=0,
                                                 offset=1)
        compat = {'class': 'NCD', 'isoforms': [], 'sample': self._test_name}
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr18_alignments(self, alignments):
        chrom = self._chromosomes[17]
        exon_2_length = len(chrom.genes[0].exons[2].sequence)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        tests.base_test.trim_alignment_start(alignment, 34)
        tests.base_test.extend_alignment_end(alignment, 10)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=10)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        mismatch_len = 20
        offset = exon_2_length - mismatch_len
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=offset,
                                             mismatch_len=mismatch_len)
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=2,
                                             offset_from_start=0,
                                             mismatch_len=mismatch_len)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=1)
        tests.base_test.extend_alignment_end(alignment, 30)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 10)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 10)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr19_alignments(self, alignments):
        chrom = self._chromosomes[18]
        exon_2_length = len(chrom.genes[0].exons[2].sequence)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        tests.base_test.trim_alignment_start(alignment, 34)
        tests.base_test.extend_alignment_end(alignment, 10)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=10)
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        mismatch_len = 20
        offset = exon_2_length - mismatch_len
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=offset,
                                             mismatch_len=mismatch_len)
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=2,
                                             offset_from_start=0,
                                             mismatch_len=mismatch_len)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=1)
        tests.base_test.extend_alignment_start(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 2)
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr20_alignments(self, alignments):
        chrom = self._chromosomes[19]
        exon_1_length = len(chrom.genes[0].exons[1].sequence)
        exon_2_length = len(chrom.genes[0].exons[2].sequence)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 5])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3])
        mismatch_len = 20
        offset = exon_1_length - mismatch_len
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=0,
                                             offset_from_start=offset,
                                             mismatch_len=mismatch_len)
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=0,
                                             mismatch_len=mismatch_len)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=0,
                                                 offset=1)
        tests.base_test.extend_alignment_start(alignment, 4)
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3])
        mismatch_len = 20
        offset = exon_2_length - mismatch_len
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=offset,
                                             mismatch_len=mismatch_len)
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=2,
                                             offset_from_start=0,
                                             mismatch_len=mismatch_len)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=1)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr21_alignments(self, alignments):
        chrom = self._chromosomes[20]
        exon_1_length = len(chrom.genes[0].exons[1].sequence)
        exon_2_length = len(chrom.genes[0].exons[2].sequence)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 5])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3])
        mismatch_len = 20
        offset = exon_1_length - mismatch_len
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=0,
                                             offset_from_start=offset,
                                             mismatch_len=mismatch_len)
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=0,
                                             mismatch_len=mismatch_len)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=0,
                                                 offset=1)
        tests.base_test.extend_alignment_start(alignment, 8)
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3])
        mismatch_len = 20
        offset = exon_2_length - mismatch_len
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=offset,
                                             mismatch_len=mismatch_len)
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=2,
                                             offset_from_start=0,
                                             mismatch_len=mismatch_len)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=1)
        tests.base_test.extend_alignment_end(alignment, 8)
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr22_alignments(self, alignments):
        chrom = self._chromosomes[21]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1])
        tests.base_test.append_copies(alignments, alignment, 1)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1])
        tests.base_test.extend_alignment_start(alignment, 20)
        tests.base_test.trim_alignment_end(alignment, 20)
        tests.base_test.add_soft_clipping_at_alignment_start(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 3)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1])
        tests.base_test.trim_alignment_start(alignment, 40)
        tests.base_test.extend_alignment_end(alignment, 40)
        tests.base_test.add_soft_clipping_at_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 5)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1])
        tests.base_test.extend_alignment_start(alignment, 40)
        tests.base_test.trim_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 7)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1])
        tests.base_test.trim_alignment_start(alignment, 60)
        tests.base_test.extend_alignment_end(alignment, 60)
        tests.base_test.append_copies(alignments, alignment, 7)

    def _add_chr23_alignments(self, alignments):
        chrom = self._chromosomes[22]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 3, 4])
        tests.base_test.trim_alignment_end(alignment, 20)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [3, 4])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 10)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 10)
        tests.base_test.append_copies(alignments, alignment, 2)

    def _add_chr24_alignments(self, alignments):
        chrom = self._chromosomes[23]
        exon_2_length = len(chrom.genes[0].exons[2].sequence)
        exon_3_length = len(chrom.genes[0].exons[3].sequence)
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4, 5])
        tests.base_test.append_copies(alignments, alignment, 1)
        tests.base_test.trim_alignment_start(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 6)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        tests.base_test.append_copies(alignments, alignment, 3)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [2, 3, 4, 5])
        tests.base_test.extend_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 5)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4, 5])
        mismatch_len = 20
        offset = exon_3_length - mismatch_len
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=2,
                                             offset_from_start=offset,
                                             mismatch_len=mismatch_len)
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=3,
                                             offset_from_start=0,
                                             mismatch_len=mismatch_len)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=2,
                                                 offset=1)
        tests.base_test.append_copies(alignments, alignment, 7)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4, 5])
        mismatch_len = 20
        offset = exon_2_length - mismatch_len
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=1,
                                             offset_from_start=offset,
                                             mismatch_len=mismatch_len)
        tests.base_test.add_mismatch_in_exon(alignment,
                                             exon_i=2,
                                             offset_from_start=0,
                                             mismatch_len=mismatch_len)
        tests.base_test.adjust_start_of_junction(alignment,
                                                 junction_i=1,
                                                 offset=1)
        tests.base_test.trim_alignment_start(alignment, 10)
        tests.base_test.extend_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 9)

    def _add_chr25_alignments(self, alignments):
        chrom = self._chromosomes[24]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 3, 4, 5])
        compat = {
            'class': 'FSM',
            'isoforms': [chrom.genes[0].isoforms[1].id],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 2)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        tests.base_test.trim_alignment_start(alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.trim_alignment_start(alignment, 10)
        tests.base_test.extend_alignment_end(alignment, 2)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.trim_alignment_start(alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 36)
        tests.base_test.append_copies(alignments, alignment, 2)
        compat = {
            'class': 'NIC',
            'isoforms':
            [tests.base_test.get_espresso_novel_id(chrom.name, 24, 0)],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 4)
        compat = {'sample': self._test_name}
        self._expected_compat_isoforms.extend([compat] * 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3])
        tests.base_test.trim_alignment_start(alignment, 12)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 4)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 4)
        tests.base_test.extend_alignment_end(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        compat = {
            'class': 'ISM',
            'isoforms':
            [tests.base_test.get_espresso_novel_id(chrom.name, 24, 0)],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 4)
        compat = {'sample': self._test_name}
        self._expected_compat_isoforms.extend([compat] * 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 3])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 4)
        tests.base_test.append_copies(alignments, alignment, 2)
        compat = {
            'class': 'ISM',
            'isoforms': [chrom.genes[0].isoforms[1].id],
            'sample': self._test_name
        }
        self._expected_compat_isoforms.extend([compat] * 4)
        compat = {'sample': self._test_name}
        self._expected_compat_isoforms.extend([compat] * 2)

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
            '--inserted_cont_cutoff', '100', '--cont_del_max', '100'
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

        # chr22 is a single exon isoform
        sj_lists_without_22 = sj_simplified_lists[:]
        chr_22_sj_simplified = sj_lists_without_22.pop(21)

        all_sjs = os.path.join(self._work_dir, 'SJ_group_all.fa')
        sam_list = os.path.join(self._work_dir, '0', 'sam.list3')
        sj_list = os.path.join(self._work_dir, '0', 'sj.list')

        self._assert_exists_and_non_empty(
            [self._samples_updated, all_sjs, sam_list, sj_list] +
            sj_lists_without_22)
        self._assert_does_not_exist_or_is_empty([chr_22_sj_simplified])

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

            if chrom.name == 'chr01':
                self.assertEqual(len(parsed_read_final), 16)
                for read_i in range(16):
                    read_id = 'read-{}'.format(read_i)
                    self.assertIn(read_id, parsed_read_final)
                    read_details = parsed_read_final[read_id]
                    sjs = read_details['SJs']
                    added_first = read_details.get('added_first')
                    added_last = read_details.get('added_last')
                    if added_first:
                        sjs = [added_first] + sjs
                    if added_last:
                        sjs = sjs + [added_last]

                    num_sjs = len(sjs)
                    num_added = len([
                        sj for sj in sjs if sj['corrected_status'] == 'added'
                    ])
                    if read_i in [0, 1]:
                        self.assertEqual(num_sjs, 2, read_id)
                        self.assertEqual(num_added, 0, read_id)
                    elif read_i in [2, 3]:
                        self.assertEqual(num_sjs, 3, read_id)
                        self.assertEqual(num_added, 0, read_id)
                    elif read_i in [4, 5, 8, 9, 12, 13]:
                        self.assertEqual(num_sjs, 2, read_id)
                        self.assertEqual(num_added, 1, read_id)
                    else:
                        self.assertEqual(num_sjs, 1, read_id)
                        self.assertEqual(num_added, 0, read_id)
            elif chrom.name == 'chr02':
                self.assertEqual(len(parsed_read_final), 16)
                for read_i in range(16):
                    # account for 16 reads from chr1 in read_id
                    read_id = 'read-{}'.format(read_i + 16)
                    self.assertIn(read_id, parsed_read_final)
                    read_details = parsed_read_final[read_id]
                    sjs = read_details['SJs']
                    added_first = read_details.get('added_first')
                    added_last = read_details.get('added_last')
                    if added_first:
                        sjs = [added_first] + sjs
                    if added_last:
                        sjs = sjs + [added_last]

                    num_sjs = len(sjs)
                    num_added = len([
                        sj for sj in sjs if sj['corrected_status'] == 'added'
                    ])
                    if read_i in [0, 1]:
                        self.assertEqual(num_sjs, 2, read_id)
                        self.assertEqual(num_added, 0, read_id)
                    elif read_i in [2, 3]:
                        self.assertEqual(num_sjs, 3, read_id)
                        self.assertEqual(num_added, 0, read_id)
                    elif read_i in [4, 5, 8, 9, 12, 13]:
                        self.assertEqual(num_sjs, 2, read_id)
                        self.assertEqual(num_added, 1, read_id)
                    else:
                        self.assertEqual(num_sjs, 1, read_id)
                        self.assertEqual(num_added, 0, read_id)
            elif chrom.name in ['chr17', 'chr18', 'chr19', 'chr20', 'chr21']:
                num_reads_with_fail = 0
                for read_details in parsed_read_final.values():
                    for sj in read_details['SJs']:
                        if sj['corrected_status'] == 'fail':
                            num_reads_with_fail += 1
                            break

                if chrom.name == 'chr17':
                    self.assertEqual(len(parsed_read_final), 12)
                    self.assertEqual(num_reads_with_fail, 8)
                if chrom.name == 'chr18':
                    self.assertEqual(len(parsed_read_final), 10)
                    self.assertEqual(num_reads_with_fail, 6)
                if chrom.name == 'chr19':
                    self.assertEqual(len(parsed_read_final), 10)
                    self.assertEqual(num_reads_with_fail, 6)
                if chrom.name == 'chr20':
                    self.assertEqual(len(parsed_read_final), 8)
                    self.assertEqual(num_reads_with_fail, 4)
                if chrom.name == 'chr21':
                    self.assertEqual(len(parsed_read_final), 8)
                    self.assertEqual(num_reads_with_fail, 4)
            else:
                self.assertGreater(len(parsed_read_final), 0)

    def _run_espresso_q(self):
        compat_basename = '{}_compatible_isoform.tsv'.format(
            self._out_file_prefix)
        self._compat_isoforms = os.path.join(self._work_dir, compat_basename)
        command = [
            'perl', self._espresso_q, '-L', self._samples_updated, '-A',
            self._gtf, '-V', self._compat_isoforms
        ]

        q_log = os.path.join(self._log_dir, 'espresso_q.log')
        tests.base_test.run_command_with_log(command, q_log)

    def _check_q_output(self):
        gtf_basename = '{}_updated.gtf'.format(self._out_file_prefix)
        self._updated_gtf = os.path.join(self._work_dir, gtf_basename)
        abundance_basename = '{}_abundance.esp'.format(self._out_file_prefix)
        self._abundance = os.path.join(self._work_dir, abundance_basename)
        self._assert_exists_and_non_empty(
            [self._updated_gtf, self._abundance, self._compat_isoforms])
        detected_isoforms = self._parse_updated_gtf(self._updated_gtf)
        abundance_rows = self._parse_abundance(self._abundance,
                                               self._sample_names)
        compat_isoforms = self._parse_compat_isoforms(self._compat_isoforms)
        self._check_compat_isoforms(self._expected_compat_isoforms,
                                    compat_isoforms)

        expected_isoforms = list()
        # chr1: [0, 1, 2]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['abundance'] = {self._test_name: 8}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[0].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[0].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 1, 2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[0].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[0].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)

        # chr2: [0, 1, 2]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[1].name
        expected_isoform['gene'] = self._chromosomes[1].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['abundance'] = {self._test_name: 8}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[1].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[1].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 1, 2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[1].name
        expected_isoform['gene'] = self._chromosomes[1].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[1].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[1].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)

        # chr3: [1, 2, 3]
        chrom_sequence = self._chromosomes[2].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[2].name
        expected_isoform['gene'] = self._chromosomes[2].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1].copy()]
        tests.base_test.trim_exon_start(expected_isoform['exons'][-1], 5)
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[3].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 30,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 18}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[2].name, 2,
                                                  0))
        expected_isoforms.append(expected_isoform)

        # chr4: [1, 2, 3]
        chrom_sequence = self._chromosomes[3].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[3].name
        expected_isoform['gene'] = self._chromosomes[3].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[3].name, 3,
                                                  2))
        expected_isoforms.append(expected_isoform)
        # [1, 2]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[3].name
        expected_isoform['gene'] = self._chromosomes[3].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1].copy()]
        tests.base_test.extend_exon_start(expected_isoform['exons'][-1], 40,
                                          chrom_sequence)
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[2].copy())
        tests.base_test.trim_exon_end(expected_isoform['exons'][-1], 5)
        expected_isoform['abundance'] = {self._test_name: 6}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[3].name, 3,
                                                  0))
        expected_isoforms.append(expected_isoform)
        # [2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[3].name
        expected_isoform['gene'] = self._chromosomes[3].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[2].copy()]
        tests.base_test.extend_exon_start(expected_isoform['exons'][-1], 10,
                                          chrom_sequence)
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[3].copy())
        tests.base_test.trim_exon_end(expected_isoform['exons'][-1], 5)
        expected_isoform['abundance'] = {self._test_name: 6}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[3].name, 3,
                                                  1))
        expected_isoforms.append(expected_isoform)

        # chr5: [1, 2, 3]
        chrom_sequence = self._chromosomes[4].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[4].name
        expected_isoform['gene'] = self._chromosomes[4].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[4].name, 4,
                                                  2))
        expected_isoforms.append(expected_isoform)
        # [1, 2]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[4].name
        expected_isoform['gene'] = self._chromosomes[4].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1].copy()]
        tests.base_test.trim_exon_start(expected_isoform['exons'][-1], 5)
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[2].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 10,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 6}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[4].name, 4,
                                                  0))
        expected_isoforms.append(expected_isoform)
        # [2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[4].name
        expected_isoform['gene'] = self._chromosomes[4].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[2].copy()]
        tests.base_test.trim_exon_start(expected_isoform['exons'][-1], 5)
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[3].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 40,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 6}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[4].name, 4,
                                                  1))
        expected_isoforms.append(expected_isoform)

        # chr6: [0, 1, 4]
        chrom_sequence = self._chromosomes[5].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[5].name
        expected_isoform['gene'] = self._chromosomes[5].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[5].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[5].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 2, 4]
        chrom_sequence = self._chromosomes[5].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[5].name
        expected_isoform['gene'] = self._chromosomes[5].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[5].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[5].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)
        # [1, 2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[5].name
        expected_isoform['gene'] = self._chromosomes[5].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 8}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[5].name, 5,
                                                  0))
        expected_isoforms.append(expected_isoform)

        # chr7: [0, 1, 4]
        chrom_sequence = self._chromosomes[6].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[6].name
        expected_isoform['gene'] = self._chromosomes[6].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[6].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[6].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 2, 4]
        chrom_sequence = self._chromosomes[6].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[6].name
        expected_isoform['gene'] = self._chromosomes[6].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[6].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[6].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)
        # [1, 2, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[6].name
        expected_isoform['gene'] = self._chromosomes[6].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[6].name, 6,
                                                  0))
        expected_isoforms.append(expected_isoform)
        # [1, 2]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[6].name
        expected_isoform['gene'] = self._chromosomes[6].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1].copy()]
        tests.base_test.extend_exon_start(expected_isoform['exons'][-1], 40,
                                          chrom_sequence)
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[2].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 8,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 6}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[6].name, 6,
                                                  1))
        expected_isoforms.append(expected_isoform)

        # chr8: [0, 1, 2, 3, 4, 5]
        chrom_sequence = self._chromosomes[7].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[7].name
        expected_isoform['gene'] = self._chromosomes[7].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[7].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[7].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [1, 3, 4]
        chrom_sequence = self._chromosomes[7].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[7].name
        expected_isoform['gene'] = self._chromosomes[7].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['abundance'] = {self._test_name: 8}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[7].name, 7,
                                                  0))
        expected_isoforms.append(expected_isoform)

        # chr9: [0, 1, 2, 3, 4, 5]
        chrom_sequence = self._chromosomes[8].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[8].name
        expected_isoform['gene'] = self._chromosomes[8].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[8].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[8].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [1, 3, 4]
        chrom_sequence = self._chromosomes[8].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[8].name
        expected_isoform['gene'] = self._chromosomes[8].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[8].name, 8,
                                                  1))
        expected_isoforms.append(expected_isoform)
        # [1, 3]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[8].name
        expected_isoform['gene'] = self._chromosomes[8].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1].copy()]
        tests.base_test.extend_exon_start(expected_isoform['exons'][-1], 40,
                                          chrom_sequence)
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[3].copy())
        tests.base_test.trim_exon_end(expected_isoform['exons'][-1], 40)
        expected_isoform['abundance'] = {self._test_name: 6}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[8].name, 8,
                                                  0))
        expected_isoforms.append(expected_isoform)

        # chr10: [0, 1, 2, 3, 4, 5]
        chrom_sequence = self._chromosomes[9].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[9].name
        expected_isoform['gene'] = self._chromosomes[9].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 4}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[9].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[9].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [1, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[9].name
        expected_isoform['gene'] = self._chromosomes[9].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['abundance'] = {self._test_name: 4}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[9].name, 9,
                                                  0))
        expected_isoforms.append(expected_isoform)

        # chr11: [0, 1, 2, 3, 4, 5]
        chrom_sequence = self._chromosomes[10].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[10].name
        expected_isoform['gene'] = self._chromosomes[10].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[10].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[10].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [1, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[10].name
        expected_isoform['gene'] = self._chromosomes[10].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[10].name,
                                                  10, 0))
        expected_isoforms.append(expected_isoform)

        # chr12: [2, 3, 5]
        chrom_sequence = self._chromosomes[11].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[11].name
        expected_isoform['gene'] = self._chromosomes[11].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[2]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 7.33}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[11].name,
                                                  11, 0))
        expected_isoforms.append(expected_isoform)
        # [1, 2, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[11].name
        expected_isoform['gene'] = self._chromosomes[11].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[4].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 4,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 14.67}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[11].name,
                                                  11, 1))
        expected_isoforms.append(expected_isoform)

        # chr13: [2, 3, 5]
        chrom_sequence = self._chromosomes[12].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[12].name
        expected_isoform['gene'] = self._chromosomes[12].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[2]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 10}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[12].name,
                                                  12, 0))
        expected_isoforms.append(expected_isoform)

        # chr14: [0, 2, 3, 4]
        chrom_sequence = self._chromosomes[13].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[13].name
        expected_isoform['gene'] = self._chromosomes[13].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[4].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 8,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 4.2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[13].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[13].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)
        # [2, 3, 5]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[13].name
        expected_isoform['gene'] = self._chromosomes[13].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[2]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 7.21}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[13].name,
                                                  13, 0))
        expected_isoforms.append(expected_isoform)
        # [1, 2, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[13].name
        expected_isoform['gene'] = self._chromosomes[13].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[4].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 4,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 12.59}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[13].name,
                                                  13, 1))
        expected_isoforms.append(expected_isoform)

        # chr15: [0, 1, 2, 3, 4, 5]
        chrom_sequence = self._chromosomes[14].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[14].name
        expected_isoform['gene'] = self._chromosomes[14].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 9.6}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[14].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[14].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [2, 3, 5]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[14].name
        expected_isoform['gene'] = self._chromosomes[14].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[2]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 7.2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[14].name,
                                                  14, 0))
        expected_isoforms.append(expected_isoform)
        # [2, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[14].name
        expected_isoform['gene'] = self._chromosomes[14].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[2]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[4].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 8,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 7.2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[14].name,
                                                  14, 1))
        expected_isoforms.append(expected_isoform)

        # chr16: [0, 1, 2, 3, 4, 5]
        chrom_sequence = self._chromosomes[15].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[15].name
        expected_isoform['gene'] = self._chromosomes[15].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 9.45}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[15].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[15].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 2, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[15].name
        expected_isoform['gene'] = self._chromosomes[15].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[4].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 8,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 9.45}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[15].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[15].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)
        # [2, 3, 5]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[15].name
        expected_isoform['gene'] = self._chromosomes[15].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[2]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 7.09}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[15].name,
                                                  15, 0))
        expected_isoforms.append(expected_isoform)

        # chr17: [0, 1, 2, 3, 4, 5]
        chrom_sequence = self._chromosomes[16].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[16].name
        expected_isoform['gene'] = self._chromosomes[16].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 5}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[16].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[16].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [1, 2, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[16].name
        expected_isoform['gene'] = self._chromosomes[16].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1].copy()]
        tests.base_test.trim_exon_start(expected_isoform['exons'][-1], 34)
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[2].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 10,
                                        chrom_sequence)
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[4].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 10,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 5}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[16].name,
                                                  16, 0))
        expected_isoforms.append(expected_isoform)

        # chr18: [0, 1, 2, 3, 4, 5]
        chrom_sequence = self._chromosomes[17].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[17].name
        expected_isoform['gene'] = self._chromosomes[17].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[17].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[17].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [1, 2, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[17].name
        expected_isoform['gene'] = self._chromosomes[17].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1].copy()]
        tests.base_test.trim_exon_start(expected_isoform['exons'][-1], 34)
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[2].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 10,
                                        chrom_sequence)
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[4].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 10,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 8}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[17].name,
                                                  17, 0))
        expected_isoforms.append(expected_isoform)

        # chr19: [0, 1, 2, 3, 4, 5]
        chrom_sequence = self._chromosomes[18].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[18].name
        expected_isoform['gene'] = self._chromosomes[18].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 8}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[18].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[18].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [1, 2, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[18].name
        expected_isoform['gene'] = self._chromosomes[18].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1].copy()]
        tests.base_test.trim_exon_start(expected_isoform['exons'][-1], 34)
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[2].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 10,
                                        chrom_sequence)
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[4].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 10,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[18].name,
                                                  18, 0))
        expected_isoforms.append(expected_isoform)

        # chr20: [0, 1, 2, 3, 4, 5]
        chrom_sequence = self._chromosomes[19].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[19].name
        expected_isoform['gene'] = self._chromosomes[19].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 4}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[19].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[19].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 1, 2, 3, 5]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[19].name
        expected_isoform['gene'] = self._chromosomes[19].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 4}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[19].name,
                                                  19, 0))
        expected_isoforms.append(expected_isoform)

        # chr21: [0, 1, 2, 3, 4, 5]
        chrom_sequence = self._chromosomes[20].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[20].name
        expected_isoform['gene'] = self._chromosomes[20].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[20].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[20].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 1, 2, 3, 5]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[20].name
        expected_isoform['gene'] = self._chromosomes[20].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[20].name,
                                                  20, 0))
        expected_isoforms.append(expected_isoform)

        # chr22: [1]
        chrom_sequence = self._chromosomes[21].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[21].name
        expected_isoform['gene'] = self._chromosomes[21].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['abundance'] = {self._test_name: 3.38}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[21].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[21].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [1] with endpoints shifted by 20nt
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[21].name
        expected_isoform['gene'] = self._chromosomes[21].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1].copy()]
        tests.base_test.trim_exon_start(expected_isoform['exons'][-1], 20)
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 20,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 5.62}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[21].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[21].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)

        # chr23: [1, 2, 3, 4]
        chrom_sequence = self._chromosomes[22].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[22].name
        expected_isoform['gene'] = self._chromosomes[22].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['abundance'] = {self._test_name: 5}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[22].name,
                                                  22, 1))
        expected_isoforms.append(expected_isoform)
        # [1, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[22].name
        expected_isoform['gene'] = self._chromosomes[22].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[4].copy())
        tests.base_test.trim_exon_end(expected_isoform['exons'][-1], 20)
        expected_isoform['abundance'] = {self._test_name: 5}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[22].name,
                                                  22, 0))
        expected_isoforms.append(expected_isoform)

        # chr24: [1, 2, 3, 4, 5]
        chrom_sequence = self._chromosomes[23].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[23].name
        expected_isoform['gene'] = self._chromosomes[23].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 7}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[23].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[23].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [1, 2, 3, 4, 5] with exon 1 trimmed by 40nt
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[23].name
        expected_isoform['gene'] = self._chromosomes[23].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1].copy()]
        tests.base_test.trim_exon_start(expected_isoform['exons'][-1], 40)
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 14}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[23].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[23].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)

        # chr25: [0, 1, 3, 4, 5]
        chrom_sequence = self._chromosomes[24].get_sequence_and_set_coords()
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[24].name
        expected_isoform['gene'] = self._chromosomes[24].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 6}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[24].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[24].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)
        # [1, 2, 3, 4] with exon 1 trimmed by 12nt and exon 4 extended by 4
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[24].name
        expected_isoform['gene'] = self._chromosomes[24].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1].copy()]
        tests.base_test.trim_exon_start(expected_isoform['exons'][-1], 12)
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[4].copy())
        tests.base_test.extend_exon_end(expected_isoform['exons'][-1], 4,
                                        chrom_sequence)
        expected_isoform['abundance'] = {self._test_name: 8}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[24].name,
                                                  24, 0))
        expected_isoforms.append(expected_isoform)

        expected_by_id = self._get_transcript_ids_from_detected(
            expected_isoforms, detected_isoforms)
        self.assertEqual(len(abundance_rows), len(expected_isoforms))
        self._check_abundance_row_values(expected_by_id, abundance_rows,
                                         self._sample_names)


class NoExternalBoundaryTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'no_external_boundary'
        self._test_dir = os.path.dirname(__file__)
        self._data_dir = os.path.join(self._test_dir,
                                      '{}_data'.format(self._test_name))
        self._out_dir = os.path.join(self._test_dir,
                                     '{}_out'.format(self._test_name))
        self._log_dir = os.path.join(self._test_dir,
                                     '{}_logs'.format(self._test_name))
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
        self._compat_isoforms = None

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
        chr_id_format = 'chr{:02d}'
        gene_id_format = 'ENSG{:02d}'
        gene_name_format = 'GENE{:02d}'
        # chr1: Always pass external boundary check
        # gene has annotated exons:
        #      [0, 1, 2, 3, 4, 5] (also with terminal exons shortened by 40nt)
        #      [1, 2, 3, 4]
        #      [4]
        # NIC: [1, 3, 4]
        # ISM: [1, 2, 3]
        #      [3, 4]
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        chr_i += 1
        gene_i += 1
        # chr2: like chr1 but with additional internal boundary failures
        gene = self._get_random_gene(num_exons=6)
        gene.id = gene_id_format.format(gene_i)
        gene.name = gene_name_format.format(gene_i)
        chromosome = self._make_chromosome_from_genes([gene])
        chromosome.name = chr_id_format.format(chr_i)
        self._chromosomes.append(chromosome)

        tests.base_test.write_fasta_from_chroms(fasta_path, self._chromosomes)

    def _create_test_gtf(self, gtf_path):
        gene = self._chromosomes[0].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '1')
        exon = gene.isoforms[-1].exons[0].copy()
        gene.isoforms[-1].exons[0] = exon
        tests.base_test.trim_exon_start(exon, 40)
        exon = gene.isoforms[-1].exons[-1].copy()
        gene.isoforms[-1].exons[-1] = exon
        tests.base_test.trim_exon_end(exon, 40)
        gene.add_isoform_for_exon_numbers([1, 2, 3, 4], '2')
        gene.add_isoform_for_exon_numbers([4], '3')

        gene = self._chromosomes[1].genes[0]
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '0')
        gene.add_isoform_for_exon_numbers([0, 1, 2, 3, 4, 5], '1')
        exon = gene.isoforms[-1].exons[0].copy()
        gene.isoforms[-1].exons[0] = exon
        tests.base_test.trim_exon_start(exon, 40)
        exon = gene.isoforms[-1].exons[-1].copy()
        gene.isoforms[-1].exons[-1] = exon
        tests.base_test.trim_exon_end(exon, 40)
        gene.add_isoform_for_exon_numbers([1, 2, 3, 4], '2')
        gene.add_isoform_for_exon_numbers([4], '3')

        tests.base_test.write_gtf_from_chroms(gtf_path, self._chromosomes)

    def _create_test_sam(self, sam_path):
        alignments = list()
        read_id_offset = 0
        chrom = self._chromosomes[0]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.append_copies(alignments, alignment, 1)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.trim_alignment_start(alignment, 40)
        tests.base_test.trim_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 3)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 3, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3])
        tests.base_test.extend_alignment_start(alignment, 30)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 30)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 30)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [3, 4])
        tests.base_test.extend_alignment_end(alignment, 30)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 30)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 30)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [4])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 40)
        tests.base_test.extend_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        chrom = self._chromosomes[1]
        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.append_copies(alignments, alignment, 1)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [0, 1, 2, 3, 4, 5])
        tests.base_test.trim_alignment_start(alignment, 40)
        tests.base_test.trim_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 3)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 3, 4])
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [1, 2, 3])
        tests.base_test.extend_alignment_end(alignment, 10)
        tests.base_test.extend_alignment_start(alignment, 30)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 30)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 30)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [3, 4])
        tests.base_test.extend_alignment_start(alignment, 10)
        tests.base_test.extend_alignment_end(alignment, 30)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 30)
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_end(alignment, 30)
        tests.base_test.append_copies(alignments, alignment, 2)

        alignment = tests.base_test.perfect_alignment_for_gene_and_exon_numbers(
            chrom.name, chrom.genes[0], [4])
        tests.base_test.append_copies(alignments, alignment, 2)
        tests.base_test.extend_alignment_start(alignment, 40)
        tests.base_test.extend_alignment_end(alignment, 40)
        tests.base_test.append_copies(alignments, alignment, 2)

        read_id_offset += tests.base_test.write_sam_from_alignments(
            sam_path, self._chromosomes, alignments, read_id_offset)
        temp_file = '{}.tmp.sam'.format(sam_path)
        sort_log = os.path.join(self._log_dir, 'sort.log')
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
            self._fasta, '-O', self._work_dir, '-A', self._gtf
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

        all_sjs = os.path.join(self._work_dir, 'SJ_group_all.fa')
        sam_list = os.path.join(self._work_dir, '0', 'sam.list3')
        sj_list = os.path.join(self._work_dir, '0', 'sj.list')

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
        compat_basename = '{}_compatible_isoform.tsv'.format(
            self._out_file_prefix)
        self._compat_isoforms = os.path.join(self._work_dir, compat_basename)
        command = [
            'perl', self._espresso_q, '-L', self._samples_updated, '-A',
            self._gtf, '-V', self._compat_isoforms,
            '--allow_longer_terminal_exons'
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
        # chr1: [0, 1, 2, 3, 4, 5]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 1}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[0].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[0].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 1, 2, 3, 4, 5] endpoints shortened by 40
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0].copy()]
        tests.base_test.trim_exon_start(expected_isoform['exons'][-1], 40)
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[5].copy())
        tests.base_test.trim_exon_end(expected_isoform['exons'][-1], 40)
        expected_isoform['abundance'] = {self._test_name: 3}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[0].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[0].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)
        # [1, 2, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['abundance'] = {self._test_name: 12.8}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[0].genes[0].isoforms[2].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[0].genes[0].isoforms[2].name)
        expected_isoforms.append(expected_isoform)
        # [1, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['abundance'] = {self._test_name: 3.2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[0].name, 0,
                                                  0))
        expected_isoforms.append(expected_isoform)
        # [4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[0].name
        expected_isoform['gene'] = self._chromosomes[0].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[4]]
        expected_isoform['abundance'] = {self._test_name: 4}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[0].genes[0].isoforms[3].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[0].genes[0].isoforms[3].name)
        expected_isoforms.append(expected_isoform)

        # chr2: [0, 1, 2, 3, 4, 5]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[1].name
        expected_isoform['gene'] = self._chromosomes[1].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[5])
        expected_isoform['abundance'] = {self._test_name: 1}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[1].genes[0].isoforms[0].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[1].genes[0].isoforms[0].name)
        expected_isoforms.append(expected_isoform)
        # [0, 1, 2, 3, 4, 5] endpoints shortened by 40
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[1].name
        expected_isoform['gene'] = self._chromosomes[1].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[0].copy()]
        tests.base_test.trim_exon_start(expected_isoform['exons'][-1], 40)
        expected_isoform['exons'].append(expected_isoform['gene'].exons[1])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['exons'].append(
            expected_isoform['gene'].exons[5].copy())
        tests.base_test.trim_exon_end(expected_isoform['exons'][-1], 40)
        expected_isoform['abundance'] = {self._test_name: 3}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[1].genes[0].isoforms[1].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[1].genes[0].isoforms[1].name)
        expected_isoforms.append(expected_isoform)
        # [1, 2, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[1].name
        expected_isoform['gene'] = self._chromosomes[1].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[2])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[1].genes[0].isoforms[2].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[1].genes[0].isoforms[2].name)
        expected_isoforms.append(expected_isoform)
        # [1, 3, 4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[1].name
        expected_isoform['gene'] = self._chromosomes[1].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[1]]
        expected_isoform['exons'].append(expected_isoform['gene'].exons[3])
        expected_isoform['exons'].append(expected_isoform['gene'].exons[4])
        expected_isoform['abundance'] = {self._test_name: 2}
        expected_isoform['is_novel'] = True
        expected_isoform['transcript_id'] = (
            tests.base_test.get_espresso_novel_id(self._chromosomes[1].name, 1,
                                                  0))
        expected_isoforms.append(expected_isoform)
        # [4]
        expected_isoform = dict()
        expected_isoform['chr'] = self._chromosomes[1].name
        expected_isoform['gene'] = self._chromosomes[1].genes[0]
        expected_isoform['exons'] = [expected_isoform['gene'].exons[4]]
        expected_isoform['abundance'] = {self._test_name: 4}
        expected_isoform['is_novel'] = False
        expected_isoform['transcript_id'] = (
            self._chromosomes[1].genes[0].isoforms[3].id)
        expected_isoform['transcript_name'] = (
            self._chromosomes[1].genes[0].isoforms[3].name)
        expected_isoforms.append(expected_isoform)

        expected_by_id = self._get_transcript_ids_from_detected(
            expected_isoforms, detected_isoforms)
        self.assertEqual(len(abundance_rows), len(expected_isoforms))
        self._check_abundance_row_values(expected_by_id, abundance_rows,
                                         self._sample_names)


if __name__ == '__main__':
    unittest.main(verbosity=2)
