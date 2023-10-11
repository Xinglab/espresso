import os
import os.path
import unittest

import tests.base_test


def _count_numbered_sub_dirs(dir_path):
    count = 0
    names = os.listdir(dir_path)
    for name in names:
        path = os.path.join(dir_path, name)
        if os.path.isdir(path) and name.isdigit():
            count += 1

    return count


class Cd44ExampleTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'cd44_example'
        self._test_dir = os.path.dirname(__file__)
        self._data_dir = os.path.join(self._test_dir, 'data')
        self._out_dir = os.path.join(self._test_dir, 'out')
        self._log_dir = os.path.join(self._test_dir, 'logs')
        self._s_work_dir = os.path.join(self._out_dir, 's_work_dir')
        self._c_work_dir = os.path.join(self._out_dir, 'c_work_dir')
        self._q_work_dir = os.path.join(self._out_dir, 'q_work_dir')
        self._out_file_prefix = (
            tests.base_test.out_file_prefix_based_on_params())
        self._cd44_fasta = None
        self._cd44_gtf = None
        self._cd44_bed = None
        self._gs689_fastqs = None
        self._gs689_sams = None
        self._gs689_bams = None
        self._pc3e_fastqs = None
        self._pc3e_sams = None
        self._pc3e_bams = None
        self._samples_tsv = None
        self._s_samples_updated = None
        self._q_samples_updated = None
        self._isoform_tsv = None
        self._updated_gtf = None
        self._abundance = None
        self._num_c_splits = 5
        self._target_c_reads = '1000'

    def test(self):
        self._initialize_dirs()
        self._unpack_test_data()
        self._create_bed_from_gtf()
        self._align_fastqs()
        self._create_sorted_bams()
        self._create_samples_tsv()
        self._run_espresso_s()
        self._check_s_output()
        self._split_s_for_c()
        self._run_espresso_c()
        self._check_c_output()
        self._combine_c_for_q()
        self._run_espresso_q()
        self._check_q_output()

    def _initialize_dirs(self):
        directories = [self._data_dir, self._out_dir, self._log_dir]
        tests.base_test.remove_then_create_directories(directories)

    def _unpack_test_data(self):
        test_tar_gz = os.path.join(self._test_data_dir,
                                   'test_data_espresso_cd44.tar.gz')
        command = ['tar', '-xvf', test_tar_gz, '--directory', self._data_dir]
        tar_log = os.path.join(self._log_dir, 'tar.log')
        tests.base_test.run_command_with_log(command, tar_log)
        extracted_base_dir = os.path.join(self._data_dir,
                                          'test_data_espresso_cd44')
        self._cd44_fasta = os.path.join(extracted_base_dir, 'cd44.fasta')
        self._cd44_gtf = os.path.join(extracted_base_dir, 'cd44.gtf')
        self._gs689_fastqs = list()
        self._pc3e_fastqs = list()
        for replicate in ['1', '2', '3']:
            gs_basename = 'GS689_{}_cd44.fastq'.format(replicate)
            pc_basename = 'PC3E_{}_cd44.fastq'.format(replicate)
            self._gs689_fastqs.append(
                os.path.join(extracted_base_dir, gs_basename))
            self._pc3e_fastqs.append(
                os.path.join(extracted_base_dir, pc_basename))

        all_files = ([self._cd44_fasta, self._cd44_gtf] + self._gs689_fastqs +
                     self._pc3e_fastqs)
        self._assert_exists_and_non_empty(all_files)

    def _create_bed_from_gtf(self):
        command = ['paftools.js', 'gff2bed', self._cd44_gtf]
        self._cd44_bed = '{}.bed'.format(self._cd44_gtf)
        bed_log = os.path.join(self._log_dir, 'create_bed_from_gtf.log')
        tests.base_test.run_command_with_output_and_log(
            command, self._cd44_bed, bed_log)
        self._assert_exists_and_non_empty([self._cd44_bed])

    def _align_fastqs(self):
        self._gs689_sams = list()
        for i, gs_fastq in enumerate(self._gs689_fastqs):
            gs_sam = gs_fastq.replace('.fastq', '.sam')
            gs_log = os.path.join(self._log_dir,
                                  'align_fastqs_gs_{}.log'.format(i))
            self._align_fastq(gs_fastq, gs_sam, gs_log)
            self._gs689_sams.append(gs_sam)

        self._pc3e_sams = list()
        for i, pc_fastq in enumerate(self._pc3e_fastqs):
            pc_sam = pc_fastq.replace('.fastq', '.sam')
            pc_log = os.path.join(self._log_dir,
                                  'align_fastqs_pc_{}.log'.format(i))
            self._align_fastq(pc_fastq, pc_sam, pc_log)
            self._pc3e_sams.append(pc_sam)

        self._assert_exists_and_non_empty(self._gs689_sams + self._pc3e_sams)

    def _align_fastq(self, fastq_path, sam_path, log_path):
        preset = 'splice'
        canonical_splice_sites = 'b'
        kmer_size = '14'
        window_size = '4'
        output_secondary_alignments = 'no'
        threads = '1'
        command = [
            'minimap2', '-a', '-x', preset,
            '-u{}'.format(canonical_splice_sites), '-k', kmer_size, '-w',
            window_size, '--junc-bed', self._cd44_bed,
            '--secondary={}'.format(output_secondary_alignments), '-t',
            threads, self._cd44_fasta, fastq_path
        ]
        tests.base_test.run_command_with_output_and_log(
            command, sam_path, log_path)

    def _create_sorted_bams(self):
        self._gs689_bams = list()
        for i, gs_sam in enumerate(self._gs689_sams):
            gs_bam = gs_sam.replace('.sam', '.bam')
            gs_log = os.path.join(self._log_dir,
                                  'create_sorted_bams_gs_{}.log'.format(i))
            self._create_sorted_bam(gs_sam, gs_bam, gs_log)
            self._gs689_bams.append(gs_bam)

        self._pc3e_bams = list()
        for i, pc_sam in enumerate(self._pc3e_sams):
            pc_bam = pc_sam.replace('.sam', '.bam')
            pc_log = os.path.join(self._log_dir,
                                  'create_sorted_bams_pc_{}.log'.format(i))
            self._create_sorted_bam(pc_sam, pc_bam, pc_log)
            self._pc3e_bams.append(pc_bam)

        self._assert_exists_and_non_empty(self._gs689_bams + self._pc3e_bams)

    def _create_sorted_bam(self, sam_path, bam_path, log_path):
        command = ['samtools', 'sort', '-o', bam_path, sam_path]
        tests.base_test.run_command_with_log(command, log_path)

    def _create_samples_tsv(self):
        self._samples_tsv = os.path.join(self._data_dir, 'samples.tsv')
        with open(self._samples_tsv, 'wt') as handle:
            for gs_bam in self._gs689_bams:
                columns = [gs_bam, 'GS689']
                tests.base_test.write_tsv_line(handle, columns)

            for pc_bam in self._pc3e_bams:
                columns = [pc_bam, 'PC3E']
                tests.base_test.write_tsv_line(handle, columns)

    def _run_espresso_s(self):
        command = [
            'perl', self._espresso_s, '-A', self._cd44_gtf, '-L',
            self._samples_tsv, '-F', self._cd44_fasta, '-O', self._s_work_dir
        ]
        s_log = os.path.join(self._log_dir, 'espresso_s.log')
        tests.base_test.run_command_with_log(command, s_log)

    def _check_s_output(self):
        self._s_samples_updated = os.path.join(self._s_work_dir,
                                               'samples.tsv.updated')
        chr11_sjs = os.path.join(self._s_work_dir, 'chr11_SJ_simplified.list')
        all_sjs = os.path.join(self._s_work_dir, 'SJ_group_all.fa')
        sam_lists = list()
        sj_lists = list()
        for i in range(6):
            sam_list = os.path.join(self._s_work_dir, str(i), 'sam.list3')
            sam_lists.append(sam_list)
            sj_list = os.path.join(self._s_work_dir, str(i), 'sj.list')
            sj_lists.append(sj_list)

        self._assert_exists_and_non_empty(
            [self._s_samples_updated, chr11_sjs, all_sjs] + sam_lists +
            sj_lists)

        with open(self._s_samples_updated, 'rt') as samples_handle:
            lines = list()
            for line in samples_handle:
                lines.append(line.strip())

        self.assertEqual(len(lines), 6)

    def _split_s_for_c(self):
        sort_buffer = '2G'
        command = [
            self._py_executable, self._split_s_py, '--orig-work-dir',
            self._s_work_dir, '--new-base-dir', self._c_work_dir,
            '--target-reads-per-c', self._target_c_reads, '--genome-fasta',
            self._cd44_fasta, '--sort-memory-buffer-size', sort_buffer
        ]
        split_log = os.path.join(self._log_dir, 'split_s_for_c.log')
        tests.base_test.run_command_with_log(command, split_log)
        num_splits = _count_numbered_sub_dirs(self._c_work_dir)
        self.assertEqual(num_splits, self._num_c_splits)

    def _run_espresso_c(self):
        threads = '2'
        for work_dir_i in range(self._num_c_splits):
            c_work_dir = os.path.join(self._c_work_dir, str(work_dir_i))
            split_fasta = os.path.join(self._c_work_dir, 'fastas',
                                       '{}.fa'.format(work_dir_i))
            command = [
                'perl', self._espresso_c, '-I', c_work_dir, '-F', split_fasta,
                '-X', '0', '-T', threads
            ]
            c_log = os.path.join(self._log_dir,
                                 'espresso_c_{}.log'.format(work_dir_i))
            tests.base_test.run_command_with_log(command, c_log)

    def _check_c_output(self):
        read_finals = list()
        for work_dir_i in range(self._num_c_splits):
            read_final = os.path.join(self._c_work_dir, str(work_dir_i), '0',
                                      'chr11_read_final.txt')
            read_finals.append(read_final)

        self._assert_exists_and_non_empty(read_finals)

    def _combine_c_for_q(self):
        command = [
            self._py_executable, self._combine_c_py, '--c-base-work-dir',
            self._c_work_dir, '--new-base-dir', self._q_work_dir
        ]
        combine_log = os.path.join(self._log_dir, 'combine_c_for_q.log')
        tests.base_test.run_command_with_log(command, combine_log)

    def _run_espresso_q(self):
        self._q_samples_updated = os.path.join(self._q_work_dir,
                                               'samples.tsv.updated')
        isoform_basename = '{}_isoform.tsv'.format(self._out_file_prefix)
        self._isoform_tsv = os.path.join(self._q_work_dir, isoform_basename)
        command = [
            'perl', self._espresso_q, '-A', self._cd44_gtf, '-L',
            self._q_samples_updated, '-V', self._isoform_tsv
        ]
        q_log = os.path.join(self._log_dir, 'espresso_q.log')
        tests.base_test.run_command_with_log(command, q_log)

    def _check_q_output(self):
        gtf_basename = '{}_updated.gtf'.format(self._out_file_prefix)
        self._updated_gtf = os.path.join(self._q_work_dir, gtf_basename)
        abundance_basename = '{}_abundance.esp'.format(self._out_file_prefix)
        self._abundance = os.path.join(self._q_work_dir, abundance_basename)
        self._assert_exists_and_non_empty(
            [self._updated_gtf, self._abundance, self._isoform_tsv])
        expected_headers = [
            'transcript_ID', 'transcript_name', 'gene_ID', 'GS689', 'PC3E'
        ]
        rows = dict()
        with open(self._abundance, 'rt') as abundance_handle:
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
                row['GS689'] = float(columns[3])
                row['PC3E'] = float(columns[4])
                rows[transcript_id] = row

        expected_transcripts = {
            'ENST00000415148.6_1': {
                'transcript_name': 'CD44-206',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 0,
                'PC3E': 110
            },
            'ENST00000526025.2_2': {
                'transcript_name': 'CD44-222',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 14,
                'PC3E': 10
            },
            'ENST00000525241.1_1': {
                'transcript_name': 'CD44-215',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 0,
                'PC3E': 1
            },
            'ENST00000434472.6_1': {
                'transcript_name': 'CD44-210',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 25,
                'PC3E': 8
            },
            'ENST00000527326.1_1': {
                'transcript_name': 'CD44-225',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 11,
                'PC3E': 11
            },
            'ENST00000433892.6_1': {
                'transcript_name': 'CD44-209',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 5,
                'PC3E': 477
            },
            'ENST00000531118.5_1': {
                'transcript_name': 'CD44-232',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 3,
                'PC3E': 12
            },
            'ENST00000534296.5_1': {
                'transcript_name': 'CD44-238',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 0,
                'PC3E': 1
            },
            'ENST00000531141.1_1': {
                'transcript_name': 'CD44-233',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 0,
                'PC3E': 1
            },
            'ENST00000263398.10_1': {
                'transcript_name': 'CD44-201',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 1168,
                'PC3E': 235
            },
            'ENST00000428726.8_3': {
                'transcript_name': 'CD44-208',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 0,
                'PC3E': 23
            },
            tests.base_test.get_espresso_novel_id('chr11', 0, 2): {
                'transcript_name': 'NA',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 0,
                'PC3E': 22
            },
            tests.base_test.get_espresso_novel_id('chr11', 0, 3): {
                'transcript_name': 'NA',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 0,
                'PC3E': 23
            },
            tests.base_test.get_espresso_novel_id('chr11', 0, 1): {
                'transcript_name': 'NA',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 1,
                'PC3E': 12
            },
            tests.base_test.get_espresso_novel_id('chr11', 0, 5): {
                'transcript_name': 'NA',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 0,
                'PC3E': 130
            },
            tests.base_test.get_espresso_novel_id('chr11', 0, 4): {
                'transcript_name': 'NA',
                'gene_id': 'ENSG00000026508.19_9',
                'GS689': 1,
                'PC3E': 75
            },
            tests.base_test.get_espresso_novel_id('chr11', 0, 0): {
                'transcript_name': 'NA',
                'gene_id': 'NA',
                'GS689': 0,
                'PC3E': 2
            }
        }
        self.assertEqual(len(rows), len(expected_transcripts))
        num_expected_novel = 0
        for expected_id in expected_transcripts:
            is_novel = expected_id.startswith('ESPRESSO')
            if is_novel:
                num_expected_novel += 1

        num_found_novel = 0
        for found_id in rows:
            is_novel = found_id.startswith('ESPRESSO')
            if is_novel:
                num_found_novel += 1

        self.assertEqual(num_found_novel, num_expected_novel)
        for transcript_id, expected in expected_transcripts.items():
            found = rows.get(transcript_id)
            self.assertTrue(found)
            self.assertEqual(found['gene_id'], expected['gene_id'])
            self.assertEqual(found['transcript_name'],
                             expected['transcript_name'])
            self._assert_within_x_percent_or_y(found['GS689'],
                                               expected['GS689'], 1, 1)
            self._assert_within_x_percent_or_y(found['PC3E'], expected['PC3E'],
                                               1, 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
