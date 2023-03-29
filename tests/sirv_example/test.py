import os
import os.path
import unittest

import tests.base_test


class SirvExampleTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_name = 'sirv_example'
        self._test_dir = os.path.dirname(__file__)
        self._data_dir = os.path.join(self._test_dir, 'data')
        self._out_dir = os.path.join(self._test_dir, 'out')
        self._vis_out_dir = os.path.join(self._out_dir, 'visualization')
        self._log_dir = os.path.join(self._test_dir, 'logs')
        self._out_file_prefix = (
            tests.base_test.out_file_prefix_based_on_params())
        self._sirv_sam = None
        self._sirv_fasta = None
        self._sirv_gtf = None
        self._samples_tsv = None
        self._samples_updated = None
        self._isoform_tsv = None
        self._updated_gtf = None
        self._abundance = None

    def test(self):
        self._initialize_dirs()
        self._unpack_test_data()
        self._create_samples_tsv()
        self._run_espresso_s()
        self._check_s_output()
        self._run_espresso_c()
        self._check_c_output()
        self._run_espresso_q()
        self._check_q_output()
        self._run_visualization()
        self._check_visualization_output()

    def _initialize_dirs(self):
        directories = [
            self._data_dir, self._out_dir, self._vis_out_dir, self._log_dir
        ]
        tests.base_test.remove_then_create_directories(directories)

    def _unpack_test_data(self):
        test_tar_gz = os.path.join(self._test_data_dir,
                                   'test_data_espresso_sirv.tar.gz')
        command = ['tar', '-xvf', test_tar_gz, '--directory', self._data_dir]
        tar_log = os.path.join(self._log_dir, 'tar.log')
        tests.base_test.run_command_with_log(command, tar_log)
        extracted_base_dir = os.path.join(self._data_dir,
                                          'test_data_espresso_sirv')
        self._sirv_sam = os.path.join(extracted_base_dir, 'SIRV2_3.sort.sam')
        self._sirv_fasta = os.path.join(extracted_base_dir, 'SIRV2.fasta')
        self._sirv_gtf = os.path.join(extracted_base_dir, 'SIRV_C.gtf')
        self._assert_exists_and_non_empty(
            [self._sirv_sam, self._sirv_fasta, self._sirv_gtf])

    def _create_samples_tsv(self):
        self._samples_tsv = os.path.join(self._data_dir, 'samples.tsv')
        with open(self._samples_tsv, 'wt') as handle:
            columns = [self._sirv_sam, self._test_name]
            tests.base_test.write_tsv_line(handle, columns)

    def _run_espresso_s(self):
        command = [
            'perl', self._espresso_s, '-A', self._sirv_gtf, '-L',
            self._samples_tsv, '-F', self._sirv_fasta, '-O', self._out_dir
        ]
        s_log = os.path.join(self._log_dir, 'espresso_s.log')
        tests.base_test.run_command_with_log(command, s_log)

    def _check_s_output(self):
        self._samples_updated = os.path.join(self._out_dir,
                                             'samples.tsv.updated')
        sirv2_sjs = os.path.join(self._out_dir, 'SIRV2_SJ_simplified.list')
        all_sjs = os.path.join(self._out_dir, 'SJ_group_all.fa')
        sam_list = os.path.join(self._out_dir, '0', 'sam.list3')
        sj_list = os.path.join(self._out_dir, '0', 'sj.list')
        self._assert_exists_and_non_empty(
            [self._samples_updated, sirv2_sjs, all_sjs, sam_list, sj_list])

        with open(self._samples_updated, 'rt') as samples_handle:
            lines = list()
            for line in samples_handle:
                lines.append(line.strip())

        self.assertEqual(len(lines), 1)
        parts = lines[0].split('\t')
        self.assertEqual(len(parts), 3)
        self.assertEqual(parts[0], self._sirv_sam)
        self.assertEqual(parts[1], self._test_name)
        self.assertEqual(parts[2], '0')

    def _run_espresso_c(self):
        threads = '2'
        command = [
            'perl', self._espresso_c, '-I', self._out_dir, '-F',
            self._sirv_fasta, '-X', '0', '-T', threads
        ]
        c_log = os.path.join(self._log_dir, 'espresso_c.log')
        tests.base_test.run_command_with_log(command, c_log)

    def _check_c_output(self):
        read_final = os.path.join(self._out_dir, '0', 'SIRV2_read_final.txt')
        self._assert_exists_and_non_empty([read_final])

    def _run_espresso_q(self):
        isoform_basename = '{}_isoform.tsv'.format(self._out_file_prefix)
        self._isoform_tsv = os.path.join(self._out_dir, isoform_basename)
        command = [
            'perl', self._espresso_q, '-A', self._sirv_gtf, '-L',
            self._samples_updated, '-V', self._isoform_tsv
        ]
        q_log = os.path.join(self._log_dir, 'espresso_q.log')
        tests.base_test.run_command_with_log(command, q_log)

    def _check_q_output(self):
        gtf_basename = '{}_updated.gtf'.format(self._out_file_prefix)
        self._updated_gtf = os.path.join(self._out_dir, gtf_basename)
        abundance_basename = '{}_abundance.esp'.format(self._out_file_prefix)
        self._abundance = os.path.join(self._out_dir, abundance_basename)
        self._assert_exists_and_non_empty(
            [self._updated_gtf, self._abundance, self._isoform_tsv])
        expected_headers = [
            'transcript_ID', 'transcript_name', 'gene_ID', self._test_name
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
                row['abundance'] = float(columns[3])
                rows[transcript_id] = row

        expected_transcripts = {
            'SIRV201': {
                'gene_id': 'SIRV2',
                'abundance': 296
            },
            'SIRV202': {
                'gene_id': 'SIRV2',
                'abundance': 2407
            },
            'SIRV203': {
                'gene_id': 'SIRV2',
                'abundance': 336
            }
        }
        self.assertEqual(sorted(rows.keys()),
                         sorted(expected_transcripts.keys()))
        for transcript_id, expected in expected_transcripts.items():
            found = rows.get(transcript_id)
            self.assertTrue(found)
            self.assertEqual(found['gene_id'], expected['gene_id'])
            self._assert_within_x_percent_or_y(found['abundance'],
                                               expected['abundance'], 1, 1)

    def _run_visualization(self):
        command = [
            self._py_executable, self._visualization_py, '--genome-fasta',
            self._sirv_fasta, '--updated-gtf', self._updated_gtf,
            '--abundance-esp', self._abundance, '--target-gene', 'SIRV2',
            '--minimum-count', '1', '--descriptive-name', 'SIRV',
            '--output-dir', self._vis_out_dir
        ]
        vis_log = os.path.join(self._log_dir, 'visualization.log')
        tests.base_test.run_command_with_log(command, vis_log)

    def _check_visualization_output(self):
        sirv_fasta_fai = '{}.fai'.format(self._sirv_fasta)
        big_wig_basename = '{}.bw'.format(self._test_name)
        big_wig = os.path.join(self._vis_out_dir, big_wig_basename)
        gene = 'SIRV2'
        bed_paths = list()
        for transcript in ['SIRV201', 'SIRV202', 'SIRV203']:
            bed_basename = '{}_{}_{}.bed'.format(self._test_name, gene,
                                                 transcript)
            bed_path = os.path.join(self._vis_out_dir, 'target_genes',
                                    bed_basename)
            bed_paths.append(bed_path)

        self._assert_exists_and_non_empty([sirv_fasta_fai, big_wig] +
                                          bed_paths)


if __name__ == '__main__':
    unittest.main(verbosity=2)
