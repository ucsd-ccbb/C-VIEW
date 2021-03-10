import os
import filecmp
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from qc.subset_csv import filter_csv

# TODO: add testing for dir_prefix case


class SubsetCsvTest(FileTestCase):
    def test_filter_csv_filter_lines(self):
        output_fp = f"{self.test_temp_qc_dir}/temp_test_filter_lines.csv"
        expected_fp = f"{self.test_data_dir}/dummy_filtered_lines.csv"
        args = ["subset_csv.py",
                f"{self.test_data_dir}/dummy_table.csv",
                "filtered_lines",
                "seq_run",
                "210109_A00953_0209_BHYHCVDRXX,PDH_79-233476243",
                output_fp]

        try:
            real_str, real_code = filter_csv(args)
            self.assertEqual("", real_str)
            self.assertEqual(0, real_code)
            self.assertTrue(os.path.isfile(output_fp))
            self.assertTrue(filecmp.cmp(output_fp, expected_fp))
        finally:
            try:
                os.remove(output_fp)
            except OSError:
                pass

    def test_filter_csv_filter_lines_accepted_cons_fnames(self):
        expected_str = "040idSEARCH-5367-SAN_L001_L002_L003_L004.trimmed.sorted.pileup.consensus.fa 041idSEARCH-5368-SAN_L001_L002_L003_L004.trimmed.sorted.pileup.consensus.fa 042idSEARCH-5369-SAN_L001_L002_L003_L004.trimmed.sorted.pileup.consensus.fa 044idSEARCH-5371-SAN_L001_L002_L003_L004.trimmed.sorted.pileup.consensus.fa"  #noqa 501
        args = ["subset_csv.py",
                f"{self.test_data_dir}/dummy_table.csv",
                "accepted_cons_fnames"]
        real_str, real_code = filter_csv(args)
        self.assertEqual(expected_str, real_str)
        self.assertEqual(0, real_code)

    def test_filter_csv_filter_lines_indel_flagged_cons_fnames(self):
        expected_str = "043idSEARCH-5370-SAN_L001_L002_L003_L004.trimmed.sorted.pileup.consensus.fa 045idSEARCH-5371-SAN_L001_L002_L003_L004.trimmed.sorted.pileup.consensus.fa"  #noqa 501
        args = ["subset_csv.py",
                f"{self.test_data_dir}/dummy_table.csv",
                "indel_flagged_cons_fnames"]
        real_str, real_code = filter_csv(args)
        self.assertEqual(expected_str, real_str)
        self.assertEqual(0, real_code)

    def test_filter_csv_filter_lines_error(self):
        expected_str = "Unrecognized filter type 'flagged_cons_fnames'"
        args = ["subset_csv.py",
                f"{self.test_data_dir}/dummy_table.csv",
                "flagged_cons_fnames"]
        real_str, real_code = filter_csv(args)
        self.assertEqual(expected_str, real_str)
        self.assertEqual(1, real_code)
