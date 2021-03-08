import os
import filecmp
from unittest import TestCase
from qc.seq_run_acceptance import make_collected_acceptances_tsv


class SeqRunAcceptanceTest(TestCase):
    test_data_path = os.path.abspath("./data")

    def test_merge_tables(self):
        qc_output_dir = f"{self.test_data_path}" \
                        f"/2021-02-08-ARTIC_quality_control/"
        samples_dir = f"{self.test_data_path}/" \
                      f"2021-02-08-ARTIC_samples"
        expected_results_fp = f"{qc_output_dir}/" \
                              f"2021-02-08-ARTIC_acceptance.tsv"
        qc_temp_dir = f"{self.test_data_path}/qc"
        output_fp = f"{qc_temp_dir}/temp_seq_run_acceptance.csv"
        arg_list = ["seq_run_acceptance.py", samples_dir, output_fp]

        try:
            make_collected_acceptances_tsv(arg_list)
            self.assertTrue(os.path.isfile(output_fp))
            self.assertTrue(filecmp.cmp(output_fp, expected_results_fp))
        finally:
            try:
                os.remove(output_fp)
            except OSError:
                pass
