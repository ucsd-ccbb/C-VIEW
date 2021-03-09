import os
import filecmp
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from qc.seq_run_acceptance import make_collected_acceptances_tsv


class SeqRunAcceptanceTest(FileTestCase):
    def test_merge_tables(self):
        expected_results_fp = f"{self.test_qc_dir}/" \
                              f"2021-02-08-ARTIC-acceptance.tsv"
        output_fp = f"{self.test_temp_qc_dir}/temp_seq_run_acceptance.csv"
        arg_list = ["seq_run_acceptance.py", self.test_samples_dir, output_fp]

        try:
            make_collected_acceptances_tsv(arg_list)
            self.assertTrue(os.path.isfile(output_fp))
            self.assertTrue(filecmp.cmp(output_fp, expected_results_fp))
        finally:
            try:
                os.remove(output_fp)
            except OSError:
                pass
