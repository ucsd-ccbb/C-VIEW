import os
import filecmp
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from src.seq_run_acceptance import make_collected_acceptances_tsv


class SeqRunAcceptanceTest(FileTestCase):
    def test_merge_tables(self):
        expected_results_fp = f"{self.dummy_dir}/" \
                              f"dummy_PDH_83-233854622-acceptance.tsv"
        output_fp = f"{self.test_temp_dir}/temp_seq_run_acceptance.tsv"
        arg_list = ["seq_run_acceptance.py",
                    f"{self.dummy_dir}/dummy_PDH_83-233854622_samples/",
                    output_fp]

        out_is_file = out_matches = False
        try:
            make_collected_acceptances_tsv(arg_list)

            out_is_file = os.path.isfile(output_fp)
            self.assertTrue(out_is_file)

            out_matches = filecmp.cmp(output_fp, expected_results_fp)
            self.assertTrue(out_matches)
        finally:
            if out_is_file and out_matches:
                try:
                    os.remove(output_fp)
                except OSError:
                    pass
