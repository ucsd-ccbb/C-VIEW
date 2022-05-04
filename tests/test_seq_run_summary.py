import os
import filecmp
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from src.seq_run_summary import merge_multiqc_and_acceptance


class SeqRunSummaryTest(FileTestCase):
    def test_merge_multiqc_and_acceptance(self):
        input_stats_fp = f"{self.dummy_qc_dir}/" \
                         f"2021-02-08-ARTIC_multiqc_data/" \
                         f"multiqc_general_stats.txt"
        input_acceptance_fp = f"{self.dummy_qc_dir}/" \
                              f"2021-02-08-ARTIC-acceptance.tsv"
        input_pi_metric_fp = f"{self.dummy_qc_dir}/" \
                             f"2021-02-08-ARTIC-pi_metric.tsv"
        input_n_metric_fp = f"{self.dummy_qc_dir}/" \
                            f"2021-02-08-ARTIC-n_metric.tsv"
        expected_results_fp = f"{self.dummy_qc_dir}/" \
                              f"2021-02-08-ARTIC-summary.csv"
        output_fp = f"{self.test_temp_dir}/temp_seq_run_summary.csv"
        arg_list = ["seq_run_summary.py", output_fp, input_stats_fp,
                    input_acceptance_fp, input_pi_metric_fp,
                    input_n_metric_fp]

        passes = False
        try:
            merge_multiqc_and_acceptance(arg_list)
            passes = os.path.isfile(output_fp)
            self.assertTrue(passes)

            file_match = filecmp.cmp(output_fp, expected_results_fp)
            passes = passes and file_match
            self.assertTrue(file_match)
        finally:
            if passes:
                try:
                    os.remove(output_fp)
                except OSError:
                    pass
