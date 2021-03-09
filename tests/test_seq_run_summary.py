import os
import filecmp
from unittest import TestCase
from qc.seq_run_summary import merge_multiqc_and_acceptance


class SeqRunSummaryTest(TestCase):
    test_data_path = os.path.abspath("./data")

    def test_merge_multiqc_and_acceptance(self):
        qc_output_dir = f"{self.test_data_path}" \
                        f"/2021-02-08-ARTIC_quality_control/"
        input_stats_fp = f"{qc_output_dir}" \
                         f"2021-02-08-ARTIC_multiqc_data/" \
                         f"multiqc_general_stats.txt"
        input_acceptance_fp = f"{qc_output_dir}/" \
                              f"2021-02-08-ARTIC-acceptance.tsv"
        expected_results_fp = f"{qc_output_dir}/" \
                              f"2021-02-08-ARTIC-summary.csv"
        qc_temp_dir = f"{self.test_data_path}/qc"
        output_fp = f"{qc_temp_dir}/temp_seq_run_summary.csv"
        arg_list = ["seq_run_summary.py", input_stats_fp,
                    input_acceptance_fp, output_fp]

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
