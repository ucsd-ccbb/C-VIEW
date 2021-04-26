import os
import filecmp
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from qc.integrate_bjorn_coverage import integrate_bjorn_coverage


class IntegrateBjornCoverageTests(FileTestCase):
    def test_integrate_bjorn_coverage(self):
        input_summary_wo_bjorn_cov_fp = f"{self.dummy_dir}/" \
                                        f"dummy_summary_wo_bjorn_coverage.csv"
        input_bjorn_cov_fp = f"{self.dummy_dir}/dummy-coverage.tsv"
        expected_results_fp = f"{self.dummy_dir}/" \
                              f"dummy_summary_w_bjorn_coverage.csv"
        output_fp = f"{self.test_temp_dir}/temp_summary_w_bjorn_coverage.csv"
        arg_list = ["integrate_bjorn_coverage.py",
                    input_summary_wo_bjorn_cov_fp,
                    input_bjorn_cov_fp, output_fp]

        passes = False
        try:
            integrate_bjorn_coverage(arg_list)
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
