import os
import filecmp
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from src.custom_reports_generation import make_bespoke_outputs


class CustomReportsGenerationTest(FileTestCase):
    def test_make_bespoke_outputs(self):
        input_full_summary_fp = f"{self.dummy_dir}/" \
                                f"dummy_full_metadata_for_custom_reports.csv"
        prefix_base = "temp_report"
        input_prefix_for_outputs = f"{self.test_temp_dir}/{prefix_base}"
        expected_outputs_fnames = ["temp_report_all.csv",
                                   "temp_report_CALM_RTL_wastewater.csv",
                                   "temp_report_CALM_RTL_not-wastewater.csv",
                                   "temp_report_Control Sample.csv",
                                   "temp_report_HELIX.csv",
                                   "temp_report_no-source.csv"]

        arg_list = ["custom_reports_generation.py", input_full_summary_fp,
                    input_prefix_for_outputs]

        output_files_info = []
        for curr_expected_output_fname in expected_outputs_fnames:
            putative_out_fp = f"{self.test_temp_dir}/" \
                              f"{curr_expected_output_fname}"

            known_good_name = curr_expected_output_fname.replace(
                prefix_base, "dummy_report")
            expected_out_fp = f"{self.dummy_dir}/{known_good_name}"
            # False = by default assume all expected outputs failed tests
            output_files_info.append([putative_out_fp, expected_out_fp, False])

        try:
            make_bespoke_outputs(arg_list)

            for curr_output_info in output_files_info:
                putative_out_is_file = os.path.isfile(curr_output_info[0])
                self.assertTrue(putative_out_is_file, msg=curr_output_info[0])

                putative_out_equal = filecmp.cmp(curr_output_info[0],
                                                 curr_output_info[1])
                self.assertTrue(putative_out_equal, msg=curr_output_info[0])

                # if we made it this far, the current output file passes
                curr_output_info[2] = True

        finally:
            for curr_output_info in output_files_info:
                if curr_output_info[2]:  # if this output passed testing
                    try:
                        os.remove(curr_output_info[0])
                    except OSError:
                        pass
