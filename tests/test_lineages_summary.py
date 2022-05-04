import os
import filecmp
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from src.lineages_summary import create_lineages_summary


class LineagesSummaryTest(FileTestCase):
    def test_create_lineages_summary(self):
        input_added_names_fp = f"{self.dummy_dir}/" \
                               f"dummy_added_fa_names.txt"
        input_lineage_fp = f"{self.dummy_dir}/" \
                           f"dummy_merged.lineage_report.csv"
        expected_qc_and_lineages_fp = f"{self.dummy_dir}/" \
                                      f"dummy_qc_and_lineages.csv"

        out_summary_fp = f"{self.test_temp_dir}/temp_qc_and_lineages.csv"
        arg_list = ["lineages_summary.py", input_added_names_fp,
                    self.dummy_dir, "-summary.csv",
                    input_lineage_fp,
                    out_summary_fp]

        qc_is_file = qc_equal = False
        try:
            create_lineages_summary(arg_list)

            qc_is_file = os.path.isfile(out_summary_fp)
            self.assertTrue(qc_is_file)

            qc_equal = filecmp.cmp(out_summary_fp, expected_qc_and_lineages_fp)
            self.assertTrue(qc_equal)
        finally:
            if qc_is_file and qc_equal:
                for fp in [out_summary_fp]:
                    try:
                        os.remove(fp)
                    except OSError:
                        pass
