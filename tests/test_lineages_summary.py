import os
import filecmp
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from qc.lineages_summary import create_lineages_summary_and_metadata


class LineagesSummaryTest(FileTestCase):
    def test_create_lineages_summary_and_metadata(self):
        input_added_names_fp = f"{self.dummy_dir}/" \
                               f"dummy_added_fa_names.txt"
        input_lineage_fp = f"{self.dummy_dir}/" \
                           f"dummy_merged.lineage_report.csv"
        expected_qc_and_lineages_fp = f"{self.dummy_dir}/" \
                                      f"dummy_qc_and_lineages.csv"
        expected_metadata_fp = f"{self.dummy_dir}/" \
                               f"dummy_metadata.tsv"

        out_summary_fp = f"{self.test_temp_dir}/temp_qc_and_lineages.csv"
        out_metadata_fp = f"{self.test_temp_dir}/temp_metadata.tsv"
        arg_list = ["lineages_summary.py", input_added_names_fp,
                    self.dummy_dir, "-summary.csv",
                    input_lineage_fp,
                    out_summary_fp, out_metadata_fp]

        qc_is_file = qc_equal = meta_is_file = meta_equal = False
        try:
            create_lineages_summary_and_metadata(arg_list)

            qc_is_file = os.path.isfile(expected_qc_and_lineages_fp)
            self.assertTrue(qc_is_file)

            qc_equal = filecmp.cmp(out_summary_fp, expected_qc_and_lineages_fp)
            self.assertTrue(qc_equal)

            meta_is_file = os.path.isfile(out_metadata_fp)
            self.assertTrue(meta_is_file)

            meta_equal = filecmp.cmp(out_metadata_fp, expected_metadata_fp)
            self.assertTrue(meta_equal)
        finally:
            if qc_is_file and qc_equal and meta_is_file and meta_equal:
                for fp in [out_summary_fp, out_metadata_fp]:
                    try:
                        os.remove(fp)
                    except OSError:
                        pass
