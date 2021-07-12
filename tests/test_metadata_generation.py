import os
import filecmp
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from qc.metadata_generation import merge_metadata


class MetadataGenerationTest(FileTestCase):
    def test_merge_metadata(self):
        input_qc_and_lineages_fp = f"{self.dummy_dir}/" \
                                      f"dummy_qc_and_lineages_w_bjorn.csv"
        input_metadata_fp = f"{self.dummy_dir}/dummy_inspect_metadata.csv"
        expected_full_metadata_fp = f"{self.dummy_dir}/dummy_full_metadata.csv"
        expected_bjorn_metadata_fp = f"{self.dummy_dir}/" \
                                     f"dummy_bjorn_metadata.csv"

        out_full_fp = f"{self.test_temp_dir}/temp_full_metadata.csv"
        out_bjorn_fp = f"{self.test_temp_dir}/temp_bjorn_metadata.csv"
        arg_list = ["metadata_generation.py", input_qc_and_lineages_fp,
                    input_metadata_fp, out_full_fp, out_bjorn_fp]

        full_is_file = bjorn_is_file = False
        full_equal = bjorn_equal = False
        try:
            merge_metadata(arg_list)

            full_is_file = os.path.isfile(out_full_fp)
            self.assertTrue(full_is_file)

            full_equal = filecmp.cmp(out_full_fp, expected_full_metadata_fp)
            self.assertTrue(full_equal)

            bjorn_is_file = os.path.isfile(out_bjorn_fp)
            self.assertTrue(bjorn_is_file)

            bjorn_equal = filecmp.cmp(out_bjorn_fp, expected_bjorn_metadata_fp)
            self.assertTrue(bjorn_equal)
        finally:
            if full_is_file and full_equal and bjorn_is_file and bjorn_equal:
                for fp in [out_full_fp, out_bjorn_fp]:
                    try:
                        os.remove(fp)
                    except OSError:
                        pass
