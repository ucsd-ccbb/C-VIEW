import os
import filecmp
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from qc.tree_prep import prep_files_for_tree_building


class TreePrepTest(FileTestCase):
    def test_tree_prep(self):
        input_qc_and_lineages_fp = f"{self.dummy_dir}/" \
                                      f"dummy_qc_and_lineages_w_bjorn_" \
                                      f"for_tree_prep.csv"
        input_metadata_fp = f"{self.dummy_dir}/" \
                            f"dummy_inspect_metadata_for_tree_prep.csv"
        in_fas_fp = f"{self.dummy_dir}/dummy.fas"

        expected_stringent_fas_fp = f"{self.dummy_dir}/" \
                                    f"dummy_stringent_only.fas"
        expected_all_empress_metadata_fp = f"{self.dummy_dir}/" \
                                           f"dummy_all_empress_metadata.tsv"
        expected_stringent_empress_metadata_fp = \
            f"{self.dummy_dir}/dummy_stringent_empress_metadata.tsv"

        out_stringent_fas_fp = f"{self.test_temp_dir}/temp_stringent_only.fas"
        out_all_empress_fp = f"{self.test_temp_dir}/" \
                             f"temp_all_empress_metadata.tsv"
        out_stringent_empress_fp = f"{self.test_temp_dir}/" \
                                   f"temp_stringent_empress_metadata.tsv"

        arg_list = ["tree_prep.py", input_qc_and_lineages_fp,
                    input_metadata_fp, in_fas_fp, out_stringent_fas_fp,
                    out_all_empress_fp, out_stringent_empress_fp]

        all_empress_is_file = False
        stringent_empress_is_file = False
        stringent_fas_is_file = False
        all_empress_equal = False
        stringent_empress_equal = False
        stringent_fas_equal = False
        try:
            prep_files_for_tree_building(arg_list)

            all_empress_is_file = os.path.isfile(out_all_empress_fp)
            self.assertTrue(all_empress_is_file)

            all_empress_equal = filecmp.cmp(
                out_all_empress_fp, expected_all_empress_metadata_fp)
            self.assertTrue(all_empress_equal)

            stringent_empress_is_file = os.path.isfile(
                out_stringent_empress_fp)
            self.assertTrue(stringent_empress_is_file)

            stringent_empress_equal = filecmp.cmp(
                out_stringent_empress_fp,
                expected_stringent_empress_metadata_fp)
            self.assertTrue(stringent_empress_equal)

            stringent_fas_is_file = os.path.isfile(
                out_stringent_fas_fp)
            self.assertTrue(stringent_fas_is_file)

            stringent_fas_equal = filecmp.cmp(
                out_stringent_fas_fp, expected_stringent_fas_fp)
            self.assertTrue(stringent_fas_equal)
        finally:
            if all_empress_is_file and all_empress_equal and \
                    stringent_empress_is_file and stringent_empress_equal and \
                    stringent_fas_is_file and stringent_fas_equal:
                for fp in [out_all_empress_fp, out_stringent_empress_fp,
                           out_stringent_fas_fp]:
                    try:
                        os.remove(fp)
                    except OSError:
                        pass
