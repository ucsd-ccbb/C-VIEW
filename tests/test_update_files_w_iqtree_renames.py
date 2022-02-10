import os
import filecmp
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from pipeline.update_files_w_iqtree_renames import make_files_w_iqtree_renames


class UpdateWithIqtreeRenamesTest(FileTestCase):
    def test_make_files_w_iqtree_renames(self):
        input_tree_build_log_fp = f"{self.dummy_dir}/dummy_treebuild.log"
        input_align_fp = f"{self.dummy_dir}/dummy_short.aln"
        input_metadata_fp = f"{self.dummy_dir}/" \
                            f"dummy_empress_metadata.tsv"

        expected_align_fp = f"{self.dummy_dir}/" \
                            f"dummy_short_w_iqtree_renames.aln"
        expected_metadata_fp = f"{self.dummy_dir}/" \
                               f"dummy_empress_metadata_w_iqtree_renames.tsv"

        out_align_fp = f"{self.test_temp_dir}/" \
                       f"temp_aln_w_iqtree_renames.aln"
        out_metadata_fp = f"{self.test_temp_dir}/" \
                          f"temp_empress_metadata_w_iqtree_renames.tsv"

        arg_list = ["update_files_w_iqtree_renames.py",
                    input_tree_build_log_fp, input_align_fp, input_metadata_fp,
                    out_align_fp, out_metadata_fp]

        align_is_file = metadata_is_file = False
        align_equal = metadata_equal = False
        try:
            make_files_w_iqtree_renames(arg_list)

            align_is_file = os.path.isfile(out_align_fp)
            self.assertTrue(align_is_file)

            align_equal = filecmp.cmp(out_align_fp, expected_align_fp)
            self.assertTrue(align_equal)

            metadata_is_file = os.path.isfile(out_metadata_fp)
            self.assertTrue(metadata_is_file)

            metadata_equal = filecmp.cmp(out_metadata_fp, expected_metadata_fp)
            self.assertTrue(metadata_equal)
        finally:
            if align_is_file and metadata_is_file and align_equal and \
                    metadata_equal:
                for fp in [out_align_fp, out_metadata_fp]:
                    try:
                        os.remove(fp)
                    except OSError:
                        pass
