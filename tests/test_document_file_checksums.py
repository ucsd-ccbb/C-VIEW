import os
import filecmp
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from pipeline.document_file_checksums import generate_checksums_file


class DocumentFileChecksumsTest(FileTestCase):
    def test_generate_checksums_file_all(self):
        self._test_generate_checksums_file(limit_fnames=False)

    def test_generate_checksums_file_specific(self):
        self._test_generate_checksums_file(limit_fnames=True)

    def _test_generate_checksums_file(self, limit_fnames=False):
        input_dir = self.dummy_qc_dir
        limited = "" if not limit_fnames else "_limited"
        expected_checksums_fp = f"{self.dummy_dir}/" \
                                f"dummy_input_checksums{limited}.csv"
        out_checksums_fp = f"{self.test_temp_dir}/temp_input_checksums.csv"
        arg_list = ["document_file_checksums.py", input_dir, out_checksums_fp]
        if limit_fnames:
            arg_list.extend(["2021-02-08-ARTIC-summary.csv",
                             "_general_stats.txt"])

        output_is_file = False
        output_equal = False
        try:
            generate_checksums_file(arg_list)

            output_is_file = os.path.isfile(out_checksums_fp)
            self.assertTrue(output_is_file)

            output_equal = filecmp.cmp(out_checksums_fp, expected_checksums_fp)
            self.assertTrue(output_equal)
        finally:
            if output_is_file and output_equal:
                for fp in [out_checksums_fp]:
                    try:
                        os.remove(fp)
                    except OSError:
                        pass
