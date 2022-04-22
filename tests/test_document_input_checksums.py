import os
import filecmp
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from pipeline.document_input_checksums import generate_checksums_file


class DocumentInputChecksumsTest(FileTestCase):
    def test_generate_checksums_file(self, make_all_empress_out=True):
        input_dir = self.dummy_qc_dir
        expected_checksums_fp = f"{self.dummy_dir}/dummy_input_checksums.csv"
        out_checksums_fp = f"{self.test_temp_dir}/temp_input_checksums.csv"
        arg_list = ["document_input_checksums.py", input_dir, out_checksums_fp]

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