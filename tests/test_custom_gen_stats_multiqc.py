import os
import yaml
from unittest import TestCase
from qc.custom_gen_stats_multiqc import write_custom_multiqc_yaml


class CustomGenStatsMultiQcTest(TestCase):
    def test_write_custom_multiqc_yaml(self):
        self.maxDiff = None
        test_data_path = os.path.abspath("./data")

        qc_data_path = f"{test_data_path}/qc"
        yaml_paths_fp = f"{qc_data_path}/qualimapReport_paths.txt"
        fastqc_data_paths_fp = f"{qc_data_path}/r1_r2_fastqc_data_paths.txt"
        output_fp = f"{qc_data_path}/temp_test_cust_gen_stats_multiqc.yaml"
        expected_result_fp = f"{qc_data_path}/multiqc_custom_gen_stats.yaml"

        arg_list = ["custom_gen_stats_multiqc.py", yaml_paths_fp,
                    fastqc_data_paths_fp, "pe", output_fp]

        with open(expected_result_fp) as f:
            expected_dict = yaml.load(f, Loader=yaml.FullLoader)

        out_dict = {}
        try:
            out_dict = write_custom_multiqc_yaml(arg_list)
            # TODO: add test that output file is in the location expected and
            #  not empty
        finally:
            try:
                os.remove(output_fp)
            except OSError:
                pass

        self.assertEqual(expected_dict, out_dict)
