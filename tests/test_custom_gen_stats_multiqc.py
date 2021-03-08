import os
from pathlib import Path
import yaml
from unittest import TestCase
from qc.custom_gen_stats_multiqc import write_custom_multiqc_yaml


class CustomGenStatsMultiQcTest(TestCase):
    test_data_path = os.path.abspath("./data")
    samples_path = f"{test_data_path}/2021-02-08-ARTIC_samples"
    qc_data_path = f"{test_data_path}/qc"

    def _make_paths_file(self, output_fname, glob_pattern):
        output_fp = f"{self.qc_data_path}/{output_fname}"
        relevant_paths = Path(self.samples_path).rglob(glob_pattern)
        relevant_fps = [f"{str(p)}\n" for p in relevant_paths]
        sorted_relevant_fps = sorted(relevant_fps)
        with open(output_fp, "w") as f:
            f.writelines(sorted_relevant_fps)
        return output_fp

    def test_write_custom_multiqc_yaml(self):
        self.maxDiff = None

        qualimap_paths_fp = self._make_paths_file(
            "qualimapReport_paths.txt", 'qualimapReport.html')
        fastqc_data_paths_fp = self._make_paths_file(
            "fastqc_data_paths.txt", 'fastqc_data.txt')
        output_fp = f"{self.qc_data_path}/" \
                    f"temp_test_cust_gen_stats_multiqc.yaml"
        expected_result_fp = f"{self.qc_data_path}/" \
                             f"multiqc_custom_gen_stats.yaml"

        arg_list = ["custom_gen_stats_multiqc.py", qualimap_paths_fp,
                    fastqc_data_paths_fp, "pe", output_fp]

        with open(expected_result_fp) as f:
            expected_dict = yaml.load(f, Loader=yaml.FullLoader)

        out_dict = {}
        try:
            out_dict = write_custom_multiqc_yaml(arg_list)
            self.assertTrue(os.path.isfile(output_fp))
        finally:
            for fp in [qualimap_paths_fp, fastqc_data_paths_fp, output_fp]:
                try:
                    os.remove(fp)
                except OSError:
                    pass

        self.assertEqual(expected_dict, out_dict)
