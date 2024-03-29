import os
from pathlib import Path
import yaml
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from src.custom_gen_stats_multiqc import write_custom_multiqc_yaml, \
    insert_q30_based_values, _calc_pct_aligned


class CustomGenStatsMultiQcTest(FileTestCase):
    def _make_paths_file(self, output_fname, glob_pattern):
        output_fp = f"{self.test_temp_dir}/{output_fname}"
        relevant_paths = Path(self.dummy_samples_dir).rglob(glob_pattern)
        relevant_fps = [f"{str(p)}\n" for p in relevant_paths]
        sorted_relevant_fps = sorted(relevant_fps)
        with open(output_fp, "w") as f:
            f.writelines(sorted_relevant_fps)
        return output_fp

    def test_insert_q30_based_values_se(self):
        r1 = '/Users/amandabirmingham/Work/Repositories/' \
             'covid_sequencing_analysis_pipeline/tests/data/dummy/' \
             'dummy_2021-02-08-ARTIC_samples/SEARCH-5329__LIBPLATE1__' \
             'A01__001002003004/SEARCH-5329__LIBPLATE1__A01__' \
             '001002003004_R1_q30_data.txt'
        s1 = [2440812, 2319178]
        input_dict = {}
        expected_dict = {'SEARCH-5329__LIBPLATE1__A01__001002003004': {
            'Pct >=Q30': 95.017, 'Uncapped Reads': 2440812}}
        insert_q30_based_values(input_dict, r1, s1)

        self.assertDictEqual(expected_dict, input_dict)

    def test_insert_q30_based_values_se_zero(self):
        r1 = '/Users/amandabirmingham/Work/Repositories/' \
             'covid_sequencing_analysis_pipeline/tests/data/dummy/' \
             'dummy_2021-02-08-ARTIC_samples/SEARCH-5329__LIBPLATE1__' \
             'A01__001002003004/SEARCH-5329__LIBPLATE1__A01__' \
             '001002003004_R1_q30_data.txt'
        s1 = [0, 0]
        input_dict = {}
        expected_dict = {'SEARCH-5329__LIBPLATE1__A01__001002003004': {
            'Pct >=Q30': 'NA', 'Uncapped Reads': 0}}
        insert_q30_based_values(input_dict, r1, s1)

        self.assertDictEqual(expected_dict, input_dict)

    def test_insert_q30_based_values_pe(self):
        r1 = '/Users/amandabirmingham/Work/Repositories/' \
             'covid_sequencing_analysis_pipeline/tests/data/dummy/' \
             'dummy_2021-02-08-ARTIC_samples/SEARCH-5329__LIBPLATE1__' \
             'A01__001002003004/SEARCH-5329__LIBPLATE1__A01__' \
             '001002003004_R1_q30_data.txt'
        r2 = '/Users/amandabirmingham/Work/Repositories/' \
             'covid_sequencing_analysis_pipeline/tests/data/dummy/' \
             'dummy_2021-02-08-ARTIC_samples/SEARCH-5329__LIBPLATE1__' \
             'A01__001002003004/SEARCH-5329__LIBPLATE1__A01__' \
             '001002003004_R2_q30_data.txt'
        s1 = [2440812, 2319178]
        s2 = [2440812, 2282795]
        input_dict = {}
        expected_dict = {'SEARCH-5329__LIBPLATE1__A01__001002003004': {
            'Pct >=Q30': 94.271, 'Uncapped Reads': 4881624}}
        insert_q30_based_values(input_dict, r1, s1, r2, s2)

        self.assertDictEqual(expected_dict, input_dict)

    def test_insert_q30_based_values_pe_zero(self):
        r1 = '/Users/amandabirmingham/Work/Repositories/' \
             'covid_sequencing_analysis_pipeline/tests/data/dummy/' \
             'dummy_2021-02-08-ARTIC_samples/SEARCH-5329__LIBPLATE1__' \
             'A01__001002003004/SEARCH-5329__LIBPLATE1__A01__' \
             '001002003004_R1_q30_data.txt'
        r2 = '/Users/amandabirmingham/Work/Repositories/' \
             'covid_sequencing_analysis_pipeline/tests/data/dummy/' \
             'dummy_2021-02-08-ARTIC_samples/SEARCH-5329__LIBPLATE1__' \
             'A01__001002003004/SEARCH-5329__LIBPLATE1__A01__' \
             '001002003004_R2_q30_data.txt'
        s1 = [0, 0]
        s2 = [2440812, 2282795]
        input_dict = {}
        expected_dict = {'SEARCH-5329__LIBPLATE1__A01__001002003004': {
            'Pct >=Q30': 93.526, 'Uncapped Reads': 2440812}}
        insert_q30_based_values(input_dict, r1, s1, r2, s2)

        self.assertDictEqual(expected_dict, input_dict)

    def test_insert_q30_based_values_pe_zero_both(self):
        r1 = '/Users/amandabirmingham/Work/Repositories/' \
             'covid_sequencing_analysis_pipeline/tests/data/dummy/' \
             'dummy_2021-02-08-ARTIC_samples/SEARCH-5329__LIBPLATE1__' \
             'A01__001002003004/SEARCH-5329__LIBPLATE1__A01__' \
             '001002003004_R1_q30_data.txt'
        r2 = '/Users/amandabirmingham/Work/Repositories/' \
             'covid_sequencing_analysis_pipeline/tests/data/dummy/' \
             'dummy_2021-02-08-ARTIC_samples/SEARCH-5329__LIBPLATE1__' \
             'A01__001002003004/SEARCH-5329__LIBPLATE1__A01__' \
             '001002003004_R2_q30_data.txt'
        s1 = [0, 0]
        s2 = [0, 0]
        input_dict = {}
        expected_dict = {'SEARCH-5329__LIBPLATE1__A01__001002003004': {
            'Pct >=Q30': 'NA', 'Uncapped Reads': 0}}
        insert_q30_based_values(input_dict, r1, s1, r2, s2)

        self.assertDictEqual(expected_dict, input_dict)

    def test__calc_pct_aligned(self):
        real_out = _calc_pct_aligned([5, 5])
        self.assertEqual(50.0, real_out)

    def test__calc_pct_aligned_zero_both(self):
        real_out = _calc_pct_aligned([0, 0])
        self.assertEqual('NA', real_out)

    def test_write_custom_multiqc_yaml_pe(self):
        self._help_test_write_custom_multiqc_yaml("pe")

    def test_write_custom_multiqc_yaml_se(self):
        self._help_test_write_custom_multiqc_yaml("se")

    def test_write_custom_multiqc_yaml_bam(self):
        self._help_test_write_custom_multiqc_yaml("bam")

    def _help_test_write_custom_multiqc_yaml(self, input_type):
        self.maxDiff = None

        qualimap_paths_fp = self._make_paths_file(
            "qualimapReport_paths.txt", 'qualimapReport.html')

        q30_glob_insert = '_R1_' if input_type == "se" else ""
        q30_data_paths_fp = self._make_paths_file(
            "q30_data_paths.txt", f'*{q30_glob_insert}q30_data.txt')
        sub_map_paths_fp = self._make_paths_file(
            "sub_map_paths.txt", '*_subsampled_mapping_stats.tsv')

        output_fp = f"{self.test_temp_dir}/" \
                    f"temp_test_cust_gen_stats_multiqc_{input_type}.yaml"
        expected_result_fp = f"{self.dummy_dir}/" \
                             f"dummy_2021-02-08-ARTIC-multiqc_" \
                             f"custom_gen_stats_{input_type}.yaml"

        arg_list = ["custom_gen_stats_multiqc.py", qualimap_paths_fp,
                    q30_data_paths_fp, sub_map_paths_fp,
                    input_type, output_fp]

        with open(expected_result_fp) as f:
            expected_dict = yaml.load(f, Loader=yaml.FullLoader)

        passes = False
        try:
            out_dict = write_custom_multiqc_yaml(arg_list)
            is_file = os.path.isfile(output_fp)
            dicts_match = expected_dict == out_dict
            passes = is_file and dicts_match

            self.assertEqual(expected_dict, out_dict)
            self.assertTrue(is_file)
        finally:
            if passes:
                for fp in [qualimap_paths_fp, q30_data_paths_fp,
                           sub_map_paths_fp, output_fp]:
                    try:
                        os.remove(fp)
                    except OSError:
                        pass
