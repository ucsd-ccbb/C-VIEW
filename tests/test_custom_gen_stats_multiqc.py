import os
from pathlib import Path
import yaml
from tests.test_sarscov2_consensus_acceptance import FileTestCase
from qc.custom_gen_stats_multiqc import write_custom_multiqc_yaml, \
    generate_q30_based_values, _calc_pct_aligned


class CustomGenStatsMultiQcTest(FileTestCase):
    def _make_paths_file(self, output_fname, glob_pattern):
        output_fp = f"{self.test_temp_dir}/{output_fname}"
        relevant_paths = Path(self.dummy_samples_dir).rglob(glob_pattern)
        relevant_fps = [f"{str(p)}\n" for p in relevant_paths]
        sorted_relevant_fps = sorted(relevant_fps)
        with open(output_fp, "w") as f:
            f.writelines(sorted_relevant_fps)
        return output_fp

    def test_generate_q30_based_values_se(self):
        r1 = '/Users/amandabirmingham/Work/Repositories/' \
             'covid_sequencing_analysis_pipeline/tests/data/dummy/' \
             'dummy_2021-02-08-ARTIC_samples/SEARCH-5329__LIBPLATE1__' \
             'A01__001002003004/SEARCH-5329__LIBPLATE1__A01__' \
             '001002003004_R1_q30_data.txt'
        s1 = [2440812, 2319178]
        expected_name = 'SEARCH-5329__LIBPLATE1__A01__001002003004'
        out_name, out_pctQ30_dict, out_uncapped_reads_dict = \
            generate_q30_based_values({}, r1, s1)

        self.assertEqual(expected_name, out_name)
        self.assertDictEqual({'Pct >=Q30': 0.95}, out_pctQ30_dict)
        self.assertDictEqual({'Uncapped Reads': 2440812}, out_uncapped_reads_dict)

    def test_generate_q30_based_values_se_zero(self):
        r1 = '/Users/amandabirmingham/Work/Repositories/' \
             'covid_sequencing_analysis_pipeline/tests/data/dummy/' \
             'dummy_2021-02-08-ARTIC_samples/SEARCH-5329__LIBPLATE1__' \
             'A01__001002003004/SEARCH-5329__LIBPLATE1__A01__' \
             '001002003004_R1_q30_data.txt'
        s1 = [0, 0]
        expected_name = 'SEARCH-5329__LIBPLATE1__A01__001002003004'
        out_name, out_pctQ30_dict, out_uncapped_reads_dict = \
            generate_q30_based_values({}, r1, s1)

        self.assertEqual(expected_name, out_name)
        self.assertDictEqual({'Pct >=Q30': 'NA'}, out_pctQ30_dict)
        self.assertDictEqual({'Uncapped Reads': 0}, out_uncapped_reads_dict)

    def test_generate_q30_based_values_pe(self):
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
        expected_name = 'SEARCH-5329__LIBPLATE1__A01__001002003004'
        out_name, out_pctQ30_dict, out_uncapped_reads_dict = \
            generate_q30_based_values({}, r1, s1, r2, s2)

        self.assertEqual(expected_name, out_name)
        self.assertDictEqual({'Pct >=Q30': 94.271}, out_pctQ30_dict)
        self.assertDictEqual({'Uncapped Reads': 4881624},
                             out_uncapped_reads_dict)

    def test_generate_q30_based_values_pe_zero(self):
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
        expected_name = 'SEARCH-5329__LIBPLATE1__A01__001002003004'
        out_name, out_pctQ30_dict, out_uncapped_reads_dict = \
            generate_q30_based_values({}, r1, s1, r2, s2)

        self.assertEqual(expected_name, out_name)
        self.assertDictEqual({'Pct >=Q30': 93.526}, out_pctQ30_dict)
        self.assertDictEqual({'Uncapped Reads': 2440812},
                             out_uncapped_reads_dict)

    def test_generate_q30_based_values_pe_zero_both(self):
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
        expected_name = 'SEARCH-5329__LIBPLATE1__A01__001002003004'
        out_name, out_pctQ30_dict, out_uncapped_reads_dict = \
            generate_q30_based_values({}, r1, s1, r2, s2)

        self.assertEqual(expected_name, out_name)
        self.assertDictEqual({'Pct >=Q30': 'NA'}, out_pctQ30_dict)
        self.assertDictEqual({'Uncapped Reads': 0},
                             out_uncapped_reads_dict)

    def test__calc_pct_aligned(self):
        real_out = _calc_pct_aligned([5, 5])
        self.assertEqual(50.0, real_out)

    def test__calc_pct_aligned_zero_both(self):
        real_out = _calc_pct_aligned([0, 0])
        self.assertEqual('NA', real_out)

    def test_write_custom_multiqc_yaml(self):
        self.maxDiff = None

        qualimap_paths_fp = self._make_paths_file(
            "qualimapReport_paths.txt", 'qualimapReport.html')
        q30_data_paths_fp = self._make_paths_file(
            "q30_data_paths.txt", '*q30_data.txt')
        sub_map_paths_fp = self._make_paths_file(
            "sub_map_paths.txt", '*_subsampled_mapping_stats.tsv')
        output_fp = f"{self.test_temp_dir}/" \
                    f"temp_test_cust_gen_stats_multiqc.yaml"
        expected_result_fp = f"{self.dummy_dir}/" \
                             f"dummy_2021-02-08-ARTIC-multiqc_custom_gen_" \
                             f"stats.yaml"

        arg_list = ["custom_gen_stats_multiqc.py", qualimap_paths_fp,
                    q30_data_paths_fp, sub_map_paths_fp, "pe", output_fp]

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
