import os
from pathlib import Path
from wand.image import Image
from unittest import TestCase
from qc.samtools_depth_plots import make_plots


class SamtoolsDepthPlotsTest(TestCase):
    def test_make_plots(self):
        test_data_path = os.path.abspath("./data")
        samples_path = f"{test_data_path}/2021-02-08-ARTIC_samples"
        results_path = f"{test_data_path}/2021-02-08-ARTIC_quality_control"
        qc_data_path = f"{test_data_path}/qc"
        expected_depth_lines_fp = f"{results_path}/2021-02-08-ARTIC-depth_lineplot.pdf"
        expected_depth_violin_fp = f"{results_path}/2021-02-08-ARTIC-depth_violin.pdf"

        depth_fps = [str(p) for p in Path(samples_path).rglob('*.depth.txt')]
        out_depth_lines_fp = f"{qc_data_path}/temp_test_depth_lines.pdf"
        out_depth_violin_fp = f"{qc_data_path}/temp_test_depth_violin.pdf"
        arg_list = ["samtools_depth_plots.py",
                    out_depth_lines_fp, out_depth_violin_fp]
        arg_list.extend(depth_fps)

        try:
            make_plots(arg_list)

            # PDFs are vectorial, so we need to set a resolution when
            # converting to an image
            actual_lines = Image(filename=out_depth_lines_fp, resolution=150)
            with Image(filename=expected_depth_lines_fp,
                       resolution=150) as expected_line:
                diff_lines = actual_lines.compare(
                    expected_line, metric='root_mean_square')
                self.assertLess(diff_lines[1], 0.01)

            actual_violin = Image(filename=out_depth_violin_fp, resolution=150)
            with Image(filename=expected_depth_violin_fp,
                       resolution=150) as expected_violin:
                diff_violin = actual_violin.compare(
                    expected_violin, metric='root_mean_square')
                self.assertLess(diff_violin[1], 0.01)
        finally:
            for out_fp in [out_depth_lines_fp, out_depth_violin_fp]:
                try:
                    os.remove(out_fp)
                except OSError:
                    pass





