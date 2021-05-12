import pandas as pd
from sys import argv

SAMPLE_KEY = "bam_name"
COVERAGE_KEY = "bjorn_coverage"
AVG_DEPTH_KEY = "bjorn_avg_depth"
MIN_KEY = "bjorn_min_coverage"
MAX_KEY = "bjorn_max_coverage"
ZERO_DEPTH = "bjorn_num_w_zero_coverage"
COVERAGE_COL_NAMES = [SAMPLE_KEY, COVERAGE_KEY, AVG_DEPTH_KEY, MIN_KEY,
                      MAX_KEY, ZERO_DEPTH]

SEQ_POOL_COMP_ID = "sequenced_pool_component_id"
BJORN_SEQ_POOL_COMP_ID = f"bjorn_{SEQ_POOL_COMP_ID}"


def integrate_bjorn_coverage(arg_list):
    in_summary_fp = arg_list[1]
    coverage_fp = arg_list[2]
    out_summary_fp = arg_list[3]

    in_summary_df = pd.read_csv(in_summary_fp, sep=",", dtype=str)
    # override the header names in the coverage file
    coverage_df = pd.read_csv(coverage_fp, sep="\t", dtype=str,
                              names=COVERAGE_COL_NAMES, header=0)
    coverage_df[BJORN_SEQ_POOL_COMP_ID] = coverage_df[SAMPLE_KEY].str.replace(
        ".trimmed.sorted.bam", "", regex=False)

    # TODO: some kind of check here to make sure the same ids are in each df?

    out_summary_df = in_summary_df.merge(
        coverage_df, left_on=SEQ_POOL_COMP_ID,
        right_on=BJORN_SEQ_POOL_COMP_ID,
        how="outer")
    out_summary_df.pop(BJORN_SEQ_POOL_COMP_ID)
    out_summary_df.to_csv(out_summary_fp, index=False)


if __name__ == '__main__':
    integrate_bjorn_coverage(argv)
