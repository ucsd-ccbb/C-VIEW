from sys import argv
import pandas as pd

MULTIQC_SEQ_POOL_COMP_ID = "Sample"
METRIC_SEQ_POOL_COMP_ID = "sequenced_pool_component_id"


# Makes a combined table from:
# 1. sequencing run's multiqc_general_stats.txt file
# 2. sequencing run's <seq_run>-acceptance.tsv file
def merge_multiqc_and_acceptance(arg_list):
    output_fp = arg_list[1]
    multiqc_stats_fp = arg_list[2]
    per_sample_metric_fps = arg_list[3:]

    combined_df = None
    for curr_per_sample_metric_fp in per_sample_metric_fps:
        curr_metric_df = pd.read_csv(curr_per_sample_metric_fp,
                                     sep="\t", dtype=str)

        if combined_df is None:
            combined_df = curr_metric_df.copy()
        else:
            # *outer* merge all fields in metric df into summary df
            combined_df = combined_df.merge(curr_metric_df,
                                            on=METRIC_SEQ_POOL_COMP_ID,
                                            how="outer")

    multiqc_df = pd.read_csv(multiqc_stats_fp, sep="\t", dtype=str)

    # remove unwanted leading strings from the multiqc column names
    unwanted_strs = ["QualiMap_mqc-generalstats-qualimap-",
                     "my_genstats_mqc-generalstats-my_genstats-"]
    for curr_str in unwanted_strs:
        multiqc_df.columns = multiqc_df.columns.str.replace(curr_str, "")

    # remove unwanted column
    multiqc_df = multiqc_df.drop("general_error_rate", axis=1)

    # *outer* merge all fields in combined metric df into src df
    output_df = multiqc_df.merge(combined_df,
                                 left_on=MULTIQC_SEQ_POOL_COMP_ID,
                                 right_on=METRIC_SEQ_POOL_COMP_ID,
                                 how="outer")

    output_df.to_csv(output_fp, index=False)


if __name__ == '__main__':
    merge_multiqc_and_acceptance(argv)
