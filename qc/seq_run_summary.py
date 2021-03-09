from sys import argv
import pandas as pd


# Makes a combined table from:
# 1. sequencing run's multiqc_general_stats.txt file
# 2. sequencing run's <seq_run>-acceptance.tsv file
def merge_multiqc_and_acceptance(arg_list):
    multiqc_stats_fp = arg_list[1]
    sum_acceptance_fp = arg_list[2]
    output_fp = arg_list[3]

    multiqc_df = pd.read_csv(multiqc_stats_fp, sep="\t")
    accept_df = pd.read_csv(sum_acceptance_fp, sep="\t")

    # remove unwanted leading strings from the multiqc column names
    unwanted_strs = ["QualiMap_mqc-generalstats-qualimap-",
                     "my_genstats_mqc-generalstats-my_genstats-"]
    for curr_str in unwanted_strs:
        multiqc_df.columns = multiqc_df.columns.str.replace(curr_str, "")

    # remove unwanted column
    multiqc_df = multiqc_df.drop("general_error_rate", axis=1)

    # *outer* merge all fields in acceptance df into multiqc df
    combined = multiqc_df.merge(accept_df, left_on="Sample",
                                right_on="fastq_id", how="outer")
    combined.to_csv(output_fp, index=False)


if __name__ == '__main__':
    merge_multiqc_and_acceptance(argv)
