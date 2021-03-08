# Gathers seq-run_summary.csv files and merges with pangolin lineage table
# Takes 2 arguments:
# 1. Path to directory containing per-sequence-run summary.csv files
# 2. Path of pangolin lineage file


import pandas as pd
from sys import argv
import glob
import os


def format_taxon(x):
    """ Parses sample name from taxon column in pangolin lineage report """
    sampleName = x.split(".")[0]
    sampleName = sampleName.replace("Consensus_", "")
    return sampleName


def create_lineages_summary(arg_list):
    run_summaries_fp = arg_list[1]
    run_summary_suffix = arg_list[2]
    lineage_fp = arg_list[3]
    output_fp = arg_list[4]

    # Merge summaries with a single header
    summaries_pattern = os.path.join(run_summaries_fp,
                                     f"*{run_summary_suffix}")
    merged_summaries_df = pd.concat(
        map(pd.read_csv, glob.glob(summaries_pattern)))

    # Load pangolin file and add a "Sample" column to it
    lineage_df = pd.read_csv(lineage_fp)
    lineage_df["Sample"] = lineage_df["taxon"].apply(format_taxon)

    # Merge tables and write to file
    output_df = merged_summaries_df.merge(
        lineage_df, left_on="Sample", right_on="Sample", how="outer")
    output_df.to_csv(output_fp)


if __name__ == '__main__':
    create_lineages_summary(argv)
