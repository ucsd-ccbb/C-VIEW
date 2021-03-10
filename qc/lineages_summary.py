# Gathers seq-run_summary.csv files and merges with pangolin lineage table
# Takes 2 arguments:
# 1. Path to directory containing per-sequence-run summary.csv files
# 2. Path of pangolin lineage file


import pandas as pd
from sys import argv
import glob
import os

SAMPLE_ID = "sample_id"


def format_taxon(x):
    """ Parses sample name from taxon column in pangolin lineage report """
    sampleName = x.split(".")[0]
    sampleName = sampleName.replace("Consensus_", "")
    return sampleName


def merge_summaries(run_summaries_fp, run_summary_suffix):
    # Merge summaries with a single header
    summaries_pattern = os.path.join(run_summaries_fp,
                                     f"*{run_summary_suffix}")
    matching_fps = glob.glob(summaries_pattern)
    matching_dfs = []
    for fp in matching_fps:
        curr_df = pd.read_csv(fp, dtype=str)
        matching_dfs.append(curr_df)
    merged_summaries_df = pd.concat(matching_dfs)
    return merged_summaries_df


def expand_with_added_fa_names(merged_summaries_df, added_fa_names_fp):
    # read this in as a tsv (even though just one column) bc some of these
    # names have commas in them so can't read as csv ...
    added_fastq_ids_df = pd.read_csv(added_fa_names_fp, sep="\t", dtype=str)
    # NB: per seq_run_summary.py, "Sample" is "fastq_id" (or, here, "fasta_id")
    added_fastq_ids_df.rename(columns={"fasta_id": "Sample"}, inplace=True)
    added_fastq_ids_df[SAMPLE_ID] = added_fastq_ids_df["Sample"]
    expanded_df = merged_summaries_df.merge(added_fastq_ids_df,
                                            left_on=["Sample", SAMPLE_ID],
                                            right_on=["Sample", SAMPLE_ID],
                                            how="outer")
    return expanded_df


def generate_metadata_df(expanded_summaries_df, lineage_df):
    # RIGHT merge expanded summaries with lineages (include ONLY lines
    # for samples that went through lineage calling)
    metadata_df = expanded_summaries_df.merge(
        lineage_df, left_on="Sample", right_on="Sample", how="right")

    # rearrange columns--want "sample_id" (which will hopefully be SEARCH
    # id where that exists) as first column
    # shift column 'Name' to first position
    first_column = metadata_df.pop(SAMPLE_ID)
    metadata_df.insert(0, SAMPLE_ID, first_column)
    return metadata_df


def create_lineages_summary_and_metadata(arg_list):
    added_fa_names_fp = arg_list[1]
    run_summaries_fp = arg_list[2]
    run_summary_suffix = arg_list[3]
    lineage_fp = arg_list[4]
    out_summary_fp = arg_list[5]
    out_metadata_fp = arg_list[6]

    merged_summaries_df = merge_summaries(run_summaries_fp, run_summary_suffix)
    expanded_summaries_df = expand_with_added_fa_names(
        merged_summaries_df, added_fa_names_fp)
    expanded_summaries_copy_df = expanded_summaries_df.copy(deep=True)

    # Load pangolin file and add a "Sample" column to it
    lineage_df = pd.read_csv(lineage_fp, dtype=str)
    lineage_df["Sample"] = lineage_df["taxon"].apply(format_taxon)

    # outer merge expanded summaries with lineages (includes lines for
    # both samples that went through lineage calling and those that didn't)
    output_df = expanded_summaries_df.merge(
        lineage_df, left_on="Sample", right_on="Sample", how="outer")
    output_df.to_csv(out_summary_fp, index=False)

    # RIGHT merge expanded summaries with lineages (include ONLY lines
    # for samples that went through lineage calling because those are also
    # the only ones that go to tree building.)
    # NB that empress metadata files must be tsv
    metadata_df = generate_metadata_df(expanded_summaries_copy_df, lineage_df)
    metadata_df.to_csv(out_metadata_fp, sep='\t', index=False)


if __name__ == '__main__':
    create_lineages_summary_and_metadata(argv)
