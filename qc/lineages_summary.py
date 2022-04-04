import pandas as pd
from sys import argv
import glob
import os

TAXON_KEY = "taxon"
PANGOLIN_STATUS_KEY = "qc_status"
PASSES_PANG_STATUS_KEY = "passed_qc"
SAMPLE_NAME = "Sample"
SEARCH_ID = "search_id"
SEQ_POOL_COMP_ID = "sequenced_pool_component_id"
CONS_NAME = "consensus_seq_name"
MOD_CONS_NAME = "modded_consensus_seq_name"
USABLE_NAME = "usable_for"
NOTHING_VAL = "nothing"
VARIANT_VAL = "variant"
VARIANT_AND_EP_VAL = "variant_and_epidemiology"
IS_HIST_OR_REF = "is_hist_or_ref"
OVERALL_FAIL_KEY = "overall_fail"


# recreate un-reversable pangolin name munge;
# code pulled and very slightly modded from
# https://github.com/cov-lineages/pangolin/blob/
#  1763ac04da0dff41bd778cfa72f41a361457d81d/pangolin/command.py#L144-L147
def _perform_pangolin_name_munge(seq_name):
    mod_seq_name = seq_name.replace(' ', '_')
    if "," in mod_seq_name:
        mod_seq_name = mod_seq_name.replace(",", "_")
    return mod_seq_name


def merge_summaries(run_summaries_fp, run_summary_suffix):
    # Merge summaries with a single header
    summaries_pattern = os.path.join(run_summaries_fp,
                                     f"*{run_summary_suffix}")
    matching_fps = glob.glob(summaries_pattern)
    matching_dfs = []
    for fp in matching_fps:
        curr_df = pd.read_csv(fp)  # , dtype=str)
        matching_dfs.append(curr_df)
    merged_summaries_df = pd.concat(matching_dfs)
    return merged_summaries_df


def expand_with_added_fa_names(merged_summaries_df, added_fa_names_fp):
    fas_col_name = "fasta_id"

    # read this in as a tsv (even though just one column) bc some of these
    # names have commas in them so can't read as csv ...
    added_fastq_ids_df = pd.read_csv(added_fa_names_fp, sep="\t", dtype=str)

    # add a column marking these as historic or reference
    added_fastq_ids_df[IS_HIST_OR_REF] = True

    # make a column to hold the sequenced pool component id;
    # for these added ids, punt to this being the same as the fas name
    added_fastq_ids_df[SEQ_POOL_COMP_ID] = added_fastq_ids_df[fas_col_name]

    # also copy it into "Sample" column for now, just so it has something there
    added_fastq_ids_df[SAMPLE_NAME] = added_fastq_ids_df[fas_col_name]

    # rename the "fasta_id" column "consensus_seq_name"
    added_fastq_ids_df.rename(
        columns={fas_col_name: CONS_NAME}, inplace=True)

    expanded_df = merged_summaries_df.merge(
        added_fastq_ids_df,
        left_on=[CONS_NAME, SAMPLE_NAME, SEQ_POOL_COMP_ID],
        right_on=[CONS_NAME, SAMPLE_NAME, SEQ_POOL_COMP_ID], how="outer")
    expanded_df.fillna({CONS_NAME: ''}, inplace=True)
    expanded_df.fillna({IS_HIST_OR_REF: False}, inplace=True)

    # add a "modded_consensus_seq_name" col
    # by modifying the consensus name column contents according to
    # pangolin's irreversible munge rules
    expanded_df[MOD_CONS_NAME] = \
        expanded_df[CONS_NAME].apply(_perform_pangolin_name_munge)

    return expanded_df


def add_final_qc_filters_inplace(qc_and_lineage_w_search_ids_df):
    qc_and_lineage_w_search_ids_df.loc[:, "mapped_reads_lt_50k"] = \
        ((qc_and_lineage_w_search_ids_df["mapped_reads"] < 50000) |
         qc_and_lineage_w_search_ids_df["mapped_reads"].isna())

    qc_and_lineage_w_search_ids_df.loc[:, "uncapped_reads_lt_100k"] = \
        ((qc_and_lineage_w_search_ids_df["Uncapped_Reads"] < 50000) |
         qc_and_lineage_w_search_ids_df["Uncapped_Reads"].isna())

    qc_and_lineage_w_search_ids_df.loc[:, "sub_map_pct_aligned_lt_50"] = \
        ((qc_and_lineage_w_search_ids_df["Sub_Map_Pct_Aligned"] < 50) |
         qc_and_lineage_w_search_ids_df["Sub_Map_Pct_Aligned"].isna())

    qc_and_lineage_w_search_ids_df.loc[:, "p25_ins_size_lt_140"] = \
        ((qc_and_lineage_w_search_ids_df["P25_Ins_size"] < 140) |
         qc_and_lineage_w_search_ids_df["P25_Ins_size"].isna())

    qc_and_lineage_w_search_ids_df.loc[:, "percent_q30_lt_90"] = \
        ((qc_and_lineage_w_search_ids_df["Pct_Q30"] < 90) |
         qc_and_lineage_w_search_ids_df["Pct_Q30"].isna())

    qc_and_lineage_w_search_ids_df.loc[:, "coverage_gte_10_reads_lt_95"] = \
        ((qc_and_lineage_w_search_ids_df["coverage_gte_10_reads"] < 0.95) |
         qc_and_lineage_w_search_ids_df["coverage_gte_10_reads"].isna())

    qc_and_lineage_w_search_ids_df.loc[:, "mean_coverage_lt_500"] = \
        ((qc_and_lineage_w_search_ids_df["mean_coverage"] < 500) |
         qc_and_lineage_w_search_ids_df["mean_coverage"].isna())

    qc_and_lineage_w_search_ids_df.loc[:, "any_fail"] = \
        (qc_and_lineage_w_search_ids_df["mapped_reads_lt_50k"] |
         qc_and_lineage_w_search_ids_df["uncapped_reads_lt_100k"] |
         qc_and_lineage_w_search_ids_df["sub_map_pct_aligned_lt_50"] |
         qc_and_lineage_w_search_ids_df["p25_ins_size_lt_140"] |
         qc_and_lineage_w_search_ids_df["percent_q30_lt_90"] |
         qc_and_lineage_w_search_ids_df["coverage_gte_10_reads_lt_95"] |
         qc_and_lineage_w_search_ids_df["mean_coverage_lt_500"])

    # some older records may not *have* an n metric;
    # these should NOT be overall fails
    qc_and_lineage_w_search_ids_df.loc[:, OVERALL_FAIL_KEY] = \
        (qc_and_lineage_w_search_ids_df["any_fail"] |
         (qc_and_lineage_w_search_ids_df["n_metric"] > 19))


def create_lineages_summary(arg_list):
    added_fa_names_fp = arg_list[1]
    run_summaries_fp = arg_list[2]
    run_summary_suffix = arg_list[3]
    lineage_fp = arg_list[4]
    out_summary_fp = arg_list[5]

    merged_summaries_df = merge_summaries(run_summaries_fp, run_summary_suffix)
    expanded_summaries_df = expand_with_added_fa_names(
        merged_summaries_df, added_fa_names_fp)

    # Load pangolin file to a dataframe and
    # copy the "taxon" column into a new col named "modded_consensus_seq_name"
    lineage_df = pd.read_csv(lineage_fp, dtype=str)
    lineage_df[MOD_CONS_NAME] = lineage_df[TAXON_KEY]

    # outer merge expanded summaries with lineages (includes lines for
    # both samples that went through lineage calling and those that didn't)
    output_df = expanded_summaries_df.merge(
        lineage_df, left_on=MOD_CONS_NAME, right_on=MOD_CONS_NAME, how="outer")

    # calculate usable_for: nothing, variant, variant_and_epidemiology
    fraction_coverage = output_df['coverage_gte_10_reads'].astype(float)
    # believe checking as below should exclude NAs ...
    passes_pangolin = output_df[PANGOLIN_STATUS_KEY] == PASSES_PANG_STATUS_KEY
    gte_70_and_passes_pangolin = passes_pangolin & (fraction_coverage >= 0.70)
    gt_95_and_passes_pangolin = passes_pangolin & (fraction_coverage > 0.95)
    output_df[USABLE_NAME] = NOTHING_VAL
    output_df.loc[gte_70_and_passes_pangolin, USABLE_NAME] = VARIANT_VAL
    output_df.loc[gt_95_and_passes_pangolin, USABLE_NAME] = VARIANT_AND_EP_VAL

    # add additional qc filter columns
    add_final_qc_filters_inplace(output_df)

    # sort to ensure deterministic output order
    output_df.sort_values(by=[SEARCH_ID, SEQ_POOL_COMP_ID], inplace=True)

    # there *shouldn't* be any rows in the lineage that aren't in the
    # expanded summaries ... if there are, something is wrong.  Raise
    # an error (but *after* writing the output file, so we have some chance of
    # figuring out what went wrong).
    output_df.to_csv(out_summary_fp, index=False)
    if len(output_df) != len(expanded_summaries_df):
        raise ValueError(f"Expected {len(expanded_summaries_df)} rows, "
                         f"got {len(output_df)}")


if __name__ == '__main__':
    create_lineages_summary(argv)
