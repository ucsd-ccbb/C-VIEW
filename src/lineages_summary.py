import pandas as pd
from sys import argv
import glob
import os

TAXON_KEY = "taxon"
PANGOLIN_STATUS_KEY = "qc_status"
CVIEW_PANG_STATUS_KEY = "status"
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
ANY_FAIL_KEY = "any_fail"


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
    # Sort here is in reverse order as a poor-man's way to put the newest
    # summary files first, as the production run names start with a date (e.g.
    # 220314 or 220401).  Whatever df is first in the list of dfs to concat
    # appears to set the order of the columns in the resulting concatted df.
    # The older summary files don't necessarily hold all the same metrics as
    # the newer summary files ... and if one of those old files is the "base"
    # df that all the others are concatted to, then the columns that don't
    # exist in the base df but do in the others get tacked on to the very end
    # of the concatted df.  Generally, we'd rather have them in whatever place
    # in the df they are in the latest summary files, since those probably
    # represent the most current use cases.  This is a way to make that likely
    # to happen without adding the fragility of explicitly forcing the column
    # order.
    sorted_matching_fps = sorted(matching_fps, reverse=True)
    matching_dfs = []
    for fp in sorted_matching_fps:
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

    # add a "modded_consensus_seq_name" col;
    # this used to be necessary because before pangolin 4+,
    # pangolin irreversibly munged consensus sequence names.
    # Pangolin 4 removes this munge, but later joins/etc are based on
    # this column, so it is easier to keep it than remove it.
    expanded_df[MOD_CONS_NAME] = expanded_df[CONS_NAME]

    return expanded_df


def add_final_qc_filters_inplace(qc_and_lineage_w_search_ids_df):
    MAPPED_READS_KEY = "mapped_reads"
    MAPPED_LT_50K_KEY = "mapped_reads_lt_50k"
    UNCAPPED_READS_KEY = "Uncapped_Reads"
    UNCAPPED_LT_100K_KEY = "uncapped_reads_lt_100k"
    SUB_MAP_KEY = "Sub_Map_Pct_Aligned"
    SUB_MAP_LT_50_KEY = "sub_map_pct_aligned_lt_50"
    P25_KEY = "P25_Ins_size"
    P25_LT_140_KEY = "p25_ins_size_lt_140"
    PCT_Q30_KEY = "Pct_Q30"
    PCT_Q30_LT_90_KEY = "percent_q30_lt_90"
    COV_GTE_10_KEY = "coverage_gte_10_reads"
    COV_GTE_10_LT_95_KEY = "coverage_gte_10_reads_lt_95"
    MEAN_COV_KEY = "mean_coverage"
    MEAN_COV_LT_500_KEY = "mean_coverage_lt_500"

    keypairs = [(MAPPED_READS_KEY, MAPPED_LT_50K_KEY, lambda a: a < 50000),
                (UNCAPPED_READS_KEY, UNCAPPED_LT_100K_KEY, lambda a: a < 100000),
                (SUB_MAP_KEY, SUB_MAP_LT_50_KEY, lambda a: a < 50),
                (P25_KEY, P25_LT_140_KEY, lambda a: a < 140),
                (PCT_Q30_KEY, PCT_Q30_LT_90_KEY, lambda a: a < 90),
                (COV_GTE_10_KEY, COV_GTE_10_LT_95_KEY, lambda a: a < 0.95),
                (MEAN_COV_KEY, MEAN_COV_LT_500_KEY, lambda a: a < 500)]

    any_fail_values = None
    for curr_tuple in keypairs:
        value_key = curr_tuple[0]
        comp_key = curr_tuple[1]
        comp_func = curr_tuple[2]

        # make a new column holding the comparison result
        qc_and_lineage_w_search_ids_df.loc[:, comp_key] = \
            (qc_and_lineage_w_search_ids_df[value_key].apply(comp_func))

        # build up the series holding values for the new any_fail column
        # by or-ing it together with each new comparison; note that for now,
        # the comparison columns are numpy bools, which treat NaN as False, so
        # if a value is missing, we assume it does NOT fail that value's check
        if any_fail_values is None:
            any_fail_values = qc_and_lineage_w_search_ids_df[comp_key]
        else:
            # NB: the "pd.Series" here is actually unnecessary (any_fail_values
            # is already a Series) but without it pycharm's linter can't tell
            # this is a pandas operation and so gives a spurious fatal error:
            # "Python version 3.9 does not allow writing union types as X | Y"
            any_fail_values = pd.Series(any_fail_values) | \
                              qc_and_lineage_w_search_ids_df[comp_key]

        # NOW convert new comparison column to "boolean"
        # NB: NOT "bool"--"bool" is a numpy type that can't be nullable,
        # whereas "boolean" is a pandas type that *can* be nullable
        qc_and_lineage_w_search_ids_df[comp_key] = \
            qc_and_lineage_w_search_ids_df[comp_key].astype("boolean")

        # if the underlying value being compared is NaN, reset the
        # (now nullable!) comparison result to None instead of False;
        # note this is done AFTER the any_fail comparison, just so we *know*
        # we had no value for these comparisons and thus assumed false for 'em
        # when calculating any_fail
        qc_and_lineage_w_search_ids_df.loc[
            qc_and_lineage_w_search_ids_df[value_key].isna(), comp_key] = None
    # next comparison

    # generally, any_fail should treat comparisons based on missing values
    # as false (not failing) as it does above. The exception is when ALL the
    # values are missing, in which case any_fail should also be missing
    # as we have ZERO info on it:

    # get subset of dataframe containing all (and only) the comparison columns
    comp_keys = [x[1] for x in keypairs]
    comparisons_df = qc_and_lineage_w_search_ids_df.loc[:, comp_keys]

    # make a mask id-ing rows for which ALL the comparison values
    # are None; these should have None for their any_fail value, too
    all_none_mask = comparisons_df.isna().all(axis=1)  # 1 means column
    any_fail_values = any_fail_values.astype("boolean")
    any_fail_values[all_none_mask] = None

    # finally, add the new any_fail column to the output dataframe
    qc_and_lineage_w_search_ids_df.loc[:, ANY_FAIL_KEY] = any_fail_values

    # now set the overall_fail value;
    # some older records may not *have* an n metric;
    # these should NOT be overall fails
    # TODO: how should overall_fail act if any_fail is None?
    # as of now, if n_metric fail is false and any_fail is None,
    # overall_fail will also be None, whereas if the n_metric check
    # fails, overall_fail will be True even if any_fail is None.
    # Note that if n_metric is missing, the n_metric fail will be False here.
    qc_and_lineage_w_search_ids_df.loc[:, OVERALL_FAIL_KEY] = \
        (qc_and_lineage_w_search_ids_df[ANY_FAIL_KEY] |
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
    # and rename the status column to a cview-specific name
    lineage_df = pd.read_csv(lineage_fp, dtype=str)
    lineage_df[MOD_CONS_NAME] = lineage_df[TAXON_KEY]
    lineage_df.rename(columns={PANGOLIN_STATUS_KEY: CVIEW_PANG_STATUS_KEY}, inplace=True)

    # outer merge expanded summaries with lineages (includes lines for
    # both samples that went through lineage calling and those that didn't)
    output_df = expanded_summaries_df.merge(
        lineage_df, left_on=MOD_CONS_NAME, right_on=MOD_CONS_NAME, how="outer")

    # calculate usable_for: nothing, variant, variant_and_epidemiology
    fraction_coverage = output_df['coverage_gte_10_reads'].astype(float)
    # believe checking as below should exclude NAs ...
    passes_pangolin = output_df[CVIEW_PANG_STATUS_KEY] == PASSES_PANG_STATUS_KEY
    gte_70_and_passes_pangolin = passes_pangolin & (fraction_coverage >= 0.70)
    gt_95_and_passes_pangolin = passes_pangolin & (fraction_coverage > 0.95)
    output_df[USABLE_NAME] = NOTHING_VAL
    output_df.loc[gte_70_and_passes_pangolin, USABLE_NAME] = VARIANT_VAL
    output_df.loc[gt_95_and_passes_pangolin, USABLE_NAME] = VARIANT_AND_EP_VAL

    # add additional src filter columns
    add_final_qc_filters_inplace(output_df)

    # sort to ensure deterministic output order
    output_df.sort_values(by=[SEARCH_ID, SEQ_POOL_COMP_ID], inplace=True)

    # there *shouldn't* be any rows in the lineage that aren't in the
    # expanded summaries ... if there are, something is wrong.  Raise
    # an error (but *after* writing the output file, so we have some chance of
    # figuring out what went wrong).  The most likely cause of this issue is
    # some records not being correctly matched up between the two data sources
    # during the outer join above.
    output_df.to_csv(out_summary_fp, index=False)
    if len(output_df) != len(expanded_summaries_df):
        raise ValueError(f"Expected {len(expanded_summaries_df)} rows, "
                         f"got {len(output_df)}")


if __name__ == '__main__':
    create_lineages_summary(argv)
