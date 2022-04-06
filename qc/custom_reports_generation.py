import pandas as pd
from sys import argv

SEARCH_ID_KEY = "search_id"
INSPECT_SAMPLE_ID_KEY = "sample_id"
EXCITE_BARCODE_KEY = "excite_barcode"
SEQ_POOL_COMP_ID = "sequenced_pool_component_id"
MULTIQC_SAMPLE_KEY = "Sample"  # dump
MULTIQC_TOT_READS_KEY = "total_reads"  # dump
PANGO_STATUS_KEY = "status"  # rename
RENAMED_STATUS_KEY = "pangolin_status"
PANGO_LINEAGE_KEY = "lineage"
SCORPIO_CALL_KEY = "scorpio_call"
PCT_30_KEY = "Pct_Q30"
UNCAPPED_READS_KEY = "Uncapped_Reads"  # rename
RENAMED_UNCAPPED_READS_KEY = "total_uncapped_reads"
SUB_MAP_PCT_ALIGN_KEY = "Sub_Map_Pct_Aligned"
MAPPED_READS_KEY = "mapped_reads"
SAMPLE_COLLECTION_KEY = "sample_collection_datetime"
SOURCE_KEY = "source"
CALM_SOURCE_VALUE = "CALM"
RTL_SOURCE_VALUE = "RTL"
SPECIMEN_TYPE_KEY = "specimen_type"
WASTEWATER_SPECIMEN_TYPE_VAL = "wastewater"
SEQ_RUN_KEY = "seq_run"
ALL_VAL = "all"
KEY_ORDER = [SEARCH_ID_KEY, INSPECT_SAMPLE_ID_KEY, EXCITE_BARCODE_KEY,
             SEQ_POOL_COMP_ID, SAMPLE_COLLECTION_KEY,
             PANGO_STATUS_KEY, PANGO_LINEAGE_KEY, SCORPIO_CALL_KEY,
             UNCAPPED_READS_KEY, PCT_30_KEY,
             MAPPED_READS_KEY, SUB_MAP_PCT_ALIGN_KEY]
RENAME_DICT = {PANGO_STATUS_KEY: RENAMED_STATUS_KEY,
               UNCAPPED_READS_KEY: RENAMED_UNCAPPED_READS_KEY}
SORT_ORDER = [SAMPLE_COLLECTION_KEY, SCORPIO_CALL_KEY,
              PANGO_LINEAGE_KEY, SEQ_POOL_COMP_ID]


def output_partial_df(curr_df, out_fp_prefix, *args):
    specifier = "_".join(args)
    curr_output_fp = f"{out_fp_prefix}_{specifier}.csv"
    curr_df.to_csv(curr_output_fp, index=False)


def make_bespoke_outputs(arg_list):
    full_fp = arg_list[1]
    out_fp_prefix = arg_list[2]

    full_df = pd.read_csv(full_fp)  # , dtype=str)

    # dump unwanted columns
    full_df.pop(MULTIQC_SAMPLE_KEY)
    full_df.pop(MULTIQC_TOT_READS_KEY)

    # rearrange (and if necessary rename) columns
    for curr_key in reversed(KEY_ORDER):
        renamed_key = RENAME_DICT.get(curr_key, curr_key)
        curr_col = full_df.pop(curr_key)
        full_df.insert(0, renamed_key, curr_col)

    # dump any records with no sequencing data
    not_sequenced_mask = full_df[SEQ_POOL_COMP_ID].isna()
    sequenced_df = full_df[~not_sequenced_mask].copy()

    # sort for end user convenience and to ensure deterministic output order
    sequenced_df.sort_values(by=SORT_ORDER, inplace=True)

    # output everything left
    all_output_fp = f"{out_fp_prefix}_{ALL_VAL}.csv"
    sequenced_df.to_csv(all_output_fp, index=False)

    sources = sequenced_df[SOURCE_KEY].unique()
    special_sources = []
    special_df = None
    for curr_source in sources:
        curr_source_name = curr_source
        curr_mask = sequenced_df[SOURCE_KEY] == curr_source
        # in pandas, can't select records on " == nan"--nans need special case
        if pd.isna(curr_source):
            curr_source_name = "no-source"
            curr_mask = sequenced_df[SOURCE_KEY].isna()
        curr_source_df = sequenced_df[curr_mask].copy()

        if curr_source == CALM_SOURCE_VALUE or curr_source == RTL_SOURCE_VALUE:
            special_sources.append(curr_source_name)
            if special_df is None:
                special_df = curr_source_df
            else:
                special_df = pd.concat([special_df, curr_source_df], ignore_index=True)
        else:
            output_partial_df(curr_source_df, out_fp_prefix, curr_source_name)
    # next

    # handle special case: rtl-relevant sources, split by wastewater or not
    if special_df is not None:
        # sort same as above (need own sort bc is mix of indiv-sorted dfs)
        special_df.sort_values(by=SORT_ORDER, inplace=True)

        wastewater_mask = \
            special_df[SPECIMEN_TYPE_KEY] == WASTEWATER_SPECIMEN_TYPE_VAL
        special_wastewater_df = special_df[wastewater_mask].copy()
        sorted_special_sources = sorted(special_sources)
        output_partial_df(special_wastewater_df, out_fp_prefix,
                          *sorted_special_sources,
                          WASTEWATER_SPECIMEN_TYPE_VAL)

        special_not_wastewater_df = special_df[~wastewater_mask].copy()
        output_partial_df(special_not_wastewater_df, out_fp_prefix,
                          *sorted_special_sources,
                          "not-" + WASTEWATER_SPECIMEN_TYPE_VAL)


if __name__ == '__main__':
    make_bespoke_outputs(argv)
