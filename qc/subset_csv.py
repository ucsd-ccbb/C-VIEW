# Filters an input csv file by stripping out any (non-header) lines that
# do not have one of the supplied allowed values in the specified column.
# Note that this only works for discrete values--
# CANNOT be used to say, e.g.,
# "keep any lines where 0 < gc_content_percent < 48.5"

from sys import argv
import pandas as pd

LINES_TO_FILE = "filtered_lines"
ACCEPTED_CONS_FNAMES = "accepted_cons_fnames"
INDEL_FLAGGED_CONS_FNAMES = "indel_flagged_cons_fnames"
NO_NANS_CONS_FNAMES = "not_na_cons_fnames"
FILTER_TYPES = [LINES_TO_FILE, ACCEPTED_CONS_FNAMES, INDEL_FLAGGED_CONS_FNAMES,
                NO_NANS_CONS_FNAMES]


def subset_csv(csv_fp, csv_col_name, allowed_vals_str):
    orig_df = _read_df(csv_fp)

    # Filter by keeping only records that have an allowed value in the
    # specified column
    allowed_vals = allowed_vals_str.split(",")
    filtered_df = orig_df[orig_df[csv_col_name].isin(allowed_vals)]
    return filtered_df


def get_consensus_fnames_w_allowed_vals(
        csv_fp, csv_col_name, allowed_vals_str, dir_prefix=None):
    filtered_df = subset_csv(csv_fp, csv_col_name, allowed_vals_str)
    return _get_consensus_fnames(filtered_df, dir_prefix)


def get_consensus_fnames_wo_nans(csv_fp, csv_col_name, dir_prefix=None):
    orig_df = _read_df(csv_fp)

    # Filter by keeping only records that aren't na in the specified column
    filtered_df = orig_df[orig_df[csv_col_name].notna()]
    return _get_consensus_fnames(filtered_df, dir_prefix)


def _read_df(csv_fp):
    # read csv file into dataframe
    return pd.read_csv(csv_fp, dtype=str)


def _get_consensus_fnames(filtered_df, dir_prefix):
    fastq_ids = filtered_df['fastq_id'].to_list()
    if dir_prefix is not None:
        dir_prefix = dir_prefix + "/"
    else:
        dir_prefix = ""
    consensus_fnames = [f"{dir_prefix}{x}.trimmed.sorted.pileup.consensus.fa"
                        for x in fastq_ids]
    return " ".join(consensus_fnames)


def filter_csv(arg_list):
    result_str = ""
    result_code = 1

    # file path of the csv file to filter
    csv_fp = arg_list[1]

    # filter type
    filter_type = arg_list[2]

    if filter_type == LINES_TO_FILE:
        csv_col_name = arg_list[3]  # column on which to filter
        # comma-separated list of allowed values (**without spaces**)
        allowed_vals_str = arg_list[4]
        output_fp = arg_list[5]

        filtered_df = subset_csv(csv_fp, csv_col_name, allowed_vals_str)
        filtered_df.to_csv(output_fp, index=False)
        result_code = 0
    else:
        dir_prefix = arg_list[3] if len(arg_list) == 4 else None

        if filter_type == ACCEPTED_CONS_FNAMES:
            result_str = get_consensus_fnames_w_allowed_vals(
                csv_fp, "is_accepted", "True", dir_prefix)
            result_code = 0
        elif filter_type == INDEL_FLAGGED_CONS_FNAMES:
            result_str = get_consensus_fnames_w_allowed_vals(
                csv_fp, "indels_flagged", "True", dir_prefix)
            result_code = 0
        elif filter_type == NO_NANS_CONS_FNAMES:
            result_str = get_consensus_fnames_wo_nans(
                csv_fp, "consensus_seq_name", dir_prefix)
            result_code = 0
        else:
            result_str = f"Unrecognized filter type '{filter_type}'"

    return result_str, result_code


if __name__ == '__main__':
    # out_loc = stderr
    out_str, out_code = filter_csv(argv)

    # if out_code == 0:
    #     out_loc = stdout

    print(out_str)  # , file=out_loc)
    exit(out_code)
