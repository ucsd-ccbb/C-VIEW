import pandas as pd
from sys import argv
import fastaparser


SEARCH_ID_KEY = "search_id"
CONS_NAME = "consensus_seq_name"
USABLE_NAME = "usable_for"
VARIANT_VAL = "variant"
VARIANT_AND_EP_VAL = "variant_and_epidemiology"
IS_HIST_OR_REF = "is_hist_or_ref"
INDEX_COL_NAME = "Unnamed: 0"


def generate_base_empress_df(raw_empress_df):
    # include ONLY samples that went through lineage calling
    has_pangolin_status = raw_empress_df["status"].notna()
    empress_df = raw_empress_df[has_pangolin_status]

    # rearrange columns--want CONS_NAME as first column to match up
    # with the fas record names, which are used as the tree node names
    # in the tree file
    # shift column 'consensus_seq_name' to first position
    first_column = empress_df.pop(CONS_NAME)
    empress_df.insert(0, CONS_NAME, first_column)
    return empress_df


def winnow_fasta(fasta_fp, base_empress_df, out_loose_fp, out_stringent_fp):
    with open(out_loose_fp, 'w') as loose_file:
        with open(out_stringent_fp, 'w') as stringent_file:
            with open(fasta_fp) as fasta_file:
                parser = fastaparser.Reader(fasta_file)
                for seq in parser:
                    # match header to consensus_seq_name column of metadata
                    # and get the usable_for value for it
                    seq_metadata_df = base_empress_df.loc[
                        base_empress_df[CONS_NAME] == seq.id]

                    if len(seq_metadata_df) == 0:
                        continue
                    elif len(seq_metadata_df) > 1:
                        raise ValueError(f"More than one metadata row with"
                                         f"consensus sequence name "
                                         f"'{seq.header}' found")
                    else:
                        seq_usable_for = seq_metadata_df.loc[
                                         :, USABLE_NAME].iat[0]
                        fasta_str = seq.formatted_fasta() + '\n'
                        if seq_usable_for == VARIANT_VAL:
                            loose_file.write(fasta_str)
                        elif seq_usable_for == VARIANT_AND_EP_VAL:
                            stringent_file.write(fasta_str)
                        # ignore any other cases


def prep_files_for_tree_building(arg_list):
    qc_and_lineage_fp = arg_list[1]
    metadata_fp = arg_list[2]  # May be None
    full_fasta_fp = arg_list[3]
    out_loose_fasta_fp = arg_list[4]
    out_stringent_fasta_fp = arg_list[5]
    out_empress_fp = arg_list[6]
    out_empress_var_fp = arg_list[7]
    out_empress_var_and_ep_fp = arg_list[8]

    qc_and_lineage_df = pd.read_csv(qc_and_lineage_fp)  # , dtype=str)

    if metadata_fp is not None:
        metadata_df = pd.read_csv(metadata_fp, dtype=str)
        # if input metadata file has an unnamed first (index) column, drop it
        if metadata_df.columns[0] == INDEX_COL_NAME:
            metadata_df.pop(INDEX_COL_NAME)

        metadata_w_search_id = metadata_df[SEARCH_ID_KEY].notna()
        metadata_w_search_ids_df = metadata_df[metadata_w_search_id].copy()
        raw_empress_df = qc_and_lineage_df.merge(
            metadata_w_search_ids_df, on=SEARCH_ID_KEY, how="outer")
    else:
        raw_empress_df = qc_and_lineage_df.copy()

    base_empress_df = generate_base_empress_df(raw_empress_df)
    # NB that empress metadata files must be tsv
    base_empress_df.to_csv(out_empress_fp, sep='\t', index=False)

    usable_is_variant = base_empress_df[USABLE_NAME] == VARIANT_VAL
    usable_is_var_and_ep = base_empress_df[USABLE_NAME] == VARIANT_AND_EP_VAL

    # these are cumulative: the var_empress_df
    # contains all records with a "usable_for" value of "variant" OR
    # "variant_and_epidemiology"--OR the reference and historical records
    # NB: Do NOT let PEP make you change this from == True to is True--
    # "is True" won't work!
    is_hist_or_ref = base_empress_df[IS_HIST_OR_REF] == True  # noqa E712
    var_empress_mask = (usable_is_variant | usable_is_var_and_ep |
                        is_hist_or_ref)
    var_empress_df = base_empress_df[var_empress_mask].copy()
    var_empress_df.to_csv(out_empress_var_fp, sep='\t', index=False)

    var_and_ep_empress_df = var_empress_df.loc[~usable_is_variant].copy()
    var_and_ep_empress_df.to_csv(out_empress_var_and_ep_fp,
                                 sep='\t', index=False)

    # the output fastas are NOT cumulative
    winnow_fasta(full_fasta_fp, base_empress_df,
                 out_loose_fasta_fp, out_stringent_fasta_fp)


if __name__ == '__main__':
    prep_files_for_tree_building(argv)
