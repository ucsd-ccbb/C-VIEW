import pandas as pd
from sys import argv
import fastaparser


SEARCH_ID_KEY = "search_id"
SEQ_POOL_COMP_ID = "sequenced_pool_component_id"
CONS_NAME = "consensus_seq_name"
USABLE_NAME = "usable_for"
NOTHING_VAL = "nothing"
VARIANT_VAL = "variant"
VARIANT_AND_EP_VAL = "variant_and_epidemiology"
CONSENSUS_S3 = "consensus_s3"
IS_HIST_OR_REF = "is_hist_or_ref"

BJORN_COL_NAMES = ["Sample ID", "SEARCH SampleID", "Ready for release?",
                   "New sequences ready for release", "Released?",
                   "Submitter", "FASTA filename", "Virus name", "Type",
                   "Passage details/history", "Collection date", "Location",
                   "Additional location information", "Host",
                   "Additional host information", "Gender", "Patient age",
                   "Patient status", "Specimen source", "Outbreak",
                   "Last vaccinated", "Treatment", "Sequencing technology",
                   "Assembly method", "Coverage", "Originating lab",
                   "Address", "Sample ID given by the sample provider",
                   "Submitting lab", "Address.1",
                   "Sample ID given by the submitting laboratory", "Authors",
                   "Comment", "Comment Icon", "Released",
                   "Sequenced Pool Component Id"]


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


def generate_full_metadata(qc_and_lineages_df, inspect_metadata_df):
    records_w_search_ids = qc_and_lineages_df[SEARCH_ID_KEY].notna()
    qc_and_lineage_w_search_ids_df = qc_and_lineages_df[
        records_w_search_ids].copy()

    add_final_qc_filters_inplace(qc_and_lineage_w_search_ids_df)

    full_df = qc_and_lineage_w_search_ids_df.merge(
        inspect_metadata_df, on=SEARCH_ID_KEY, how="outer")
    return full_df


def filter_metadata_for_bjorn(full_df):
    has_no_search_id = full_df[SEARCH_ID_KEY].isna()
    has_no_specimen_type = full_df["specimen_type"].isna()
    has_no_consensus_seq = full_df["consensus_seq_name"].isna()
    is_control = full_df["source"] == "Control Sample"

    is_human = full_df["subject_species"] == "Human"
    has_gender = full_df["subject_gender"].isin(["Male", "Female", "Unknown"])
    has_positive_age = ~full_df["subject_age"].isna() & full_df[
        "subject_age"] > 0
    has_demographics = has_gender & has_positive_age
    is_human_wo_demographics = is_human & ~has_demographics

    is_excluded = (has_no_search_id | has_no_specimen_type |
                   has_no_consensus_seq | is_control |
                   is_human_wo_demographics)
    filtered_df = full_df[~is_excluded].copy()
    return filtered_df


def generate_bjorn_df(filtered_df):
    filtered_df[["sample_collection_date", "sample_collection_time"]] = \
        filtered_df["sample_collection_datetime"].str.split(expand=True)
    filtered_df[["sample_collection_year",
                 "sample_collection_month",
                 "sample_collection_day"]] = filtered_df[
        "sample_collection_date"].str.split("-", expand=True)

    output_df = pd.DataFrame()
    output_df.loc[:, 'sample_id'] = filtered_df['sample_id']
    output_df.loc[:, 'search_id'] = filtered_df['search_id']

    release_mask = filtered_df["any_fail"] == False  # noqa 712
    output_df.loc[:, "ready_for_release"] = "No"
    output_df.loc[release_mask, "ready_for_release"] = "Yes"

    output_df.loc[:, "new_seqs_ready_for_release"] = "Yes"
    output_df.loc[:, "released"] = ""
    output_df.loc[:, "submitter"] = ""
    output_df.loc[:, "fasta_fname"] = ""
    output_df.loc[:, "virus_name"] = "hCoV-19/USA/" + output_df["search_id"] \
                                     + "/" \
                                     + filtered_df["sample_collection_year"]
    output_df.loc[:, "type"] = "betacoronavirus"
    output_df.loc[:, "passage_details"] = filtered_df["passage_details"]
    output_df.loc[:, "sample_collection_date"] = filtered_df[
        "sample_collection_date"]
    output_df.loc[:, "sample_collection_location"] = filtered_df[
        "sample_collection_location"]
    output_df.loc[:, "zip"] = filtered_df["zip"]
    output_df.loc[:, "subject_species"] = filtered_df["subject_species"]
    output_df.loc[:, "addl_host_info"] = "N/A"
    output_df.loc[:, "subject_gender"] = filtered_df["subject_gender"]
    output_df.loc[:, "subject_age"] = filtered_df["subject_age"]
    output_df.loc[:, "patient_status"] = "N/A"
    output_df.loc[:, "specimen_type"] = filtered_df["specimen_type"]
    output_df.loc[:, "outbreak"] = "N/A"
    output_df.loc[:, "last_vaccinated"] = "N/A"
    output_df.loc[:, "treatment"] = "N/A"
    output_df.loc[:, "sequencing_tech"] = "Illumina"
    output_df.loc[:, "assembly_method"] = filtered_df[
        "assembly_method"].str.replace(" version", "")
    output_df.loc[:, "bjorn_coverage"] = filtered_df["bjorn_coverage"]
    output_df.loc[:, "originating_lab"] = filtered_df["originating_lab"]
    output_df.loc[:, "originating_lab_address"] = filtered_df[
        "originating_lab_address"]
    output_df.loc[:, "sample_provider_sample_id"] = filtered_df["sample_id"]
    output_df.loc[:, "submitting_lab"] = "Andersen lab at Scripps Research"
    output_df.loc[:, "submitting_lab_address"] = \
        "10550 North Torrey Pines Road, La Jolla, CA 92037"
    output_df.loc[:, "submitting_lab_search_id"] = filtered_df["search_id"]
    output_df.loc[:, "project_authors"] = filtered_df["project_authors"]
    output_df.loc[:, "project_name"] = filtered_df["project_name"]
    output_df.loc[:, "comment_icon"] = ""
    output_df.loc[:, "released_2"] = ""
    output_df.loc[:, SEQ_POOL_COMP_ID] = filtered_df[SEQ_POOL_COMP_ID]

    return output_df


def generate_base_empress_df(filtered_df):
    # include ONLY samples that went through lineage calling
    has_pangolin_status = filtered_df["status"].notna()
    empress_df = filtered_df[has_pangolin_status]

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
                        seq_usable_for = seq_metadata_df.loc[:,
                                         USABLE_NAME].iat[0]
                        fasta_str = seq.formatted_fasta() + '\n'
                        if seq_usable_for == VARIANT_VAL:
                            loose_file.write(fasta_str)
                        elif seq_usable_for == VARIANT_AND_EP_VAL:
                            stringent_file.write(fasta_str)
                        # ignore any other cases


def merge_metadata(arg_list):
    qc_and_lineage_fp = arg_list[1]
    metadata_fp = arg_list[2]
    out_full_fp = arg_list[3]
    out_bjorn_fp = arg_list[4]
    out_empress_fp = arg_list[5]
    out_empress_var_fp = arg_list[6]
    out_empress_var_and_ep_fp = arg_list[7]
    full_fasta_fp = arg_list[8]
    out_loose_fasta_fp = arg_list[9]
    out_stringent_fasta_fp = arg_list[10]

    qc_and_lineage_df = pd.read_csv(qc_and_lineage_fp)  # , dtype=str)
    metadata_df = pd.read_csv(metadata_fp, dtype=str)

    full_df = generate_full_metadata(qc_and_lineage_df, metadata_df)
    full_df.to_csv(out_full_fp, index=False)

    filtered_df = filter_metadata_for_bjorn(full_df)
    bjorn_metadata_df = generate_bjorn_df(filtered_df)
    bjorn_metadata_df.columns = BJORN_COL_NAMES
    bjorn_metadata_df.to_csv(out_bjorn_fp, index=False)

    # NB that empress metadata files must be tsv
    metadata_w_search_id = metadata_df[SEARCH_ID_KEY].notna()
    metadata_w_search_ids_df = metadata_df[metadata_w_search_id].copy()
    raw_empress_df = qc_and_lineage_df.merge(
        metadata_w_search_ids_df, on=SEARCH_ID_KEY, how="outer")
    base_empress_df = generate_base_empress_df(raw_empress_df)
    base_empress_df.to_csv(out_empress_fp, sep='\t', index=False)

    usable_is_variant = base_empress_df[USABLE_NAME] == VARIANT_VAL
    usable_is_var_and_ep = base_empress_df[USABLE_NAME] == VARIANT_AND_EP_VAL

    # these are cumulative: the var_empress_df
    # contains all records with a "usable_for" value of "variant" OR
    # "variant_and_epidemiology"--OR the reference and historical records
    is_hist_or_ref = base_empress_df[IS_HIST_OR_REF] == True
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
    merge_metadata(argv)
