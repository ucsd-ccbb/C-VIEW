import pandas as pd
from sys import argv


SEARCH_ID_KEY = "search_id"
CONS_NAME = "consensus_seq_name"

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
                   "Comment", "Comment Icon", "Released"]


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
    filtered_df = full_df[~is_excluded]
    return filtered_df


def generate_bjorn_df(filtered_df):
    filtered_df[["sample_collection_date", "sample_collection_time"]] = \
        filtered_df["sample_collection_datetime"].str.split(expand=True)
    filtered_df[["sample_collection_year",
                 "sample_collection_month",
                 "sample_collection_day"]] = filtered_df[
        "sample_collection_date"].str.split("-", expand=True)

    output_df = pd.DataFrame()
    output_df['sample_id'] = filtered_df['sample_id']
    output_df['search_id'] = filtered_df['search_id']
    output_df["ready_for_release"] = "Yes"
    output_df["new_seqs_ready_for_release"] = "Yes"
    output_df["released"] = ""
    output_df["submitter"] = ""
    output_df["fasta_fname"] = ""
    output_df["virus_name"] = "hCoV-19/USA/" + output_df["search_id"] + "/" + \
                              filtered_df["sample_collection_year"]
    output_df["type"] = "betacoronavirus"
    output_df["passage_details"] = filtered_df["passage_details"]
    output_df["sample_collection_date"] = filtered_df["sample_collection_date"]
    output_df["sample_collection_location"] = filtered_df[
        "sample_collection_location"]
    output_df["zip"] = filtered_df["zip"]
    output_df["subject_species"] = filtered_df["subject_species"]
    output_df["addl_host_info"] = "N/A"
    output_df["subject_gender"] = filtered_df["subject_gender"]
    output_df["subject_age"] = filtered_df["subject_age"]
    output_df["patient_status"] = "N/A"
    output_df["specimen_type"] = filtered_df["specimen_type"]
    output_df["outbreak"] = "N/A"
    output_df["last_vaccinated"] = "N/A"
    output_df["treatment"] = "N/A"
    output_df["sequencing_tech"] = "Illumina"
    output_df["assembly_method"] = filtered_df["assembly_method"].str.replace(
        " version", "")
    output_df["bjorn_coverage"] = filtered_df["bjorn_coverage"]
    output_df["originating_lab"] = filtered_df["originating_lab"]
    output_df["originating_lab_address"] = filtered_df[
        "originating_lab_address"]
    output_df["sample_provider_sample_id"] = filtered_df["sample_id"]
    output_df["submitting_lab"] = "Andersen lab at Scripps Research"
    output_df["submitting_lab_address"] = \
        "10550 North Torrey Pines Road, La Jolla, CA 92037"
    output_df["submitting_lab_search_id"] = filtered_df["search_id"]
    output_df["project_authors"] = filtered_df["project_authors"]
    output_df["project_name"] = filtered_df["project_name"]
    output_df["comment_icon"] = ""
    output_df["released_2"] = ""

    return output_df


def generate_empress_df(filtered_df):
    # include ONLY samples that went through lineage calling)
    has_pangolin_status = filtered_df["status"].notna()
    empress_df = filtered_df[has_pangolin_status]

    # rearrange columns--want CONS_NAME as first column to match up
    # with the fas record names, which are used as the tree node names
    # in the tree file
    # shift column 'consensus_seq_name' to first position
    first_column = empress_df.pop(CONS_NAME)
    empress_df.insert(0, CONS_NAME, first_column)
    return empress_df


def merge_metadata(arg_list):
    qc_and_lineage_fp = arg_list[1]
    metadata_fp = arg_list[2]
    out_full_fp = arg_list[3]
    out_bjorn_fp = arg_list[4]
    out_empress_fp = arg_list[5]

    qc_and_lineage_df = pd.read_csv(qc_and_lineage_fp, dtype=str)
    records_w_search_ids = qc_and_lineage_df[SEARCH_ID_KEY].notna()
    qc_and_lineage_w_search_ids_df = qc_and_lineage_df[records_w_search_ids]
    metadata_df = pd.read_csv(metadata_fp, dtype=str)

    full_df = qc_and_lineage_w_search_ids_df.merge(
        metadata_df, on=SEARCH_ID_KEY, how="outer")
    full_df.to_csv(out_full_fp, index=False)

    filtered_df = filter_metadata_for_bjorn(full_df)
    bjorn_metadata_df = generate_bjorn_df(filtered_df)
    bjorn_metadata_df.columns = BJORN_COL_NAMES
    bjorn_metadata_df.to_csv(out_bjorn_fp, index=False)

    # NB that empress metadata files must be tsv
    metadata_w_search_id = metadata_df[SEARCH_ID_KEY].notna()
    metadata_w_search_ids_df = metadata_df[metadata_w_search_id]
    raw_empress_df = qc_and_lineage_df.merge(
        metadata_w_search_ids_df, on=SEARCH_ID_KEY, how="outer")
    empress_df = generate_empress_df(raw_empress_df)
    empress_df.to_csv(out_empress_fp, sep='\t', index=False)


if __name__ == '__main__':
    merge_metadata(argv)
