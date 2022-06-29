import pandas as pd
from sys import argv

SEARCH_ID_KEY = "search_id"
SEQ_POOL_COMP_ID = "sequenced_pool_component_id"
VARIANT_S3 = "variants_s3"
CONSENSUS_S3 = "consensus_s3"
BAM_S3 = "trimmed_bam_s3"
SOURCE = "source"
SEQ_RUN = "seq_run"
OVERALL_FAIL = "overall_fail"
SUBMIT_TO_GISAID = "submit_to_gisaid"
INDEX_COL_NAME = "Unnamed: 0"
SUBJ_AGE = "subject_age"
METADATA_CLEARED = "metadata_cleared"
SEQUENCING_TECH_KEY = "sequencing_tech"

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
                   "Sequenced Pool Component Id", "Variant File S3 URL",
                   "Consensus File S3 URL", "BAM File S3 URL",
                   "Source", "Sequencing Run", "Overall Fail",
                   "Inspect Submit-to-GISAID", METADATA_CLEARED]


def filter_metadata_for_bjorn(full_df):
    has_no_search_id = full_df[SEARCH_ID_KEY].isna()
    has_no_specimen_type = full_df["specimen_type"].isna()
    has_no_consensus_seq = full_df["consensus_seq_name"].isna()
    is_control = full_df["specimen_type"] == "control"

    is_human = full_df["subject_species"] == "Human"
    has_gender = full_df["subject_gender"].isin(["Male", "Female", "Unknown"])
    # NB: the age column is a mix of strings and NaNs; we have to convert it to
    # a number type (temporarily) to do the ">0" comparison, and although all
    # the real values are integers, NaN won't cast to an integer--only a float
    has_positive_age = (~full_df[SUBJ_AGE].isna()) & \
                       (full_df[SUBJ_AGE].astype(float) > 0)
    has_demographics = has_gender & has_positive_age
    is_scripps_sequencing = full_df[SOURCE] == "Scripps-Sequencing"
    is_human_wo_demographics_not_scripps = \
        is_human & ~has_demographics & ~is_scripps_sequencing

    is_excluded = (has_no_search_id | has_no_specimen_type |
                   has_no_consensus_seq | is_control |
                   is_human_wo_demographics_not_scripps)
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

    # The metadata reads in as strings (hence 'True') while the src and lineage
    # data reads in parsed (hence False, not 'False')
    inspect_approval = filtered_df[SUBMIT_TO_GISAID] == 'True'  # noqa 712
    # NB: if sample's metadata_cleared is NA, sample IS allowed for release.
    # This is because we have no knowledge about whether that metadata is bad.
    metadata_bad = filtered_df[METADATA_CLEARED] == 'False'  # noqa 712
    # NB: this check is false for overall_fail values of None, as it should be:
    # don't release items with overall_fail == None since they are unknown
    no_overall_fail = filtered_df[OVERALL_FAIL] == False  # noqa 712

    release_mask = inspect_approval & no_overall_fail & ~metadata_bad
    output_df.loc[:, "ready_for_release"] = "No"
    output_df.loc[release_mask, "ready_for_release"] = "Yes"

    output_df.loc[:, "new_seqs_ready_for_release"] = "Yes"
    output_df.loc[:, "released"] = ""
    output_df.loc[:, "submitter"] = ""
    output_df.loc[:, "fasta_fname"] = ""

    # Note: if any of the virus name input fields is non-existent, then the
    # resulting virus name is non-existent, which seems reasonable.
    # I believe this justifies hardcoding the "USA": this is appropriate iff
    # the sample has a state code, and the string will not be generated if the
    # sample does NOT have a state code.
    output_df.loc[:, "virus_name"] = "hCoV-19/USA/" \
                                     + filtered_df["state_code"] \
                                     + "-" + filtered_df["search_id"] \
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
    output_df.loc[:, SEQUENCING_TECH_KEY] = filtered_df[SEQUENCING_TECH_KEY]
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
    output_df.loc[:, VARIANT_S3] = filtered_df[VARIANT_S3]
    output_df.loc[:, CONSENSUS_S3] = filtered_df[CONSENSUS_S3]
    output_df.loc[:, BAM_S3] = filtered_df[BAM_S3]
    output_df.loc[:, SOURCE] = filtered_df[SOURCE]
    output_df.loc[:, SEQ_RUN] = filtered_df[SEQ_RUN]
    output_df.loc[:, OVERALL_FAIL] = filtered_df[OVERALL_FAIL]
    output_df.loc[:, SUBMIT_TO_GISAID] = filtered_df[SUBMIT_TO_GISAID]
    output_df.loc[:, METADATA_CLEARED] = filtered_df[METADATA_CLEARED]

    return output_df


def merge_metadata(arg_list):
    qc_and_lineage_fp = arg_list[1]
    metadata_fp = arg_list[2]
    out_full_fp = arg_list[3]
    out_bjorn_fp = arg_list[4]

    qc_and_lineage_df = pd.read_csv(qc_and_lineage_fp)  # , dtype=str)
    metadata_df = pd.read_csv(metadata_fp, dtype=str)
    # if the input metadata file has an unnamed first (index) column, drop it
    if metadata_df.columns[0] == INDEX_COL_NAME:
        metadata_df.pop(INDEX_COL_NAME)

    records_w_search_ids = qc_and_lineage_df[SEARCH_ID_KEY].notna()
    qc_and_lineage_w_search_ids_df = qc_and_lineage_df[
        records_w_search_ids].copy()
    full_df = qc_and_lineage_w_search_ids_df.merge(
        metadata_df, on=SEARCH_ID_KEY, how="outer")
    full_df.to_csv(out_full_fp, index=False)

    filtered_df = filter_metadata_for_bjorn(full_df)
    if len(filtered_df) > 0:
        bjorn_metadata_df = generate_bjorn_df(filtered_df)
        bjorn_metadata_df.columns = BJORN_COL_NAMES
        bjorn_metadata_df.to_csv(out_bjorn_fp, index=False)
    else:
        # create a bjorn file that has only a header line
        with open(out_bjorn_fp, 'w') as f:
            f.write(",".join(BJORN_COL_NAMES))


if __name__ == '__main__':
    merge_metadata(argv)
