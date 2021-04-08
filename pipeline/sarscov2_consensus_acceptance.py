from sys import argv
import os
import json
import nwalign3 as nw

ID_DELIMITER = "__"
REF_FIRST_ORF_START_1BASED = 266
REF_LAST_ORF_END_1BASED = 29674
DEPTH_THRESH = 10
FRACTION_THRESH = 0.95

COL_NA = "NA"
SEQ_POOL_COMP_ID = "sequenced_pool_component_id"
IS_ACCEPTED = "is_accepted"
INDELS_FLAGGED = "indels_flagged"
PASS_DEPTH_FRACTION = "coverage_gte_10_reads"
PASS_BASE_IDENTITY_FRACTION = "fraction_acgt_bases"
NUM_INSERTS = "num_inserts_in_consensus"
NUM_DELS = "num_deletions_in_consensus"
SEARCH_ID = "search_id"
SEQ_RUN = "seq_run"
TIMESTAMP = "timestamp"
SE_OR_PE = "se_or_pe"
IVAR_VER = "assembly_method"
CONS_SEQ_NAME = "consensus_seq_name"
REF_SEQ_NAME = "reference_seq_name"
CONSENSUS_S3 = "consensus_s3"
TRIMMED_BAM_S3 = "trimmed_bam_s3"
VARIANT_S3 = "variants_s3"

FIELD_ORDER = [SEQ_POOL_COMP_ID, SEARCH_ID,
               PASS_DEPTH_FRACTION, PASS_BASE_IDENTITY_FRACTION,
               NUM_INSERTS, NUM_DELS,
               CONS_SEQ_NAME, IVAR_VER,
               TIMESTAMP, SE_OR_PE, SEQ_RUN,
               CONSENSUS_S3, TRIMMED_BAM_S3, VARIANT_S3]
REF_ALIGNMENT = "ref_alignment"
CONS_ALIGNMENT = "cons_alignment"
CONS_FIRST_ORF_START_0B = "cons_first_orf_start_0based"
CONS_LAST_ORF_END_0B = "cons_last_orf_end_0based"


def check_acceptance(consensus_seq, consensus_depths,
                     align_results,
                     depth_threshold, fraction_threshold):
    # range[,)
    result = {}
    depth_pass_list = []
    identity_pass_list = []
    pass_list = []
    for i in range(align_results[CONS_LAST_ORF_END_0B] + 1):
        if i >= align_results[CONS_FIRST_ORF_START_0B]:
            curr_base = consensus_seq[i]
            curr_depth = consensus_depths[i]

            base_pass = curr_base in ['A', 'C', 'G', 'T']
            identity_pass_list.append(base_pass)

            depth_pass = curr_depth >= depth_threshold
            depth_pass_list.append(int(depth_pass))

            pass_list.append(int(base_pass and depth_pass))

    depth_pass_fraction = sum(depth_pass_list)/len(depth_pass_list)
    identity_pass_fraction = sum(identity_pass_list)/len(identity_pass_list)
    pass_fraction = sum(pass_list)/len(pass_list)
    fraction_passes = pass_fraction >= fraction_threshold

    # if the number of ref gaps (consensus insertions) and/or
    # the number of consensus gaps (consensus deletions) is NOT
    # a perfect multiple of 3, flag this consensus sequence as
    # possibly suspicious
    indels_flagged = \
        align_results[NUM_INSERTS] % 3 != 0 or \
        align_results[NUM_DELS] % 3 != 0

    is_accepted = fraction_passes  # and (not indels_flagged)

    result[PASS_DEPTH_FRACTION] = depth_pass_fraction
    result[PASS_BASE_IDENTITY_FRACTION] = identity_pass_fraction
    result[INDELS_FLAGGED] = indels_flagged
    result[IS_ACCEPTED] = is_accepted
    return result


def get_search_id(pipeline_sample_name):
    # split on double underscores; the first entry is the search id
    name_pieces = pipeline_sample_name.split(ID_DELIMITER)
    return name_pieces[0]


def _read_input_fps(consensus_fa_fp, depth_txt_fp, ref_genome_fas_fp):
    result = None

    if os.path.isfile(consensus_fa_fp) and os.path.isfile(depth_txt_fp):
        with open(consensus_fa_fp) as consensus_file:
            with open(depth_txt_fp) as depth_file:
                with open(ref_genome_fas_fp) as ref_genome_file:
                    result = _read_input_files(
                        consensus_file, depth_file, ref_genome_file)

    return result


def _read_input_files(consensus_file, depth_file, ref_genome_file):
    result = None
    consensus_lines = _read_consensus_file(consensus_file)
    consensus_depths = _read_depths(depth_file)
    ref_genome_lines = _read_ref_genome_file(ref_genome_file)

    if consensus_lines is not None and len(consensus_depths) > 0 \
            and ref_genome_lines is not None:
        result = (consensus_lines, consensus_depths,
                  ref_genome_lines)

    return result


def _read_consensus_file(consensus_file):
    consensus_lines = consensus_file.readlines()
    return _read_fa_lines(consensus_lines)


def _read_ref_genome_file(ref_genome_fas_file):
    ref_lines = _read_ref_genome_lines(ref_genome_fas_file.readlines())
    return _read_fa_lines(ref_lines)


def _read_fa_lines(lines):
    result = None

    # a .fa must have at 2 lines in order to contain
    # a real fasta record
    if len(lines) == 2:
        seq_name = lines[0].strip().replace(">", "")
        seq = lines[1].strip()

        if seq != "":
            result = (seq_name, seq)

    return result


def _read_ref_genome_lines(lines):
    revised_lines = []
    pieces = []
    for curr_line in lines:
        stripped_line = curr_line.strip()
        if stripped_line.startswith(">"):
            revised_lines.append(stripped_line)
        else:
            pieces.append(stripped_line)

    revised_lines.append("".join(pieces))
    return revised_lines


def _read_depths(depth_file):
    result = []
    for curr_line in depth_file.readlines():
        curr_pieces = curr_line.split("\t")
        result.append(int(curr_pieces[2]))
    return result


def _pairwise_align(consensus_info, ref_genome_info,
                    ref_first_orf_start_1based, ref_last_orf_end_1based):
    ref_first_orf_start_0based = ref_first_orf_start_1based - 1
    ref_last_orf_end_0based = ref_last_orf_end_1based - 1

    ref_gapped_seq, consensus_gapped_seq = nw.global_align(
        ref_genome_info[1], consensus_info[1])

    num_ref_gaps_in_orf_region = num_cons_gaps_in_orf_region = 0
    curr_ref_index = curr_cons_index = -1
    cons_first_orf_start_0based = cons_last_orf_end_0based = None
    for gapped_index in range(len(ref_gapped_seq)):
        curr_cons_base = consensus_gapped_seq[gapped_index]

        if curr_cons_base != "-":
            curr_cons_index += 1
        else:
            # if we are within the orf-containing region, keep
            # track of how many gap bases we see
            if cons_first_orf_start_0based is not None and \
                    cons_last_orf_end_0based is None:
                num_cons_gaps_in_orf_region += 1

        curr_ref_base = ref_gapped_seq[gapped_index]
        if curr_ref_base != "-":
            curr_ref_index += 1
            if curr_ref_index == ref_first_orf_start_0based:
                cons_first_orf_start_0based = curr_cons_index
            elif curr_ref_index == ref_last_orf_end_0based:
                cons_last_orf_end_0based = curr_cons_index
                break
        else:
            if cons_first_orf_start_0based is not None and \
                    cons_last_orf_end_0based is None:
                num_ref_gaps_in_orf_region += 1

    result = {}
    result[CONS_SEQ_NAME] = consensus_info[0]
    result[REF_SEQ_NAME] = ref_genome_info[0]
    result[REF_ALIGNMENT] = ref_gapped_seq
    result[CONS_ALIGNMENT] = consensus_gapped_seq
    result[CONS_FIRST_ORF_START_0B] = cons_first_orf_start_0based
    result[CONS_LAST_ORF_END_0B] = cons_last_orf_end_0based
    result[NUM_INSERTS] = num_ref_gaps_in_orf_region
    result[NUM_DELS] = num_cons_gaps_in_orf_region

    return result


def _acceptance_check_inputs(file_inputs_tuple):
    align_result = accept_result = None
    if file_inputs_tuple is not None:
        i_consensus_info = file_inputs_tuple[0]
        i_depths = file_inputs_tuple[1]
        i_ref_info = file_inputs_tuple[2]

        align_result = _pairwise_align(
            i_consensus_info, i_ref_info,
            REF_FIRST_ORF_START_1BASED, REF_LAST_ORF_END_1BASED)

        accept_result = check_acceptance(
            i_consensus_info[1], i_depths, align_result,
            DEPTH_THRESH, FRACTION_THRESH)

    return align_result, accept_result


def _generate_header_and_data_lines(output_dict):
    data_strs = []
    for curr_field_name in FIELD_ORDER:
        data_strs.append(str(output_dict.get(curr_field_name, COL_NA)))

    header_line = "\t".join(FIELD_ORDER) + "\n"
    data_line = "\t".join(data_strs) + "\n"
    return [header_line, data_line]


def generate_acceptance_tsv(arg_list):
    seq_run = arg_list[1]
    timestamp = arg_list[2]
    se_or_pe = arg_list[3]
    ivar_ver_string = arg_list[4]
    seq_pool_comp_id = arg_list[5]
    input_consensus_fa_fp = arg_list[6]
    input_depth_txt_fp = arg_list[7]
    input_ref_genome_fas_fp = arg_list[8]
    trimmed_bam_fname = arg_list[9]
    variants_tsv_fname = arg_list[10]
    s3_dir = arg_list[11]
    output_table_fp = arg_list[12]
    output_align_fp = arg_list[13]

    ivar_version = ivar_ver_string.splitlines()[0].strip()
    search_id = get_search_id(seq_pool_comp_id)
    output_fields = {SEQ_POOL_COMP_ID: seq_pool_comp_id,
                     SEARCH_ID: search_id,
                     SEQ_RUN: seq_run,
                     TIMESTAMP: timestamp,
                     SE_OR_PE: se_or_pe,
                     IVAR_VER: ivar_version,
                     IS_ACCEPTED: False}

    # add the S3 URLs
    consensus_fname = os.path.basename(input_consensus_fa_fp)
    output_fields[CONSENSUS_S3] = os.path.join(s3_dir, consensus_fname)
    output_fields[TRIMMED_BAM_S3] = os.path.join(s3_dir, trimmed_bam_fname)
    output_fields[VARIANT_S3] = os.path.join(s3_dir, variants_tsv_fname)

    contents_tuple = _read_input_fps(
        input_consensus_fa_fp, input_depth_txt_fp, input_ref_genome_fas_fp)

    align_result, accept_result = _acceptance_check_inputs(contents_tuple)
    align_result = {} if align_result is None else align_result
    with open(output_align_fp, 'w') as fp:
        json.dump(align_result, fp)

    if accept_result is not None:
        # The relevant values in the align and accept dicts
        # will overwrite the default ones
        output_fields.update(align_result)
        output_fields.update(accept_result)

    output_lines = _generate_header_and_data_lines(output_fields)
    with open(output_table_fp, 'w') as output_f:
        output_f.writelines(output_lines)


if __name__ == '__main__':
    generate_acceptance_tsv(argv)
