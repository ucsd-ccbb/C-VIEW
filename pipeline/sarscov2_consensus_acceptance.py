from sys import argv, stderr
import os
import nwalign3 as nw

SEARCH_ID_PREFIX = "SEARCH"
REF_FIRST_ORF_START_1BASED = 266
REF_LAST_ORF_END_1BASED = 29674
DEPTH_THRESH = 10
FRACTION_THRESH = 0.95

COL_NA = "NA"
SAMP_NAME = "sample_name"
IS_ACCEPTED = "is_accepted"
COVERAGE = "coverage_gte_10_reads"
NUM_INSERTS = "num_inserts_in_consensus"
NUM_DELS = "num_deletions_in_consensus"
PUTATIVE_SAMP_ID = "putative_sample_id"
SEQ_RUN = "seq_run"
TIMESTAMP = "timestamp"
SE_OR_PE = "se_or_pe"
IVAR_VER = "assembly_method"
CONS_SEQ_NAME = "consensus_seq_name"
FIELD_ORDER = [SAMP_NAME, IS_ACCEPTED, COVERAGE, NUM_INSERTS, NUM_DELS,
               PUTATIVE_SAMP_ID, CONS_SEQ_NAME, IVAR_VER, TIMESTAMP,
               SE_OR_PE, SEQ_RUN]


def check_consensus_acceptance(consensus_seq, consensus_depths, ref_genome_seq,
                               ref_first_orf_start_1based,
                               ref_last_orf_end_1based,
                               depth_threshold, fraction_threshold):

    ref_first_orf_start_0based = ref_first_orf_start_1based - 1
    ref_last_orf_end_0based = ref_last_orf_end_1based - 1

    cons_first_orf_start_0based, cons_last_orf_end_0based, \
    num_ref_gaps_in_orf_region, num_cons_gaps_in_orf_region = \
        _get_consensus_orfs_start_end(consensus_seq, ref_genome_seq,
                                      ref_first_orf_start_0based,
                                      ref_last_orf_end_0based)

    is_accepted, depth_pass_fraction = _verify_fraction_acceptable_bases(
        consensus_seq, consensus_depths,
        cons_first_orf_start_0based, cons_last_orf_end_0based,
        depth_threshold, fraction_threshold)

    result = {}
    result[NUM_INSERTS] = num_ref_gaps_in_orf_region
    result[NUM_DELS] = num_cons_gaps_in_orf_region
    result[IS_ACCEPTED] = is_accepted
    result[COVERAGE] = depth_pass_fraction

    return result


# this sucker is a heuristic ... don't expect too much
def _extract_putative_sample_id(putative_sample_name):
    result = putative_sample_name

    lane_split = putative_sample_name.split("_L00")
    if len(lane_split) < 2:
        return result

    result = _attempt_extract_search_id(lane_split[0])
    return result


def _attempt_extract_search_id(putative_sample_id):
    # expect sample_id to look something like
    # 002idSEARCH-5329-SAN

    result = putative_sample_id
    search_split = putative_sample_id.split(SEARCH_ID_PREFIX)
    if len(search_split) != 2:
        return result

    result = SEARCH_ID_PREFIX + search_split[1]
    return result


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
    ref_genome_seq = _read_ref_genome_fas(ref_genome_file)

    if consensus_lines is not None and \
            len(consensus_depths) > 0 and ref_genome_seq != "":
        result = (consensus_lines, consensus_depths,
                  ref_genome_seq)

    return result


def _read_consensus_file(consensus_file):
    result = None
    consensus_lines = consensus_file.readlines()

    # consensus.fa must have at least 2 lines in order to contain
    # a real fasta record
    if len(consensus_lines) == 2:
        consensus_seq_name = consensus_lines[0].strip().replace(">", "")
        consensus_seq = consensus_lines[1].strip()

        if consensus_seq != "":
            result = (consensus_seq_name, consensus_seq)

    return result


def _read_depths(depth_file):
    result = []
    for curr_line in depth_file.readlines():
        curr_pieces = curr_line.split("\t")
        result.append(int(curr_pieces[2]))
    return result


def _read_ref_genome_fas(ref_genome_fas_file):
    pieces = []
    for curr_line in ref_genome_fas_file.readlines():
        if not curr_line.startswith(">"):
            pieces.append(curr_line.strip())
    result = "".join(pieces)
    return result


def _get_consensus_orfs_start_end(consensus_seq, ref_genome_seq,
                                  ref_first_orf_start_0based,
                                  ref_last_orf_end_0based):
    """Get first orf start and last orf end indexes in consensus sequence

    Note that since we start with the start position (inclusive) of the first
    orf and the end position (inclusive) of last orf *on the reference*, and
    since the consensus may have insertions or deletions relative
    to the reference, this requires performing a *global* alignment of the
    consensus against the reference and using the results to determine what
    positions on the consensus correspond to the start of the first orf and
    the end of the last orf on the reference.
    """

    ref_gapped_seq, consensus_gapped_seq = nw.global_align(
        ref_genome_seq, consensus_seq)

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

    return cons_first_orf_start_0based, cons_last_orf_end_0based, \
           num_ref_gaps_in_orf_region, num_cons_gaps_in_orf_region


def _verify_fraction_acceptable_bases(consensus_seq, consensus_depths,
                                      cons_first_orf_start_0based,
                                      cons_last_orf_end_0based,
                                      depth_threshold, fraction_threshold):
    # range[,)
    depth_pass_list = []
    pass_list = []
    for i in range(cons_last_orf_end_0based + 1):
        if i >= cons_first_orf_start_0based:
            curr_base = consensus_seq[i]
            curr_depth = consensus_depths[i]

            base_pass = curr_base in ['A', 'C', 'G', 'T']
            depth_pass = curr_depth >= depth_threshold
            depth_pass_list.append(int(depth_pass))
            pass_list.append(int(base_pass and depth_pass))

    depth_pass_fraction = sum(depth_pass_list)/len(depth_pass_list)
    pass_fraction = sum(pass_list)/len(pass_list)
    result = pass_fraction >= fraction_threshold
    return result, depth_pass_fraction


def _acceptance_check_inputs(file_inputs_tuple):
    result = None
    if file_inputs_tuple is not None:
        i_consensus_seq = file_inputs_tuple[0][1]
        i_depths = file_inputs_tuple[1]
        i_ref_seq = file_inputs_tuple[2]

        result = check_consensus_acceptance(
            i_consensus_seq, i_depths, i_ref_seq,
            REF_FIRST_ORF_START_1BASED, REF_LAST_ORF_END_1BASED,
            DEPTH_THRESH, FRACTION_THRESH)

    return result


def _generate_header_and_data_lines(output_dict):
    data_strs = []
    for curr_field_name in FIELD_ORDER:
        data_strs.append(str(output_dict[curr_field_name]))

    header_line = "\t".join(FIELD_ORDER) + "\n"
    data_line = "\t".join(data_strs) + "\n"
    return [header_line, data_line]


if __name__ == '__main__':
    USAGE = "USAGE: %s <sequencing run name> <timestamp> <se or pe> " \
            "<ivar version> <sample name> <consensus.fa file path> " \
            "<depth.txt file path> <reference genome.fas file path> " \
            % argv[0]

    argv = ["python sarscov2_consensus_acceptance.py",
            "2021-02-08-ARTIC",
            "2021-02-26_19-40-24",
            "pe",
            "iVar version 1.3.1\nPlease raise issues and bug reports at "
            "https://github.com/andersen-lab/ivar/",
            "039idSEARCH-5366-SAN_L001_L002_L003_L004",
            #"/Users/amandabirmingham/Downloads/SU002_S13_L001.trimmed."
            #"sorted.pileup.consensus.fa",
            "/Users/amandabirmingham/Downloads/039idSEARCH-5366-SAN_L001_L002"
            "_L003_L004.trimmed.sorted.pileup.consensus.fa",
            "/Users/amandabirmingham/Downloads/039idSEARCH-5366-SAN_L001_L002"
            "_L003_L004.trimmed.sorted.depth.txt",
            "/Users/amandabirmingham/Work/Repositories/covid_sequencing_"
            "analysis_pipeline/reference_files/NC_045512.2.fas"]

    # check command line for validity
    if len(argv) != 9:
        print(USAGE, file=stderr)
        exit(1)

    seq_run = timestamp = se_or_pe = ivar_ver_string = sample_name = None
    input_consensus_fa_fp = input_depth_txt_fp = input_ref_genome_fas_fp = None

    try:
        seq_run = argv[1]
        timestamp = argv[2]
        se_or_pe = argv[3]
        ivar_ver_string = argv[4]
        sample_name = argv[5]
        input_consensus_fa_fp = argv[6]
        input_depth_txt_fp = argv[7]
        input_ref_genome_fas_fp = argv[8]
    except Exception as e:
        print(f"Error parsing arguments: {e}")
        exit(1)

    ivar_version = ivar_ver_string.splitlines()[0].strip()
    putative_sample_id = _extract_putative_sample_id(sample_name)
    output_fields = {SAMP_NAME: sample_name,
                     PUTATIVE_SAMP_ID: putative_sample_id,
                     SEQ_RUN: seq_run,
                     TIMESTAMP: timestamp,
                     SE_OR_PE: se_or_pe,
                     IVAR_VER: ivar_version,
                     IS_ACCEPTED: False,
                     CONS_SEQ_NAME: COL_NA,
                     COVERAGE: COL_NA,
                     NUM_INSERTS: COL_NA,
                     NUM_DELS: COL_NA}

    contents_tuple = _read_input_fps(
        input_consensus_fa_fp, input_depth_txt_fp, input_ref_genome_fas_fp)

    conclusions_dict = _acceptance_check_inputs(contents_tuple)
    if conclusions_dict is not None:
        # overwrite default (dummy) consensus sequence name
        output_fields[CONS_SEQ_NAME] = contents_tuple[0][0]
        # The accepted, coverage, num inserts, and num dels values in
        # the conclusions_dict will overwrite the default ones
        output_fields.update(conclusions_dict)

    output_lines = _generate_header_and_data_lines(output_fields)

    dir_fp = os.path.dirname(input_consensus_fa_fp)
    output_fp = os.path.join(dir_fp, f"{sample_name}.acceptance.tsv")
    with open(output_fp, 'w') as output_f:
        output_f.writelines(output_lines)
