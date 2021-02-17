from sys import argv, stderr
import os
import errno
import nwalign3 as nw

SEARCH_ID_PREFIX = "SEARCH"
REF_FIRST_ORF_START_1BASED = 266
REF_LAST_ORF_END_1BASED = 29674
DEPTH_THRESH = 10
FRACTION_THRESH = 0.95
USAGE = "USAGE: %s <consensus.fa file path> <depth.txt file path> " \
        "<reference genome.fas file path>" % argv[0]
NO_DEPTH_PASS_INFO = "NA"


def check_consensus_acceptance_by_fps(consensus_fa_fp, depth_txt_fp,
                                      ref_genome_fas_fp):
    if not ((os.path.isfile(consensus_fa_fp)
             and os.path.isfile(depth_txt_fp))):
        return False, NO_DEPTH_PASS_INFO

    with open(consensus_fa_fp) as consensus_file:
        with open(depth_txt_fp) as depth_file:
            with open(ref_genome_fas_fp) as ref_genome_file:
                result, depth_pass_fraction = check_consensus_acceptance(
                    consensus_file, depth_file, ref_genome_file,
                    REF_FIRST_ORF_START_1BASED, REF_LAST_ORF_END_1BASED,
                    DEPTH_THRESH, FRACTION_THRESH)

    return result, depth_pass_fraction


def check_consensus_acceptance(consensus_file, depth_file, ref_genome_file,
                               ref_first_orf_start_1based,
                               ref_last_orf_end_1based,
                               depth_threshold, fraction_threshold):

    consensus_lines = consensus_file.readlines()
    consensus_depths = _read_depths(depth_file.readlines())

    # consensus.fa must have at least 2 lines in order to contain
    # a real fasta record
    if len(consensus_lines) < 2 or len(consensus_depths) == 0:
        return False, NO_DEPTH_PASS_INFO

    consensus_seq = consensus_lines[1].strip()
    if consensus_seq == "":
        return False, NO_DEPTH_PASS_INFO

    ref_genome_seq = _read_ref_genome_fas(ref_genome_file)

    ref_first_orf_start_0based = ref_first_orf_start_1based - 1
    ref_last_orf_end_0based = ref_last_orf_end_1based - 1
    cons_first_orf_start_0based, cons_last_orf_end_0based = \
        _get_consensus_orfs_start_end(consensus_seq, ref_genome_seq,
                                      ref_first_orf_start_0based,
                                      ref_last_orf_end_0based)

    result, depth_pass_fraction = _verify_fraction_acceptable_bases(
        consensus_seq, consensus_depths,
        cons_first_orf_start_0based, cons_last_orf_end_0based,
        depth_threshold, fraction_threshold)

    return result, depth_pass_fraction


def _read_depths(depth_lines):
    result = []
    for curr_line in depth_lines:
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

    curr_ref_index = curr_cons_index = -1
    cons_first_orf_start_0based = cons_last_orf_end_0based = None
    for gapped_index in range(len(ref_gapped_seq)):
        curr_cons_base = consensus_gapped_seq[gapped_index]
        curr_cons_index += 1 if curr_cons_base != "-" else 0

        curr_ref_base = ref_gapped_seq[gapped_index]
        if curr_ref_base != "-":
            curr_ref_index += 1
            if curr_ref_index == ref_first_orf_start_0based:
                cons_first_orf_start_0based = curr_cons_index
            elif curr_ref_index == ref_last_orf_end_0based:
                cons_last_orf_end_0based = curr_cons_index
                break

    return cons_first_orf_start_0based, cons_last_orf_end_0based


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


def _write_acceptance_check_to_file(consensus_fa_fp, is_accepted,
                                    depth_pass_fract):
    dir_fp = os.path.dirname(consensus_fa_fp)
    base_name = os.path.splitext(os.path.basename(consensus_fa_fp))[0]
    putative_sample_name = base_name.replace(
        ".trimmed.sorted.pileup.consensus", "")
    output_fp = os.path.join(dir_fp, f"{putative_sample_name}.acceptance.tsv")
    search_id = _extract_search_id(putative_sample_name)

    with open(output_fp, 'w') as output_f:
        header_line = "fastq_id\tis_accepted\t" \
                      "coverage_gte_10_reads\tsearch_id\n"
        data_line = f"{putative_sample_name}\t{is_accepted}" \
                    f"\t{depth_pass_fract}\t{search_id}\n"
        output_f.writelines([header_line, data_line])

    return output_fp


def _extract_search_id(sample_id):
    # expect sample_id to look something like
    # 002idSEARCH-5329-SAN_L001_L002_L003_L004

    result = ""
    search_split = sample_id.split(SEARCH_ID_PREFIX)
    if len(search_split) != 2:
        return result

    underscore_split = search_split[1].split("_")
    if len(underscore_split) < 2:
        return result

    result = SEARCH_ID_PREFIX + underscore_split[0]
    return result


if __name__ == '__main__':
    # check command line for validity
    if len(argv) != 4:
        print(USAGE, file=stderr)
        exit(1)

    input_consensus_fa_fp = input_depth_txt_fp = input_ref_genome_fas_fp = None
    try:
        input_consensus_fa_fp = argv[1]
        input_depth_txt_fp = argv[2]
        input_ref_genome_fas_fp = argv[3]
    except Exception as e:
        print(f"Error parsing arguments: {e}")
        exit(1)

    input_is_accepted, i_depth_pass_fract = check_consensus_acceptance_by_fps(
        input_consensus_fa_fp, input_depth_txt_fp, input_ref_genome_fas_fp)
    _write_acceptance_check_to_file(input_consensus_fa_fp, input_is_accepted,
                                    i_depth_pass_fract)
