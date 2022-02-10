import pandas as pd
from sys import argv

# The "doReadFasta" method in IQ-TREE 2 (see
# https://github.com/iqtree/iqtree2/blob/
# 3e6df5205fab8cbf3b8f3c01086cd7938fd26554/alignment/alignment.cpp#L2098-L2191
# ) actively *changes* the sequence labels provided in the input fasta-format
# alignment, chewing back the name at each whitespace and checking if it is
# still unique within the alignment until it either (a) gets a name that has
# some whitespace in it but is the smallest unique fraction of the original or
# gets a unique name with no whitespace in it. This affects, for example, the
# reference sars-cov2 sequence from NCBI, which has fasta label
# NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-
# Hu-1, complete genome
# and gets munged into
# NC_045512.2_Severe
#
# This issue leads to a mismatch between the (munged) ids in IQ-TREE's output
# newark tree file and the ids in the alignment file (and all our other files,
# like the empress metadata file).  This could, in theory, affect ANY fasta
# label that includes any kind of whitespace.

# One option would be to launder all the sequence ids into something without
# whitespace/etc before we put them through our pipeline, but (a) that makes us
# guilty of the very same sin of irrevocably *changing* ids, and (b) this
# approach could break if IQ-TREE changes its precise munging approach.  Though
# still unpleasant, it seems safer and more flexible to read the IQ-TREE output
# to find out what it changed, via notations like
# NOTE: Change sequence name 'NC_045512.2 Severe acute respiratory syndrome cor
# onavirus 2 isolate Wuhan-Hu-1, complete genome' -> NC_045512.2_Severe
# This should catch ANY affected label, and should be robust to IQ-TREE
# changing its criteria and/or munging rules.

CHANGE_NAME_STR = "NOTE: Change sequence name "
NAMES_DIVIDER = "->"
CONS_SEQ_NAME_COL = "consensus_seq_name"


def _find_name_changes(tree_log_fp):
    changed_names = []
    with open(tree_log_fp) as tree_log_file:
        for curr_line in tree_log_file.readlines():
            if curr_line.startswith(CHANGE_NAME_STR):
                change_pieces = curr_line.replace(
                    CHANGE_NAME_STR, "").split(NAMES_DIVIDER)
                # trim off whitespace and remove single quotes
                clean_pieces = [x.strip().replace("'", "") for x
                                in change_pieces]
                changed_names.append(clean_pieces)

    return changed_names


def _write_alignment_w_iqtree_renames(changed_names, align_orig_fp,
                                      align_out_fp):
    with open(align_orig_fp) as align_orig_file:
        with open(align_out_fp, 'w') as align_out_file:
            for curr_line in align_orig_file.readlines():
                if curr_line.startswith(">"):
                    curr_label = curr_line.replace(">", "").strip()

                    for curr_name_change in changed_names:
                        if curr_label == curr_name_change[0]:
                            curr_line = f">{curr_name_change[1]}\n"
                            break

                align_out_file.write(curr_line)


def _write_metadata_w_iqtree_renames(changed_names, metadata_orig_fp,
                                     metadata_output_fp):
    metadata_df = pd.read_csv(metadata_orig_fp, sep="\t", dtype=str)
    old_names = [x[0] for x in changed_names]
    new_names = [x[1] for x in changed_names]
    metadata_df[CONS_SEQ_NAME_COL] = metadata_df[CONS_SEQ_NAME_COL].replace(
        old_names, new_names)
    metadata_df.to_csv(metadata_output_fp, sep='\t', index=False)


def make_files_w_iqtree_renames(arg_list):
    tree_log_fp = arg_list[1]
    in_align_fp = arg_list[2]
    in_metadata_fp = arg_list[3]
    out_align_fp = arg_list[4]
    out_metadata_fp = arg_list[5]

    name_changes = _find_name_changes(tree_log_fp)
    _write_alignment_w_iqtree_renames(name_changes, in_align_fp, out_align_fp)
    _write_metadata_w_iqtree_renames(name_changes, in_metadata_fp,
                                     out_metadata_fp)


if __name__ == '__main__':
    make_files_w_iqtree_renames(argv)
