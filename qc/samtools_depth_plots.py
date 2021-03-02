# Modified from code by Niema Moshiri 2021

import os
from sys import argv, stderr
from os.path import isdir, isfile, basename

from seaborn import violinplot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# headless matplotlib
import matplotlib as mpl
mpl.use('Agg')

FOLDER_NOT_FILE = "ERROR: Argument is a folder, not a file"
FILE_NOT_FOUND = "ERROR: File not found"
NOT_PILEUP = "ERROR: File is not a samtools pileup output"
DUP_FILE = "ERROR: Duplicate file in arguments"
MULTIPLE_REFS = "ERROR: Multiple reference IDs were found " \
                "(needs to be exactly 1)"


def _validate_argvs(argv):
    input_fnames = dict()

    for curr_fp in argv[1:]:
        if not isfile(curr_fp):
            if isdir(curr_fp):
                return "%s: %s" % (FOLDER_NOT_FILE, curr_fp)
            else:
                return "%s: %s" % (FILE_NOT_FOUND, curr_fp)

        curr_fname = basename(curr_fp)
        if curr_fname in input_fnames:
            return "%s: %s" % (DUP_FILE, curr_fname)

        input_fnames[curr_fname] = True
    return ""


# TODO: Reunite this with its copy-paste ancestor in
#  sarscov2_consensus_acceptance.py
def _read_depths(depth_lines):
    ref_genome_ids = set()
    result = []
    for curr_line in depth_lines:
        curr_pieces = curr_line.split("\t")
        ref_genome_ids.add(curr_pieces[0])
        result.append(int(curr_pieces[2]))

        if len(ref_genome_ids) != 1:
            raise ValueError(f"Multiple reference genome ids found: "
                             f"{', '.join(ref_genome_ids)}")
    return result


def _load_depths(fps_list):
    depths_by_fname = {}

    for curr_fp in fps_list:
        with open(curr_fp) as curr_f:
            curr_fname = os.path.basename(curr_fp)
            curr_depths = _read_depths(curr_f.readlines())
            depths_by_fname[curr_fname] = curr_depths

    return depths_by_fname


def _save_line_plot(depths_by_fname, output_fname):
    filtered_depths_by_fname = {}
    for curr_name, curr_depths in depths_by_fname.items():
        if len(curr_depths) > 0:
            filtered_depths_by_fname[curr_name] = curr_depths

    YMAX = max([max(curr_depths) for curr_depths
                in filtered_depths_by_fname.values()])
    with PdfPages(output_fname) as pdf_pages:
        for i, curr_fname \
                in enumerate(sorted(filtered_depths_by_fname.keys())):
            curr_depths = filtered_depths_by_fname[curr_fname]
            num_positions = len(curr_depths)

            fig = plt.figure(i)
            plt.plot([1 + i for i in range(num_positions)], curr_depths)
            plt.xlim((1, num_positions + 1))
            plt.ylim((1, YMAX))
            plt.yscale('log')
            plt.xlabel('Consensus Sequence Position')
            plt.ylabel('Mapping Depth')
            plt.title(curr_fname, fontsize=8)
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            pdf_pages.savefig(fig)
            plt.close(fig)


def _make_parallel_lists(depths_by_fname):
    flattened_list_of_depths = []
    flattened_list_of_fnames_for_each_depth = []
    for curr_fname in sorted(depths_by_fname.keys()):
        curr_depths = depths_by_fname[curr_fname]
        flattened_list_of_depths.extend(curr_depths)
        flattened_list_of_fnames_for_each_depth.extend(
            [curr_fname] * len(curr_depths))

    return flattened_list_of_fnames_for_each_depth, flattened_list_of_depths


def _save_violin_plot(x, y, sample_names, output_fname):
    # create violin plot

    INCH_PER_SAMPLE = 0.25

    # TODO: Is headless used implicitly somewhere in plotting? else axe
    try:
        fig, ax = plt.subplots()
        headless = False
    except:
        # TODO: does this duplicate what is done at start of script?
        import matplotlib
        matplotlib.use('Agg')
        fig, ax = plt.subplots()
        headless = True
    fig.set_size_inches(INCH_PER_SAMPLE * len(sample_names), 4.8)
    violinplot(x=x, y=y)
    plt.ylim(ymin=1)
    ax.set_yscale('log')
    plt.title("Distribution of Mapping Depth")
    plt.xlabel("Sample")
    plt.ylabel("Mapping Depth (per site)")
    plt.xticks(rotation=90)
    fig.savefig(output_fname, format='pdf', bbox_inches='tight')


if __name__ == '__main__':
    USAGE = "USAGE: %s <file1.depth.txt> [file2.depth.txt] [file3.depth.txt]" \
            " ..." % argv[0]
    if len(argv) == 1 or argv[1].lower() in {'-h', '--help', '-help'}:
        print(USAGE, file=stderr)
        exit(1)

    # check command line for validity
    inputs_error_msg = _validate_argvs(argv)
    if inputs_error_msg != "":
        print(inputs_error_msg, file=stderr)
        exit(1)

    loaded_depths_by_fname = _load_depths(argv[1:])
    x_vals, y_vals = _make_parallel_lists(loaded_depths_by_fname)

    _save_line_plot(loaded_depths_by_fname, 'depth_lineplot.pdf')
    _save_violin_plot(x_vals, y_vals, loaded_depths_by_fname.keys(),
                      'depth_violin.pdf')
