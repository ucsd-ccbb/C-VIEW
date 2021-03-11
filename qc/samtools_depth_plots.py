# Modified from code by Niema Moshiri 2021

import os
from sys import argv, stderr

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


def _save_line_plot(depths_by_fname, output_fp):
    filtered_depths_by_fname = {}
    for curr_name, curr_depths in depths_by_fname.items():
        if len(curr_depths) > 0:
            filtered_depths_by_fname[curr_name] = curr_depths

    YMAX = max([max(curr_depths) for curr_depths
                in filtered_depths_by_fname.values()])
    with PdfPages(output_fp) as pdf_pages:
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


def _save_violin_plot(x, y, sample_names, output_fp):
    # create violin plot

    INCH_PER_SAMPLE = 0.25

    # TODO: Is headless used implicitly somewhere in plotting? else axe
    try:
        fig, ax = plt.subplots()
        # headless = False
    except:  # noqa E722
        # TODO: does this duplicate what is done at start of script?
        import matplotlib
        matplotlib.use('Agg')
        fig, ax = plt.subplots()
        # headless = True
    fig.set_size_inches(INCH_PER_SAMPLE * len(sample_names), 4.8)
    violinplot(x=x, y=y)
    plt.ylim(ymin=1)
    ax.set_yscale('log')
    plt.title("Distribution of Mapping Depth")
    plt.xlabel("Sample")
    plt.ylabel("Mapping Depth (per site)")
    plt.xticks(rotation=90)
    fig.savefig(output_fp, format='pdf', bbox_inches='tight')


def make_plots(args_list):
    lineplot_out_fp = args_list[1]
    violin_out_fp = args_list[2]

    loaded_depths_by_fname = _load_depths(args_list[3:])
    x_vals, y_vals = _make_parallel_lists(loaded_depths_by_fname)

    _save_line_plot(loaded_depths_by_fname, lineplot_out_fp)
    _save_violin_plot(x_vals, y_vals, loaded_depths_by_fname.keys(),
                      violin_out_fp)


if __name__ == '__main__':
    make_plots(argv)
