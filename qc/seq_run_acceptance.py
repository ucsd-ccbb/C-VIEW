from pathlib import Path
from sys import argv


def make_collected_acceptances_tsv(arg_list):
    workspace_dir = arg_list[1]
    output_fp = arg_list[2]

    relevant_paths = Path(workspace_dir).rglob("*.acceptance.tsv")
    relevant_fps = [str(p) for p in relevant_paths]

    header_line = None
    output_lines = []
    for curr_fp in sorted(relevant_fps):
        with open(curr_fp) as curr_acceptance_f:
            curr_acceptance_lines = curr_acceptance_f.readlines()
            curr_header_line = curr_acceptance_lines[0]

            if header_line is None:
                header_line = curr_header_line
                output_lines.append(header_line)
            elif header_line != curr_header_line:
                raise ValueError(f"Error: Non-matching header lines found: "
                                 f"'{header_line}' and '{curr_header_line}'")

            output_lines.append(curr_acceptance_lines[1])

    if len(output_lines) == 0:
        raise ValueError("Error: No acceptance summary info found")

    with open(output_fp, "w") as output_f:
        output_f.writelines(output_lines)


if __name__ == '__main__':
    make_collected_acceptances_tsv(argv)
