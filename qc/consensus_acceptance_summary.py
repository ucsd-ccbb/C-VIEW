import os
from sys import argv, stderr


ACCEPTANCE_SUFFIX = ".acceptance.tsv"
header_line = None
output_lines = []
output_fname = f"summary{ACCEPTANCE_SUFFIX}"
USAGE = "USAGE: %s <workspace directory>" % argv[0]

# check command line for validity
if len(argv) != 2:
    print(USAGE, file=stderr)
    exit(1)

workspace_dir = argv[1]
curr_dir_contents = os.listdir(workspace_dir)
for curr_item_name in sorted(curr_dir_contents):
    # ignore any previous output of this script and ignore any file that
    # does not look like a per-sample acceptance file
    if curr_item_name == output_fname or \
            not curr_item_name.endswith(".acceptance.tsv"):
        continue

    curr_item_fp = os.path.join(workspace_dir, curr_item_name)
    with open(curr_item_fp) as curr_acceptance_f:
        curr_acceptance_lines = curr_acceptance_f.readlines()
        curr_header_line = curr_acceptance_lines[0]

        if header_line is None:
            header_line = curr_header_line
            output_lines.append(header_line)
        elif header_line != curr_header_line:
            print(f"Error: Non-matching header lines found: "
                  f"'{header_line}' and '{curr_header_line}'")
            exit(1)

        output_lines.append(curr_acceptance_lines[1])

if len(output_lines) == 0:
    print("Error: No acceptance summary info found")
    exit(1)

output_fp = os.path.join(workspace_dir, output_fname)
with open(output_fp, "w") as output_f:
    output_f.writelines(output_lines)
