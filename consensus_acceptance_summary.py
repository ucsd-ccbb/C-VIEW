import os
from sys import exit


ACCEPTANCE_SUFFIX = ".acceptance.tsv"
header_line = None
output_lines = []
curr_dir_contents = os.listdir(".")
for curr_item_name in sorted(curr_dir_contents):
    if ".acceptance.tsv" not in curr_item_name:
        continue

    with open(curr_item_name) as curr_acceptance_f:
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

with open(f"summary{ACCEPTANCE_SUFFIX}", "w") as output_f:
    output_f.writelines(output_lines)
