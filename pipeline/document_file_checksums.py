import os
import re
from hashlib import md5
from sys import argv


def generate_input_checksums(workspace_fp, fsuffixes_to_include=None):
    result = {}

    for base_fp, _, fnames in os.walk(workspace_fp):
        for curr_fname in fnames:
            if curr_fname.startswith("."):
                continue
            elif fsuffixes_to_include:
                matches = [x for x in fsuffixes_to_include if x in curr_fname]
                if len(matches) == 0:
                    continue

            fpath = os.path.join(base_fp, curr_fname)
            with open(fpath, "rb") as f:
                filebytes = f.read()
                filehash = md5(filebytes)
                filehash_hex = filehash.hexdigest()
                if curr_fname in result:
                    # TODO: remove test print
                    print(curr_fname)
                    raise ValueError(f"file name {curr_fname} found >1 time")
                result[curr_fname] = filehash_hex

    return result


def generate_checksums_file(args_list):
    input_dir = args_list[1]
    output_fp = args_list[2]
    fnames_to_include = None if len(args_list) < 4 else args_list[3:]

    fname_and_checksum_dict = generate_input_checksums(
        input_dir, fnames_to_include)

    output_lines = ["file_name,md5_hex_checksum\n"]
    output_entries = [f"{k},{fname_and_checksum_dict[k]}\n" for k
                      in sorted(fname_and_checksum_dict)]
    output_lines.extend(output_entries)

    with open(output_fp, "w") as f:
        f.writelines(output_lines)


if __name__ == '__main__':
    generate_checksums_file(argv)
