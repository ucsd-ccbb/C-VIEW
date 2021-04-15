import sys
import os
import subprocess


def main(args):
    seq_run = args[1]
    s3_path = args[2]

    if s3_path[-1] != '/':  # needs to end with /
        s3_path = s3_path + '/'
    print(f"s3 path: {s3_path}")

    # capture the output of aws s3 ls, parse out the old sequence names
    s3_ls_bytestr = subprocess.check_output("aws s3 ls " + s3_path, shell=True)
    s3_ls_str = s3_ls_bytestr.decode("utf-8")

    # the last entry is blank, don't include it
    s3_ls_lines = s3_ls_str.split("\n")[:-1]

    # skip subdirectories (look for '/')
    relevant_lines = [s for s in s3_ls_lines if s.find('/') == -1]

    # the 3rd entry is the old filename.
    # split with no arguments splits on whitespace.
    old_names = [s.split()[3] for s in relevant_lines]

    # loop over all old_names and rename the files
    for curr_old_name in old_names:
        # only rename if it is a fastq
        if curr_old_name.find('.fastq.gz') > -1:
            # only add seq_run if not already in name
            if curr_old_name.find(seq_run) != -1:
                print(f"Name {curr_old_name} already has sequence name "
                      f"inserted and cannot be rewritten")
            else:
                curr_new_name = add_seq_run(curr_old_name, seq_run)

                # if the file was renameable, then move
                if curr_new_name is not None:
                    print(f"{curr_old_name} : {curr_new_name}")

                    # rename the file
                    cmd = 'aws s3 mv ' + s3_path + curr_old_name + ' ' + s3_path + curr_new_name
                    #print(cmd)

                    # TODO: uncomment this when script is fully debugged
                    #os.system(cmd)


def add_seq_run(old_name, seq_run):

    # example old_name: 'SEARCH-16659__D101866__O02__S6_L001_R1_001.fastq.gz'
    # example seq_run: '210411_A00953_0274_BH3MJMDSX2'
    DELIMITER = "__"
    new_name = None
    old_name_split = old_name.split(DELIMITER)
    if len(old_name_split) != 4:
        print(f"Name '{old_name}' does not follow naming convention and "
              f"cannot be rewritten")
    else:
        first_piece = DELIMITER.join([old_name_split[i] for i in [0, 1, 2]])
        new_pieces = [first_piece, seq_run, old_name_split[3]]
        new_name = DELIMITER.join(new_pieces)

    return new_name


if __name__ == "__main__":
    # example_inputs = ["add_seq_run_to_fastq_names.py",
    #                   '210411_A00953_0274_BH3MJMDSX2',
    #                   's3://helix-all/210411_A00953_0274_BH3MJMDSX2/'
    #                   '210411_A00953_0274_BH3MJMDSX2_fastq/']

    main(sys.argv)
