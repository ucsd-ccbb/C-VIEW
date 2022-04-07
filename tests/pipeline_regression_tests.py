from sys import argv
import os
import boto3
#from wand.image import Image
import pandas as pd
import filecmp

# remember to run aws configure from the command line before running this
# script to record the AWS access for your user

_S3_PREFIX = "s3://"


def _split_s3_url(s3_results_url):
    results_url = s3_results_url.replace(_S3_PREFIX, "")
    bucket_name = results_url.split("/")[0]
    results_path = results_url.replace(f"{bucket_name}/", "")
    if not results_path.endswith("/"):
        results_path += "/"
    return bucket_name, results_path


def compare_pdfs(expected_pdf_fp, real_pdf_fp):
    raise ValueError("Comparing pdfs is not currently supported")
    # real_img = Image(filename=real_pdf_fp, resolution=150)
    # with Image(filename=expected_pdf_fp,
    #            resolution=150) as expected_img:
    #     imgs_diff = real_img.compare(
    #         expected_img, metric='root_mean_square')
    #     return imgs_diff[1]


def regression_test(gs_fp, temp_dir, s3_results_url):
    s3_bucket, s3_results_dir = _split_s3_url(s3_results_url)
    s3 = boto3.client('s3')

    # take the input dir of the known-good files
    for base_fp, sub_dirs, fnames in os.walk(gs_fp):
        for curr_fname in fnames:
            if curr_fname.startswith("."):
                continue

            print(curr_fname)
            curr_fp = os.path.join(base_fp, curr_fname)
            rel_fp = os.path.relpath(curr_fp, start=gs_fp)
            s3_fp = os.path.join(s3_results_dir, rel_fp)
            temp_fp = os.path.join(temp_dir, curr_fname)

            # generate its path in aws and download the aws version to cwd
            s3.download_file(s3_bucket, s3_fp, temp_fp)

            _, curr_fext = os.path.splitext(curr_fp)
            curr_files_equal = False
            try:
                if curr_fext == ".pdf":
                    img_diff = compare_pdfs(curr_fp, temp_fp)
                    curr_files_equal = img_diff < 0.01
                elif curr_fext == ".csv" or curr_fext == ".tsv":
                    sep = "," if curr_fext == ".csv" else "\t"
                    curr_df = pd.read_csv(curr_fp, sep=sep, dtype=str)
                    temp_df = pd.read_csv(temp_fp, sep=sep, dtype=str)

                    timestamp_col = "timestamp"
                    if timestamp_col in curr_df.columns:
                        curr_df.pop(timestamp_col)
                        temp_df.pop(timestamp_col)

                    curr_files_equal = curr_df.equals(temp_df)
                else:
                    curr_files_equal = filecmp.cmp(curr_fp, temp_fp)
            finally:
                if curr_files_equal:
                    try:
                        os.remove(temp_fp)
                    except OSError:
                        pass
                else:
                    print("F\n")


if __name__ == '__main__':
    # TODO: remove test arguments!
    argv = ["pipeline_regression_tests.py", "/Users/abirmingham/Work/Repositories/cview/tests/data/gold_standard", "/Users/abirmingham/Desktop/regression_temp", "s3://ucsd-rtl-test/2021-02-08-ARTIC_partial/2021-02-08-ARTIC_partial_results/2022-04-07_19-53-44_pe/"]
    regression_test(argv[1], argv[2], argv[3])
