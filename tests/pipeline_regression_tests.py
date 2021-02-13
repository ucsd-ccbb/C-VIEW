import os
from zipfile import ZipFile
import filecmp
import boto3

# remember to run aws configure from the command line before running this
# script to record the AWS access for your user

AWS_RESULTS_S3_BUCKET = "ucsd-ccbb-projects"
AWS_RESULTS_S3_PATH = "2021/20210208_COVID_sequencing/2021-02-08-ARTIC/results_20210212/"

_GOLD_STANDARD_DIR = "./tests/data/gold_standard_2021-02-08-ARTIC"
_FNAME_CRUFT = ".trimmed.sorted.pileup."
_CONSENSUS_FNAME = "consensus.fa"
_VARIANTS_FNAME = "variants.tsv"

_AWS_CONSENSUS_DIRNAME = "consensus"
_AWS_VARIANT_DIRNAME = "variants"

_TEST_FNAME_PREFIX = "test_"
_ZIP_EXT = ".zip"
_GSC = "gs_consensus"
_GSV = "gs_variants"
_TC = "test_consensus"
_TV = "test_variants"


def temp():
    curr_match = None
    num_match = num_mismatch = 0
    result_details = {}

    gold_standard_filenames = os.listdir(_GOLD_STANDARD_DIR)
    for curr_gs_fname in sorted(gold_standard_filenames):
        if curr_gs_fname.endswith(_ZIP_EXT):
            #  compare consensus sequence and variant calls ONLY
            curr_gs_basename = curr_gs_fname.replace(_ZIP_EXT, "")
            fnames = _formulate_curr_fnames(curr_gs_basename)

            try:
                _load_files_to_cwd(curr_gs_basename, fnames)

                curr_match, result_details[curr_gs_basename] = \
                    _compare_files(curr_gs_basename, fnames)
            finally:
                # remove any and all of the local compared files
                for curr_fname in fnames.values():
                    try:
                        os.remove(curr_fname)
                    except OSError:
                        pass

            if curr_match:
                num_match += 1
                print(".", end='')
            else:
                num_mismatch +=1
                print("F", end='')

    print("")
    print(f"matches: {num_match}")
    print(f"non-matches: {num_mismatch}")
    if num_mismatch > 0:
        for curr_key in sorted(result_details):
            print(result_details[curr_key])


def _formulate_curr_fnames(curr_gs_basename):
    fnames = {}
    fnames[_GSC] = f"{curr_gs_basename}" \
                   f"{_FNAME_CRUFT}{_CONSENSUS_FNAME}"
    fnames[_GSV] = f"{curr_gs_basename}{_FNAME_CRUFT}" \
                   f"{_VARIANTS_FNAME}"
    fnames[_TC] = f"{_TEST_FNAME_PREFIX}{fnames[_GSC]}"
    fnames[_TV] = f"{_TEST_FNAME_PREFIX}{fnames[_GSV]}"
    return fnames


def _formulate_aws_fp(fnames, fname_key):
    dirname = gs_key = None
    if fname_key == _TC:
        dirname = _AWS_CONSENSUS_DIRNAME
        gs_key = _GSC
    elif fname_key == _TV:
        dirname = _AWS_VARIANT_DIRNAME
        gs_key = _GSV
    return os.path.join(AWS_RESULTS_S3_PATH, dirname, fnames[gs_key])


def _load_files_to_cwd(curr_gs_basename, fnames):
    curr_zip_fp = os.path.join(_GOLD_STANDARD_DIR,
                               f"{curr_gs_basename}{_ZIP_EXT}")
    with ZipFile(curr_zip_fp) as curr_zip:
        curr_zip.extract(fnames[_GSC])
        curr_zip.extract(fnames[_GSV])

    s3 = boto3.client('s3')
    aws_results_curr_cons_fp = _formulate_aws_fp(fnames, _TC)
    aws_results_curr_var_fp = _formulate_aws_fp(fnames, _TV)
    s3.download_file(AWS_RESULTS_S3_BUCKET, aws_results_curr_cons_fp,
                     fnames[_TC])
    s3.download_file(AWS_RESULTS_S3_BUCKET, aws_results_curr_var_fp,
                     fnames[_TV])


def _compare_files(curr_gs_basename, fnames):
    curr_cons_equal = filecmp.cmp(fnames[_GSC], fnames[_TC])
    curr_var_equal = filecmp.cmp(fnames[_GSV], fnames[_TV])
    curr_match = curr_cons_equal and curr_var_equal
    return curr_match, [curr_gs_basename, curr_cons_equal, curr_var_equal]


if __name__ == '__main__':
    temp()
