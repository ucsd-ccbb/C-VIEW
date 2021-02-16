from sys import argv, stderr
import os
import filecmp
import boto3

# remember to run aws configure from the command line before running this
# script to record the AWS access for your user

_GOLD_STANDARD_DIR = "./tests/data/gold_standard_2021-02-08-ARTIC"
_FNAME_CRUFT = ".trimmed.sorted.pileup."
_CONSENSUS_FNAME = "consensus.fa"
_VARIANTS_FNAME = "variants.tsv"

_S3_PREFIX = "s3://"
_AWS_CONSENSUS_DIRNAME = "consensus"
_AWS_VARIANT_DIRNAME = "variants"

_TEST_FNAME_PREFIX = "test_"
_GSC = "gs_consensus"
_GSV = "gs_variants"
_TC = "test_consensus"
_TV = "test_variants"


def regression_test(s3_results_url):
    curr_match = None
    num_match = num_mismatch = 0
    result_details = {}

    gold_standard_basenames = os.listdir(_GOLD_STANDARD_DIR)
    for curr_gs_basename in sorted(gold_standard_basenames):
        if not "SEARCH-" in curr_gs_basename:
            continue

        # compare consensus sequence and variant calls ONLY
        fnames = _formulate_curr_fnames(curr_gs_basename)

        try:
            _download_files_to_cwd(s3_results_url, fnames)

            curr_match, result_details[curr_gs_basename] = \
                _compare_files(curr_gs_basename, fnames)
        finally:
            # remove any and all of the to-test files downloaded
            for curr_key in [_TC, _TV]:
                try:
                    os.remove(fnames[curr_key])
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


def _split_s3_url(s3_results_url):
    results_url = s3_results_url.replace(_S3_PREFIX, "")
    bucket_name = results_url.split("/")[0]
    results_path = results_url.replace(f"{bucket_name}/", "")
    if not results_path.endswith("/"):
        results_path += "/"
    return bucket_name, results_path


def _formulate_curr_fnames(curr_gs_basename):
    fnames = {}
    fnames[_GSC] = f"{curr_gs_basename}" \
                   f"{_FNAME_CRUFT}{_CONSENSUS_FNAME}"
    fnames[_GSV] = f"{curr_gs_basename}{_FNAME_CRUFT}" \
                   f"{_VARIANTS_FNAME}"
    fnames[_TC] = f"{_TEST_FNAME_PREFIX}{fnames[_GSC]}"
    fnames[_TV] = f"{_TEST_FNAME_PREFIX}{fnames[_GSV]}"
    return fnames


def _formulate_aws_fp(s3_results_path, fnames, fname_key):
    dirname = gs_key = None
    if fname_key == _TC:
        dirname = _AWS_CONSENSUS_DIRNAME
        gs_key = _GSC
    elif fname_key == _TV:
        dirname = _AWS_VARIANT_DIRNAME
        gs_key = _GSV
    return os.path.join(s3_results_path, dirname, fnames[gs_key])


def _download_files_to_cwd(s3_results_url, fnames):
    s3_bucket, s3_results_path = _split_s3_url(s3_results_url)
    aws_results_curr_cons_fp = _formulate_aws_fp(s3_results_path, fnames, _TC)
    aws_results_curr_var_fp = _formulate_aws_fp(s3_results_path, fnames, _TV)

    s3 = boto3.client('s3')
    s3.download_file(s3_bucket, aws_results_curr_cons_fp, fnames[_TC])
    s3.download_file(s3_bucket, aws_results_curr_var_fp, fnames[_TV])


def _compare_files(curr_gs_basename, fnames):
    gs_cons_fp = os.path.join(_GOLD_STANDARD_DIR, curr_gs_basename,
                              fnames[_GSC])
    curr_cons_equal = filecmp.cmp(gs_cons_fp, fnames[_TC])

    gs_var_fp = os.path.join(_GOLD_STANDARD_DIR, curr_gs_basename,
                             fnames[_GSV])
    curr_var_equal = filecmp.cmp(gs_var_fp, fnames[_TV])

    curr_match = curr_cons_equal and curr_var_equal
    return curr_match, [curr_gs_basename, curr_cons_equal, curr_var_equal]


if __name__ == '__main__':
    if len(argv) != 2 or not argv[1].startswith(_S3_PREFIX):
        print(f"USAGE: {argv[0]} <{_S3_PREFIX} url of results to test>",
              file=stderr)
        exit(1)

    regression_test(argv[1])
