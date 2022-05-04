#!/bin/bash

export PATH=$PATH:/shared/workspace/software/IQTree/iqtree-2.1.2-Linux/bin:/shared/workspace/software/viralMSA:/shared/workspace/software/MinVar-Rooting-master:/shared/workspace/software/anaconda3/envs/cview/bin
CVIEWDIR=/shared/workspace/software/cview
ANACONDADIR=/shared/workspace/software/anaconda3/bin
S3HELIX=s3://helix-all
S3UCSD=s3://ucsd-all
S3TEST=s3://ucsd-rtl-test
S3INSPECT=s3://ucsd-inspect
INSPECT_METADATA_PATTERN="all_samples_search_ids_"
THREADS=8
rm -rf $WORKSPACE
mkdir -p $WORKSPACE
mkdir -p $WORKSPACE/stringent

# Based on inputs, decide what to download from where
EXCLUDEMASK=""
INCLUDEMASK="*.*"
if [[ "$SEQ_RUN" != all ]]; then
  EXCLUDEMASK="*"
  INCLUDEMASK="$SEQ_RUN""$INCLUDEMASK"
fi

DOWNLOAD_BUCKETS=()
if [[ "$ISTEST" == false ]]; then
  # if this is a real run, always download helix data
  DOWNLOAD_BUCKETS+=($S3HELIX)
  if [[ "$ORGANIZATION" == ucsd ]]; then
    # if we're allowed, also download ucsd
    DOWNLOAD_BUCKETS+=($S3UCSD)
    S3UPLOAD=$S3UCSD
  else
    S3UPLOAD=$S3HELIX
  fi
else
  # for test run, only download test data
  DOWNLOAD_BUCKETS+=($S3TEST)
  S3UPLOAD=$S3TEST
  S3INSPECT=$S3TEST
fi

runLineages () {
  echo "$VERSION_INFO" >> $WORKSPACE/"$PROCESSINGID".version.log

  # Actually do the downloads
  for CURR_BUCKET in "${DOWNLOAD_BUCKETS[@]}"
  do
      aws s3 cp $CURR_BUCKET/phylogeny/cumulative_data/consensus/ $WORKSPACE/  --recursive --quiet --exclude "$EXCLUDEMASK" --include "$INCLUDEMASK"
      aws s3 cp $CURR_BUCKET/phylogeny/cumulative_data/historic/ $WORKSPACE/ --recursive --quiet
  done

  # find the most recently created file on the inspect bucket
  # that matches the inspect metadata file naming convention and download it
  # Note that this is NOT recursive: it is looking in the top directory only;
  # this is because (a) that's where it has always been put and (b) doing recursive
  # screws up when running in test mode where the download and inspect buckets
  # are the same--it picks the copy of the metadata in whatever the latest
  # phylogeny run is :|
  INSPECT_METADATA_FNAME=$(aws s3 ls $S3INSPECT/ |  grep $INSPECT_METADATA_PATTERN| sort | tail -n 1 | awk '{print $NF}')
  # echo $S3INSPECT/$INSPECT_METADATA_FNAME >> $WORKSPACE/"$PROCESSINGID"-lineages.debug.log
  aws s3 cp "$S3INSPECT"/"$INSPECT_METADATA_FNAME" $WORKSPACE/$INSPECT_METADATA_FNAME

	# start with the reference sequence
	cat $CVIEWDIR/reference_files/NC_045512.2.fas > $WORKSPACE/"$PROCESSINGID"_refs_hist.fas

	# add a bat coronavirus sequence that is used as the outgroup for tree rooting
	cat $CVIEWDIR/reference_files/RmYN02.fas >> $WORKSPACE/"$PROCESSINGID"_refs_hist.fas

	# add B.1.1.7 sequence
	cat $CVIEWDIR/reference_files/hCoV-19_USA_CA-SEARCH-5574_2020.fasta >> $WORKSPACE/"$PROCESSINGID"_refs_hist.fas

  # add the historic fas sequences
	cat $WORKSPACE/*historic.fas >> $WORKSPACE/"$PROCESSINGID"_refs_hist.fas

  # find every fasta header line in the "$PROCESSINGID"_refs_hist.fas,
  # cut off its first char (the >),
  # then put it into the added_fa_names.txt file
  # (we'll use this later for the src and lineages file creation)
  echo "fasta_id" > $WORKSPACE/added_fa_names.txt
  grep "^>" $WORKSPACE/"$PROCESSINGID"_refs_hist.fas | cut -c 2- >> $WORKSPACE/added_fa_names.txt

  # add the ref and historic fast to the passing fas from all the sequencing runs
  cat $WORKSPACE/*passQC.fas >> $WORKSPACE/"$PROCESSINGID"_passQC.fas
  cat $WORKSPACE/"$PROCESSINGID"_passQC.fas $WORKSPACE/"$PROCESSINGID"_refs_hist.fas >> $WORKSPACE/"$PROCESSINGID"_passQC_refs_hist.fas

	# pangolin
	source $ANACONDADIR/activate pangolin
	pangolin --update
	pangolin -t $THREADS --analysis-mode fast --outfile $WORKSPACE/"$PROCESSINGID".lineage_report.csv $WORKSPACE/"$PROCESSINGID"_passQC_refs_hist.fas
  echo -e "pangolin exit code: $?" >> $WORKSPACE/"$PROCESSINGID"-lineages.exit.log

  # deactivate the pangolin environment and re-activate the main pipeline environment
  source $ANACONDADIR/deactivate pangolin
  source $ANACONDADIR/activate cview

  # generate file of input checksums, for record-keeping
  python $CVIEWDIR/qc/document_file_checksums.py \
    $WORKSPACE $WORKSPACE/"$PROCESSINGID"_input_checksums.csv \
    "-passQC.fas" "-summary.csv"
  echo -e "document_file_checksums.py exit code: $?" >> $WORKSPACE/"$PROCESSINGID"-lineages.exit.log

  # produce merged qc_and_lineages.csv
  python $CVIEWDIR/qc/lineages_summary.py $WORKSPACE/added_fa_names.txt $WORKSPACE "-summary.csv" $WORKSPACE/"$PROCESSINGID".lineage_report.csv $WORKSPACE/"$PROCESSINGID".qc_and_lineages.csv
  echo -e "lineages_summary.py exit code: $?" >> $WORKSPACE/"$PROCESSINGID"-lineages.exit.log

  # merge with inspect metadata to produce full summary and bjorn summary
  python $CVIEWDIR/qc/metadata_generation.py \
    $WORKSPACE/"$PROCESSINGID".qc_and_lineages.csv \
    $WORKSPACE/$INSPECT_METADATA_FNAME \
    $WORKSPACE/"$PROCESSINGID".full_summary.csv \
    $WORKSPACE/"$PROCESSINGID".bjorn_summary.csv

  echo -e "metadata_generation.py exit code: $?" >> $WORKSPACE/"$PROCESSINGID"-lineages.exit.log

  # generate customized summary file slices for RTL constituents
  python $CVIEWDIR/qc/custom_reports_generation.py \
    $WORKSPACE/"$PROCESSINGID".full_summary.csv \
    $WORKSPACE/"$PROCESSINGID"_summary-report

  echo -e "custom_reports_generation.py exit code: $?" >> $WORKSPACE/"$PROCESSINGID"-lineages.exit.log

  # generate empress metadata and winnowed fastas in preparation for
  # (possible) tree building
  python $CVIEWDIR/qc/tree_prep.py \
    $WORKSPACE/"$PROCESSINGID".qc_and_lineages.csv \
    $WORKSPACE/$INSPECT_METADATA_FNAME \
    $WORKSPACE/"$PROCESSINGID"_passQC.fas \
    $WORKSPACE/stringent/"$PROCESSINGID"_stringent_only.fas \
    $WORKSPACE/stringent/"$PROCESSINGID"_stringent_refs_hist_empress_metadata.tsv

  echo -e "tree_prep.py exit code: $?" >> $WORKSPACE/"$PROCESSINGID"-lineages.exit.log

  # add the refs_hist.fas to the stringent_only.fas
  cat $WORKSPACE/"$PROCESSINGID"_refs_hist.fas $WORKSPACE/stringent/"$PROCESSINGID"_stringent_only.fas >> $WORKSPACE/stringent/"$PROCESSINGID"_stringent_refs_hist.fas

  aws s3 cp $WORKSPACE/"$PROCESSINGID".full_summary.csv $S3INSPECT/"$PROCESSINGID".full_summary.csv
}

{ time ( runLineages ) ; } > $WORKSPACE/"$PROCESSINGID"-lineages.log 2>&1
aws s3 cp $WORKSPACE/"$PROCESSINGID"-lineages.log $S3UPLOAD/phylogeny/$PROCESSINGID/\

grep -v "exit code: 0" $WORKSPACE/"$PROCESSINGID"-lineages.exit.log | head -n 1 >> $WORKSPACE/"$PROCESSINGID"-lineages.error.log
aws s3 cp $WORKSPACE/ $S3UPLOAD/phylogeny/$PROCESSINGID/ --recursive --quiet --include "*"  --exclude "*.fas" --exclude "*-summary.csv" --include "*_passQC_refs_hist.fas" --include "*_stringent_refs_hist.fas"
