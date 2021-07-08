#!/bin/bash

export PATH=$PATH:/shared/workspace/software/IQTree/iqtree-2.1.2-Linux/bin:/shared/workspace/software/viralMSA:/shared/workspace/software/MinVar-Rooting-master:/shared/workspace/software/anaconda3/envs/covid1.2/bin
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
ANACONDADIR=/shared/workspace/software/anaconda3/bin
S3HELIX=s3://helix-all
S3UCSD=s3://ucsd-all
S3TEST=s3://ucsd-rtl-test
S3INSPECT=s3://ucsd-inspect
INSPECT_METADATA_PATTERN="all_samples_search_ids_"
THREADS=8
rm -rf $WORKSPACE
mkdir -p $WORKSPACE
mkdir -p $WORKSPACE/passQC
mkdir -p $WORKSPACE/loose_stringent
mkdir -p $WORKSPACE/stringent

echo "$VERSION_INFO" >> $WORKSPACE/"$PROCESSINGID".version.log

# Based on inputs, decide what to download from where
EXCLUDEMASK=""
INCLUDEMASK="*.*"
if [[ "$SEQ_RUN" != "NA" ]]; then
  EXCLUDEMASK="*"
  INCLUDEMASK="$SEQ_RUN""$INCLUDEMASK"
fi

# echo $EXCLUDEMASK >> $WORKSPACE/"$PROCESSINGID"-phylogeny.exit.log
# echo $INCLUDEMASK >> $WORKSPACE/"$PROCESSINGID"-phylogeny.exit.log

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

# echo $S3UPLOAD >> $WORKSPACE/"$PROCESSINGID"-phylogeny.exit.log
# echo $S3INSPECT >> $WORKSPACE/"$PROCESSINGID"-phylogeny.exit.log
# echo "${DOWNLOAD_BUCKETS[@]}" >> $WORKSPACE/"$PROCESSINGID"-phylogeny.exit.log

# Actually do the downloads
for CURR_BUCKET in "${DOWNLOAD_BUCKETS[@]}"
do
    # echo $CURR_BUCKET >> $WORKSPACE/"$PROCESSINGID"-phylogeny.exit.log
    aws s3 cp $CURR_BUCKET/phylogeny/cumulative_data/consensus/ $WORKSPACE/  --recursive --quiet --exclude "$EXCLUDEMASK" --include "$INCLUDEMASK"
    aws s3 cp $CURR_BUCKET/phylogeny/cumulative_data/historic/ $WORKSPACE/ --recursive --quiet
done

# find the most recently created file on the inspect bucket
# that matches the inspect metadata file naming convention and download it
INSPECT_METADATA_FNAME=$(aws s3 ls s3://ucsd-inspect/ --recursive |  grep $INSPECT_METADATA_PATTERN| sort | tail -n 1 | awk '{print $NF}')
aws s3 cp $S3INSPECT/$INSPECT_METADATA_FNAME $WORKSPACE/$INSPECT_METADATA_FNAME

runPangolin () {

	# start with the reference sequence
	cat $PIPELINEDIR/reference_files/NC_045512.2.fas > $WORKSPACE/"$PROCESSINGID"_refs_hist.fas

	# add a bat coronavirus sequence that is used as the outgroup for tree rooting
	cat $PIPELINEDIR/reference_files/RmYN02.fas >> $WORKSPACE/"$PROCESSINGID"_refs_hist.fas

	# add B.1.1.7 sequence
	cat $PIPELINEDIR/reference_files/hCoV-19_USA_CA-SEARCH-5574_2020.fasta >> $WORKSPACE/"$PROCESSINGID"_refs_hist.fas

  # add the historic fas sequences
	cat $WORKSPACE/*historic.fas >> $WORKSPACE/"$PROCESSINGID"_refs_hist.fas

  # find every fasta header line in the "$PROCESSINGID"_refs_hist.fas,
  # cut off its first char (the >),
  # then put it into the added_fa_names.txt file
  # (we'll use this later for the qc and lineages file creation)
  echo "fasta_id" > $WORKSPACE/added_fa_names.txt
  grep "^>" $WORKSPACE/"$PROCESSINGID"_refs_hist.fas | cut -c 2- >> $WORKSPACE/added_fa_names.txt

  # add the ref and historic fast to the passing fas from all the sequencing runs
  cat $WORKSPACE/*passQC.fas >> $WORKSPACE/"$PROCESSINGID"_passQC.fas
  cat $WORKSPACE/"$PROCESSINGID"_passQC.fas $WORKSPACE/"$PROCESSINGID"_refs_hist.fas >> $WORKSPACE/passQC/"$PROCESSINGID"_passQC_refs_hist.fas

	# pangolin
	source $ANACONDADIR/activate pangolin2
	pangolin --update
	pangolin -t $THREADS --outfile $WORKSPACE/"$PROCESSINGID".lineage_report.csv $WORKSPACE/passQC/"$PROCESSINGID"_passQC_refs_hist.fas
  echo -e "pangolin exit code: $?" >> $WORKSPACE/"$PROCESSINGID"-phylogeny.exit.log

  # produce merged_qc_and_lineages.csv
  python $PIPELINEDIR/qc/lineages_summary.py $WORKSPACE/added_fa_names.txt $WORKSPACE "-summary.csv" $WORKSPACE/"$PROCESSINGID".lineage_report.csv $WORKSPACE/"$PROCESSINGID".qc_and_lineages.csv
  echo -e "lineages_summary.py exit code: $?" >> $WORKSPACE/"$PROCESSINGID"-phylogeny.exit.log

  # merge with inspect metadata to produce full summary, bjorn summary, empress metadata, and winnowed fastas
  python $PIPELINEDIR/qc/metadata_generation.py \
    $WORKSPACE/"$PROCESSINGID".qc_and_lineages.csv \
    $WORKSPACE/$INSPECT_METADATA_FNAME \
    $WORKSPACE/"$PROCESSINGID".full_summary.csv \
    $WORKSPACE/"$PROCESSINGID".bjorn_summary.csv

  echo -e "metadata_generation.py exit code: $?" >> $WORKSPACE/"$PROCESSINGID"-phylogeny.exit.log

  # generate empress metadata and winnowed fastas
  python $PIPELINEDIR/qc/tree_prep.py \
    $WORKSPACE/"$PROCESSINGID".qc_and_lineages.csv \
    $WORKSPACE/$INSPECT_METADATA_FNAME \
    $WORKSPACE/"$PROCESSINGID"_passQC.fas \
    $WORKSPACE/loose_stringent/"$PROCESSINGID"_loose_only.fas \
    $WORKSPACE/stringent/"$PROCESSINGID"_stringent_only.fas \
    $WORKSPACE/passQC/"$PROCESSINGID"_passQC_refs_hist_empress_metadata.tsv \
    $WORKSPACE/loose_stringent/"$PROCESSINGID"_loose_stringent_refs_hist_empress_metadata.tsv \
    $WORKSPACE/stringent/"$PROCESSINGID"_stringent_refs_hist_empress_metadata.tsv

  echo -e "tree_prep.py exit code: $?" >> $WORKSPACE/"$PROCESSINGID"-phylogeny.exit.log

  # add the refs_hist.fas to the stringent_only.fas
  cat $WORKSPACE/"$PROCESSINGID"_refs_hist.fas $WORKSPACE/stringent/"$PROCESSINGID"_stringent_only.fas >> $WORKSPACE/stringent/"$PROCESSINGID"_stringent_refs_hist.fas

  # add the loose_only.fas to the stringent_refs_hist.fas
  cat $WORKSPACE/stringent/"$PROCESSINGID"_stringent_refs_hist.fas $WORKSPACE/loose_stringent/"$PROCESSINGID"_loose_only.fas >> $WORKSPACE/loose_stringent/"$PROCESSINGID"_loose_stringent_refs_hist.fas

  aws s3 cp $WORKSPACE/"$PROCESSINGID".full_summary.csv $S3INSPECT/"$PROCESSINGID".full_summary.csv
}

{ time ( runPangolin ) ; } > $WORKSPACE/"$PROCESSINGID"-pangolin.log 2>&1
aws s3 cp $WORKSPACE/"$PROCESSINGID"-pangolin.log $S3UPLOAD/phylogeny/$PROCESSINGID/\

grep -v "exit code: 0" $WORKSPACE/"$PROCESSINGID"-phylogeny.exit.log | head -n 1 >> $WORKSPACE/"$PROCESSINGID"-phylogeny.error.log
aws s3 cp $WORKSPACE/ $S3UPLOAD/phylogeny/$PROCESSINGID/ --recursive --quiet
