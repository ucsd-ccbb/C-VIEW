#!/bin/bash

export PATH=$PATH:/shared/workspace/software/IQTree/iqtree-2.1.2-Linux/bin:/shared/workspace/software/viralMSA:/shared/workspace/software/MinVar-Rooting-master:/shared/workspace/software/anaconda3/envs/covid1.2/bin
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
ANACONDADIR=/shared/workspace/software/anaconda3/bin
S3HELIX=s3://helix-all
S3UCSD=s3://ucsd-all
S3TEST=s3://ucsd-rtl-test
S3INSPECT=s3://ucsd-inspect
INSPECT_INPUT_FNAME=metadata.csv
THREADS=8
rm -rf $WORKSPACE
mkdir -p $WORKSPACE
mkdir -p $WORKSPACE/passQC
mkdir -p $WORKSPACE/loose_stringent
mkdir -p $WORKSPACE/stringent

echo "$VERSION_INFO" >> $WORKSPACE/"$TIMESTAMP".version.log

if [[ "$ISTEST" == false ]]; then
  if [[ "$ORGANIZATION" == ucsd ]]; then
    # only if NOT is_helix
    aws s3 cp $S3UCSD/phylogeny/cumulative_data/consensus/ $WORKSPACE/  --recursive --quiet
    aws s3 cp $S3UCSD/phylogeny/cumulative_data/historic/ $WORKSPACE/ --recursive --quiet
    S3UPLOAD=$S3UCSD
  else
    S3UPLOAD=$S3HELIX
  fi

  # always download helix data
  aws s3 cp $S3HELIX/phylogeny/cumulative_data/consensus/ $WORKSPACE/  --recursive --quiet
  aws s3 cp $S3HELIX/phylogeny/cumulative_data/historic/ $WORKSPACE/ --recursive --quiet
else
  aws s3 cp $S3TEST/phylogeny/cumulative_data/consensus/ $WORKSPACE/  --recursive --quiet
  aws s3 cp $S3TEST/phylogeny/cumulative_data/historic/ $WORKSPACE/ --recursive --quiet
  S3UPLOAD=$S3TEST
  S3INSPECT=$S3TEST
fi

# TODO: Need real file name
aws s3 cp $S3INSPECT/$INSPECT_INPUT_FNAME $WORKSPACE/$INSPECT_INPUT_FNAME

runPangolin () {

	# start with the reference sequence
	cat $PIPELINEDIR/reference_files/NC_045512.2.fas > $WORKSPACE/"$TIMESTAMP"_refs_hist.fas

	# add a bat coronavirus sequence that is used as the outgroup for tree rooting
	cat $PIPELINEDIR/reference_files/RmYN02.fas >> $WORKSPACE/"$TIMESTAMP"_refs_hist.fas

	# add B.1.1.7 sequence
	cat $PIPELINEDIR/reference_files/hCoV-19_USA_CA-SEARCH-5574_2020.fasta >> $WORKSPACE/"$TIMESTAMP"_refs_hist.fas

  # add the historic fas sequences
	cat $WORKSPACE/*historic.fas >> $WORKSPACE/"$TIMESTAMP"_refs_hist.fas

  # find every fasta header line in the "$TIMESTAMP"_refs_hist.fas,
  # cut off its first char (the >),
  # then put it into the added_fa_names.txt file
  # (we'll use this later for the qc and lineages file creation)
  echo "fasta_id" > $WORKSPACE/added_fa_names.txt
  grep "^>" $WORKSPACE/"$TIMESTAMP"_refs_hist.fas | cut -c 2- >> $WORKSPACE/added_fa_names.txt

  # add the ref and historic fast to the passing fas from all the sequencing runs
  cat $WORKSPACE/*passQC.fas >> $WORKSPACE/"$TIMESTAMP"_passQC.fas
  cat $WORKSPACE/"$TIMESTAMP"_passQC.fas $WORKSPACE/"$TIMESTAMP"_refs_hist.fas >> $WORKSPACE/passQC/"$TIMESTAMP"_passQC_refs_hist.fas

	# pangolin
	source $ANACONDADIR/activate pangolin2
	pangolin --update
	pangolin -t $THREADS --outfile $WORKSPACE/"$TIMESTAMP".lineage_report.csv $WORKSPACE/passQC/"$TIMESTAMP"_passQC_refs_hist.fas
  echo -e "pangolin exit code: $?" >> $WORKSPACE/"$TIMESTAMP"-phylogeny.exit.log

  # produce merged_qc_and_lineages.csv
  python $PIPELINEDIR/qc/lineages_summary.py $WORKSPACE/added_fa_names.txt $WORKSPACE "-summary.csv" $WORKSPACE/"$TIMESTAMP".lineage_report.csv $WORKSPACE/"$TIMESTAMP".qc_and_lineages.csv
  echo -e "lineages_summary.py exit code: $?" >> $WORKSPACE/"$TIMESTAMP"-phylogeny.exit.log

  # merge with inspect metadata to produce full summary, bjorn summary, empress metadata, and winnowed fastas
  python $PIPELINEDIR/qc/metadata_generation.py \
    $WORKSPACE/"$TIMESTAMP".qc_and_lineages.csv \
    $WORKSPACE/$INSPECT_INPUT_FNAME \
    $WORKSPACE/"$TIMESTAMP".full_summary.csv \
    $WORKSPACE/"$TIMESTAMP".bjorn_summary.csv \
    $WORKSPACE/passQC/"$TIMESTAMP"_passQC_refs_hist_empress_metadata.tsv \
    $WORKSPACE/loose_stringent/"$TIMESTAMP"_loose_stringent_refs_hist_empress_metadata.tsv \
    $WORKSPACE/stringent/"$TIMESTAMP"_stringent_refs_hist_empress_metadata.tsv \
    $WORKSPACE/"$TIMESTAMP"_passQC.fas \
    $WORKSPACE/loose_stringent/"$TIMESTAMP"_loose_only.fas \
    $WORKSPACE/stringent/"$TIMESTAMP"_stringent_only.fas

  echo -e "metadata_generation.py exit code: $?" >> $WORKSPACE/"$TIMESTAMP"-phylogeny.exit.log

  # add the refs_hist.fas to the stringent_only.fas
  cat $WORKSPACE/"$TIMESTAMP"_refs_hist.fas $WORKSPACE/stringent/"$TIMESTAMP"_stringent_only.fas >> $WORKSPACE/stringent/"$TIMESTAMP"_stringent_refs_hist.fas

  # add the loose_only.fas to the stringent_refs_hist.fas
  cat $WORKSPACE/stringent/"$TIMESTAMP"_stringent_refs_hist.fas $WORKSPACE/loose_stringent/"$TIMESTAMP"_loose_only.fas >> $WORKSPACE/loose_stringent/"$TIMESTAMP"_loose_stringent_refs_hist.fas

  aws s3 cp $WORKSPACE/"$TIMESTAMP".full_summary.csv $S3INSPECT/full_summary.csv
}

{ time ( runPangolin ) ; } > $WORKSPACE/"$TIMESTAMP"-pangolin.log 2>&1
aws s3 cp $WORKSPACE/"$TIMESTAMP"-pangolin.log $S3UPLOAD/phylogeny/$TIMESTAMP/\

grep -v "exit code: 0" $WORKSPACE/"$TIMESTAMP"-phylogeny.exit.log | head -n 1 >> $WORKSPACE/"$TIMESTAMP"-phylogeny.error.log
aws s3 cp $WORKSPACE/ $S3UPLOAD/phylogeny/$TIMESTAMP/ --recursive --quiet
