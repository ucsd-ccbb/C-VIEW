#!/bin/bash

export PATH=$PATH:/shared/workspace/software/IQTree/iqtree-2.1.2-Linux/bin:/shared/workspace/software/viralMSA:/shared/workspace/software/MinVar-Rooting-master:/shared/workspace/software/anaconda3/envs/covid1.2/bin
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
ANACONDADIR=/shared/workspace/software/anaconda3/bin
S3HELIX=s3://helix-all
S3UCSD=s3://ucsd-all
S3TEST=s3://ucsd-rtl-test
THREADS=8
rm -rf $WORKSPACE
mkdir -p $WORKSPACE

if [[ "$ISTEST" = false ]]; then
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
fi

# TODO: Download inspect metadata file

	
runPangolin () {

	# start with the reference sequence
	cat $PIPELINEDIR/reference_files/NC_045512.2.fas > $WORKSPACE/"$TIMESTAMP".fas

	# add a bat coronavirus sequence that is used as the outgroup for tree rooting
	cat $PIPELINEDIR/reference_files/RmYN02.fas >> $WORKSPACE/"$TIMESTAMP".fas

	# add B.1.1.7 sequence
	cat $PIPELINEDIR/reference_files/hCoV-19_USA_CA-SEARCH-5574_2020.fasta >> $WORKSPACE/"$TIMESTAMP".fas

  # add the historic fas sequences
	cat $WORKSPACE/*historic.fas >> $WORKSPACE/"$TIMESTAMP".fas

  # BEFORE adding all the passQC fas files from the runs, stop and
  # find every fasta header line in the "$TIMESTAMP".fas,
  # cut off its first char (the >),
  # then put it into the added_fa_names.txt file
  # (we'll use this later for the qc and lineages file creation)
  echo "fasta_id" > $WORKSPACE/added_fa_names.txt
  grep "^>" $WORKSPACE/"$TIMESTAMP".fas | cut -c 2- >> $WORKSPACE/added_fa_names.txt

    # add in the passing fas from all the sequencing runs
	cat $WORKSPACE/*passQC.fas >> $WORKSPACE/"$TIMESTAMP".fas

	# pangolin
	source $ANACONDADIR/activate pangolin
	pangolin --update
	pangolin -t $THREADS --outfile $WORKSPACE/"$TIMESTAMP".lineage_report.csv $WORKSPACE/"$TIMESTAMP".fas
  echo -e "pangolin exit code: $?" >> $WORKSPACE/"$TIMESTAMP"-phylogeny.exit.log

  # TODO: expand script to take in inspect metadata and qc+lineages file,
  #  output merged file and Andersen lab sample_sheet_metadata.csv

  # produce merged_qc_and_lineages.csv and "$TIMESTAMP".empress_metadata.tsv
  python $PIPELINEDIR/qc/lineages_summary.py $WORKSPACE/added_fa_names.txt $WORKSPACE "-summary.csv" $WORKSPACE/"$TIMESTAMP".lineage_report.csv $WORKSPACE/"$TIMESTAMP".qc_and_lineages.csv $WORKSPACE/"$TIMESTAMP".empress_metadata.tsv
  echo -e "lineages_summary.py exit code: $?" >> $WORKSPACE/"$TIMESTAMP"-phylogeny.exit.log

  # TODO: upload new merged file to inspect bucket
}

buildTree () {
	# Must use biopy env due to numpy conflicts
	source $ANACONDADIR/activate biopy
	ViralMSA.py -s $WORKSPACE/"$TIMESTAMP".fas -r SARS-CoV-2 -o $WORKSPACE/viralmsa_out -t $THREADS -e aws-CCBB@health.ucsd.edu
    echo -e "ViralMSA.py exit code: $?" >> $WORKSPACE/"$TIMESTAMP"-phylogeny.exit.log

	python $PIPELINEDIR/pipeline/trim_msa.py -i $WORKSPACE/viralmsa_out/"$TIMESTAMP".fas.aln -s 100 -e 50 -o $WORKSPACE/"$TIMESTAMP".trimmed.aln
    echo -e "trim_msa.py exit code: $?" >> $WORKSPACE/"$TIMESTAMP"-phylogeny.exit.log

	iqtree2 -T $THREADS -m GTR+F+G4 --polytomy -blmin 1e-9 -s $WORKSPACE/"$TIMESTAMP".trimmed.aln
    echo -e "iqtree2 exit code: $?" >> $WORKSPACE/"$TIMESTAMP"-phylogeny.exit.log

	python /shared/workspace/software/MinVar-Rooting-master/FastRoot.py -i $WORKSPACE/"$TIMESTAMP".trimmed.aln.treefile -o $WORKSPACE/"$TIMESTAMP".trimmed.aln.rooted.treefile -m OG -g "hCoV-19/bat/Yunnan/RmYN02/2019|EPI_ISL_412977|2019-06-25"
    echo -e "iFastRoot.py exit code: $?" >> $WORKSPACE/"$TIMESTAMP"-phylogeny.exit.log

	# tree building 
	source $ANACONDADIR/activate qiime2-2020.11

	empress tree-plot --tree $WORKSPACE/"$TIMESTAMP".trimmed.aln.rooted.treefile --feature-metadata $WORKSPACE/"$TIMESTAMP".empress_metadata.tsv --output-dir $WORKSPACE/tree-viz
    echo -e "empress tree-plot exit code: $?" >> $WORKSPACE/"$TIMESTAMP"-phylogeny.exit.log
}

# Always run Pangolin
{ time ( runPangolin ) ; } > $WORKSPACE/"$TIMESTAMP"-pangolin.log 2>&1
aws s3 cp $WORKSPACE/"$TIMESTAMP"-pangolin.log $S3UPLOAD/phylogeny/$TIMESTAMP/

# Run tree building if TREE_BUILD == true
if [[ "$TREE_BUILD" == true ]]; then
	{ time ( buildTree ) ; } > $WORKSPACE/"$TIMESTAMP"-treebuild.log 2>&1
	aws s3 cp $WORKSPACE/"$TIMESTAMP"-treebuild.log $S3UPLOAD/phylogeny/$TIMESTAMP/
fi

CURRDIR=$(pwd)
cd $PIPELINEDIR
bash $PIPELINEDIR/show_version.sh >> $WORKSPACE/"$TIMESTAMP".version.log
echo -e "show_version.sh exit code: $?" >> $WORKSPACE/"$TIMESTAMP"-phylogeny.exit.log
cd $CURRDIR

grep -v "exit code: 0" $WORKSPACE/"$TIMESTAMP"-phylogeny.exit.log | head -n 1 >> $WORKSPACE/"$TIMESTAMP"-phylogeny.error.log
aws s3 cp $WORKSPACE/ $S3UPLOAD/phylogeny/$TIMESTAMP/ --recursive --quiet
