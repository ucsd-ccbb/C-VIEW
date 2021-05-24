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

echo "$VERSION_INFO" >> $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist.version.log

if [[ "$ISTEST" == true ]]; then
  S3DOWNLOAD=$S3TEST
  S3UPLOAD=$S3TEST
else
  if [[ "$ORGANIZATION" == ucsd ]]; then
    S3DOWNLOAD=$S3UCSD
    S3UPLOAD=$S3UCSD
  else
    S3DOWNLOAD=$S3HELIX
    S3UPLOAD=$S3HELIX
  fi
fi

echo "ORGANIZATION is $ORGANIZATION"
echo "S3DOWNLOAD is $S3DOWNLOAD"

aws s3 cp $S3DOWNLOAD/phylogeny/$TIMESTAMP/"$DATASET"/"$TIMESTAMP"_"$DATASET"_refs_hist.fas $WORKSPACE/
aws s3 cp $S3DOWNLOAD/phylogeny/$TIMESTAMP/"$DATASET"/"$TIMESTAMP"_"$DATASET"_refs_hist_empress_metadata.tsv $WORKSPACE/

buildTree () {
	# Must use biopy env due to numpy conflicts
	source $ANACONDADIR/activate biopy
	ViralMSA.py -s $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist.fas -r SARS-CoV-2 -o $WORKSPACE/viralmsa_out -t $THREADS -e aws-CCBB@health.ucsd.edu
    echo -e "ViralMSA.py exit code: $?" >> $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist-phylogeny.exit.log

	python $PIPELINEDIR/pipeline/trim_msa.py -i $WORKSPACE/viralmsa_out/"$TIMESTAMP"_"$DATASET"_refs_hist.fas.aln -s 100 -e 50 -o $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist.trimmed.aln
    echo -e "trim_msa.py exit code: $?" >> $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist-phylogeny.exit.log

	iqtree2 -T $THREADS -m GTR+F+G4 --polytomy -blmin 1e-9 -s $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist.trimmed.aln
    echo -e "iqtree2 exit code: $?" >> $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist-phylogeny.exit.log

	python /shared/workspace/software/MinVar-Rooting-master/FastRoot.py -i $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist.trimmed.aln.treefile -o $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist.trimmed.aln.rooted.treefile -m OG -g "hCoV-19/bat/Yunnan/RmYN02/2019|EPI_ISL_412977|2019-06-25"
    echo -e "iFastRoot.py exit code: $?" >> $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist-phylogeny.exit.log

	# tree building 
	source $ANACONDADIR/activate qiime2-2020.11

	empress tree-plot --tree $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist.trimmed.aln.rooted.treefile --feature-metadata $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist_empress_metadata.tsv --output-dir $WORKSPACE/"$TIMESTAMP"_"$DATASET"_tree-viz
    echo -e "empress tree-plot exit code: $?" >> $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist-phylogeny.exit.log
  mv $WORKSPACE/tree-viz/empress.html $WORKSPACE/tree-viz/"$TIMESTAMP"_"$DATASET"_empress.html
}

{ time ( buildTree ) ; } > $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist-treebuild.log 2>&1
aws s3 cp $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist-treebuild.log $S3UPLOAD/phylogeny/$TIMESTAMP/$DATASET/


grep -v "exit code: 0" $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist-phylogeny.exit.log | head -n 1 >> $WORKSPACE/"$TIMESTAMP"_"$DATASET"_refs_hist-phylogeny.error.log
aws s3 cp $WORKSPACE/ $S3UPLOAD/phylogeny/$TIMESTAMP/$DATASET/ --recursive --quiet
