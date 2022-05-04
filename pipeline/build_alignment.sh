#!/bin/bash

export PATH=$PATH:/shared/workspace/software/viralMSA:/shared/workspace/software/MinVar-Rooting-master:/shared/workspace/software/anaconda3/envs/cview/bin
CVIEWDIR=/shared/workspace/software/cview
ANACONDADIR=/shared/workspace/software/anaconda3/bin
S3HELIX=s3://helix-all
S3UCSD=s3://ucsd-all
S3TEST=s3://ucsd-rtl-test
S3INSPECT=s3://ucsd-inspect

THREADS=8
rm -rf $WORKSPACE
mkdir -p $WORKSPACE

echo "$VERSION_INFO" >> $WORKSPACE/"$PROCESSINGID"_"$DATASET"_refs_hist.version.log

if [[ "$ISTEST" == true ]]; then
  S3DOWNLOAD=$S3TEST
  S3UPLOAD=$S3TEST
  S3INSPECT=$S3TEST
else
  if [[ "$ORGANIZATION" == ucsd ]]; then
    S3DOWNLOAD=$S3UCSD
    S3UPLOAD=$S3UCSD
  else
    S3DOWNLOAD=$S3HELIX
    S3UPLOAD=$S3HELIX
  fi
fi

aws s3 cp $S3DOWNLOAD/phylogeny/$PROCESSINGID/"$DATASET"/"$PROCESSINGID"_"$DATASET"_refs_hist.fas $WORKSPACE/
#aws s3 cp $S3DOWNLOAD/phylogeny/$PROCESSINGID/"$DATASET"/"$PROCESSINGID"_"$DATASET"_refs_hist_empress_metadata.tsv $WORKSPACE/

buildTree () {
	# Must use biopy env due to numpy conflicts
	source $ANACONDADIR/activate biopy
	ViralMSA.py -s $WORKSPACE/"$PROCESSINGID"_"$DATASET"_refs_hist.fas -r SARS-CoV-2 -o $WORKSPACE/viralmsa_out -t $THREADS -e aws-CCBB@health.ucsd.edu
    echo -e "ViralMSA.py exit code: $?" >> $WORKSPACE/"$PROCESSINGID"_"$DATASET"_refs_hist-phylogeny.exit.log

	python $CVIEWDIR/src/trim_msa.py -i $WORKSPACE/viralmsa_out/"$PROCESSINGID"_"$DATASET"_refs_hist.fas.aln -s 100 -e 50 -o $WORKSPACE/"$PROCESSINGID"_"$DATASET"_refs_hist.trimmed.aln
    echo -e "trim_msa.py exit code: $?" >> $WORKSPACE/"$PROCESSINGID"_"$DATASET"_refs_hist-phylogeny.exit.log

	source $ANACONDADIR/deactivate
}

{ time ( buildTree ) ; } > $WORKSPACE/"$PROCESSINGID"_"$DATASET"_refs_hist-treebuild.log 2>&1
aws s3 cp $WORKSPACE/"$PROCESSINGID"_"$DATASET"_refs_hist-treebuild.log $S3UPLOAD/phylogeny/$PROCESSINGID/$DATASET/


grep -v "exit code: 0" $WORKSPACE/"$PROCESSINGID"_"$DATASET"_refs_hist-phylogeny.exit.log | head -n 1 >> $WORKSPACE/"$PROCESSINGID"_"$DATASET"_refs_hist-phylogeny.error.log
aws s3 cp $WORKSPACE/ $S3UPLOAD/phylogeny/$PROCESSINGID/$DATASET/ --recursive --quiet
