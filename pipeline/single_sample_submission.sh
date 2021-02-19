#!/bin/bash

INPUT=$1 # Sample Sheet with header - batch,sample,s3download,s3upload,primers,reads

[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
sed 1d $INPUT | while IFS=',' read BATCH S3DOWNLOAD S3UPLOAD PRIMER_SET FQ
do
	echo "Batch: $BATCH"
	echo "Sample: $SAMPLE"
	echo "S3 Fastq path: $S3DOWNLOAD"
	echo "S3 Results Path: $S3UPLOAD"
	echo "Primers: $PRIMER_SET"
	echo "Fastq Reads: $FQ"


	TIMESTAMP=$(date +'%Y-%m-%d_%H-%M-%S')
	PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline

	if [[ ! "$PRIMER_SET" =~ ^(artic|swift_v2)$ ]]; then
		echo "Parameter 2 - PRIMER_SET must be one of artic or swift_v2"
		exit 1
	fi

	if [[ ! "$FQ" =~ ^(se|pe)$ ]]; then
	  echo "FQ must be one of se or pe"
	  exit 1
	fi


	qsub -v SAMPLE=$SAMPLE \
		-v S3DOWNLOAD=$S3DOWNLOAD \
		-v S3UPLOAD=$S3UPLOAD \
		-v PRIMER_SET=$PRIMER_SET \
		-v FQ=$FQ \
		-v TIMESTAMP=$TIMESTAMP \
		-N Covid19_"$SAMPLE" \
		-wd /shared/workspace/projects/covid/logs \
		-pe smp 1 \
		-S /bin/bash \
		$PIPELINEDIR/pipeline/sarscov2_consensus_pipeline.sh

done