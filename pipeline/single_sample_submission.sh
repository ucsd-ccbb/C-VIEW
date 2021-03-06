#!/bin/bash

INPUT=$1 # Sample Sheet with header - seq_run,sample,s3download,s3upload,primers,reads
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline

[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
sed 1d $INPUT | while IFS=',' read SAMPLE SEQ_RUN S3DOWNLOAD PRIMER_SET FQ MERGE_LANES TIMESTAMP
do

	echo "Sample: $SAMPLE"
	echo "Seq_Run: $SEQ_RUN"
	echo "S3 Fastq path: $S3DOWNLOAD"
	echo "Primers: $PRIMER_SET"
	echo "Fastq Reads: $FQ"
	echo "Merge: $MERGE_LANES"
	echo "TIMESTAMP: $TIMESTAMP"

	if [[ ! "$PRIMER_SET" =~ ^(artic|swift_v2)$ ]]; then
		echo "Parameter 2 - PRIMER_SET must be one of artic or swift_v2"
		exit 1
	fi

	if [[ ! "$FQ" =~ ^(se|pe)$ ]]; then
	  echo "FQ must be one of se or pe"
	  exit 1
	fi

	qsub $QSUBSAMPLEPARAMS \
		-v SEQ_RUN="$SEQ_RUN" \
		-v SAMPLE=$SAMPLE \
		-v S3DOWNLOAD=$S3DOWNLOAD \
		-v PRIMER_SET=$PRIMER_SET \
		-v MERGE_LANES=$MERGE_LANES \
		-v FQ=$FQ \
		-v TIMESTAMP=$TIMESTAMP \
		-N Covid19_"$SEQ_RUN"_"$SAMPLE" \
		-wd /shared/workspace/projects/covid/logs \
		-pe smp 2 \
		-S /bin/bash \
		$PIPELINEDIR/pipeline/sarscov2_consensus_pipeline.sh

done