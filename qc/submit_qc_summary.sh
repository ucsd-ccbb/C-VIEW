#!/bin/bash


INPUT=$1 # Sample Sheet with header - seq_run,s3download,s3upload,primers,reads
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline

[ ! -f $INPUT ] && { echo "Error: $INPUT file not found"; exit 99; }
sed 1d $INPUT | while IFS=',' read SEQ_RUN S3DOWNLOAD TIMESTAMP FQ
do
	echo "Seq_Run: $SEQ_RUN"
	echo "QC results path: $S3DOWNLOAD/"$SEQ_RUN"_fastq"
	echo "Timestamp of run results: TIMESTAMP"
	echo "Fastq Reads: $FQ"
	WORKSPACE=/scratch/$SEQ_RUN/$TIMESTAMP
	echo "Workspace : $WORKSPACE"
	
	
	# Perform QC and summary on seq_run when all samples have completed
	qsub \
		-v SEQ_RUN=$SEQ_RUN \
		-v S3DOWNLOAD=$S3DOWNLOAD \
		-v WORKSPACE=/scratch/$SEQ_RUN/"$TIMESTAMP"_"$FQ" \
		-v FQ=$FQ \
		-v TIMESTAMP=$TIMESTAMP \
		-N QC_summary_"$SEQ_RUN" \
		-wd /shared/workspace/projects/covid/logs \
		-pe smp 32 \
		-S /bin/bash \
    	$PIPELINEDIR/qc/qc_summary.sh

done