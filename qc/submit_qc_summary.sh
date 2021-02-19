#!/bin/bash


INPUT=$1 # Sample Sheet with header - batch,s3download,s3upload,primers,reads
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
WORKSPACE=/scratch/$(date +"%T" | sed 's/:/-/g')

[ ! -f $INPUT ] && { echo "Error: $INPUT file not found"; exit 99; }
sed 1d $INPUT | while IFS=',' read BATCH S3DOWNLOAD S3UPLOAD PRIMER_SET FQ
do
	echo "Batch: $BATCH"
	echo "S3 Fastq path: $S3DOWNLOAD"
	echo "S3 Results Path: $S3UPLOAD"
	echo "Primers: $PRIMER_SET"
	echo "Fastq Reads: $FQ"


	if [[ ! "$PRIMER_SET" =~ ^(artic|swift_v2)$ ]]; then
		echo "Error: Parameter PRIMER_SET must be one of artic or swift_v2"
		exit 1
	fi

	if [[ ! "$FQ" =~ ^(se|pe)$ ]]; then
	  echo "Error: FQ must be one of se or pe"
	  exit 1
	fi

	
	qsub -v S3DOWNLOAD=$S3DOWNLOAD \
		 -v S3UPLOAD=$S3UPLOAD \
		 -v BATCH=$BATCH \
		 -v WORKSPACE=$WORKSPACE \
		 -wd /shared/workspace/projects/covid/logs \
		 -pe smp 1 \
		 -S /bin/bash \
	     $PIPELINEDIR/qc/qc_summary.sh

done