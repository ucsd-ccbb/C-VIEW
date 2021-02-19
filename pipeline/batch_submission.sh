#!/bin/bash

S3DOWNLOAD=$1 # e.g. s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/
PRIMER_SET=$2 # artic or swift_v2 - the primer set used to create these reads
FQ=$3 # se or pe - single end or paired end reads
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


SAMPLE_LIST=$(aws s3 ls $S3DOWNLOAD/ | grep fastq | awk '{print $NF}' | awk -F '.fastq' '{print $1}' | sort | uniq | grep -v Undetermined)

for SAMPLE in $SAMPLE_LIST; do
	qsub -v SAMPLE=$SAMPLE \
		-v FQ=$FQ \
		-v PRIMER_SET=$PRIMER_SET \
		-v TIMESTAMP=$TIMESTAMP \
		-v S3DOWNLOAD=$S3DOWNLOAD \
		-N Covid19_"$SAMPLE" \
		-wd /shared/workspace/projects/covid/logs \
		-pe smp 1 \
		-S /bin/bash \
		$PIPELINEDIR/pipeline/sarscov2_consensus_pipeline.sh
done
