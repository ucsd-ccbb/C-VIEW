#!/bin/bash

S3DOWNLOAD=$1 # e.g. s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/
FQ=$2 # se or pe - single end or paired end reads
MERGED=$3 # merged or unmerged - merge lanes
IS_ARTIC=$4 # true or false - is this a pipeline test-running on ARTIC samples?

# Make sure MERGED/FQ parameters work
if [[ ! "$FQ" =~ ^(se|pe)$ ]]; then
	exit 1
fi

if [[ ! "$MERGED" =~ ^(merged|unmerged)$ ]]; then
	exit 1
fi

if [[ ! "$IS_ARTIC" =~ ^(true|false)$ ]]; then
	exit 1
fi

if $MERGED == merged; then
	SAMPLE_LIST=$(aws s3 ls $S3DOWNLOAD/ | grep fastq | awk '{print $NF}' | awk -F '_L00' '{print $1}' | sort | uniq | grep -v Undetermined)
else
	SAMPLE_LIST=$(aws s3 ls $S3DOWNLOAD/ | grep fastq | awk '{print $NF}' | awk -F '_R' '{print $1}' | sort | uniq | grep -v Undetermined)
fi

for SAMPLE in $SAMPLE_LIST; do
	qsub -v SAMPLE=$SAMPLE \
		 -v FQ=$FQ \
		 -v MERGED=$MERGED \
		 -v IS_ARTIC=$IS_ARTIC \
		 -v S3DOWNLOAD=$S3DOWNLOAD \
		 -N $SAMPLE \
		 -wd /shared/workspace/projects/covid/logs \
		 -pe smp 1 \
		 -S /bin/bash \
		 /shared/workspace/software/covid_sequencing_analysis_pipeline/covid.sh
done

# When all samples have finished
# qsub -v S3DOWNLOAD=$S3DOWNLOAD \
# 	 -v FQ=$FQ \
# 	 -v MERGED=$MERGED \
#      /shared/workspace/software/covid_sequencing_analysis_pipeline/qc_summary.sh

     