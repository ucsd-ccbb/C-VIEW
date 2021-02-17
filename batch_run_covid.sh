#!/bin/bash

S3DOWNLOAD=$1 # e.g. s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/
IS_ARTIC=$2 # true or false - is this a pipeline test-running on ARTIC samples?
FQ=$3 # se or pe - single end or paired end reads
MERGED=$4 # merged or unmerged - merge lanes

if [[ ! "$IS_ARTIC" =~ ^(true|false)$ ]]; then
	echo "Parameter 2 - IS_ARTIC must be one of true or false"
	exit 1
fi

if [[ "$IS_ARTIC" == true ]]; then
	FQ=NA
	MERGED=NA

	SAMPLE_LIST=$(aws s3 ls $S3DOWNLOAD/ | grep fastq | awk '{print $NF}' | awk -F '_R' '{print $1}' | sort | uniq | grep -v Undetermined)

else
	
	if [[ "$IS_ARTIC" == false ]]; then

	# Make sure IS_ARTIC/MERGED/FQ parameters work
		if [[ ! "$FQ" =~ ^(se|pe)$ ]]; then
			echo "FQ must be one of se or pe"
			exit 1
		fi
		if [[ ! "$MERGED" =~ ^(merged|unmerged)$ ]]; then
			echo "MERGED must be one of merged or unmerged"
			exit 1
		fi

	if [[ "$MERGED" == merged ]]; then
		SAMPLE_LIST=$(aws s3 ls $S3DOWNLOAD/ | grep fastq | awk '{print $NF}' | awk -F '_L00' '{print $1}' | sort | uniq | grep -v Undetermined)
	else
		SAMPLE_LIST=$(aws s3 ls $S3DOWNLOAD/ | grep fastq | awk '{print $NF}' | awk -F '_R' '{print $1}' | sort | uniq | grep -v Undetermined)
	fi



	fi
fi


for SAMPLE in $SAMPLE_LIST; do
	qsub -v SAMPLE=$SAMPLE \
		-v FQ=$FQ \
		-v MERGED=$MERGED \
		-v IS_ARTIC=$IS_ARTIC \
		-v S3DOWNLOAD=$S3DOWNLOAD \
		-N Covid19_"$SAMPLE" \
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
