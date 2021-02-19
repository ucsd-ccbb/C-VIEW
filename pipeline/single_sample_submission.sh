#!/bin/bash

S3DOWNLOAD=$1 # e.g. s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/
IS_ARTIC=$2 # true or false - is this a pipeline test-running on ARTIC samples?
FQ=$3 # se or pe - single end or paired end reads
MERGED=$4 # merged or unmerged - merge lanes
SAMPLE=$5
RESULTSDATE=$6 # Manually input if rerunning a single sample
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline

if [[ ! "$IS_ARTIC" =~ ^(true|false)$ ]]; then
	echo "Parameter 2 - IS_ARTIC must be one of true or false"
	exit 1
fi

if [[ "$IS_ARTIC" == true ]]; then
	FQ=NA
	MERGED=NA
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

	fi
fi


qsub -v SAMPLE=$SAMPLE \
	-v S3DOWNLOAD=$S3DOWNLOAD \
	-v IS_ARTIC=$IS_ARTIC \
	-v FQ=$FQ \
	-v MERGED=$MERGED \
	-v RESULTSDATE=$RESULTSDATE \
	-N Covid19_"$SAMPLE" \
	-wd /shared/workspace/projects/covid/logs \
	-pe smp 1 \
	-S /bin/bash \
	$PIPELINEDIR/pipeline/sarscov2_consensus_pipeline.sh
