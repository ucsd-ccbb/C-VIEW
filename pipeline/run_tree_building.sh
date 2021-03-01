#!/bin/bash

# Only run tree building

TIMESTAMP=$(date +'%Y-%m-%d_%H-%M-%S')
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
S3_TREEBUILD=s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/tree_building

qsub \
	-v S3DOWNLOAD=$S3_TREEBUILD \
	-v TIMESTAMP=$TIMESTAMP \
	-v WORKSPACE=/scratch/tree_build \
	-N tree_building \
	-wd /shared/workspace/projects/covid/logs \
	-pe smp 96 \
	-S /bin/bash \
	$PIPELINEDIR/pipeline/tree_building_merged.sh
