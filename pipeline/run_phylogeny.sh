#!/bin/bash

# Only run tree building

TIMESTAMP=$(date +'%Y-%m-%d_%H-%M-%S')
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
S3_TREEBUILD=s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/tree_building
TREE_BUILD=$1

if [[ ! "$TREE_BUILD" =~ ^(true|false)$ ]]; then
  echo "Error: TREE_BUILD must be one of true or false"
  exit 1
fi

qsub \
	-hold_jid 'QC_summary_'$SEQ_RUN'' \
	-v S3DOWNLOAD=$S3_TREEBUILD \
	-v TREE_BUILD=$TREE_BUILD \
	-v TIMESTAMP=$TIMESTAMP \
	-v WORKSPACE=/scratch/lineage/$TIMESTAMP \
	-N tree_building \
	-wd /shared/workspace/projects/covid/logs \
	-pe smp 96 \
	-S /bin/bash \
	$PIPELINEDIR/pipeline/phylogeny.sh
