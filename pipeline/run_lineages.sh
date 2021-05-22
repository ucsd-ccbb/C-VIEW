#!/bin/bash

INPUT=$1 # Sample Sheet with header
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
S3HELIX=s3://helix-all
S3UCSD=s3://ucsd-all
TIMESTAMP=$(date +'%Y-%m-%d_%H-%M-%S')
VERSION_INFO=$(bash $PIPELINEDIR/pipeline/show_version.sh)
QSUBPARAMS=''

[ ! -f $INPUT ] && { echo "Error: $INPUT file not found"; exit 99; }
sed 1d $INPUT | while IFS=',' read ORGANIZATION LINEAGE TREE_BUILD ISTEST
do

	if [[ ! "$ORGANIZATION" =~ ^(ucsd|helix)$ ]]; then
		echo "Error: Parameter ORGANIZATION must be one of ucsd or helix"
		exit 1
	fi

	if [[ "$ORGANIZATION" == ucsd ]]; then
		S3DOWNLOAD=$S3UCSD
	elif [[ "$ORGANIZATION" == helix ]]; then
		S3DOWNLOAD=$S3HELIX
	fi

	if [[ ! "$LINEAGE" =~ ^(true|false)$ ]]; then
	  echo "Error: LINEAGE must be one of true or false"
	  exit 1
	fi

	if [[ ! "$TREE_BUILD" =~ ^(true|false)$ ]]; then
	  echo "Error: TREE_BUILD must be one of true or false"
	  exit 1
	fi

	echo "Organization: $ORGANIZATION"
	echo "Lineage with Pangolin: $LINEAGE"
	echo "Run tree building: $TREE_BUILD"
	echo "Is test run: $ISTEST"

	if [[ "$LINEAGE" == true ]]; then
	    qsub \
			-v ORGANIZATION=$ORGANIZATION \
			-v TREE_BUILD=$TREE_BUILD \
			-v TIMESTAMP=$TIMESTAMP \
			-v WORKSPACE=/scratch/phylogeny/$TIMESTAMP \
			-v ISTEST=$ISTEST \
			-v VERSION_INFO="$VERSION_INFO" \
			-N lineage \
			-wd /shared/workspace/projects/covid/logs \
			-pe smp 16 \
			-S /bin/bash \
	    	$PIPELINEDIR/pipeline/lineages.sh

		QSUBPARAMS=' -hold_jid lineage'
	fi

    if [[ "$TREE_BUILD" == true ]]; then
	    for DATASET in stringent loose_stringent passQC; do
		    qsub \
				$QSUBPARAMS \
				-v ORGANIZATION=$ORGANIZATION \
				-v TIMESTAMP=$TIMESTAMP \
				-v DATASET=$DATASET \
				-v ISTEST=$ISTEST \
				-v WORKSPACE=/scratch/treebuilding/$TIMESTAMP/$DATASET \
				-v VERSION_INFO="$VERSION_INFO" \
				-N treebuilding_"$DATASET" \
				-wd /shared/workspace/projects/covid/logs \
				-pe smp 16 \
				-S /bin/bash \
		    	$PIPELINEDIR/pipeline/treebuild.sh
		done
	fi

done
