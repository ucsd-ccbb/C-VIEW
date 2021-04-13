#!/bin/bash

INPUT=$1 # Sample Sheet with header
#TIMESTAMP=$(date +'%Y-%m-%d_%H-%M-%S')
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
S3HELIX=s3://helix-all
S3UCSD=s3://ucsd-all
TIMESTAMP=$(date +'%Y-%m-%d_%H-%M-%S')

[ ! -f $INPUT ] && { echo "Error: $INPUT file not found"; exit 99; }
sed 1d $INPUT | while IFS=',' read ORGANIZATION SEQ_RUN PRIMER_SET FQ MERGE_LANES VARIANTS QC LINEAGE TREE_BUILD ISTEST
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

	if [[ ! "$PRIMER_SET" =~ ^(artic|swift_v2)$ ]]; then
		echo "Error: Parameter PRIMER_SET must be one of artic or swift_v2"
		exit 1
	fi

	if [[ ! "$FQ" =~ ^(se|pe)$ ]]; then
	  echo "Error: FQ must be one of se or pe"
	  exit 1
	fi

	if [[ ! "$MERGE_LANES" =~ ^(true|false)$ ]]; then
	  echo "Error: MERGE_LANES must be one of true or false"
	  exit 1
	fi

	if [[ ! "$VARIANTS" =~ ^(true|false)$ ]]; then
	  echo "Error: VARIANTS must be one of true or false"
	  exit 1
	fi

	if [[ ! "$QC" =~ ^(true|false)$ ]]; then
	  echo "Error: QC must be one of true or false"
	  exit 1
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
	echo "Seq_Run: $SEQ_RUN"
	echo "S3 Fastq path: $S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_fastq"
	echo "Primers: $PRIMER_SET"
	echo "Fastq Reads: $FQ"
	echo "Merge Lanes: $MERGE_LANES"
	echo "Call Variants: $VARIANTS"
	echo "Run QC: $QC"
	echo "Lineage with Pangolin: $LINEAGE"
	echo "Run tree building: $TREE_BUILD"
	echo "Is test run: $ISTEST"

    qsub \
		-v ORGANIZATION=$ORGANIZATION \
		-v TREE_BUILD=$TREE_BUILD \
		-v TIMESTAMP=$TIMESTAMP \
		-v WORKSPACE=/scratch/phylogeny/$TIMESTAMP \
		-v ISTEST=$ISTEST \
		-N tree_building \
		-wd /shared/workspace/projects/covid/logs \
		-pe smp 96 \
		-S /bin/bash \
    	$PIPELINEDIR/pipeline/phylogeny.sh

done
