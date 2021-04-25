#!/bin/bash

INPUT=$1
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
S3HELIX=s3://helix-all
S3UCSD=s3://ucsd-all

VERSION_INFO=$(bash $PIPELINEDIR/pipeline/show_version.sh)

[ ! -f $INPUT ] && { echo "Error: $INPUT file not found"; exit 99; }
sed 1d $INPUT | while IFS=',' read SAMPLE ORGANIZATION SEQ_RUN PRIMER_SET FQ MERGE_LANES VARIANTS QC LINEAGE TREE_BUILD READ_CAP TIMESTAMP
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

	if [[ "$MERGE_LANES" == true ]]; then
	  echo "Error: MERGE_LANES is currently not supported. Exiting."
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

	re='^[0-9]+$'
	if [[ ! $READ_CAP =~ ^($re|all)$ ]] ; then
	   echo "Error: READ_CAP must be an integer or 'all'"
	   exit 1
	fi

	echo "Sample: $SAMPLE"
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
	echo "Timestamp: $TIMESTAMP"


	qsub \
		-v SEQ_RUN="$SEQ_RUN" \
		-v SAMPLE=$SAMPLE \
		-v S3DOWNLOAD=$S3DOWNLOAD \
		-v PRIMER_SET=$PRIMER_SET \
		-v MERGE_LANES=$MERGE_LANES \
		-v FQ=$FQ \
		-v TIMESTAMP=$TIMESTAMP \
		-v READ_CAP=$READ_CAP \
		-v VERSION_INFO="$VERSION_INFO" \
		-N Covid19_"$SEQ_RUN"_"$SAMPLE" \
		-wd /shared/workspace/projects/covid/logs \
		-pe smp 2 \
		-S /bin/bash \
		$PIPELINEDIR/pipeline/sarscov2_consensus_pipeline.sh

done
