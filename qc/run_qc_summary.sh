#!/bin/bash

INPUT=$1 # Sample Sheet with header
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
S3HELIX=s3://helix-all
S3UCSD=s3://ucsd-all

[ ! -f $INPUT ] && { echo "Error: $INPUT file not found"; exit 99; }
sed 1d $INPUT | while IFS=',' read ORGANIZATION SEQ_RUN PRIMER_SET FQ MERGE_LANES VARIANTS QC LINEAGE TREE_BUILD TIMESTAMP ISTEST
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

	if [[ "$TIMESTAMP" == "" ]]; then
	  echo "Error: Must input run timestamp"
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
	echo "Previous run timestamp: $TIMESTAMP"

	RUN_EXISTS=$(aws s3 ls $S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/"$SEQ_RUN"_samples/)
	if [[ "$RUN_EXISTS" == "" ]]; then
		echo "Run $TIMESTAMP does not exist"
		exit 1
	fi

	qsub \
		-v SEQ_RUN=$SEQ_RUN \
		-v S3DOWNLOAD=$S3DOWNLOAD \
		-v WORKSPACE=/scratch/$SEQ_RUN/$TIMESTAMP \
		-v FQ=$FQ \
		-v TIMESTAMP=$TIMESTAMP \
		-v ISTEST=$ISTEST \
		-N QC_summary_"$SEQ_RUN" \
		-wd /shared/workspace/projects/covid/logs \
		-pe smp 32 \
		-S /bin/bash \
    	$PIPELINEDIR/qc/qc_summary.sh


	echo "organization,seq_run,primers,reads,merge,variants,qc,lineage,tree_build,is_test" > "$SEQ_RUN"-"$TIMESTAMP".csv
	echo "$ORGANIZATION,$SEQ_RUN,$PRIMER_SET,$FQ,$MERGE_LANES,$VARIANTS,$QC,$LINEAGE,$TREE_BUILD,$ISTEST" >> "$SEQ_RUN"-"$TIMESTAMP".csv
	aws s3 cp "$SEQ_RUN"-"$TIMESTAMP".csv $S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/
	rm "$SEQ_RUN"-"$TIMESTAMP".csv
done
