#!/bin/bash

INPUT=$1 # Sample Sheet with header - seq_run,s3download,s3upload,primers,reads
TIMESTAMP=$(date +'%Y-%m-%d_%H-%M-%S')
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
S3_TREEBUILD=s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/tree_building
QSUBSAMPLEPARAMS=''

[ ! -f $INPUT ] && { echo "Error: $INPUT file not found"; exit 99; }
sed 1d $INPUT | while IFS=',' read SEQ_RUN S3DOWNLOAD PRIMER_SET FQ MERGE_LANES TREE_BUILD
do
	echo "Seq_Run: $SEQ_RUN"
	echo "S3 Fastq path: $S3DOWNLOAD/"$SEQ_RUN"_fastq"
	echo "Primers: $PRIMER_SET"
	echo "Fastq Reads: $FQ"
	echo "Merge Lanes: $MERGE_LANES"
	echo "Run tree building: $TREE_BUILD" 

	# Append Results URL
	RESULTS="$TIMESTAMP"_"$FQ"

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

	if [[ ! "$TREE_BUILD" =~ ^(true|false)$ ]]; then
	  echo "Error: TREE_BUILD must be one of true or false"
	  exit 1
	fi

	# Merge fastq files from multiple lanes
	if [[ "$MERGE_LANES" == true ]]; then
		qsub -v SEQ_RUN=$SEQ_RUN \
			 -v WORKSPACE=/scratch/$SEQ_RUN/$TIMESTAMP \
			 -v S3DOWNLOAD=$S3DOWNLOAD/"$SEQ_RUN"_fastq \
			 -wd /shared/workspace/projects/covid/logs \
			 -N merge_fq_lanes_"$SEQ_RUN" \
			 -pe smp 16 \
			 -S /bin/bash \
			 $PIPELINEDIR/pipeline/merge_lanes.sh
		DELIMITER=_L00
		QSUBSAMPLEPARAMS=' -hold_jid merge_fq_lanes_'$SEQ_RUN''
	else
		DELIMITER=_R
	fi

	SAMPLE_LIST=$(aws s3 ls $S3DOWNLOAD/"$SEQ_RUN"_fastq/ | grep fastq.gz | awk '{print $NF}' | awk -F $DELIMITER '{print $1}' | sort | uniq | grep -v Undetermined)

	# Run pipeline on each sample
	for SAMPLE in $SAMPLE_LIST; do
		qsub $QSUBSAMPLEPARAMS \
			-v SEQ_RUN="$SEQ_RUN" \
			-v SAMPLE=$SAMPLE \
			-v S3DOWNLOAD=$S3DOWNLOAD \
			-v PRIMER_SET=$PRIMER_SET \
			-v MERGE_LANES=$MERGE_LANES \
			-v FQ=$FQ \
			-v TIMESTAMP=$TIMESTAMP \
			-N Covid19_"$SEQ_RUN"_"$SAMPLE" \
			-wd /shared/workspace/projects/covid/logs \
			-pe smp 1 \
			-S /bin/bash \
			$PIPELINEDIR/pipeline/sarscov2_consensus_pipeline.sh
	done


	# Perform QC and summary on seq_run when all samples have completed
	qsub \
		-hold_jid 'Covid19_'$SEQ_RUN'_*' \
		-v SEQ_RUN=$SEQ_RUN \
		-v S3DOWNLOAD=$S3DOWNLOAD \
		-v WORKSPACE=/scratch/$SEQ_RUN/$TIMESTAMP \
		-v FQ=$FQ \
		-v TIMESTAMP=$TIMESTAMP \
		-N QC_summary_"$SEQ_RUN" \
		-wd /shared/workspace/projects/covid/logs \
		-pe smp 32 \
		-S /bin/bash \
    	$PIPELINEDIR/qc/qc_summary.sh

    # Tree building
    if [[ "$TREE_BUILD" == true ]]; then
	    qsub \
			-hold_jid 'QC_summary_'$SEQ_RUN'' \
			-v S3DOWNLOAD=$S3_TREEBUILD \
			-v TIMESTAMP=$TIMESTAMP \
			-v WORKSPACE=/scratch/tree_build \
			-N tree_building \
			-wd /shared/workspace/projects/covid/logs \
			-pe smp 96 \
			-S /bin/bash \
	    	$PIPELINEDIR/pipeline/tree_building_merged.sh

    fi

	echo "seq_run,s3download,primers,reads" > "$SEQ_RUN"-"$TIMESTAMP".csv
	echo "$SEQ_RUN,$S3DOWNLOAD,$PRIMER_SET,$FQ" >> "$SEQ_RUN"-"$TIMESTAMP".csv
	aws s3 cp "$SEQ_RUN"-"$TIMESTAMP".csv $S3DOWNLOAD/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/
	rm "$SEQ_RUN"-"$TIMESTAMP".csv
done
