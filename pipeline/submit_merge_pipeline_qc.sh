#!/bin/bash

INPUT=$1 # Sample Sheet with header - seq_run,s3download,s3upload,primers,reads
TIMESTAMP=$(date +'%Y-%m-%d_%H-%M-%S')
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
QSUBSAMPLEPARAMS=''

[ ! -f $INPUT ] && { echo "Error: $INPUT file not found"; exit 99; }
sed 1d $INPUT | while IFS=',' read SEQ_RUN S3DOWNLOAD S3UPLOAD PRIMER_SET FQ MERGE_LANES
do
	echo "Seq_Run: $SEQ_RUN"
	echo "S3 Fastq path: $S3DOWNLOAD"
	echo "S3 Results Path: $S3UPLOAD"
	echo "Primers: $PRIMER_SET"
	echo "Fastq Reads: $FQ"
	echo "Merge Lanes: $MERGE_LANES"

	WORKSPACE=/scratch/$SEQ_RUN/$TIMESTAMP
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

	# Merge fastq files from multiple lanes
	if [[ "$MERGE_LANES" == true ]]; then
		qsub -v SEQ_RUN=$SEQ_RUN \
			 -v WORKSPACE=$WORKSPACE \
			 -v S3DOWNLOAD=$S3DOWNLOAD \
			 -wd /shared/workspace/projects/covid/logs \
			 -N merge_lanes_"$SEQ_RUN" \
			 -pe smp 16 \
			 $PIPELINEDIR/pipeline/merge_lanes.sh
		S3DOWNLOAD=$S3DOWNLOAD/merged_lanes
		S3UPLOAD=$S3UPLOAD/merged_lanes
		QSUBSAMPLEPARAMS=' -hold_jid merge_lanes_'$SEQ_RUN''
	fi

	SAMPLE_LIST=$(aws s3 ls $S3DOWNLOAD/ | grep fastq | awk '{print $NF}' | awk -F '_R' '{print $1}' | sort | uniq | grep -v Undetermined)

	# Run pipeline on each sample
	for SAMPLE in $SAMPLE_LIST; do
		qsub $QSUBSAMPLEPARAMS \
			-v SAMPLE=$SAMPLE \
			-v S3DOWNLOAD=$S3DOWNLOAD \
			-v S3UPLOAD=$S3UPLOAD/$RESULTS \
			-v PRIMER_SET=$PRIMER_SET \
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
		-v S3DOWNLOAD=$S3UPLOAD/$RESULTS \
		-v S3UPLOAD=$S3UPLOAD/$RESULTS \
		-v SEQ_RUN=$SEQ_RUN \
		-v WORKSPACE=$WORKSPACE \
		-N QC_summary_"$SEQ_RUN" \
		-wd /shared/workspace/projects/covid/logs \
		-pe smp 32 \
		-S /bin/bash \
    	$PIPELINEDIR/qc/qc_summary.sh

	echo "seq_run,s3download,s3upload,primers,reads" > "$SEQ_RUN"-"$TIMESTAMP".csv
	echo "$SEQ_RUN,$S3DOWNLOAD,$S3UPLOAD,$PRIMER_SET,$FQ" >> "$SEQ_RUN"-"$TIMESTAMP".csv
	aws s3 cp "$SEQ_RUN"-"$TIMESTAMP".csv $S3UPLOAD/"$TIMESTAMP"_"$FQ"/
	rm "$SEQ_RUN"-"$TIMESTAMP".csv
done
