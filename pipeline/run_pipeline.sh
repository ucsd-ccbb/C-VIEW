#!/bin/bash

INPUT=$1 # Sample Sheet with header
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
S3HELIX=s3://helix-all
S3UCSD=s3://ucsd-all
S3TEST=s3://ucsd-rtl-test
ANACONDADIR=/shared/workspace/software/anaconda3/bin
source $ANACONDADIR/activate covid1.2

VERSION_INFO=$(bash $PIPELINEDIR/pipeline/show_version.sh)

[ ! -f $INPUT ] && { echo "Error: $INPUT file not found"; exit 99; }
INPUT_FIELDS="VARIANTS,QC,LINEAGE,TREE_BUILD,ORGANIZATION,SEQ_RUN,MERGE_LANES,PRIMER_SET,FQ,READ_CAP,SAMPLE,TIMESTAMP,ISTEST"
sed 1d $INPUT | while IFS=',' read VARIANTS QC LINEAGE TREE_BUILD ORGANIZATION SEQ_RUN MERGE_LANES PRIMER_SET FQ READ_CAP SAMPLE TIMESTAMP ISTEST

# for each row in the input file
do
  # set per-row variables
  INPUT_VALS="$VARIANTS,$QC,$LINEAGE,$TREE_BUILD,$ORGANIZATION,$SEQ_RUN,$MERGE_LANES,$PRIMER_SET,$FQ,$READ_CAP,$SAMPLE,$TIMESTAMP,$ISTEST"
  QSUBSAMPLEPARAMS=''
  QSUBLINEAGEPARAMS=''

	# echo the inputs to the screen
	echo $INPUT_FIELDS
	echo $INPUT_VALS

  # validate the inputs
  # --------------------
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

	if [[ ! "$ORGANIZATION" =~ ^(ucsd|helix)$ ]]; then
		echo "Error: Parameter ORGANIZATION must be one of ucsd or helix"
		exit 1
	fi

	if [[ ! "$MERGE_LANES" =~ ^(true|false)$ ]]; then
	  echo "Error: MERGE_LANES must be one of true or false"
	  exit 1
	fi

	if [[ ! "$PRIMER_SET" =~ ^(artic|swift_v2)$ ]]; then
		echo "Error: Parameter PRIMER_SET must be one of artic or swift_v2"
		exit 1
	fi

	if [[ ! "$FQ" =~ ^(se|pe)$ ]]; then
	  echo "Error: FQ must be one of se or pe"
	  exit 1
	fi

	re='^[0-9]+$'
	if [[ ! $READ_CAP =~ ^($re|all)$ ]] ; then
	   echo "Error: READ_CAP must be an integer or 'all'"
	   #exit 1
	fi

	if [[ ! "$ISTEST" =~ ^(true|false)$ ]]; then
	  echo "Error: ISTEST must be one of true or false"
	  exit 1
	fi

  # prevent currently-unsupported uses of SAMPLE
	if [[ "$SAMPLE" != NA ]]; then
	  if [[ "$VARIANTS" != true ]]; then
      echo "Error: SAMPLE is supported only if VARIANTS is true"
      exit 1
    else
      if [ "$QC" == true ] || [ "$LINEAGE" == true ] || [ "$TREE_BUILD" == true ]; then
        echo "Error: SAMPLE is supported only if non-VARIANTS steps are false"
        exit 1
      fi
    fi
  fi

  # prevent currently-unsupported uses of TIMESTAMP
	if [[ "$TIMESTAMP" != NA ]]; then
	  if [ "$VARIANTS" != true ] && [ "$QC" != true ]; then
      echo "Error: TIMESTAMP is supported only if VARIANTS and/or QC is true"
      exit 1
    else
      if [ "$LINEAGE" == true ] || [ "$TREE_BUILD" == true ]; then
        echo "Error: TIMESTAMP is supported only if LINEAGE and TREE_BUILD steps are false"
        exit 1
      fi
    fi

    if [ "$VARIANTS" == true ]; then
      if [ "$SAMPLE" == NA ]; then
        echo "Error: TIMESTAMP is supported for VARIANTS only if SAMPLE is not NA"
        exit 1
      fi
    fi
  else
    if [ "$VARIANTS" != true ] && [ "$QC" == true ]; then
      echo "Error: TIMESTAMP must be provided if QC is true and VARIANTS is false"
      exit 1
    fi
  fi

  # prevent currently-unsupported use of seq_run
  if [ "$SEQ_RUN" == all ]; then
    if [ "$VARIANTS" != false ] || [ "$QC" != false ]; then
      echo "Error: a specific SEQ_RUN must be provided for VARIANTS and/or QC step"
      exit 1
    fi
  fi
  # --------------------
  # end validating inputs

	# identify where to get data from
	if [[ "$ISTEST" == false ]]; then
	  if [[ "$ORGANIZATION" == ucsd ]]; then
      S3DOWNLOAD=$S3UCSD
    elif [[ "$ORGANIZATION" == helix ]]; then
      S3DOWNLOAD=$S3HELIX
    fi
  else
    S3DOWNLOAD=$S3TEST
  fi

	# set timestamp if not already specified
	if [[ $TIMESTAMP == NA ]]; then
	  TIMESTAMP=$(date +'%Y-%m-%d_%H-%M-%S')
	fi

	# if we are doing per-sample processing
	if [[ "$VARIANTS" == true ]]; then
	  FASTQS_PATH=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_fastq

	  # if we are processing just one sample
	  if [[ "$SAMPLE" != NA ]]; then
	    SAMPLE_LIST=($SAMPLE)
	  else
	    # set up to find all the relevant fastqs
	    DELIMITER=_R1_001.fastq.gz
      BOTH_FASTQS=$(aws s3 ls $FASTQS_PATH/ |  grep ".fastq.gz" | sort -k3 -n | awk '{print $NF}' | sort | uniq | grep -v Undetermined)
	    R1_FASTQS=$(echo "$BOTH_FASTQS" |  grep $DELIMITER )

      # If necessary, merge fastq files from multiple lanes
      if [[ "$MERGE_LANES" == true ]]; then
        # get the (non-unique) list of sample identifiers without lane/read info (but with sample number)
        INSPECT_DELIMITER=__
        SAMPLES_WO_LANES_LIST=()
        for R1_FASTQ in $(echo $R1_FASTQS); do
          SEQUENCING_INFO=$(echo $R1_FASTQ | awk -F $INSPECT_DELIMITER '{print $NF}')
          SAMPLE_NUM=$(echo $SEQUENCING_INFO | awk -F _ '{print $1}')
          A_SAMPLE=$(echo $R1_FASTQ | sed "s/$SEQUENCING_INFO/$SAMPLE_NUM/g")
          SAMPLES_WO_LANES_LIST+=($A_SAMPLE)
        done

        # reduce above list to only unique values and loop over them
        FINAL_R1_FASTQS=()
        for SAMPLE in $(printf '%s\n' "${SAMPLES_WO_LANES_LIST[@]}" | sort | uniq ); do
          LANES=$(echo "$BOTH_FASTQS" | grep "$SAMPLE" | awk -F $INSPECT_DELIMITER '{print $NF}'| awk -F '_L|_R' '{print $2}' | sort | uniq | grep 00)
          LANES_COMBINED=$(echo $LANES | sed 's/ //g')

          FINAL_R1_FASTQS+=("$SAMPLE"_"$LANES_COMBINED""$DELIMITER")
        done

        qsub -v SEQ_RUN=$SEQ_RUN \
           -v WORKSPACE=/scratch/$SEQ_RUN/$TIMESTAMP \
           -v S3DOWNLOAD=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_fastq \
           -wd /shared/workspace/projects/covid/logs \
           -N m_"$SEQ_RUN" \
           -pe smp 16 \
           -S /bin/bash \
           $PIPELINEDIR/pipeline/merge_lanes.sh

        R1_FASTQS=$(printf '%s\n' "${FINAL_R1_FASTQS[@]}")
        QSUBSAMPLEPARAMS=' -hold_jid m_'$SEQ_RUN''
      fi # end if we are merging lanes

	    SAMPLE_LIST=$(echo "$R1_FASTQS" | awk -F $DELIMITER '{print $1}' | sort | uniq)
	  fi # end if we are handling all samples, not just one

	  # sanity-check we have at least one sample to run
    if [[ "$SAMPLE_LIST" == "" ]]; then
      echo "Error: There are no samples to run in $FASTQS_PATH"
      exit 1
    fi

		# Run VARIANTS step on each sample
		for SAMPLE in $SAMPLE_LIST; do
			qsub $QSUBSAMPLEPARAMS \
				-v SEQ_RUN="$SEQ_RUN" \
				-v SAMPLE=$SAMPLE \
				-v S3DOWNLOAD=$S3DOWNLOAD \
				-v PRIMER_SET=$PRIMER_SET \
				-v MERGE_LANES=$MERGE_LANES \
				-v FQ=$FQ \
				-v TIMESTAMP=$TIMESTAMP \
				-v VERSION_INFO="$VERSION_INFO" \
				-v READ_CAP=$READ_CAP \
				-N v_"$SEQ_RUN"_"$TIMESTAMP"_"$SAMPLE" \
				-wd /shared/workspace/projects/covid/logs \
				-pe smp 2 \
				-S /bin/bash \
				$PIPELINEDIR/pipeline/sarscov2_consensus_pipeline.sh
		done
	fi  # end if we are calling variants

	if [[ "$QC" == true ]]; then
    qsub \
      -hold_jid 'v_'$SEQ_RUN'_'$TIMESTAMP'_*' \
      -v SEQ_RUN=$SEQ_RUN \
      -v S3DOWNLOAD=$S3DOWNLOAD \
      -v WORKSPACE=/scratch/$SEQ_RUN/$TIMESTAMP \
      -v FQ=$FQ \
      -v TIMESTAMP=$TIMESTAMP \
      -v ISTEST=$ISTEST \
      -v VERSION_INFO="$VERSION_INFO" \
      -N q_"$SEQ_RUN" \
      -wd /shared/workspace/projects/covid/logs \
      -pe smp 32 \
      -S /bin/bash \
        $PIPELINEDIR/qc/qc_summary.sh
	fi # end if we are running qc

  # lineage and tree output is stored to a different location than
  # the outputs of the earlier steps, so it needs a slightly different
  # identifier for its process
  PROCESSINGID="$TIMESTAMP"-"$SEQ_RUN"

	if [[ "$LINEAGE" == true ]]; then
	    qsub \
      -hold_jid 'q_'$SEQ_RUN'' \
      -v ORGANIZATION=$ORGANIZATION \
			-v PROCESSINGID=$PROCESSINGID \
			-v SEQ_RUN=$SEQ_RUN \
			-v ISTEST=$ISTEST \
			-v WORKSPACE=/scratch/phylogeny/$PROCESSINGID \
			-v VERSION_INFO="$VERSION_INFO" \
			-N l_"$PROCESSINGID" \
			-wd /shared/workspace/projects/covid/logs \
			-pe smp 16 \
			-S /bin/bash \
	    	$PIPELINEDIR/pipeline/lineages.sh

		QSUBLINEAGEPARAMS=' -hold_jid l_'$PROCESSINGID
	fi # end if we are calling lineages

  if [[ "$TREE_BUILD" == true ]]; then
    for DATASET in stringent loose_stringent passQC; do
      qsub \
        $QSUBLINEAGEPARAMS \
        -v ORGANIZATION=$ORGANIZATION \
        -v PROCESSINGID=$PROCESSINGID \
        -v DATASET=$DATASET \
        -v ISTEST=$ISTEST \
        -v WORKSPACE=/scratch/treebuilding/$PROCESSINGID/$DATASET \
        -v VERSION_INFO="$VERSION_INFO" \
        -N t_"$DATASET"_"$PROCESSINGID" \
        -wd /shared/workspace/projects/covid/logs \
        -pe smp 16 \
        -S /bin/bash \
          $PIPELINEDIR/pipeline/treebuild.sh
    done
  fi # end if we are building trees

  # if we did variant or qc processing, write a file of the
  # settings used into the results directory for this seq_run/timestamp
  if [ "$VARIANTS" == true ] || [ "$QC" == true ]; then
    SETTINGS_FNAME="$SEQ_RUN"-"$TIMESTAMP".csv
    echo $INPUT_FIELDS > $SETTINGS_FNAME
    echo $INPUT_VALS >> $SETTINGS_FNAME
    aws s3 cp $SETTINGS_FNAME $S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/
	  rm $SETTINGS_FNAME
  fi

  # ensure that, if we are processing additional input rows,
  # they will have a distinct timestamp
  sleep 1
done
