#!/bin/bash

INPUT=$1 # Sample Sheet with header
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
S3HELIX=s3://helix-all
S3UCSD=s3://ucsd-all
S3TEST=s3://ucsd-rtl-test
ANACONDADIR=/shared/workspace/software/anaconda3/bin
source $ANACONDADIR/activate covid1.2

VERSION_INFO=$(bash $PIPELINEDIR/pipeline/show_version.sh)

# TODO: start requiring that first line exactly match expected headers so
# can ensure we don't mis-interpret old file w diff number/order of fields
[ ! -f $INPUT ] && { echo "Error: $INPUT file not found"; exit 99; }
INPUT_FIELDS="FUNCTION,ORGANIZATION,SEQ_RUN,MERGE_LANES,PRIMER_SET,FQ,READ_CAP,SAMPLE,TIMESTAMP,ISTEST"
sed 1d $INPUT | while IFS=',' read FUNCTION ORGANIZATION SEQ_RUN MERGE_LANES PRIMER_SET FQ READ_CAP SAMPLE TIMESTAMP ISTEST

# for each row in the input file
do
  # set per-row variables
  INPUT_VALS="$FUNCTION,$ORGANIZATION,$SEQ_RUN,$MERGE_LANES,$PRIMER_SET,$FQ,$READ_CAP,$SAMPLE,$TIMESTAMP,$ISTEST"
  QSUBSAMPLEPARAMS=''
  QSUBLINEAGEPARAMS=''

	# echo the inputs to the screen
	echo " "
	echo $INPUT_FIELDS
	echo $INPUT_VALS

  # initialize variables that will be reset to run-specific
  # values based on the inputs
  VARIANTS=false
  QC=false
  LINEAGE=false
  TREE_BUILD=false
  PHYLO_SEQ_RUN=$SEQ_RUN

	declare -A FIELD_IGNORED
	FIELD_IGNORED[SEQ_RUN]=$SEQ_RUN
	FIELD_IGNORED[MERGE_LANES]=$MERGE_LANES
	FIELD_IGNORED[PRIMER_SET]=$PRIMER_SET
	FIELD_IGNORED[FQ]=$FQ
	FIELD_IGNORED[READ_CAP]=$READ_CAP
	FIELD_IGNORED[SAMPLE]=$SAMPLE
	FIELD_IGNORED[TIMESTAMP]=$TIMESTAMP

  # validate the inputs
  # --------------------

  if [ "$FUNCTION" != pipeline ] && [ "$FUNCTION" != cumulative_pipeline ] && [ "$FUNCTION" != variants ] &&  [ "$FUNCTION" != sample ] &&  [ "$FUNCTION" != qc ] &&  [ "$FUNCTION" != lineages ] &&  [ "$FUNCTION" != phylogeny ] &&  [ "$FUNCTION" != cumulative_lineages ] && [ "$FUNCTION" != cumulative_phylogeny ]; then
		echo "Error: Parameter FUNCTION must be one of pipeline, cumulative_pipeline, variants, sample, qc, lineages, phylogeny, cummulative_lineages, or cumulative_phylogeny"
		exit 1
	fi

  if [[ ! "$ORGANIZATION" =~ ^(ucsd|helix)$ ]]; then
		echo "Error: Parameter ORGANIZATION must be one of ucsd or helix"
		exit 1
	fi

	if [[ ! "$ISTEST" =~ ^(true|false)$ ]]; then
	  echo "Error: ISTEST must be one of true or false"
	  exit 1
	fi

  if [ "$FUNCTION" == pipeline ] || [ "$FUNCTION" == cumulative_pipeline ] || [ "$FUNCTION" == variants ] || [ "$FUNCTION" == sample ]; then
    VARIANTS=true

    if [[ ! "$MERGE_LANES" =~ ^(true|false)$ ]]; then
      echo "Error: MERGE_LANES must be one of true or false"
      exit 1
    fi
	  unset FIELD_IGNORED[MERGE_LANES]

    if [[ ! "$PRIMER_SET" =~ ^(artic|swift_v2|mini_artic)$ ]]; then
      echo "Error: Parameter PRIMER_SET must be one of artic, swift_v2, or mini_artic"
      exit 1
    fi

    if [[ "$PRIMER_SET" == artic ]]; then
      PRIMER_BED_FNAME="nCoV-2019.primer.bed"
    elif [[ "$PRIMER_SET" == swift_v2 ]]; then
      PRIMER_BED_FNAME="sarscov2_v2_primers.bed"
    elif [[ "$PRIMER_SET" == mini_artic ]]; then
      PRIMER_BED_FNAME="SARS2_short_primers_V3_no_adapter_sequences_mini_artic.bed"
    fi

    unset FIELD_IGNORED[PRIMER_SET]

    re='^[0-9]+$'
    if [[ ! $READ_CAP =~ ^($re|all)$ ]] ; then
       echo "Error: READ_CAP must be an integer or 'all'"
       #exit 1
    fi

    # TODO: add check if this is NOT using bams (if it is, this field stays ignored)
    unset FIELD_IGNORED[READ_CAP]


    if [[ $READ_CAP == all ]] ; then
      # C/C++ Unsigned long max = 4294967295
      READ_CAP=4294967295
    fi
  fi

  if [ "$FUNCTION" == pipeline ] || [ "$FUNCTION" == cumulative_pipeline ] || [ "$FUNCTION" == variants ] || [ "$FUNCTION" == sample ] || [ "$FUNCTION" == qc ]; then
    if [[ ! "$FQ" =~ ^(se|pe)$ ]]; then
      echo "Error: FQ must be one of se or pe"
      exit 1
    fi
	  unset FIELD_IGNORED[FQ]
  fi

  if [ "$FUNCTION" == pipeline ] || [ "$FUNCTION" == cumulative_pipeline ] || [ "$FUNCTION" == qc ]; then
    QC=true
  fi

  if [ "$FUNCTION" == pipeline ] || [ "$FUNCTION" == cumulative_pipeline ] || [ "$FUNCTION" == lineages ] || [ "$FUNCTION" == phylogeny ] || [ "$FUNCTION" == cumulative_lineages ] || [ "$FUNCTION" == cumulative_phylogeny ]; then
    LINEAGE=true
  fi

  if [ "$FUNCTION" == pipeline ] || [ "$FUNCTION" == cumulative_pipeline ] || [ "$FUNCTION" == phylogeny ] || [ "$FUNCTION" == cumulative_phylogeny ]; then
    TREE_BUILD=true
  fi

  if [ "$FUNCTION" == cumulative_pipeline ] || [ "$FUNCTION" == cumulative_lineages ] || [ "$FUNCTION" == cumulative_phylogeny ] ; then
    PHYLO_SEQ_RUN=all
  fi

  if [ "$FUNCTION" == sample ] || [ "$FUNCTION" == qc ]; then
    unset FIELD_IGNORED[TIMESTAMP]
  else
	  TIMESTAMP=$(date +'%Y-%m-%d_%H-%M-%S')
  fi

  if [ "$FUNCTION" == sample ]; then
    unset FIELD_IGNORED[SAMPLE]
  else
    SAMPLE=NA
  fi

  if [ "$FUNCTION" != cumulative_lineages ] && [ "$FUNCTION" != cumulative_phylogeny ] ; then
    unset FIELD_IGNORED[SEQ_RUN]
  else
    SEQ_RUN=all
  fi

  if [ "${#FIELD_IGNORED[@]}" -gt 0 ]; then
    echo " "
    echo "For input FUNCTION '$FUNCTION', the following inputs will be *ignored*: "
    for i in "${!FIELD_IGNORED[@]}"; do
      echo "$i=${FIELD_IGNORED[$i]}"
    done
  fi
  echo " "
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

	# if we are doing per-sample processing
	if [[ "$VARIANTS" == true ]]; then
	  RAW_FILE_TYPE="_bam"
	  # TODO add check for if we are doing genexus bams, in which case this
	  # ends with "_bam" instead of "_fastq"
	  FASTQS_PATH=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN""$RAW_FILE_TYPE"

	  # if we are processing just one sample
	  if [[ "$SAMPLE" != NA ]]; then
	    SAMPLE_LIST=($SAMPLE)
	  else
	    # TODO: For now, for POC, just hardcode the sample names :-|
	    # If decide to do this for real, will need to decide whether to
	    # require bams come in with names that are already correct or if
	    # we try to intuit them like we do with fastqs ...

#	    # set up to find all the relevant fastqs
#	    DELIMITER=_R1_001.fastq.gz
#      BOTH_FASTQS=$(aws s3 ls $FASTQS_PATH/ |  grep ".fastq.gz" | sort -k3 -n | awk '{print $NF}' | sort | uniq | grep -v Undetermined)
#	    R1_FASTQS=$(echo "$BOTH_FASTQS" |  grep $DELIMITER )
#
#      # If necessary, merge fastq files from multiple lanes
#      if [[ "$MERGE_LANES" == true ]]; then
#        # get the (non-unique) list of sample identifiers without lane/read info (but with sample number)
#        INSPECT_DELIMITER=__
#        SAMPLES_WO_LANES_LIST=()
#        for R1_FASTQ in $(echo $R1_FASTQS); do
#          SEQUENCING_INFO=$(echo $R1_FASTQ | awk -F $INSPECT_DELIMITER '{print $NF}')
#          SAMPLE_NUM=$(echo $SEQUENCING_INFO | awk -F _ '{print $1}')
#          A_SAMPLE=$(echo $R1_FASTQ | sed "s/$SEQUENCING_INFO/$SAMPLE_NUM/g")
#          SAMPLES_WO_LANES_LIST+=($A_SAMPLE)
#        done
#
#        # reduce above list to only unique values and loop over them
#        FINAL_R1_FASTQS=()
#        for SAMPLE in $(printf '%s\n' "${SAMPLES_WO_LANES_LIST[@]}" | sort | uniq ); do
#          LANES=$(echo "$BOTH_FASTQS" | grep "$SAMPLE" | awk -F $INSPECT_DELIMITER '{print $NF}'| awk -F '_L|_R' '{print $2}' | sort | uniq | grep 00)
#          LANES_COMBINED=$(echo $LANES | sed 's/ //g')
#
#          FINAL_R1_FASTQS+=("$SAMPLE"_"$LANES_COMBINED""$DELIMITER")
#        done
#
#        qsub -v SEQ_RUN=$SEQ_RUN \
#           -v WORKSPACE=/scratch/$SEQ_RUN/$TIMESTAMP \
#           -v S3DOWNLOAD=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_fastq \
#           -wd /shared/workspace/projects/covid/logs \
#           -N m_"$SEQ_RUN" \
#           -pe smp 16 \
#           -S /bin/bash \
#           $PIPELINEDIR/pipeline/merge_lanes.sh
#
#        R1_FASTQS=$(printf '%s\n' "${FINAL_R1_FASTQS[@]}")
#        QSUBSAMPLEPARAMS=' -hold_jid m_'$SEQ_RUN''
#      fi # end if we are merging lanes

      # TODO: fill this in with real POC sample names (no file extensions)
      SAMPLE_LIST=(AssayDev_EXC_MW4_208844__E0002080__L07__SARS-CoV-2InsightResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_EXC_MW4_350739__E0002080__M01__SARS-CoV-2InsightResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_EXC_MW4_478321__E0002080__I06__SARS-CoV-2InsightResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_EXC_MW4_478505__E0002080__H19__SARS-CoV-2InsightResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_EXC_MW4_616272__E0002080__I13__SARS-CoV-2InsightResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_EXC_MW4_631160__E0002080__C18__SARS-CoV-2InsightResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_EXC_MW5_369577__E0002080__F04__SARS-CoV-2InsightResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_EXC_MW5_461830__E0002080__L14__SARS-CoV-2InsightResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_PLJAN0922__1__000__SARS-CoV-2InsightWastewaterResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_PLJAN0922__2__000__SARS-CoV-2InsightWastewaterResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_PLJAN1022__1__000__SARS-CoV-2InsightWastewaterResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_PLJAN1022__2__000__SARS-CoV-2InsightWastewaterResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_PLJAN322__1__000__SARS-CoV-2InsightWastewaterResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_PLJAN322__2__000__SARS-CoV-2InsightWastewaterResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_PLJAN522__1__000__SARS-CoV-2InsightWastewaterResearchAssay_Lineage_7XUOM_7XUOM__001 AssayDev_PLJAN522__2__000__SARS-CoV-2InsightWastewaterResearchAssay_Lineage_7XUOM_7XUOM__001)
	    # SAMPLE_LIST=$(echo "$R1_FASTQS" | awk -F $DELIMITER '{print $1}' | sort | uniq)
	  fi # end if we are handling all samples, not just one

	  # sanity-check we have at least one sample to run
    if [[ "$SAMPLE_LIST" == "" ]]; then
      echo "Error: There are no samples to run in $FASTQS_PATH"
      exit 1
    fi

    # TODO: here genexus_bam is hardcoded to true; change if
		# Run VARIANTS step on each sample
		for SAMPLE in $SAMPLE_LIST; do
			qsub $QSUBSAMPLEPARAMS \
				-v SEQ_RUN="$SEQ_RUN" \
				-v SAMPLE=$SAMPLE \
				-v GENEXUS_BAM=true \
				-v S3DOWNLOAD=$S3DOWNLOAD \
        -v PRIMER_BED_FNAME=$PRIMER_BED_FNAME \
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
  PROCESSINGID="$TIMESTAMP"-"$PHYLO_SEQ_RUN"

	if [[ "$LINEAGE" == true ]]; then
	    qsub \
      -hold_jid 'q_'$SEQ_RUN'' \
      -v ORGANIZATION=$ORGANIZATION \
			-v PROCESSINGID=$PROCESSINGID \
			-v SEQ_RUN=$PHYLO_SEQ_RUN \
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
    # TODO: this is a temporary fix to build only stringent trees
    #  In the future, we will want to make this more tunable
    for DATASET in stringent; do # loose_stringent passQC; do
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
    SETTINGS_FNAME="$SEQ_RUN"-"$TIMESTAMP"-"$FUNCTION"-$(date +'%Y-%m-%d_%H-%M-%S').csv
    echo $INPUT_FIELDS > $SETTINGS_FNAME
    echo $INPUT_VALS >> $SETTINGS_FNAME
    aws s3 cp $SETTINGS_FNAME $S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/
	  rm $SETTINGS_FNAME
  fi

  # ensure that, if we are processing additional input rows,
  # they will have a distinct timestamp
  sleep 1
done
