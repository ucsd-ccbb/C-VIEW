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
INPUT_FIELDS="FUNCTION,ORGANIZATION,SEQ_RUN,MERGE_LANES,PRIMER_SET,FQ,READ_CAP,SAMPLE,TIMESTAMP,ISTEST"
sed 1d $INPUT | while IFS=',' read FUNCTION ORGANIZATION SEQ_RUN MERGE_LANES PRIMER_SET FQ READ_CAP SAMPLE TIMESTAMP ISTEST

# for each row in the input file
do
  # set per-row variables
  INPUT_VALS="$FUNCTION,$ORGANIZATION,$SEQ_RUN,$MERGE_LANES,$PRIMER_SET,$FQ,$READ_CAP,$SAMPLE,$TIMESTAMP,$ISTEST"
  SBATCHSAMPLEPARAMS=''
  SBATCHLINEAGEPARAMS=''

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

  if [ "$FUNCTION" != pipeline ] && [ "$FUNCTION" != cumulative_pipeline ] && [ "$FUNCTION" != variants_qc ] && [ "$FUNCTION" != variants ] &&  [ "$FUNCTION" != sample ] &&  [ "$FUNCTION" != qc ] &&  [ "$FUNCTION" != lineages ] &&  [ "$FUNCTION" != phylogeny ] &&  [ "$FUNCTION" != cumulative_lineages ] && [ "$FUNCTION" != cumulative_phylogeny ]; then
		echo "Error: Parameter FUNCTION must be one of pipeline, cumulative_pipeline, variants_qc, variants, sample, qc, lineages, phylogeny, cummulative_lineages, or cumulative_phylogeny"
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

	if [ "$FUNCTION" == pipeline ] || [ "$FUNCTION" == cumulative_pipeline ] || [ "$FUNCTION" == variants_qc ] || [ "$FUNCTION" == variants ] || [ "$FUNCTION" == sample ] || [ "$FUNCTION" == qc ]; then
    if [[ ! "$FQ" =~ ^(se|pe|bam)$ ]]; then
      echo "Error: FQ must be one of se, pe, or bam"
      exit 1
    fi
	  unset FIELD_IGNORED[FQ]
  fi

  if [ "$FUNCTION" == pipeline ] || [ "$FUNCTION" == cumulative_pipeline ] || [ "$FUNCTION" == variants_qc ] || [ "$FUNCTION" == variants ] || [ "$FUNCTION" == sample ]; then
    VARIANTS=true

    if [[ "$FQ" =~ ^(se|pe)$ ]]; then
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
         exit 1
      fi
      unset FIELD_IGNORED[READ_CAP]

      if [[ $READ_CAP == all ]] ; then
        # C/C++ Unsigned long max = 4294967295
        READ_CAP=4294967295
      fi

      INPUT_TYPE="fastq"
      INPUT_SUFFIX=.fastq.gz
      INPUT_DELIMITER=_R1_001.fastq.gz
    else
      MERGE_LANES=false
      INPUT_TYPE="bam"
      INPUT_SUFFIX=.trimmed.sorted.unfiltered.bam
      INPUT_DELIMITER=$INPUT_SUFFIX
    fi
  fi

  if [ "$FUNCTION" == pipeline ] || [ "$FUNCTION" == cumulative_pipeline ] || [ "$FUNCTION" == variants_qc ] || [ "$FUNCTION" == qc ]; then
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
    echo "For provided settings, the following inputs will be *ignored*: "
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
	  INPUTS_PATH=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_"$INPUT_TYPE"

	  # if we are processing just one sample
	  if [[ "$SAMPLE" != NA ]]; then
	    SAMPLE_LIST=($SAMPLE)
	  else
	    # set up to find all the relevant input files
      INPUT_BASE_NAMES=$(aws s3 ls $INPUTS_PATH/ |  grep $INPUT_SUFFIX | sort -k3 -n | awk '{print $NF}' | sort | uniq | grep -v Undetermined)

      # if dealing with fastq files
      if [[ "$FQ" =~ (se|pe)$ ]]; then
        INPUT_NAMES=$(echo "$INPUT_BASE_NAMES" |  grep $INPUT_DELIMITER )

        # If necessary, merge fastq files from multiple lanes (functionality not available for bams)
        if [[ "$MERGE_LANES" == true ]]; then
          # get the (non-unique) list of sample identifiers without lane/read info (but with sample number)
          INSPECT_DELIMITER=__
          SAMPLES_WO_LANES_LIST=()
          for R1_FASTQ in $(echo $INPUT_NAMES); do
            SEQUENCING_INFO=$(echo $R1_FASTQ | awk -F $INSPECT_DELIMITER '{print $NF}')
            SAMPLE_NUM=$(echo $SEQUENCING_INFO | awk -F _ '{print $1}')
            A_SAMPLE=$(echo $R1_FASTQ | sed "s/$SEQUENCING_INFO/$SAMPLE_NUM/g")
            SAMPLES_WO_LANES_LIST+=($A_SAMPLE)
          done

          # reduce above list to only unique values and loop over them
          FINAL_INPUT_NAMES=()
          for SAMPLE in $(printf '%s\n' "${SAMPLES_WO_LANES_LIST[@]}" | sort | uniq ); do
            LANES=$(echo "$INPUT_BASE_NAMES" | grep "$SAMPLE" | awk -F $INSPECT_DELIMITER '{print $NF}'| awk -F '_L|_R' '{print $2}' | sort | uniq | grep 00)
            LANES_COMBINED=$(echo $LANES | sed 's/ //g')

            FINAL_INPUT_NAMES+=("$SAMPLE"_"$LANES_COMBINED""$INPUT_DELIMITER")
          done

          M_SLURM_JOB_ID=$(sbatch --export=SEQ_RUN=$SEQ_RUN,WORKSPACE=/scratch/$SEQ_RUN/$TIMESTAMP,S3DOWNLOAD=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_fastq \
            -D /shared/workspace/projects/covid/logs \
            -J m_"$SEQ_RUN" \
            -c 16 \
            $PIPELINEDIR/pipeline/merge_lanes.sh)

          INPUT_NAMES=$(printf '%s\n' "${FINAL_INPUT_NAMES[@]}")
          SBATCHSAMPLEPARAMS=' --dependency=afterok:${M_SLURM_JOB_ID##* }'
        fi # end if we are merging lanes
      else
        INPUT_NAMES=$INPUT_BASE_NAMES
      fi # end if input is/isn't fastq

	    SAMPLE_LIST=$(echo "$INPUT_NAMES" | awk -F $INPUT_DELIMITER '{print $1}' | sort | uniq)
	  fi # end if we are handling all samples, not just one

	  # sanity-check we have at least one sample to run
    if [[ "$SAMPLE_LIST" == "" ]]; then
      echo "Error: There are no samples to run in $INPUTS_PATH"
      exit 1
    fi

		# Run VARIANTS step on each sample
		V_SLURM_JOB_IDS=""
		for SAMPLE in $SAMPLE_LIST; do
      V_SLURM_JOB_IDS=$V_SLURM_JOB_IDS:$(sbatch $SBATCHSAMPLEPARAMS \
        --export=$(echo "SEQ_RUN=$SEQ_RUN,\
                  SAMPLE=$SAMPLE,\
                  S3DOWNLOAD=$S3DOWNLOAD,\
                  PRIMER_BED_FNAME=$PRIMER_BED_FNAME,\
                  MERGE_LANES=$MERGE_LANES,\
                  FQ=$FQ,\
                  TIMESTAMP=$TIMESTAMP,\
                  VERSION_INFO="$VERSION_INFO",\
                  READ_CAP=$READ_CAP \
                  INPUT_TYPE=$INPUT_TYPE \
                  INPUT_SUFFIX=$INPUT_SUFFIX" | sed 's/ //g') \
        -J v_"$SEQ_RUN"_"$TIMESTAMP"_"$SAMPLE" \
        -D /shared/workspace/projects/covid/logs \
        -c 2 \
        $PIPELINEDIR/pipeline/sarscov2_consensus_pipeline.sh)
		done
	fi  # end if we are calling variants

  V_SLURM_JOB_IDS=$(echo $V_SLURM_JOB_IDS | sed 's/Submitted batch job //g')

  Q_SLURM_JOB_ID=""
	if [[ "$QC" == true ]]; then
	        Q_SLURM_JOB_ID=$Q_SLURM_JOB_ID:$(sbatch \
        --dependency=afterok$V_SLURM_JOB_IDS \
        --export=$(echo "SEQ_RUN=$SEQ_RUN,\
                 S3DOWNLOAD=$S3DOWNLOAD,\
                 WORKSPACE=/scratch/$SEQ_RUN/$TIMESTAMP,\
                 FQ=$FQ,\
                 TIMESTAMP=$TIMESTAMP,\
                 VERSION_INFO="$VERSION_INFO",\
                 ISTEST=$ISTEST" | sed 's/ //g') \
        -J q_$SEQ_RUN \
        -D /shared/workspace/projects/covid/logs \
        -c 32 \
        $PIPELINEDIR/qc/qc_summary.sh)
	fi # end if we are running qc
	Q_SLURM_JOB_ID=${Q_SLURM_JOB_ID##* }

  # lineage and tree output is stored to a different location than
  # the outputs of the earlier steps, so it needs a slightly different
  # identifier for its process
  PROCESSINGID="$TIMESTAMP"-"$PHYLO_SEQ_RUN"

	if [[ "$LINEAGE" == true ]]; then
	    L_SLURM_JOB_ID=$(sbatch \
      --dependency=afterok:$Q_SLURM_JOB_ID \
      --export=$(echo "ORGANIZATION=$ORGANIZATION,\
                      PROCESSINGID=$PROCESSINGID,\
                      SEQ_RUN=$PHYLO_SEQ_RUN,\
                      ISTEST=$ISTEST,\
                      VERSION_INFO="$VERSION_INFO", \
                      WORKSPACE=/scratch/phylogeny/$PROCESSINGID" | sed 's/ //g') \
			-J l_"$PROCESSINGID" \
			-D /shared/workspace/projects/covid/logs \
			-c 16 \
	    $PIPELINEDIR/pipeline/lineages.sh)

		SBATCHLINEAGEPARAMS=" --dependency=afterok:${L_SLURM_JOB_ID##* }"
	fi # end if we are calling lineages

  if [[ "$TREE_BUILD" == true ]]; then
    # TODO: this is a temporary fix to build only stringent trees
    #  In the future, we will want to make this more tunable
    for DATASET in stringent; do # loose_stringent passQC; do
      sbatch \
        $SBATCHLINEAGEPARAMS \
        --export=$(echo "ORGANIZATION=$ORGANIZATION,\
                        PROCESSINGID=$PROCESSINGID,\
                        DATASET=$DATASET,\
                        ISTEST=$ISTEST,\
                        VERSION_INFO="$VERSION_INFO", \
                        WORKSPACE=/scratch/treebuilding/$PROCESSINGID/$DATASET" | sed 's/ //g') \
        -J t_"$DATASET"_"$PROCESSINGID" \
        -D /shared/workspace/projects/covid/logs \
        -c 16 \
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
