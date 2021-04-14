#!/bin/bash

mkdir -p $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq

aws s3 cp $S3DOWNLOAD/ $WORKSPACE/ --recursive --exclude "*" --include "*fastq.gz"

DELIMITER=_R1_001.fastq.gz
BOTH_FASTQS=$(ls $WORKSPACE |  grep ".fastq.gz" | sort | uniq | grep -v Undetermined)
R1_FASTQS=$(echo "$BOTH_FASTQS" |  grep $DELIMITER )

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
for SAMPLE in $(printf '%s\n' "${SAMPLES_WO_LANES_LIST[@]}" | sort | uniq ); do
  LANES=$(echo "$BOTH_FASTQS" | grep "$SAMPLE" | awk -F $INSPECT_DELIMITER '{print $NF}'| awk -F '_L|_R' '{print $2}' | sort | uniq | grep 00)
  LANES_COMBINED=$(echo $LANES | sed 's/ //g')
  SAMPLE_W_LANES_COMBINED="$SAMPLE"_"$LANES_COMBINED"

  cat $WORKSPACE/"$SAMPLE"*R1_001.fastq.gz > $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq/"$SAMPLE_W_LANES_COMBINED"_R1_001.fastq.gz
  cat $WORKSPACE/"$SAMPLE"*R2_001.fastq.gz > $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq/"$SAMPLE_W_LANES_COMBINED"_R2_001.fastq.gz
done

aws s3 cp $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq/ $S3DOWNLOAD/"$SEQ_RUN"_lane_merged_fastq/ --recursive --exclude "*" --include "*fastq.gz" --include "*merge.log"
rm -rf $WORKSPACE