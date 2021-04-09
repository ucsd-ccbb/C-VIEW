#!/bin/bash

mkdir -p $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq

#aws s3 cp $S3DOWNLOAD/ $WORKSPACE/ --recursive --exclude "*" --include "*fastq.gz"
IFS=' ' read -ra SAMPLES_WO_LANE_INFO <<< "$SAMPLES_WO_LANE_INFO_STR"
IFS=' ' read -ra SAMPLES_W_LANES_COMBINED <<< "$SAMPLES_W_LANES_COMBINED_STR"

# sanity check that arrays are same length
if [[ ! ${#SAMPLES_WO_LANE_INFO[@]} == ${#SAMPLES_W_LANES_COMBINED[@]} ]] ; then
   echo "Error: must have same number of sample-without-lane and sample-with-lane-combined entries" >> $WORKSPACE/merge.log
   aws s3 cp $WORKSPACE/ $S3DOWNLOAD/"$SEQ_RUN"_lane_merged_fastq/ --recursive --exclude "*" --include "*.log"
   exit 1
fi

END=${#SAMPLES_WO_LANE_INFO[@]}
for (( c=1; c<=$END; c++ ))
do
  SAMPLE=${SAMPLES_WO_LANE_INFO[$i]}
  SAMPLE_W_LANES_COMBINED=${SAMPLES_W_LANES_COMBINED[$i]}
  echo $SAMPLE >> merge.log
  echo $SAMPLE_W_LANES_COMBINED >> $WORKSPACE/merge.log
	#cat $WORKSPACE/"$SAMPLE"*R1_001.fastq.gz > $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq/"$SAMPLE_W_LANES_COMBINED"_R1_001.fastq.gz
	#cat $WORKSPACE/"$SAMPLE"*R2_001.fastq.gz > $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq/"$SAMPLE_W_LANES_COMBINED"_R2_001.fastq.gz
done

aws s3 cp $WORKSPACE/ $S3DOWNLOAD/"$SEQ_RUN"_lane_merged_fastq/ --recursive --exclude "*" --include "*.log"

#aws s3 cp $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq/ $S3DOWNLOAD/"$SEQ_RUN"_lane_merged_fastq/ --recursive --exclude "*" --include "*fastq.gz"
rm -rf $WORKSPACE