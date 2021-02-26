#!/bin/bash

mkdir -p $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq

aws s3 cp $S3DOWNLOAD/ $WORKSPACE/ --recursive --exclude "*" --include "*fastq.gz"
for SAMPLE in $(ls $WORKSPACE | grep fastq.gz | awk -F '_L00' '{print $1}' | sort | uniq | grep -v Undetermined); do
	cat $WORKSPACE/"$SAMPLE"*R1*fastq.gz > $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq/"$SAMPLE"_R1.fastq.gz
	cat $WORKSPACE/"$SAMPLE"*R2*fastq.gz > $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq/"$SAMPLE"_R2.fastq.gz
done

aws s3 cp $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq/ $S3DOWNLOAD/"$SEQ_RUN"_lane_merged_fastq/ --recursive --exclude "*" --include "*fastq.gz"

rm -rf $WORKSPACE