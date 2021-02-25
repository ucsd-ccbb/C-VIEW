#!/bin/bash

mkdir -p $WORKSPACE/merged_lanes

aws s3 cp $S3DOWNLOAD/ $WORKSPACE/ --recursive --exclude "*" --include "*fastq.gz"
for SAMPLE in $(ls $WORKSPACE | grep fastq | awk -F '_L00' '{print $1}' | sort | uniq | grep -v Undetermined); do
	cat $WORKSPACE/"$SAMPLE"*R1*fastq.gz > $WORKSPACE/merged_lanes/"$SAMPLE"_R1.fastq.gz
	cat $WORKSPACE/"$SAMPLE"*R2*fastq.gz > $WORKSPACE/merged_lanes/"$SAMPLE"_R2.fastq.gz
done

aws s3 cp $WORKSPACE/merged_lanes/ $S3DOWNLOAD/merged_lanes/ --recursive --exclude "*" --include "*fastq.gz"

rm -rf $WORKSPACE
