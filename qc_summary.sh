#!/bin/bash

results_date=$(date +%F)
WORKSPACE=/scratch
mkdir -p $WORKSPACE/variants
mkdir -p $WORKSPACE/consensus
mkdir -p $WORKSPACE/depth
mkdir -p $WORKSPACE/summary

aws s3 cp $S3DOWNLOAD/variants/ $WORKSPACE/variants/ --recursive --exclude "*" --include "*.variants.tsv"
aws s3 cp $S3DOWNLOAD/consensus/ $WORKSPACE/consensus/ --recursive --exclude "*" --include "*.consensus.fa"
aws s3 cp $S3DOWNLOAD/depth/ $WORKSPACE/depth/ --recursive --exclude "*" --include "*.depth.txt"

cd $WORKSPACE && zip -9 variants/"$results_date"-variants.zip variants/*.variants.tsv
cd $WORKSPACE && zip -9 consensus/"$results_date"-consensus.zip consensus/*.consensus.fa
cd $WORKSPACE && zip -9 depth/"$results_date"-depth.zip depth/*.depth.txt

python /shared/workspace/software/SD-COVID-Sequencing/samtools_depth_violinplot.py $WORKSPACE/depth/*.depth.txt && mv depth_violin.pdf $WORKSPACE/summary
python /shared/workspace/software/SD-COVID-Sequencing/samtools_depth_lineplot.py $WORKSPACE/depth/*.depth.txt && mv depth_lineplot.pdf $WORKSPACE/summary
python /shared/workspace/software/SD-COVID-Sequencing/samtools_depth_concat.py $WORKSPACE/depth/*.depth.txt > $WORKSPACE/summary/depth.tsv
python /shared/workspace/software/SD-COVID-Sequencing/samtools_depth_low.py 266 29674 10 $WORKSPACE/depth/*.depth.txt > $WORKSPACE/summary/depth_below_10.tsv

aws s3 cp $WORKSPACE/variants/"$results_date"-variants.zip $S3DOWNLOAD/
aws s3 cp $WORKSPACE/consensus/"$results_date"-consensus.zip $S3DOWNLOAD/
aws s3 cp $WORKSPACE/depth/"$results_date"-depth.zip $S3DOWNLOAD/
aws s3 cp $WORKSPACE/summary/ $S3DOWNLOAD/ --recursive