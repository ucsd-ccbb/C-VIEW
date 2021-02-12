#!/bin/bash

results_date=$(date +%F)
WORKSPACE=/scratch
STATSDIR=/shared/workspace/software/SD-COVID-Sequencing
MULTIQCCONFIG=/shared/workspace/software/covid_sequencing_analysis_pipeline/covid_custom_config.yaml
mkdir -p $WORKSPACE/variants
mkdir -p $WORKSPACE/consensus
mkdir -p $WORKSPACE/depth
mkdir -p $WORKSPACE/stats
mkdir -p $WORKSPACE/summary
mkdir -p $WORKSPACE/qualimap
mkdir -p $WORKSPACE/fastqc

aws s3 cp $S3DOWNLOAD/variants/ $WORKSPACE/variants/ --recursive --exclude "*" --include "*.variants.tsv"
aws s3 cp $S3DOWNLOAD/consensus/ $WORKSPACE/consensus/ --recursive --exclude "*" --include "*.consensus.fa"
aws s3 cp $S3DOWNLOAD/depth/ $WORKSPACE/depth/ --recursive --exclude "*" --include "*.depth.txt"
aws s3 cp $S3DOWNLOAD/qualimap/ $WORKSPACE/qualimap/ --recursive --include "*.sorted.stats.tar.gz" --exclude "*.trimmed.sorted.stats.tar.gz" 
aws s3 cp $S3DOWNLOAD/fastqc/ $WORKSPACE/fastqc/ --recursive

cd $WORKSPACE && zip -9 variants/"$results_date"-variants.zip variants/*.variants.tsv
cd $WORKSPACE && zip -9 consensus/"$results_date"-consensus.zip consensus/*.consensus.fa
cd $WORKSPACE && zip -9 depth/"$results_date"-depth.zip depth/*.depth.txt
find $WORKSPACE/qualimap/ -name  '*.sorted.stats.tar.gz' -exec tar -xf '{}' --strip-components=2 -C $WORKSPACE/qualimap \; 

# summary figures and stats
python $STATSDIR/samtools_depth_violinplot.py $WORKSPACE/depth/*.depth.txt && mv depth_violin.pdf $WORKSPACE/summary
python $STATSDIR/samtools_depth_lineplot.py $WORKSPACE/depth/*.depth.txt && mv depth_lineplot.pdf $WORKSPACE/summary
python $STATSDIR/samtools_depth_concat.py $WORKSPACE/depth/*.depth.txt > $WORKSPACE/summary/depth.tsv
python $STATSDIR/samtools_depth_low.py 266 29674 10 $WORKSPACE/depth/*.depth.txt > $WORKSPACE/summary/"$results_date"-depth_below_10.tsv
# Multiqc
multiqc --config $MULTIQCCONFIG $WORKSPACE

aws s3 cp $WORKSPACE/variants/"$results_date"-variants.zip $S3DOWNLOAD/
aws s3 cp $WORKSPACE/consensus/"$results_date"-consensus.zip $S3DOWNLOAD/
aws s3 cp $WORKSPACE/depth/"$results_date"-depth.zip $S3DOWNLOAD/
aws s3 cp $WORKSPACE/summary/ $S3DOWNLOAD/ --recursive
aws s3 cp $WORKSPACE/multiqc_data/ $S3DOWNLOAD/multiqc_data/ --recursive
aws s3 cp $WORKSPACE/multiqc_report.html $S3DOWNLOAD/
