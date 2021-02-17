#!/bin/bash

WORKSPACE=/scratch
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
mkdir -p $WORKSPACE/qc

aws s3 cp $S3DOWNLOAD/ $WORKSPACE/ \
	--recursive \
	--exclude "*" \
	--include "*.variants.tsv" \
	--include "*.consensus.fa" \
	--include "*.depth.txt" \
	--include "*fastqc.zip" \
	--include "*.sorted.stats" \
	--include "*.acceptance.tsv"

cd $WORKSPACE && zip -9 variants.zip */*.variants.tsv
cd $WORKSPACE && zip -9 consensus.zip */*.consensus.fa
cd $WORKSPACE && zip -9 depth.zip */*.depth.txt

# summary figures and stats
python $PIPELINEDIR/qc/samtools_depth_violinplot.py $WORKSPACE/*/*.depth.txt && mv depth_violin.pdf $WORKSPACE/qc
python $PIPELINEDIR/qc/samtools_depth_lineplot.py $WORKSPACE/*/*.depth.txt && mv depth_lineplot.pdf $WORKSPACE/qc
python $PIPELINEDIR/qc/samtools_depth_concat.py $WORKSPACE/*/*.depth.txt > $WORKSPACE/qc/depth.tsv
python $PIPELINEDIR/qc/samtools_depth_low.py 266 29674 10 $WORKSPACE/*/*.depth.txt > $WORKSPACE/qc/depth_below_10.tsv
python $PIPELINEDIR/qc/consensus_acceptance_summary.py 

# Multiqc
multiqc --config $PIPELINEDIR/qc/covid_custom_config.yaml $WORKSPACE

aws s3 cp $WORKSPACE/variants.zip $S3DOWNLOAD/
aws s3 cp $WORKSPACE/consensus.zip $S3DOWNLOAD/
aws s3 cp $WORKSPACE/depth.zip $S3DOWNLOAD/
aws s3 cp $WORKSPACE/multiqc_data/ $S3DOWNLOAD/qc/ --recursive
aws s3 cp $WORKSPACE/multiqc_report.html $S3DOWNLOAD/qc/
aws s3 cp $WORKSPACE/qc/ $S3DOWNLOAD/qc/ --recursive
aws s3 cp $WORKSPACE/summary.acceptance.tsv $S3DOWNLOAD/qc/
