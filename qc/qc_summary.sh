#!/bin/bash

export PATH=/shared/workspace/software/anaconda3/envs/covid1.1/bin:$PATH

PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
mkdir -p $WORKSPACE/qc


runQC () {

	aws s3 cp $S3DOWNLOAD/ $WORKSPACE/ \
		--quiet \
		--recursive \
		--exclude "*" \
		--include "*.variants.tsv" \
		--include "*.consensus.fa" \
		--include "*.depth.txt" \
		--include "*fastqc.zip" \
		--include "*.sorted.stats" \
		--include "*.acceptance.tsv"

	mv $WORKSPACE/*/*.variants.tsv $WORKSPACE/*/*.consensus.fa $WORKSPACE/*/*.depth.txt $WORKSPACE/*/*.acceptance.tsv $WORKSPACE
	cat $WORKSPACE/*.consensus.fa > $WORKSPACE/consensus.fas
	cd $WORKSPACE && zip -9 variants.zip *.variants.tsv && zip -9 consensus.zip *.consensus.fa && zip -9 depth.zip *.depth.txt

	# summary figures and stats
	echo "Generating a violin plot of mapping depth across all samples."
	python $PIPELINEDIR/qc/samtools_depth_violinplot.py $WORKSPACE/*.depth.txt && mv depth_violin.pdf $WORKSPACE/qc
	echo "Generating line plots of mapping depth across all samples."
	python $PIPELINEDIR/qc/samtools_depth_lineplot.py $WORKSPACE/*.depth.txt && mv depth_lineplot.pdf $WORKSPACE/qc
	echo "Concatenating samtools depth output across all samples into a single TSV."
	python $PIPELINEDIR/qc/samtools_depth_concat.py $WORKSPACE/*.depth.txt > $WORKSPACE/qc/depth.tsv
	echo "Listing the positions with a depth below 10 reads."
	python $PIPELINEDIR/qc/samtools_depth_low.py 266 29674 10 $WORKSPACE/*.depth.txt > $WORKSPACE/qc/depth_below_10.tsv
	echo "Summarizing consensus QC."
	python $PIPELINEDIR/qc/consensus_acceptance_summary.py $WORKSPACE
	
	# Multiqc
	# Get qualimapReport.html paths
	find $WORKSPACE -name "qualimapReport.html" | sort -n > $PIPELINEDIR/qc/qualimapReport_paths.txt
	# Unzip fastqc.zip and get fastqc_data.txt paths
	unzip $WORKSPACE"/qc/fastqc/*_fastqc.zip" -d $WORKSPACE/qc/fastqc
	find $WORKSPACE -name "fastqc_data.txt" | sort -n > $PIPELINEDIR/qc/fastqc_data_paths.txt
	# Run script to generate custom stats in multiqc_custom_gen_stats.yaml
	python $PIPELINEDIR/qc/custom_gen_stats_multiqc.py $PIPELINEDIR/qc/qualimapReport_paths.txt $PIPELINEDIR/qc/fastqc_data_paths.txt
	cat $PIPELINEDIR/qc/covid_custom_config.yaml multiqc_custom_gen_stats.yaml > $PIPELINEDIR/qc/custom_gen_stats_config.yaml
	
	multiqc --config $PIPELINEDIR/qc/custom_gen_stats_config.yaml --ignore *fastqc $WORKSPACE

	echo "Uploading QC results."
	aws s3 cp $WORKSPACE/variants.zip $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/consensus.zip $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/consensus.fas $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/depth.zip $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/multiqc_data/ $S3DOWNLOAD/qc/multiqc_data/ --recursive
	aws s3 cp $WORKSPACE/multiqc_report.html $S3DOWNLOAD/qc/
	aws s3 cp $WORKSPACE/qc/ $S3DOWNLOAD/qc/ --recursive
	aws s3 cp $WORKSPACE/summary.acceptance.tsv $S3DOWNLOAD/qc/

}

{ time ( runQC ) ; } > $WORKSPACE/qc/qc_summary.log 2>&1

aws s3 cp $WORKSPACE/qc/qc_summary.log $S3DOWNLOAD/qc/

