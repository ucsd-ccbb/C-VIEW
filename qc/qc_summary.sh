#!/bin/bash

export PATH=/shared/workspace/software/anaconda3/envs/covid1.1/bin:$PATH

PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
mkdir -p $WORKSPACE/qc/fastqc

runQC () {

	aws s3 cp $S3DOWNLOAD/ $WORKSPACE/ \
		--quiet \
		--recursive \
		--exclude "*" \
		--include "*.variants.tsv" \
		--include "*.consensus.fa" \
		--include "*.depth.txt" \
		--include "*fastqc.zip" \
		--include "*sorted.stats*" \
		--include "*.acceptance.tsv"

	# Zip files
	mv $WORKSPACE/*/*.variants.tsv $WORKSPACE/*/*.consensus.fa $WORKSPACE/*/*.depth.txt $WORKSPACE/*/*.acceptance.tsv $WORKSPACE
	cd $WORKSPACE && zip -9 "$SEQ_RUN"-variants.zip *.variants.tsv && zip -9 "$SEQ_RUN"-consensus.zip *.consensus.fa && zip -9 "$SEQ_RUN"-depth.zip *.depth.txt

	# summary figures and stats
	echo "Generating a violin plot of mapping depth across all samples and line plots of mapping depth per sample."
	python $PIPELINEDIR/qc/samtools_depth_plots.py $WORKSPACE/*.depth.txt
	mv depth_violin.pdf $WORKSPACE/qc/"$SEQ_RUN"-depth_violin.pdf
	mv depth_lineplot.pdf $WORKSPACE/qc/"$SEQ_RUN"-depth_lineplot.pdf
	echo "Summarizing consensus QC."
	python $PIPELINEDIR/qc/consensus_acceptance_summary.py $WORKSPACE
	mv $WORKSPACE/summary.acceptance.tsv $WORKSPACE/"$SEQ_RUN"-summary.acceptance.tsv

	# Multiqc
	echo "Configuring Multiqc"
	find $WORKSPACE -name "qualimapReport.html" | sort -n > $WORKSPACE/qc/qualimapReport_paths.txt
	for z in $WORKSPACE/*/fastqc/*fastqc.zip; do unzip -q $z -d $WORKSPACE/qc/fastqc; done
	find $WORKSPACE -name "fastqc_data.txt" | sort -n > $WORKSPACE/qc/fastqc_data_paths.txt
	python $PIPELINEDIR/qc/custom_gen_stats_multiqc.py $WORKSPACE/qc/qualimapReport_paths.txt $WORKSPACE/qc/fastqc_data_paths.txt
	cat $PIPELINEDIR/qc/covid_custom_config.yaml $WORKSPACE/multiqc_custom_gen_stats.yaml > $WORKSPACE/qc/"$SEQ_RUN"-custom_gen_stats_config.yaml
	multiqc --config $WORKSPACE/qc/"$SEQ_RUN"-custom_gen_stats_config.yaml --module qualimap $WORKSPACE

	# Pangolin
	cat $WORKSPACE/*.consensus.fa > $WORKSPACE/"$SEQ_RUN".fas
	bash $PIPELINEDIR/qc/pangolin.sh $WORKSPACE/$SEQ_RUN

	echo "Uploading QC results."
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-variants.zip $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-consensus.zip $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-depth.zip $S3DOWNLOAD/

	aws s3 cp $WORKSPACE/multiqc_data/ $S3DOWNLOAD/"$SEQ_RUN"-qc/multiqc_data/ --recursive --quiet
	aws s3 cp $WORKSPACE/multiqc_report.html $S3DOWNLOAD/"$SEQ_RUN"-qc/
	aws s3 cp $WORKSPACE/qc/ $S3DOWNLOAD/"$SEQ_RUN"-qc/ --recursive --quiet
	
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-summary.acceptance.tsv $S3DOWNLOAD/"$SEQ_RUN"-qc/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-passQC.samples.tsv $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-passQC.fas $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/"$SEQ_RUN".fas $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/"$SEQ_RUN".lineage_report.csv $S3DOWNLOAD/
}

{ time ( runQC ) ; } > $WORKSPACE/qc/"$SEQ_RUN"-qc_summary.log 2>&1

aws s3 cp $WORKSPACE/qc/"$SEQ_RUN"-qc_summary.log $S3UPLOAD/"$SEQ_RUN"-qc/

