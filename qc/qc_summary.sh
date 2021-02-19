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
		--include "*qualimapReport.html" \
		--include "*.acceptance.tsv"

	# Zip files
	mv $WORKSPACE/*/*.variants.tsv $WORKSPACE/*/*.consensus.fa $WORKSPACE/*/*.depth.txt $WORKSPACE/*/*.acceptance.tsv $WORKSPACE
	cd $WORKSPACE && zip -9 "$BATCH"-variants.zip *.variants.tsv && zip -9 "$BATCH"-consensus.zip *.consensus.fa && zip -9 "$BATCH"-depth.zip *.depth.txt

	# summary figures and stats
	echo "Generating a violin plot of mapping depth across all samples."
	python $PIPELINEDIR/qc/samtools_depth_violinplot.py $WORKSPACE/*.depth.txt && mv depth_violin.pdf $WORKSPACE/qc/"$BATCH"-depth_violin.pdf
	echo "Generating line plots of mapping depth across all samples."
	python $PIPELINEDIR/qc/samtools_depth_lineplot.py $WORKSPACE/*.depth.txt && mv depth_lineplot.pdf $WORKSPACE/qc/"$BATCH"-depth_violin.pdf
	echo "Summarizing consensus QC."
	python $PIPELINEDIR/qc/consensus_acceptance_summary.py $WORKSPACE
	mv $WORKSPACE/summary.acceptance.tsv $WORKSPACE/"$BATCH"-summary.acceptance.tsv

	# Multiqc
	echo "Configuring Multiqc"
	find $WORKSPACE -name "qualimapReport.html" | sort -n > $WORKSPACE/qc/qualimapReport_paths.txt
	for z in $WORKSPACE/*/fastqc/*fastqc.zip; do unzip -q $z -d $WORKSPACE/qc/fastqc; done
	find $WORKSPACE -name "fastqc_data.txt" | sort -n > $WORKSPACE/qc/fastqc_data_paths.txt
	python $PIPELINEDIR/qc/custom_gen_stats_multiqc.py $WORKSPACE/qc/qualimapReport_paths.txt $WORKSPACE/qc/fastqc_data_paths.txt
	cat $PIPELINEDIR/qc/covid_custom_config.yaml $WORKSPACE/multiqc_custom_gen_stats.yaml > $WORKSPACE/qc/"$BATCH"-custom_gen_stats_config.yaml
	multiqc --config $WORKSPACE/qc/"$BATCH"-custom_gen_stats_config.yaml --ignore *fastqc $WORKSPACE

	# Pangolin
	cat $WORKSPACE/*.consensus.fa > $WORKSPACE/"$BATCH".fas
	bash $PIPELINEDIR/qc/pangolin.sh $WORKSPACE/$BATCH

	echo "Uploading QC results."
	aws s3 cp $WORKSPACE/"$BATCH"-variants.zip $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/"$BATCH"-consensus.zip $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/"$BATCH"-depth.zip $S3DOWNLOAD/

	aws s3 cp $WORKSPACE/multiqc_data/ $S3DOWNLOAD/"$BATCH"-qc/multiqc_data/ --recursive --quiet
	aws s3 cp $WORKSPACE/multiqc_report.html $S3DOWNLOAD/"$BATCH"-qc/
	aws s3 cp $WORKSPACE/qc/ $S3DOWNLOAD/"$BATCH"-qc/ --recursive --quiet
	
	aws s3 cp $WORKSPACE/"$BATCH"-summary.acceptance.tsv $S3DOWNLOAD/"$BATCH"-qc/
	aws s3 cp $WORKSPACE/"$BATCH"-passQC.samples.tsv $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/"$BATCH"-passQC.fas $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/"$BATCH".fas $S3DOWNLOAD/
	aws s3 cp $WORKSPACE/"$BATCH".lineage_report.csv $S3DOWNLOAD/
}

{ time ( runQC ) ; } > $WORKSPACE/qc/"$BATCH"-qc_summary.log 2>&1

aws s3 cp $WORKSPACE/qc/"$BATCH"-qc_summary.log $S3UPLOAD/"$BATCH"-qc/

