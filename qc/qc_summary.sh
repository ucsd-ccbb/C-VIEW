#!/bin/bash

PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
PHYLORESULTS=$S3DOWNLOAD/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/"$SEQ_RUN"_phylogenetic_results
QCRESULTS=$S3DOWNLOAD/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/"$SEQ_RUN"_quality_control
# Activate conda env covid1.2
ANACONDADIR=/shared/workspace/software/anaconda3/bin
source $ANACONDADIR/activate covid1.2
# clear workspace if node is being reused
rm -rf $WORKSPACE/*
mkdir -p $WORKSPACE/qc/fastqc

runQC () {

	aws s3 cp $S3DOWNLOAD/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/"$SEQ_RUN"_samples/ $WORKSPACE/ \
		--quiet \
		--recursive \
		--exclude "*" \
		--include "*.variants.tsv" \
		--include "*.consensus.fa" \
		--include "*.depth.txt" \
		--include "*fastqc.zip" \
		--include "*.sorted.stats*" \
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
	python $PIPELINEDIR/qc/custom_gen_stats_multiqc.py $WORKSPACE/qc/qualimapReport_paths.txt $WORKSPACE/qc/fastqc_data_paths.txt $FQ
	cat $PIPELINEDIR/qc/covid_custom_config.yaml $WORKSPACE/multiqc_custom_gen_stats.yaml > $WORKSPACE/qc/"$SEQ_RUN"-custom_gen_stats_config.yaml
	multiqc --config $WORKSPACE/qc/"$SEQ_RUN"-custom_gen_stats_config.yaml --module qualimap $WORKSPACE

	# # Make merged consensus of passing samples
	cat $WORKSPACE/*.consensus.fa > $WORKSPACE/"$SEQ_RUN".fas
	# filter the true samples
	awk '{ if ($2 == "True") { print } }' $WORKSPACE/"$SEQ_RUN"-summary.acceptance.tsv > $WORKSPACE/"$SEQ_RUN"-summary.acceptance.true.tsv
	awk '{print $1}' $WORKSPACE/"$SEQ_RUN"-summary.acceptance.true.tsv > $WORKSPACE/"$SEQ_RUN"-passQC.samples.tsv

	# loop over individual .fa files, keep the ones which are in passQC.samples.tsv
	touch $WORKSPACE/"$SEQ_RUN"-passQC.fas # initialize the file
	for f in *.fa; do
	    fshort="$(cut -d'.' -f1 <<<$f)"
	    echo $fshort
	    if (grep -qF $fshort $WORKSPACE/"$SEQ_RUN"-passQC.samples.tsv); then
	       cat $f >> $WORKSPACE/"$SEQ_RUN"-passQC.fas
	    fi
	done

	# Pangolin
	#bash $PIPELINEDIR/qc/pangolin.sh $WORKSPACE/$SEQ_RUN

	# Make QC table
	python $PIPELINEDIR/qc/makeQCSummaryTable.py $WORKSPACE/multiqc_data/multiqc_general_stats.txt $WORKSPACE/"$SEQ_RUN"-summary.acceptance.tsv
	mv $WORKSPACE/QCSummaryTable.csv $WORKSPACE/"$SEQ_RUN"-QCSummaryTable.csv

	# Upload Results
	echo "Uploading QC and summary results."
	# phylogenetic results folder
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-variants.zip $PHYLORESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-consensus.zip $PHYLORESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-depth.zip $PHYLORESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-passQC.samples.tsv $PHYLORESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-passQC.fas $PHYLORESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN".fas $PHYLORESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-summary.acceptance.tsv $PHYLORESULTS/
	#aws s3 cp $WORKSPACE/"$SEQ_RUN".lineage_report.csv $PHYLORESULTS/

	# quality control folder
	aws s3 cp $WORKSPACE/multiqc_data/ $QCRESULTS/"$SEQ_RUN"_multiqc_data/ --recursive --quiet
	aws s3 cp $WORKSPACE/multiqc_report.html $QCRESULTS/"$SEQ_RUN"_multiqc_report.html
	aws s3 cp $WORKSPACE/qc/ $QCRESULTS/ --recursive --quiet
  aws s3 cp $WORKSPACE/"$SEQ_RUN"-QCSummaryTable.csv $QCRESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-summary.acceptance.tsv $QCRESULTS/

	# Tree building data
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-passQC.fas s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/tree_building/consensus/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-summary.acceptance.tsv s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/tree_building/acceptance/
  aws s3 cp $WORKSPACE/"$SEQ_RUN"-QCSummaryTable.csv s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/tree_building/qc_summary/
}

{ time ( runQC ) ; } > $WORKSPACE/qc/"$SEQ_RUN"-qc_summary.log 2>&1

aws s3 cp $WORKSPACE/qc/"$SEQ_RUN"-qc_summary.log $QCRESULTS/

