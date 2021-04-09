#!/bin/bash

PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
QCRESULTS=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/"$SEQ_RUN"_summary_files
S3TEST=s3://ucsd-rtl-test

# Activate conda env covid1.2
ANACONDADIR=/shared/workspace/software/anaconda3/bin
source $ANACONDADIR/activate covid1.2
# clear workspace if node is being reused
rm -rf $WORKSPACE/*
mkdir -p $WORKSPACE/qc

runQC () {

	aws s3 cp $S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/"$SEQ_RUN"_samples/ $WORKSPACE/ \
		--quiet \
		--recursive \
		--exclude "*" \
		--include "*.variants.tsv" \
		--include "*.consensus.fa" \
		--include "*.depth.txt" \
		--include "*q30_reads.txt" \
		--include "*.sorted.stats*" \
		--include "*.acceptance.tsv" \
		--include "*coverage.tsv" \
		--include "*error.log"

	# Zip files
	mv $WORKSPACE/*/*.variants.tsv $WORKSPACE/*/*.consensus.fa $WORKSPACE/*/*.depth.txt $WORKSPACE/*/*.acceptance.tsv $WORKSPACE
	cd $WORKSPACE && zip -9 "$SEQ_RUN"-variants.zip *.variants.tsv && zip -9 "$SEQ_RUN"-consensus.zip *.consensus.fa && zip -9 "$SEQ_RUN"-depth.zip *.depth.txt

	# summary figures and stats
	# echo "Generating a violin plot of mapping depth across all samples and line plots of mapping depth per sample."
	# python $PIPELINEDIR/qc/samtools_depth_plots.py $WORKSPACE/qc/"$SEQ_RUN"-depth_lineplot.pdf $WORKSPACE/qc/"$SEQ_RUN"-depth_violin.pdf $WORKSPACE/*.depth.txt
	# mv depth_violin.pdf $WORKSPACE/qc/"$SEQ_RUN"-depth_violin.pdf
	# mv depth_lineplot.pdf $WORKSPACE/qc/"$SEQ_RUN"-depth_lineplot.pdf
	echo "Summarizing consensus QC."
	python $PIPELINEDIR/qc/seq_run_acceptance.py $WORKSPACE $WORKSPACE/"$SEQ_RUN"-acceptance.tsv
  echo -e "seq_run_acceptance.py exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

	# Multiqc
	echo "Configuring Multiqc"
	find $WORKSPACE -name "qualimapReport.html" | sort -n > $WORKSPACE/qc/qualimapReport_paths.txt
	# for z in $WORKSPACE/*/fastqc/*fastqc.zip; do unzip -q $z -d $WORKSPACE/qc/fastqc; done
	# find $WORKSPACE -name "fastqc_data.txt" | sort -n > $WORKSPACE/qc/fastqc_data_paths.txt
	find $WORKSPACE -name "*q30_reads.txt" | sort -n > $WORKSPACE/qc/q30_reads_paths.txt
	# python $PIPELINEDIR/qc/custom_gen_stats_multiqc.py $WORKSPACE/qc/qualimapReport_paths.txt $WORKSPACE/qc/fastqc_data_paths.txt $FQ $WORKSPACE/multiqc_custom_gen_stats.yaml
	python $PIPELINEDIR/qc/custom_gen_stats_multiqc.py $WORKSPACE/qc/qualimapReport_paths.txt $WORKSPACE/qc/q30_reads_paths.txt $FQ $WORKSPACE/multiqc_custom_gen_stats.yaml
	cat $PIPELINEDIR/qc/covid_custom_config.yaml $WORKSPACE/multiqc_custom_gen_stats.yaml > $WORKSPACE/qc/"$SEQ_RUN"-custom_gen_stats_config.yaml
	multiqc --config $WORKSPACE/qc/"$SEQ_RUN"-custom_gen_stats_config.yaml --module qualimap --module custom_content $WORKSPACE
  echo -e "multiqc exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

    # Concatenate coverage files
    echo "Concatenating coverage files"
    echo -e "SAMPLE\tCOVERAGE\tAVG_DEPTH\tMIN\tMAX\tZERO_DEPTH" > $WORKSPACE/"$SEQ_RUN"-coverage.tsv
    cat $WORKSPACE/*/*coverage.tsv | sort -n -k 2 >> $WORKSPACE/"$SEQ_RUN"-coverage.tsv
    echo -e "coverage cat exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

    # TODO: This is where the download of the full InspectSeq metadata should go

	# Make summary table
	echo "Making run summary table."
	# TODO: extend seq_run_summary.py to include merging in (only relevant records from) InspectSeq metadata
	python $PIPELINEDIR/qc/seq_run_summary.py $WORKSPACE/multiqc_data/multiqc_general_stats.txt $WORKSPACE/"$SEQ_RUN"-acceptance.tsv $WORKSPACE/"$SEQ_RUN"-summary.csv
  echo -e "seq_run_summary.py exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

	# Concatenate all consensus files to a .fas file
	cat $WORKSPACE/*.consensus.fa > $WORKSPACE/"$SEQ_RUN".fas

  # Id only passing consensus files and write them to a *-passQC.fas file
  echo "Merging passing consensus files."
  PASSING_CONS_FNAMES=$(python $PIPELINEDIR/qc/subset_csv.py $WORKSPACE/"$SEQ_RUN"-summary.csv not_na_cons_fnames $WORKSPACE)
  echo -e "subset_csv.py exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log
  cat $PASSING_CONS_FNAMES > $WORKSPACE/"$SEQ_RUN"-passQC.fas

	# Exit codes
	echo "Gathering per-sample exit codes."
	cat $WORKSPACE/*/*error.log > $WORKSPACE/"$SEQ_RUN".error.log
	grep -v "exit code: 0" $WORKSPACE/"$SEQ_RUN"-qc.exit.log | head -n 1 >> $WORKSPACE/"$SEQ_RUN".error.log

	# Upload Results
	echo "Uploading QC and summary results."
	# summary files folder
	aws s3 cp $WORKSPACE/multiqc_data/ $QCRESULTS/"$SEQ_RUN"_multiqc_data/ --recursive --quiet
	aws s3 cp $WORKSPACE/multiqc_report.html $QCRESULTS/"$SEQ_RUN"_multiqc_report.html
	aws s3 cp $WORKSPACE/qc/ $QCRESULTS/ --recursive --quiet
  aws s3 cp $WORKSPACE/"$SEQ_RUN"-summary.csv $QCRESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-acceptance.tsv $QCRESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-coverage.tsv $QCRESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN".error.log $QCRESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-variants.zip $QCRESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-consensus.zip $QCRESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-depth.zip $QCRESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-passQC.fas $QCRESULTS/
	aws s3 cp $WORKSPACE/"$SEQ_RUN".fas $QCRESULTS/

	# cumulative data folder
  if [[ "$ISTEST" = false ]]; then
    S3CUMULATIVE=$S3DOWNLOAD
  else
    S3CUMULATIVE=$S3TEST
  fi
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-passQC.fas $S3CUMULATIVE/phylogeny/cumulative_data/consensus/
	aws s3 cp $WORKSPACE/"$SEQ_RUN".fas $S3CUMULATIVE/phylogeny/cumulative_data/consensus/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-summary.csv $S3CUMULATIVE/phylogeny/cumulative_data/consensus/
}

{ time ( runQC ) ; } > $WORKSPACE/qc/"$SEQ_RUN"-qc_summary.log 2>&1

aws s3 cp $WORKSPACE/qc/"$SEQ_RUN"-qc_summary.log $QCRESULTS/

