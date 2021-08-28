#!/bin/bash

PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
QCRESULTS=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/"$SEQ_RUN"_summary_files

# Activate conda env covid1.2
ANACONDADIR=/shared/workspace/software/anaconda3/bin
source $ANACONDADIR/activate covid1.2
# clear workspace if node is being reused
rm -rf $WORKSPACE/*
mkdir -p $WORKSPACE/qc

runQC () {
  echo "$VERSION_INFO" >> $WORKSPACE/"$SEQ_RUN".version.log

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
		# TODO: would be nice if these three metrics were defined in just one place in this script rather than 3
		--include "*pi-metric.tsv" \
		--include "*n-metric.tsv" \
		--include "*trimmed_bam_read_count.tsv"
		# end TODO
		--include "*_subsampled_mapping_stats.tsv" \
		--include "*error.log"

	# Exit codes
	echo "Gathering per-sample exit codes."
	cat $WORKSPACE/*/*error.log > $WORKSPACE/"$SEQ_RUN".error.log
	grep -v "exit code: 0" $WORKSPACE/"$SEQ_RUN"-qc.exit.log | head -n 1 >> $WORKSPACE/"$SEQ_RUN".error.log

  # if the error log is NOT empty
  # create a complete_w_failure.txt sentinel file holding these errors
  # upload that to the output dir
  # exit

	## Zip files
	# mv $WORKSPACE/*/*.variants.tsv $WORKSPACE/*/*.consensus.fa $WORKSPACE/*/*.depth.txt $WORKSPACE/*/*.acceptance.tsv $WORKSPACE
	# cd $WORKSPACE && zip -9 "$SEQ_RUN"-variants.zip *.variants.tsv && zip -9 "$SEQ_RUN"-consensus.zip *.consensus.fa && zip -9 "$SEQ_RUN"-depth.zip *.depth.txt

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
	find $WORKSPACE -name "*q30_reads.txt" | sort -n > $WORKSPACE/qc/q30_reads_paths.txt
	find $WORKSPACE -name "*_subsampled_mapping_stats.tsv" | sort -n > $WORKSPACE/qc/subsampled_mapping_stats_paths.txt

	python $PIPELINEDIR/qc/custom_gen_stats_multiqc.py $WORKSPACE/qc/qualimapReport_paths.txt $WORKSPACE/qc/q30_reads_paths.txt $WORKSPACE/qc/subsampled_mapping_stats_paths.txt $FQ $WORKSPACE/multiqc_custom_gen_stats.yaml
  echo -e "custom_gen_stats_multiqc.py exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

	cat $PIPELINEDIR/qc/covid_custom_config.yaml $WORKSPACE/multiqc_custom_gen_stats.yaml > $WORKSPACE/qc/"$SEQ_RUN"-custom_gen_stats_config.yaml
	multiqc --config $WORKSPACE/qc/"$SEQ_RUN"-custom_gen_stats_config.yaml --module qualimap --module custom_content $WORKSPACE
    echo -e "multiqc exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

    # Concatenate coverage files
    echo "Concatenating coverage files"
    echo -e "SAMPLE\tCOVERAGE\tAVG_DEPTH\tMIN\tMAX\tZERO_DEPTH" > $WORKSPACE/"$SEQ_RUN"-coverage.tsv
    cat $WORKSPACE/*/*coverage.tsv | sort -n -k 2 >> $WORKSPACE/"$SEQ_RUN"-coverage.tsv
    echo -e "coverage cat exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

  # Concatenate per-sample files for per-sample metrics
  for PER_SAMPLE_METRIC in trimmed_bam_read_count n-metric pi-metric ; do
    echo "Concatenating $PER_SAMPLE_METRIC files"
    echo -e "sequenced_pool_component_id\t$PER_SAMPLE_METRIC" > $WORKSPACE/"$SEQ_RUN"-"$PER_SAMPLE_METRIC".tsv
    cat $WORKSPACE/*/*"$PER_SAMPLE_METRIC".tsv | sort -n -k 2 >> $WORKSPACE/"$SEQ_RUN"-"$PER_SAMPLE_METRIC".tsv
    echo -e "$PER_SAMPLE_METRIC cat exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log
  done

	# Make summary table
	echo "Making run summary table."
	# TODO: would be nice to use the same list of per-sample metrics here as above to avoid having to remember to add new ones in two places
	python $PIPELINEDIR/qc/seq_run_summary.py $WORKSPACE/"$SEQ_RUN"-temp-summary.csv $WORKSPACE/multiqc_data/multiqc_general_stats.txt $WORKSPACE/"$SEQ_RUN"-acceptance.tsv $WORKSPACE/"$SEQ_RUN"-trimmed_bam_read_count.tsv $WORKSPACE/"$SEQ_RUN"-pi-metric.tsv $WORKSPACE/"$SEQ_RUN"-n-metric.tsv
    echo -e "seq_run_summary.py exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log
	python $PIPELINEDIR/qc/integrate_bjorn_coverage.py $WORKSPACE/"$SEQ_RUN"-temp-summary.csv $WORKSPACE/"$SEQ_RUN"-coverage.tsv $WORKSPACE/"$SEQ_RUN"-summary.csv
    echo -e "integrate_bjorn_coverage.py exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

	# Concatenate all consensus files to a .fas file
	cat $WORKSPACE/*.consensus.fa > $WORKSPACE/"$SEQ_RUN".fas

    # Id only passing consensus files and write them to a *-passQC.fas file
    echo "Merging passing consensus files."
    PASSING_CONS_FNAMES=$(python $PIPELINEDIR/qc/subset_csv.py $WORKSPACE/"$SEQ_RUN"-summary.csv not_na_cons_fnames $WORKSPACE)
    echo -e "subset_csv.py exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log
    cat $PASSING_CONS_FNAMES > $WORKSPACE/"$SEQ_RUN"-passQC.fas

  # CURRDIR=$(pwd)
  # cd $PIPELINEDIR
  # bash $PIPELINEDIR/show_version.sh >> $WORKSPACE/"$SEQ_RUN".version.log
  # echo -e "show_version.sh exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log
  # cd $CURRDIR

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
#    aws s3 cp $WORKSPACE/"$SEQ_RUN"-summary.csv $QCRESULTS/
#	aws s3 cp $WORKSPACE/"$SEQ_RUN"-acceptance.tsv $QCRESULTS/
#	aws s3 cp $WORKSPACE/"$SEQ_RUN"-coverage.tsv $QCRESULTS/
#	aws s3 cp $WORKSPACE/"$SEQ_RUN".error.log $QCRESULTS/
#	aws s3 cp $WORKSPACE/"$SEQ_RUN"-variants.zip $QCRESULTS/
#	aws s3 cp $WORKSPACE/"$SEQ_RUN"-consensus.zip $QCRESULTS/
#	aws s3 cp $WORKSPACE/"$SEQ_RUN"-depth.zip $QCRESULTS/
#	aws s3 cp $WORKSPACE/"$SEQ_RUN"-passQC.fas $QCRESULTS/
#	aws s3 cp $WORKSPACE/"$SEQ_RUN".fas $QCRESULTS/
#	aws s3 cp $WORKSPACE/"$SEQ_RUN".version.log $QCRESULTS/

	# cumulative data folder
	S3CUMULATIVE=$S3DOWNLOAD
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-passQC.fas $S3CUMULATIVE/phylogeny/cumulative_data/consensus/
	aws s3 cp $WORKSPACE/"$SEQ_RUN".fas $S3CUMULATIVE/phylogeny/cumulative_data/consensus/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-summary.csv $S3CUMULATIVE/phylogeny/cumulative_data/consensus/
}

{ time ( runQC ) ; } > $WORKSPACE/qc/"$SEQ_RUN"-qc_summary.log 2>&1

# copy only top-level results
aws s3 cp $WORKSPACE $QCRESULTS/ --recursive --include "*.*" --exclude "*/*.*"
aws s3 cp $WORKSPACE/qc/"$SEQ_RUN"-qc_summary.log $QCRESULTS/

