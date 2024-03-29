#!/bin/bash

CVIEWDIR=/shared/workspace/software/cview
QCRESULTS=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/"$SEQ_RUN"_summary_files

# Activate conda env cview
ANACONDADIR=/shared/workspace/software/anaconda3/bin
source $ANACONDADIR/activate cview
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
		--include "*pi_metric.tsv" \
		--include "*n_metric.tsv" \
		--include "*trimmed_bam_read_count.tsv" \
		--include "*_subsampled_mapping_stats.tsv" \
		--include "*error.log"
		# TODO: would be nice if the three per-sample metrics were defined in just one place in this script rather than 3

	## Zip files
	mv $WORKSPACE/*/*.depth.txt $WORKSPACE
	mv $WORKSPACE/*/*.consensus.fa $WORKSPACE/*/*.acceptance.tsv $WORKSPACE
	cd $WORKSPACE
  zip -9 "$SEQ_RUN"-depth.zip *.depth.txt

	# summary figures and stats
	# echo "Generating a violin plot of mapping depth across all samples and line plots of mapping depth per sample."
	# python $CVIEWDIR/src/samtools_depth_plots.py $WORKSPACE/qc/"$SEQ_RUN"-depth_lineplot.pdf $WORKSPACE/qc/"$SEQ_RUN"-depth_violin.pdf $WORKSPACE/*.depth.txt
	# mv depth_violin.pdf $WORKSPACE/qc/"$SEQ_RUN"-depth_violin.pdf
	# mv depth_lineplot.pdf $WORKSPACE/qc/"$SEQ_RUN"-depth_lineplot.pdf
	echo "Summarizing consensus QC."
	python $CVIEWDIR/src/seq_run_acceptance.py $WORKSPACE $WORKSPACE/"$SEQ_RUN"-acceptance.tsv
    echo -e "seq_run_acceptance.py exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

	# Multiqc
	echo "Configuring Multiqc"
	find $WORKSPACE -name "qualimapReport.html" | sort -n > $WORKSPACE/qc/qualimapReport_paths.txt
	find $WORKSPACE -name "*q30_reads.txt" | sort -n > $WORKSPACE/qc/q30_reads_paths.txt
	find $WORKSPACE -name "*_subsampled_mapping_stats.tsv" | sort -n > $WORKSPACE/qc/subsampled_mapping_stats_paths.txt

	python $CVIEWDIR/src/custom_gen_stats_multiqc.py $WORKSPACE/qc/qualimapReport_paths.txt $WORKSPACE/qc/q30_reads_paths.txt $WORKSPACE/qc/subsampled_mapping_stats_paths.txt $FQ $WORKSPACE/multiqc_custom_gen_stats.yaml
  echo -e "custom_gen_stats_multiqc.py exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

	cat $CVIEWDIR/src/cview_multiqc_config.yaml $WORKSPACE/multiqc_custom_gen_stats.yaml > $WORKSPACE/qc/"$SEQ_RUN"-custom_gen_stats_config.yaml
	multiqc --config $WORKSPACE/qc/"$SEQ_RUN"-custom_gen_stats_config.yaml --module qualimap --module custom_content $WORKSPACE
    echo -e "multiqc exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

    # Concatenate coverage files
    echo "Concatenating coverage files"
    echo -e "SAMPLE\tCOVERAGE\tAVG_DEPTH\tMIN\tMAX\tZERO_DEPTH" > $WORKSPACE/"$SEQ_RUN"-coverage.tsv
    cat $WORKSPACE/*/*coverage.tsv | sort -n -k 2 >> $WORKSPACE/"$SEQ_RUN"-coverage.tsv
    echo -e "coverage cat exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

  # Concatenate per-sample files for per-sample metrics
  for PER_SAMPLE_METRIC in trimmed_bam_read_count n_metric pi_metric ; do
    echo "Concatenating $PER_SAMPLE_METRIC files"
    echo -e "sequenced_pool_component_id\t$PER_SAMPLE_METRIC" > $WORKSPACE/"$SEQ_RUN"-"$PER_SAMPLE_METRIC".tsv
    cat $WORKSPACE/*/*"$PER_SAMPLE_METRIC".tsv | sort -n -k 2 >> $WORKSPACE/"$SEQ_RUN"-"$PER_SAMPLE_METRIC".tsv
    echo -e "$PER_SAMPLE_METRIC cat exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log
  done

	# Make summary table
	echo "Making run summary table."
	# TODO: would be nice to use the same list of per-sample metrics here as above to avoid having to remember to add new ones in two places
	python $CVIEWDIR/src/seq_run_summary.py $WORKSPACE/"$SEQ_RUN"-temp-summary.csv $WORKSPACE/multiqc_data/multiqc_general_stats.txt $WORKSPACE/"$SEQ_RUN"-acceptance.tsv $WORKSPACE/"$SEQ_RUN"-trimmed_bam_read_count.tsv $WORKSPACE/"$SEQ_RUN"-pi_metric.tsv $WORKSPACE/"$SEQ_RUN"-n_metric.tsv
    echo -e "seq_run_summary.py exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log
	python $CVIEWDIR/src/integrate_bjorn_coverage.py $WORKSPACE/"$SEQ_RUN"-temp-summary.csv $WORKSPACE/"$SEQ_RUN"-coverage.tsv $WORKSPACE/"$SEQ_RUN"-summary.csv
    echo -e "integrate_bjorn_coverage.py exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

	# Concatenate all consensus files to a .fas file
	cat $WORKSPACE/*.consensus.fa > $WORKSPACE/"$SEQ_RUN".fas

    # Id only passing consensus files and write them to a *-passQC.fas file
    echo "Merging passing consensus files."
    PASSING_CONS_FNAMES=$(python $CVIEWDIR/src/subset_csv.py $WORKSPACE/"$SEQ_RUN"-summary.csv not_na_cons_fnames $WORKSPACE)
    echo -e "subset_csv.py exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log
    cat $PASSING_CONS_FNAMES > $WORKSPACE/"$SEQ_RUN"-passQC.fas

  # generate file of summary and passQC fas checksums, for record-keeping
  python $CVIEWDIR/src/document_file_checksums.py \
    $WORKSPACE $WORKSPACE/"$SEQ_RUN"_artifact_checksums.csv \
    "$SEQ_RUN-summary.csv" "$SEQ_RUN-passQC.fas"
  echo -e "document_file_checksums.py exit code: $?" >> $WORKSPACE/"$SEQ_RUN"-qc.exit.log

	# Exit codes
	echo "Gathering per-sample exit codes."
	cat $WORKSPACE/*/*error.log > $WORKSPACE/"$SEQ_RUN".error.log
	grep -v "exit code: 0" $WORKSPACE/"$SEQ_RUN"-qc.exit.log | head -n 1 >> $WORKSPACE/"$SEQ_RUN".error.log

	# Upload Results
	echo "Uploading multiqc and qc-subfolder results."
	# summary files folder
	aws s3 cp $WORKSPACE/multiqc_data/ $QCRESULTS/"$SEQ_RUN"_multiqc_data/ --recursive --quiet
	aws s3 cp $WORKSPACE/multiqc_report.html $QCRESULTS/"$SEQ_RUN"_multiqc_report.html
	aws s3 cp $WORKSPACE/qc/ $QCRESULTS/ --recursive --quiet

	# cumulative data folder
	S3CUMULATIVE=$S3DOWNLOAD
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-passQC.fas $S3CUMULATIVE/phylogeny/cumulative_data/consensus/
	#aws s3 cp $WORKSPACE/"$SEQ_RUN".fas $S3CUMULATIVE/phylogeny/cumulative_data/consensus/
	aws s3 cp $WORKSPACE/"$SEQ_RUN"-summary.csv $S3CUMULATIVE/phylogeny/cumulative_data/consensus/
}

{ time ( runQC ) ; } > $WORKSPACE/qc/"$SEQ_RUN"-qc_summary.log 2>&1

# upload only top-level results here
aws s3 cp $WORKSPACE $QCRESULTS/ --recursive --include "*.*" --exclude "*/*.*" --exclude "*.consensus.fa" --exclude "*.acceptance.tsv" --exclude "*.depth.txt"
aws s3 cp $WORKSPACE/qc/"$SEQ_RUN"-qc_summary.log $QCRESULTS/

