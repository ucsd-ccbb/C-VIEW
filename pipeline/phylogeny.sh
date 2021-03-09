#!/bin/bash

export PATH=$PATH:/shared/workspace/software/IQTree/iqtree-2.1.2-Linux/bin:/shared/workspace/software/viralMSA:/shared/workspace/software/MinVar-Rooting-master:/shared/workspace/software/anaconda3/envs/covid1.2/bin
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
ANACONDADIR=/shared/workspace/software/anaconda3/bin
THREADS=16
rm -rf $WORKSPACE
mkdir -p $WORKSPACE

runPangolin () {

	aws s3 cp $S3DOWNLOAD/consensus/ $WORKSPACE/ --recursive --quiet
	aws s3 cp $S3DOWNLOAD/qc_summary/ $WORKSPACE/ --recursive --quiet
	# TODO: add download of historical fas and metadata files if helix

	cat $WORKSPACE/*passQC.fas > $WORKSPACE/merged.fas
	sed -i -e 's/Consensus_//g' -e 's/.trimmed.sorted.pileup.consensus_threshold_0.5_quality_20//g' $WORKSPACE/merged.fas

	# TODO: merge all historical fas files into merged.fas if helix

	# add the reference sequence
	cat $PIPELINEDIR/reference_files/RmYN02.fas >> $WORKSPACE/merged.fas

	# add B.1.1.7 sequence
	cat $PIPELINEDIR/reference_files/hCoV-19_USA_CA-SEARCH-5574_2020.fasta >> $WORKSPACE/merged.fas

	# pangolin
	source $ANACONDADIR/activate pangolin
	pangolin --update
	pangolin -t $THREADS --outfile $WORKSPACE/merged.lineage_report.csv $WORKSPACE/merged.fas

  # produce merged_qc_and_lineages.csv
  # TODO: must add merge of historical sample names BEFORE lineage merge
  $PIPELINEDIR/qc/lineages_summary.py $WORKSPACE _summary.csv $WORKSPACE/merged.lineage_report.csv $WORKSPACE/merged.qc_and_lineages.csv
  # TODO: lineage merge should be done twice--once as an outer,
  #  once as a right (where only things with lineages stay); latter becomes metadata for
  #  empress
  # TODO: must output empress metadata file as TSV, not csv
  # output as $WORKSPACE/merged.empress_metadata.tsv

	rename 's/merged/'$TIMESTAMP'/' $WORKSPACE/merged.*
	aws s3 cp $WORKSPACE/ $S3UPLOAD/ --recursive --quiet
}

buildTree () {
	# Must use biopy env due to numpy conflicts
	source $ANACONDADIR/activate biopy
	ViralMSA.py -s $WORKSPACE/merged.fas -r SARS-CoV-2 -o $WORKSPACE/viralmsa_out -t $THREADS -e aws-CCBB@health.ucsd.edu

	python $PIPELINEDIR/pipeline/trim_msa.py -i $WORKSPACE/viralmsa_out/merged.fas.aln -s 100 -e 50 -o $WORKSPACE/merged.trimmed.aln

	iqtree2 -T $THREADS -m GTR+F+G4 --polytomy -blmin 1e-9 -s $WORKSPACE/merged.trimmed.aln

	python /shared/workspace/software/MinVar-Rooting-master/FastRoot.py -i $WORKSPACE/merged.trimmed.aln.treefile -o $WORKSPACE/merged.trimmed.aln.rooted.treefile -m OG -g "hCoV-19/bat/Yunnan/RmYN02/2019|EPI_ISL_412977|2019-06-25"

	# tree building 
	source $ANACONDADIR/activate qiime2-2020.11

	empress tree-plot --tree $WORKSPACE/merged.trimmed.aln.rooted.treefile --feature-metadata $WORKSPACE/merged.empress_metadata.tsv --output-dir $WORKSPACE/tree-viz

	rename 's/merged/'$TIMESTAMP'/' $WORKSPACE/merged.*
	rename 's/merged/'$TIMESTAMP'/' $WORKSPACE/viralmsa_out/merged.*
	aws s3 cp $WORKSPACE/ $S3DOWNLOAD/trees/$TIMESTAMP/ --recursive --quiet

}

# Always run Pangolin
{ time ( runPangolin ) ; } > $WORKSPACE/"$TIMESTAMP"-pangolin.log 2>&1
aws s3 cp $WORKSPACE/"$TIMESTAMP"-pangolin.log $S3UPLOAD/

# Run tree building if TREE_BUILD == true
if [[ "$TREE_BUILD" == true ]]; then
	{ time ( buildTree ) ; } > $WORKSPACE/"$TIMESTAMP"-treebuild.log 2>&1
	aws s3 cp $WORKSPACE/"$TIMESTAMP"-treebuild.log $S3DOWNLOAD/trees/$TIMESTAMP/
fi

