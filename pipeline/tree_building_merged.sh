
export PATH=$PATH:/shared/workspace/software/IQTree/iqtree-2.1.2-Linux/bin:/shared/workspace/software/viralMSA:/shared/workspace/software/MinVar-Rooting-master:/shared/workspace/software/anaconda3/envs/pangolin/bin
mkdir -p $WORKSPACE
# S3DOWNLOAD=s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/tree_building
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
ANACONDADIR=/shared/workspace/software/anaconda3/bin
THREADS=96

buildTree () {

	aws s3 cp $S3DOWNLOAD/consensus/ $WORKSPACE/ --recursive
	aws s3 cp $S3DOWNLOAD/acceptance/ $WORKSPACE/ --recursive

	cat $WORKSPACE/*passQC.fas > $WORKSPACE/merged.fas

	# add the reference sequence
	cat $PIPELINEDIR/reference_files/RmYN02.fas >> $WORKSPACE/merged.fas

	# add B.1.1.7 sequence
	cat $PIPELINEDIR/reference_files/hCoV-19_USA_CA-SEARCH-5574_2020.fasta >> $WORKSPACE/merged.fas

	source $ANACONDADIR/activate covid1.1
	ViralMSA.py -s $WORKSPACE/merged.fas -r SARS-CoV-2 -o viralmsa_out -t $THREADS -e aws-CCBB@health.ucsd.edu

	python $PIPELINEDIR/pipeline/trim_msa.py -i $WORKSPACE/viralmsa_out/merged.fas.aln -s 100 -e 50 -o $WORKSPACE/merged.trimmed.aln

	iqtree2 -T $THREADS -m GTR+F+G4 --polytomy -blmin 1e-9 -s $WORKSPACE/merged.trimmed.aln

	python FastRoot.py -i $WORKSPACE/merged.trimmed.aln.treefile -o $WORKSPACE/merged.trimmed.aln.rooted.treefile -m OG -g "hCoV-19/bat/Yunnan/RmYN02/2019|EPI_ISL_412977|2019-06-25"

	# --- metadata MAY CHANGE ----
	pangolin -t $THREADS --outfile $WORKSPACE/merged.lineage_report.csv $WORKSPACE/merged.fas

	awk -F "," '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' $WORKSPACE/merged.lineage_report.csv > $WORKSPACE/merged.metadata.txt

	sed -i '1 a q2:types\tcategorical\tcategorical\tcategorical\tcategorical' $WORKSPACE/merged.metadata.txt

	# note currently adding batch info to metadata on local jupyter notebook... need to improve
	# -------------------------

	# eval "$(conda shell.bash hook)" # necessary to enter conda env from bash script
	# the other way, found by amanda
	source $ANACONDADIR/activate qiime2-2020.11

	empress tree-plot --tree $WORKSPACE/merged.trimmed.aln.rooted.treefile --feature-metadata $WORKSPACE/merged.metadata.txt --output-dir $WORKSPACE/tree-viz

	aws s3 cp $WORKSPACE/ $S3DOWNLOAD/trees/$TIMESTAMP/ --recursive

}

{ time ( buildTree ) ; } > $WORKSPACE/"$TIMESTAMP"-treebuild.log 2>&1

aws s3 cp $WORKSPACE/"$TIMESTAMP"-treebuild.log $S3DOWNLOAD/trees/$TIMESTAMP/