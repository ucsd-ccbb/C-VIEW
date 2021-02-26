
export PATH=$PATH:/shared/workspace/software/IQTree/iqtree-2.1.2-Linux/bin:/shared/workspace/software/viralMSA:/shared/workspace/software/MinVar-Rooting-master:/shared/workspace/software/anaconda3/envs/pangolin/bin
mkdir -p $WORKSPACE
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
ANACONDADIR=/shared/workspace/software/anaconda3/bin
THREADS=96

buildTree () {

	aws s3 cp $S3DOWNLOAD/consensus/ $WORKSPACE/ --recursive
	aws s3 cp $S3DOWNLOAD/acceptance/ $WORKSPACE/ --recursive

	cat $WORKSPACE/*passQC.fas > $WORKSPACE/merged.fas
	sed -i -e 's/Consensus_//g' -e 's/.trimmed.sorted.pileup.consensus_threshold_0.5_quality_20//g' $WORKSPACE/merged.fas
	# add the reference sequence
	cat $PIPELINEDIR/reference_files/RmYN02.fas >> $WORKSPACE/merged.fas

	# add B.1.1.7 sequence
	cat $PIPELINEDIR/reference_files/hCoV-19_USA_CA-SEARCH-5574_2020.fasta >> $WORKSPACE/merged.fas

	source $ANACONDADIR/activate covid1.1
	ViralMSA.py -s $WORKSPACE/merged.fas -r SARS-CoV-2 -o viralmsa_out -t $THREADS -e aws-CCBB@health.ucsd.edu

	python $PIPELINEDIR/pipeline/trim_msa.py -i $WORKSPACE/viralmsa_out/merged.fas.aln -s 100 -e 50 -o $WORKSPACE/merged.trimmed.aln

	iqtree2 -T $THREADS -m GTR+F+G4 --polytomy -blmin 1e-9 -s $WORKSPACE/merged.trimmed.aln

	python FastRoot.py -i $WORKSPACE/merged.trimmed.aln.treefile -o $WORKSPACE/merged.trimmed.aln.rooted.treefile -m OG -g "hCoV-19/bat/Yunnan/RmYN02/2019|EPI_ISL_412977|2019-06-25"

	# pangolin
	source $ANACONDADIR/activate pangolin
	pangolin --update
	pangolin -t $THREADS --outfile $WORKSPACE/merged.lineage_report.csv $WORKSPACE/merged.fas

	# Metadata -------------------------
	# Get SEQ_RUN from acceptance files to merge with metadata
	for seq_run in $(ls $WORKSPACE/*acceptance.tsv | awk -F '/' '{print $NF}' | awk -F '-summary' '{print $1}'); do 
		grep -v "^fastq_id" $WORKSPACE/"$seq_run"-summary.acceptance.tsv \
		| awk -v seq_run=$seq_run \
		'OFS="\t"{print $1, seq_run}' \
		> $WORKSPACE/"$seq_run"-metadata.txt
	done
	cat $WORKSPACE/*-metadata.txt > $WORKSPACE/tmp.merged.metadata.txt
	echo -e "hCoV-19/bat/Yunnan/RmYN02/2019|EPI_ISL_412977|2019-06-25\thCoV-19/bat/Yunnan/RmYN02/2019|EPI_ISL_412977|2019-06-25" >> $WORKSPACE/tmp.merged.metadata.txt
	echo -e "hCoV-19/USA/CA-SEARCH-5574/2020|EPI_ISL_751801|2020-12-29\thCoV-19/USA/CA-SEARCH-5574/2020|EPI_ISL_751801|2020-12-29" >> $WORKSPACE/tmp.merged.metadata.txt

	# Merge lineage report with sample/seqrun file to make final metadata
	join <(awk 'BEGIN {FS=",";OFS="\t"} {print $1,$2,$3,$4,$5}' $WORKSPACE/merged.lineage_report.csv | sort -k1) <(sort -k1 $WORKSPACE/tmp.merged.metadata.txt) -t $'\t'  > $WORKSPACE/final.metadata.txt
	sed -i '1i taxon\tlineage\tprobability\tpangoLEARN_version\run_name' $WORKSPACE/final.metadata.txt
	sed -i '1 a q2:types\tcategorical\tcategorical\tcategorical\tcategorical\tcategorical' $WORKSPACE/final.metadata.txt
	# -------------------------

	source $ANACONDADIR/activate qiime2-2020.11

	empress tree-plot --tree $WORKSPACE/merged.trimmed.aln.rooted.treefile --feature-metadata $WORKSPACE/final.metadata.txt --output-dir $WORKSPACE/tree-viz

	aws s3 cp $WORKSPACE/ $S3DOWNLOAD/trees/$TIMESTAMP/ --recursive

}

{ time ( buildTree ) ; } > $WORKSPACE/"$TIMESTAMP"-treebuild.log 2>&1

aws s3 cp $WORKSPACE/"$TIMESTAMP"-treebuild.log $S3DOWNLOAD/trees/$TIMESTAMP/

