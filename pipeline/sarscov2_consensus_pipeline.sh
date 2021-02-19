#!/bin/bash

export PATH=/shared/workspace/software/ivar/bin:/shared/workspace/software/anaconda3/envs/covid1.1/bin:$PATH

# Set variables
THREADS=1
WORKSPACE=/scratch/$SAMPLE
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
REF_FAS="/scratch/reference/NC_045512.2.fas"
REF_MMI="/scratch/reference/NC_045512.2.fas.mmi"
REF_GFF="/scratch/reference/NC_045512.2.gff3"
# Clear fastq directory if node is being reused
rm -rf $WORKSPACE/*
mkdir -p $WORKSPACE/fastq
mkdir -p $WORKSPACE/fastqc

# Move reference files to compute node once
if [[ ! -f "$REF_FAS" ]]; then
	mkdir -p /scratch/reference/
    cp $PIPELINEDIR/reference_files/NC_045512.2.fas $REF_FAS
    cp $PIPELINEDIR/reference_files/NC_045512.2.fas.mmi $REF_MMI
    cp $PIPELINEDIR/reference_files/NC_045512.2.gff3 $REF_GFF
fi

# Determine which primer bed file to use and ensure it is downloaded
if [[ "$PRIMER_SET" == artic ]]; then
	PRIMER_BED_FNAME="nCoV-2019.primer.bed"
fi

if [[ "$PRIMER_SET" == swift_v2 ]]; then
	PRIMER_BED_FNAME="sarscov2_v2_primers.bed"
fi

SCRATCH_PRIMER_FP=/scratch/reference/$PRIMER_BED_FNAME
if [[ ! -f "$SCRATCH_PRIMER_FP" ]]; then
  cp $PIPELINEDIR/reference_files/$PRIMER_BED_FNAME $SCRATCH_PRIMER_FP
fi

# Step 0: Download fastq
RESULTS="$TIMESTAMP"_"$FQ"
if [[ "$FQ" == se ]]; then
  # leave out any files for this sample that
  # contain in their filename the string "R2" anywhere
  # after the sample name but before the "fastq.gz"
	aws s3 cp $S3DOWNLOAD/ $WORKSPACE/fastq/ --recursive --exclude "*" \
	  --include "$SAMPLE*fastq.gz" --exclude "$SAMPLE*R2*fastq.gz"
fi

if [[ "$FQ" == pe ]]; then
	aws s3 cp $S3DOWNLOAD/"$SAMPLE*R1*fastq.gz" $WORKSPACE/fastq/
	aws s3 cp $S3DOWNLOAD/"$SAMPLE*R2*fastq.gz" $WORKSPACE/fastq/
fi

# Fastqc
fastqc $WORKSPACE/fastq/"$SAMPLE"*fastq.gz -o $WORKSPACE/fastqc

# Step 1: Map Reads + Sort
{ time ( minimap2 -t $THREADS -a -x sr $REF_MMI $WORKSPACE/fastq/"$SAMPLE"*.fastq.gz | samtools sort --threads $THREADS -o $WORKSPACE/"$SAMPLE".sorted.bam ) ; } 2> $WORKSPACE/"$SAMPLE".log.1.map.log

# Step 2: Trim Sorted BAM
{ time ( ivar trim -x 5 -e -i $WORKSPACE/"$SAMPLE".sorted.bam -b $PRIMER_BED -p $WORKSPACE/"$SAMPLE".trimmed ) ; } > $WORKSPACE/"$SAMPLE".log.2.trim.log 2>&1

# Step 3: Sort Trimmed BAM
{ time ( samtools sort --threads $THREADS -o $WORKSPACE/"$SAMPLE".trimmed.sorted.bam $WORKSPACE/"$SAMPLE".trimmed.bam && rm $WORKSPACE/"$SAMPLE".trimmed.bam ) ; } 2> $WORKSPACE/"$SAMPLE".log.3.sorttrimmed.log

# Step 4: Generate Pile-Up
{ time ( samtools mpileup -A -aa -d 0 -Q 0 --reference $REF_FAS $WORKSPACE/"$SAMPLE".trimmed.sorted.bam ) ; } > $WORKSPACE/"$SAMPLE".trimmed.sorted.pileup.txt 2> $WORKSPACE/"$SAMPLE".log.4.pileup.log

# Step 5: Call Variants
{ time ( cat $WORKSPACE/"$SAMPLE".trimmed.sorted.pileup.txt | ivar variants -r $REF_FAS -g $REF_GFF -p $WORKSPACE/"$SAMPLE".trimmed.sorted.pileup.variants.tsv -m 10 ) ; } 2> $WORKSPACE/"$SAMPLE".log.5.variants.log

# Step 6: Call Consensus
{ time ( cat $WORKSPACE/"$SAMPLE".trimmed.sorted.pileup.txt | ivar consensus -p $WORKSPACE/"$SAMPLE".trimmed.sorted.pileup.consensus -m 10 -n N -t 0.5 ) ; } > $WORKSPACE/"$SAMPLE".log.6.consensus.log 2>&1

# Step 7: Call Depth
{ time ( samtools depth -d 0 -Q 0 -q 0 -aa $WORKSPACE/"$SAMPLE".trimmed.sorted.bam ) ; } > $WORKSPACE/"$SAMPLE".trimmed.sorted.depth.txt 2> $WORKSPACE/"$SAMPLE".log.7.depth.log

# # Step 8: Qualimap
{ time ( qualimap bamqc -bam $WORKSPACE/"$SAMPLE".sorted.bam -nt $THREADS --java-mem-size=4G -outdir $WORKSPACE/"$SAMPLE".sorted.stats ) ; } > $WORKSPACE/"$SAMPLE".log.8.qualimap.sorted.log 2>&1

# QC
python $PIPELINEDIR/pipeline/sarscov2_consensus_acceptance.py $WORKSPACE/"$SAMPLE".trimmed.sorted.pileup.consensus.fa $WORKSPACE/"$SAMPLE".trimmed.sorted.depth.txt $REF_FAS

aws s3 cp $WORKSPACE/ $S3DOWNLOAD/$RESULTS/$SAMPLE/ --recursive --include "*" --exclude "*fastq.gz"
