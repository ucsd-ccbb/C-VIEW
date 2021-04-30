#!/bin/bash

export PATH=/shared/workspace/software/ivar/bin:/shared/workspace/software/q30:/shared/workspace/software/SD-COVID-Sequencing/samhead/:$PATH

# Activate conda env covid1.2
ANACONDADIR=/shared/workspace/software/anaconda3/bin
source $ANACONDADIR/activate covid1.2

# Set variables
THREADS=2
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
REF_FAS="/scratch/reference/NC_045512.2.fas"
REF_MMI="/scratch/reference/NC_045512.2.fas.mmi"
REF_GFF="/scratch/reference/NC_045512.2.gff3"
INSPECT_DELIMITER=__
INTERNAL_DELIMITER=_

FASTQBASE=$SAMPLE
SEQUENCING_INFO=$(echo $SAMPLE | awk -F $INSPECT_DELIMITER '{print $5}')
LANE_INFO=$(echo $SEQUENCING_INFO | awk -F $INTERNAL_DELIMITER '{print $2}' | sed "s/L//g")
SAMPLEID=$(echo $SAMPLE | sed "s/$SEQUENCING_INFO/$LANE_INFO/g")
WORKSPACE=/scratch/$SAMPLEID/$TIMESTAMP
RESULTS=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/"$SEQ_RUN"_samples/$SAMPLEID

# Clear fastq directory if node is being reused
rm -rf $WORKSPACE/*
mkdir -p $WORKSPACE/fastq

echo "$VERSION_INFO" >> $WORKSPACE/"$SAMPLEID".version.log

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

if [[ "$MERGE_LANES" == true ]]; then
  S3DOWNLOAD=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_fastq/"$SEQ_RUN"_lane_merged_fastq
else
  S3DOWNLOAD=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_fastq
fi

# Step 0: Download fastq
# always download read 1
aws s3 cp $S3DOWNLOAD/ $WORKSPACE/fastq/ --recursive --exclude "*" --include "$FASTQBASE*R1_001.fastq.gz"

{ time ( q30.py $WORKSPACE/fastq/"$FASTQBASE"*R1_001.fastq.gz $WORKSPACE/"$SAMPLEID"_R1_q30_reads.txt ) ; } 2> $WORKSPACE/"$SAMPLEID"_R1.log.0.q30.log
echo -e "$SAMPLEID\tq30 R1 exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

if [[ "$FQ" == pe ]]; then
	aws s3 cp $S3DOWNLOAD/ $WORKSPACE/fastq/ --recursive --exclude "*" --include "$FASTQBASE*R2_001.fastq.gz"

  { time ( q30.py $WORKSPACE/fastq/"$FASTQBASE"*R2_001.fastq.gz $WORKSPACE/"$SAMPLEID"_R2_q30_reads.txt ) ; } 2> $WORKSPACE/"$SAMPLEID"_R2.log.0.q30.log
  echo -e "$SAMPLEID\tq30 R2 exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log
fi

# Step 1: Map Reads + Sort
# C/C++ Unsigned long max = 4294967295
if [[ "$READ_CAP" == all ]]; then
  { time ( minimap2 -t $THREADS -a -x sr $REF_MMI $WORKSPACE/fastq/"$FASTQBASE"*.fastq.gz | samtools view -h | samhead 4294967295 successful 2> $WORKSPACE/"$SAMPLEID"_subsampled_mapping_stats.tsv | samtools sort --threads $THREADS -o $WORKSPACE/"$SAMPLEID".sorted.bam ) ; } 2> $WORKSPACE/"$SAMPLEID".log.1.map.log
else
  { time ( minimap2 -t $THREADS -a -x sr $REF_MMI $WORKSPACE/fastq/"$FASTQBASE"*.fastq.gz | samtools view -h | samhead $READ_CAP successful 2> $WORKSPACE/"$SAMPLEID"_subsampled_mapping_stats.tsv | samtools sort --threads $THREADS -o $WORKSPACE/"$SAMPLEID".sorted.bam ) ; } 2> $WORKSPACE/"$SAMPLEID".log.1.map.log
fi
echo -e "$SAMPLEID\tminimap2 exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Step 2: Trim Sorted BAM
{ time ( ivar trim -x 5 -e -i $WORKSPACE/"$SAMPLEID".sorted.bam -b $SCRATCH_PRIMER_FP -p $WORKSPACE/"$SAMPLEID".trimmed ) ; } > $WORKSPACE/"$SAMPLEID".log.2.trim.log 2>&1
echo -e "$SAMPLEID\tivar trim exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Step 3: Sort Trimmed BAM
{ time ( samtools sort --threads $THREADS -o $WORKSPACE/"$SAMPLEID".trimmed.sorted.bam $WORKSPACE/"$SAMPLEID".trimmed.bam && samtools index $WORKSPACE/"$SAMPLEID".trimmed.sorted.bam && rm $WORKSPACE/"$SAMPLEID".trimmed.bam) ; } 2> $WORKSPACE/"$SAMPLEID".log.3.sorttrimmed.log
echo -e "$SAMPLEID\tsamtools sort exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Step 4: Generate Pile-Up
{ time ( samtools mpileup -B -A -aa -d 0 -Q 0 --reference $REF_FAS $WORKSPACE/"$SAMPLEID".trimmed.sorted.bam ) ; } > $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.txt 2> $WORKSPACE/"$SAMPLEID".log.4.pileup.log
echo -e "$SAMPLEID\tsamtools mpileup exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Step 5: Call Variants
{ time ( cat $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.txt | ivar variants -r $REF_FAS -g $REF_GFF -p $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.variants.tsv -m 10 ) ; } 2> $WORKSPACE/"$SAMPLEID".log.5.variants.log
echo -e "$SAMPLEID\tivar variants exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Step 6: Call Consensus
{ time ( cat $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.txt | ivar consensus -p $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.consensus -m 10 -n N -t 0.5 ) ; } > $WORKSPACE/"$SAMPLEID".log.6.consensus.log 2>&1
echo -e "$SAMPLEID\tivar consensus exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Step 7: Call Depth
{ time ( samtools depth -d 0 -Q 0 -q 0 -aa $WORKSPACE/"$SAMPLEID".trimmed.sorted.bam ) ; } > $WORKSPACE/"$SAMPLEID".trimmed.sorted.depth.txt 2> $WORKSPACE/"$SAMPLEID".log.7.depth.log
echo -e "$SAMPLEID\tsamtools depth exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# # Step 8: Qualimap
{ time ( qualimap bamqc -bam $WORKSPACE/"$SAMPLEID".sorted.bam -nt $THREADS --java-mem-size=4G -outdir $WORKSPACE/"$SAMPLEID".sorted.stats ) ; } > $WORKSPACE/"$SAMPLEID".log.8.qualimap.sorted.log 2>&1
echo -e "$SAMPLEID\tqualimap exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Step 9: Acceptance
IVAR_VER=$(ivar version)
{ time ( python $PIPELINEDIR/pipeline/sarscov2_consensus_acceptance.py $SEQ_RUN $TIMESTAMP $FQ "$IVAR_VER" $SAMPLEID $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.consensus.fa $WORKSPACE/"$SAMPLEID".trimmed.sorted.depth.txt $REF_FAS "$SAMPLEID".trimmed.sorted.bam "$SAMPLEID".trimmed.sorted.pileup.variants.tsv $RESULTS $WORKSPACE/"$SAMPLEID".acceptance.tsv $WORKSPACE/"$SAMPLEID".align.json ) ; } 2> $WORKSPACE/"$SAMPLEID".log.9.acceptance.log
echo -e "$SAMPLEID\tacceptance.py exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Step 10: Coverage
{ time ( samtools depth -r NC_045512.2:266-29674 -d 0 -a $WORKSPACE/"$SAMPLEID".trimmed.sorted.bam | \
  awk -v b="$SAMPLEID.trimmed.sorted.bam" 'BEGIN{MIN=10000000000;MAX=0;NUC=0;COV=0;DEPTH=0;NUCZERO=0;}{if(MIN > $3){MIN=$3;};if(MAX < $3){MAX=$3;};if($3==0){NUCZERO+=1};if($3 >= 10){COV+=1;}NUC+=1;DEPTH+=$3;}END{if(NUC>0){print b"\t"(COV/NUC)*100"\t"DEPTH/NUC"\t"MIN"\t"MAX"\t"NUCZERO}else{print b"\t"0"\t"0"\t"MIN"\t"MAX"\t"NUCZERO}}' ) ; } >> $WORKSPACE/"$SAMPLEID".coverage.tsv 2> $WORKSPACE/"$SAMPLEID".log.10.coverage.log
echo -e "$SAMPLEID\tcoverage exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

#QC
grep -v "exit code: 0" $WORKSPACE/"$SAMPLEID".exit.log | head -n 1 > $WORKSPACE/"$SAMPLEID".error.log

aws s3 cp $WORKSPACE/ $RESULTS/ --recursive --include "*" --exclude "*fastq.gz"
