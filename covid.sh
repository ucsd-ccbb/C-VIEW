#!/bin/bash

export PATH=/shared/workspace/software/SD-COVID-Sequencing:/shared/workspace/software/ivar/bin:/shared/workspace/software/anaconda3/bin:$PATH

THREADS=1
REF_FAS="/scratch/reference/NC_045512.2.fas"
REF_MMI="/scratch/reference/NC_045512.2.fas.mmi"
REF_GFF="/scratch/reference/NC_045512.2.gff3"
PRIMER_BED="/scratch/reference/sarscov2_v2_primers.bed"
WORKSPACE=/scratch/$SAMPLE
mkdir -p $WORKSPACE

if [[ ! -f "$REF_FAS" ]]; then
	mkdir -p /scratch/reference/
    cp /shared/workspace/software/SD-COVID-Sequencing/reference_genome/NC_045512.2.fas $REF_FAS
    cp /shared/workspace/software/SD-COVID-Sequencing/reference_genome/NC_045512.2.fas.mmi $REF_MMI
    cp /shared/workspace/software/SD-COVID-Sequencing/reference_genome/NC_045512.2.gff3 $REF_GFF
    cp /shared/workspace/software/SD-COVID-Sequencing/primers/swift/sarscov2_v2_primers.bed $PRIMER_BED
fi

# Step 0: Download fastq
aws s3 cp $S3DOWNLOAD/"$SAMPLE"_R1_001.fastq.gz $WORKSPACE/
aws s3 cp $S3DOWNLOAD/"$SAMPLE"_R2_001.fastq.gz $WORKSPACE/


# Step 1: Map Reads + Sort
{ time ( minimap2 -t $THREADS -a -x sr $REF_MMI $WORKSPACE/"$SAMPLE"*.fastq.gz | samtools sort --threads $THREADS -o $WORKSPACE/"$SAMPLE".sorted.bam ) ; } 2> $WORKSPACE/"$SAMPLE".log.1.map.log

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
for x in sorted trimmed.sorted ; do
    { time ( qualimap bamqc -bam $WORKSPACE/"$SAMPLE".$x.bam -nt $THREADS --java-mem-size=4G -outdir $WORKSPACE/"$SAMPLE".$x.stats && tar c $WORKSPACE/"$SAMPLE".$x.stats | pigz -9 -p $THREADS > $WORKSPACE/"$SAMPLE".$x.stats.tar.gz && rm -rf $WORKSPACE/"$SAMPLE".$x.stats ) ; } > $WORKSPACE/"$SAMPLE".log.8.qualimap.$x.log 2>&1
done

# Step 9: Zip
cd $WORKSPACE && zip -9 "$SAMPLE".zip "$SAMPLE"*

aws s3 cp $WORKSPACE/"$SAMPLE".zip $S3DOWNLOAD/results/zip/
aws s3 cp $WORKSPACE/"$SAMPLE".trimmed.sorted.pileup.variants.tsv $S3DOWNLOAD/results/variants/
aws s3 cp $WORKSPACE/"$SAMPLE".trimmed.sorted.pileup.consensus.fa $S3DOWNLOAD/results/consensus/
aws s3 cp $WORKSPACE/"$SAMPLE".trimmed.sorted.depth.txt $S3DOWNLOAD/results/depth/
aws s3 cp $WORKSPACE/"$SAMPLE".sorted.stats.tar.gz $S3DOWNLOAD/results/stats/
aws s3 cp $WORKSPACE/"$SAMPLE".trimmed.sorted.stats.tar.gz $S3DOWNLOAD/results/stats/
