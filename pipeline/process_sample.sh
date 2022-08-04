#!/bin/bash

PATH=/shared/workspace/software:/shared/workspace/software/q30:/shared/workspace/software/samhead/:$PATH

# Activate conda env cview
ANACONDADIR=/shared/workspace/software/anaconda3/bin
source $ANACONDADIR/activate cview

# Set variables
THREADS=2
CVIEWDIR=/shared/workspace/software/cview
INSPECT_DELIMITER=__
INTERNAL_DELIMITER=_

if [[ "$INPUT_TYPE" == fastq ]]; then
  SEQUENCING_INFO=$(echo $SAMPLE | awk -F $INSPECT_DELIMITER '{print $5}')
  LANE_INFO=$(echo $SEQUENCING_INFO | awk -F $INTERNAL_DELIMITER '{print $2}' | sed "s/L//g")
  SAMPLEID=$(echo $SAMPLE | sed "s/$SEQUENCING_INFO/$LANE_INFO/g")
  REF_NAME="NC_045512.2"
  REF_ORF_LIMITS="266-29674"
fi

if [[ "$INPUT_TYPE" == bam ]]; then
  SAMPLEID=$SAMPLE
  # Thermo rep said: "The reference file used for mapping is ... the same sequence as the original
  # MN908947.3 SARS sequence, we named it as 2019-nCoV."
  # (Note that 2019-nCoV.fas is the exact same sequence as NC_045512.2.fas ...)
  REF_NAME="2019-nCoV"
  REF_ORF_LIMITS="266-29674"  # not a typo, really identical to the limits for se/pe--genexus uses a different label for its reference sequence but the sequence itself is identical
fi

REF_FAS="/scratch/reference/"$REF_NAME".fas"
REF_MMI="/scratch/reference/$REF_NAME.mmi"
REF_GFF="/scratch/reference/$REF_NAME.gff3"

WORKSPACE=/scratch/$SAMPLEID/$TIMESTAMP
RESULTS=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_results/"$TIMESTAMP"_"$FQ"/"$SEQ_RUN"_samples/$SAMPLEID
# NB: do NOT put this line before the previous RESULTS= line, since this line is *redefining* S3DOWNLOAD
S3DOWNLOAD=$S3DOWNLOAD/$SEQ_RUN/"$SEQ_RUN"_$INPUT_TYPE

# Clear input data directory if node is being reused
rm -rf $WORKSPACE/*

echo "$VERSION_INFO" >> $WORKSPACE/"$SAMPLEID".version.log

# Move reference files to compute node once
if [[ ! -f "$REF_FAS" ]]; then
	mkdir -p /scratch/reference/
    cp $CVIEWDIR/reference_files/$REF_NAME.fas $REF_FAS
    cp $CVIEWDIR/reference_files/$REF_NAME.gff3 $REF_GFF
fi

if [[ "$INPUT_TYPE" == fastq ]]; then
  # make working directory to hold the fastqs
  mkdir -p $WORKSPACE/fastq

  # will also need an MMI file for the reference; get it now
  cp $CVIEWDIR/reference_files/$REF_NAME.fas.mmi $REF_MMI

  # ensure that primer file is downloaded
  SCRATCH_PRIMER_FP=/scratch/reference/$PRIMER_BED_FNAME
  if [[ ! -f "$SCRATCH_PRIMER_FP" ]]; then
    cp $CVIEWDIR/reference_files/$PRIMER_BED_FNAME $SCRATCH_PRIMER_FP
  fi

  if [[ "$MERGE_LANES" == true ]]; then
    # NB: Once again renaming S3DOWNLOAD in this situation
    S3DOWNLOAD=$S3DOWNLOAD/"$SEQ_RUN"_lane_merged_fastq
  fi

  # Step 0: Download input data
  { time ( aws s3 cp $S3DOWNLOAD/ $WORKSPACE/fastq/ --recursive --exclude "*" --include "$SAMPLE*R1_001$INPUT_SUFFIX" ) ; } 2> $WORKSPACE/"$SAMPLEID"_R1.log.0.download.log
  echo -e "$SAMPLEID\tdownload R1 exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

  { time ( q30.py $WORKSPACE/fastq/"$SAMPLE"*R1_001$INPUT_SUFFIX $WORKSPACE/"$SAMPLEID"_R1_q30_reads.txt ) ; } 2>> $WORKSPACE/"$SAMPLEID"_R1.log.0a.q30.log
  echo -e "$SAMPLEID\tq30 R1 exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

  if [[ "$FQ" == pe ]]; then
    { time ( aws s3 cp $S3DOWNLOAD/ $WORKSPACE/fastq/ --recursive --exclude "*" --include "$SAMPLE*R2_001$INPUT_SUFFIX" ) ; } 2> $WORKSPACE/"$SAMPLEID"_R2.log.0.download.log
    echo -e "$SAMPLEID\tdownload R2 exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

    { time ( q30.py $WORKSPACE/fastq/"$SAMPLE"*R2_001$INPUT_SUFFIX $WORKSPACE/"$SAMPLEID"_R2_q30_reads.txt ) ; } 2>> $WORKSPACE/"$SAMPLEID"_R2.log.0a.q30.log
    echo -e "$SAMPLEID\tq30 R2 exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log
  fi

  # Step 1: Map Reads + Sort
  { time ( minimap2 -t $THREADS -a -x sr $REF_MMI $WORKSPACE/fastq/"$SAMPLE"*$INPUT_SUFFIX | samtools view -h | samhead $READ_CAP successful 2> $WORKSPACE/"$SAMPLEID"_subsampled_mapping_stats.tsv | samtools sort --threads $THREADS -o $WORKSPACE/"$SAMPLEID".sorted.bam) ; } 2> $WORKSPACE/"$SAMPLEID".log.1.map.log
  echo -e "$SAMPLEID\tminimap2 exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

  source $ANACONDADIR/deactivate
  source $ANACONDADIR/activate ivar

  # Step 2: Trim Sorted BAM
  { time ( ivar trim -x 5 -e -i $WORKSPACE/"$SAMPLEID".sorted.bam -b $SCRATCH_PRIMER_FP -p $WORKSPACE/"$SAMPLEID".trimmed ) ; } > $WORKSPACE/"$SAMPLEID".log.2.trim.log 2>&1
  echo -e "$SAMPLEID\tivar trim exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log
  QUALIMAP_BAM=$WORKSPACE/"$SAMPLEID".sorted.bam

  source $ANACONDADIR/deactivate
  source $ANACONDADIR/activate cview

  # Step 3: Sort Trimmed BAM
  { time ( samtools sort --threads $THREADS -o $WORKSPACE/"$SAMPLEID".trimmed.sorted.bam $WORKSPACE/"$SAMPLEID".trimmed.bam && samtools index $WORKSPACE/"$SAMPLEID".trimmed.sorted.bam && rm $WORKSPACE/"$SAMPLEID".trimmed.bam) ; } 2> $WORKSPACE/"$SAMPLEID".log.3.sorttrimmed.log
  echo -e "$SAMPLEID\tsamtools sort exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log
fi # end if the input is fastq

if [[ "$INPUT_TYPE" == bam ]]; then
  # if genexus bam, download just one file and rename it to have the same naming convention as
  # trimmed sorted bam produced by the above step
  TBAM=$WORKSPACE/"$SAMPLE$INPUT_SUFFIX"
  { time ( aws s3 cp $S3DOWNLOAD/"$SAMPLE$INPUT_SUFFIX" $TBAM ) ; } 2> $WORKSPACE/"$SAMPLEID".log.0.download.log
  echo -e "$SAMPLEID\tdownload bam exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log


  # filter out bogus empty records in the genexus bam for other "chromosomes"--other reference sequences
  # that they align everything against
  { time ( samtools reheader <(samtools view -H $TBAM | grep -P "^@HD\t" && samtools view -H $TBAM | grep -P "^@SQ\tSN:$REF_NAME\t") $TBAM > $WORKSPACE/"$SAMPLEID.trimmed.sorted.bam") ; } 2> $WORKSPACE/"$SAMPLEID".log.3.sorttrimmed.log
  # Note >> here--adding to log made above, not making new one
  { time ( samtools index $WORKSPACE/"$SAMPLEID.trimmed.sorted.bam") ; } 2>> $WORKSPACE/"$SAMPLEID".log.3.sorttrimmed.log
  QUALIMAP_BAM=$WORKSPACE/"$SAMPLEID.trimmed.sorted.bam"
fi # end if the input is a genexus pre-processed bam

# Step 3a: Count number of trimmed bam reads
{ time ( echo -e "$SAMPLEID\t"$(samtools view -c -F 260 $WORKSPACE/"$SAMPLEID".trimmed.sorted.bam ) > $WORKSPACE/"$SAMPLEID".trimmed_bam_read_count.tsv) ; } >> $WORKSPACE/"$SAMPLEID".log.2.trim.log 2>&1
echo -e "$SAMPLEID\ttrimmed bam read count exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Step 4: Generate Pile-Up
{ time ( samtools mpileup -B -A -aa -d 0 -Q 0 --reference $REF_FAS $WORKSPACE/"$SAMPLEID".trimmed.sorted.bam ) ; } > $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.txt 2> $WORKSPACE/"$SAMPLEID".log.4.pileup.log
echo -e "$SAMPLEID\tsamtools mpileup exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

source $ANACONDADIR/deactivate
source $ANACONDADIR/activate ivar
IVAR_VER=$(ivar version)

# Step 5: Call Variants
{ time ( cat $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.txt | ivar variants -r $REF_FAS -g $REF_GFF -p $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.variants.tsv -m 10 ) ; } 2> $WORKSPACE/"$SAMPLEID".log.5.variants.log
echo -e "$SAMPLEID\tivar variants exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Step 6: Call Consensus
{ time ( cat $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.txt | ivar consensus -p $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.consensus -m 10 -n N -t 0.5 ) ; } > $WORKSPACE/"$SAMPLEID".log.6.consensus.log 2>&1
echo -e "$SAMPLEID\tivar consensus exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

source $ANACONDADIR/deactivate
source $ANACONDADIR/activate cview

# Step 7: Call Depth
{ time ( samtools depth -J -d 0 -Q 0 -q 0 -aa $WORKSPACE/"$SAMPLEID".trimmed.sorted.bam ) ; } > $WORKSPACE/"$SAMPLEID".trimmed.sorted.depth.txt 2> $WORKSPACE/"$SAMPLEID".log.7.depth.log
echo -e "$SAMPLEID\tsamtools depth exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# # Step 8: Qualimap
{ time ( qualimap bamqc -bam $QUALIMAP_BAM -nt $THREADS --java-mem-size=4G -outdir $WORKSPACE/"$SAMPLEID".sorted.stats ) ; } > $WORKSPACE/"$SAMPLEID".log.8.qualimap.sorted.log 2>&1
echo -e "$SAMPLEID\tqualimap exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Step 9: Acceptance
{ time ( python $CVIEWDIR/src/sarscov2_consensus_acceptance.py $SEQ_RUN $TIMESTAMP $FQ "$IVAR_VER" $SAMPLEID $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.consensus.fa $WORKSPACE/"$SAMPLEID".trimmed.sorted.depth.txt $REF_FAS "$SAMPLEID".trimmed.sorted.bam "$SAMPLEID".trimmed.sorted.pileup.variants.tsv $RESULTS $WORKSPACE/"$SAMPLEID".acceptance.tsv $WORKSPACE/"$SAMPLEID".align.json ) ; } 2> $WORKSPACE/"$SAMPLEID".log.9.acceptance.log
echo -e "$SAMPLEID\tacceptance.py exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Step 10: Coverage
{ time ( samtools depth -r $REF_NAME:$REF_ORF_LIMITS -d 0 -a $WORKSPACE/"$SAMPLEID".trimmed.sorted.bam | \
  awk -v b="$SAMPLEID.trimmed.sorted.bam" 'BEGIN{MIN=10000000000;MAX=0;NUC=0;COV=0;DEPTH=0;NUCZERO=0;}{if(MIN > $3){MIN=$3;};if(MAX < $3){MAX=$3;};if($3==0){NUCZERO+=1};if($3 >= 10){COV+=1;}NUC+=1;DEPTH+=$3;}END{if(NUC>0){print b"\t"(COV/NUC)*100"\t"DEPTH/NUC"\t"MIN"\t"MAX"\t"NUCZERO}else{print b"\t"0"\t"0"\t"MIN"\t"MAX"\t"NUCZERO}}' ) ; } >> $WORKSPACE/"$SAMPLEID".coverage.tsv 2> $WORKSPACE/"$SAMPLEID".log.10.coverage.log
echo -e "$SAMPLEID\tcoverage exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Step 11: Heterogeneity scores
{ time (
pi_from_pileup $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.txt > $WORKSPACE/"$SAMPLEID".pi_from_pileup.tsv
echo -e "$SAMPLEID\t"$(cat $WORKSPACE/"$SAMPLEID".pi_from_pileup.tsv | tail -1 | cut -f3) > $WORKSPACE/"$SAMPLEID".pi_metric.tsv
echo -e "$SAMPLEID\tpi metric exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log
echo -e "$SAMPLEID\t"$(cat $WORKSPACE/"$SAMPLEID".trimmed.sorted.pileup.variants.tsv | awk '$11 != "ALT_FREQ" && $11 >= 0.05 && $11 <= 0.95' | grep -v "ALT_FREQ" | cut -f2,4 | sort | uniq | wc -l) >> $WORKSPACE/"$SAMPLEID".n_metric.tsv ) ;
} > $WORKSPACE/"$SAMPLEID".log.11.diversity.log 2>&1
echo -e "$SAMPLEID\tn metric exit code: $?" >> $WORKSPACE/"$SAMPLEID".exit.log

# Collect failure exit codes to error.log
# Note that the 255 exit code from qualimap is ignored: it is raised when the
# bam file has zero reads in it, which is a real case that is handled gracefully
# by the rest of the pipeline and is therefore not cause for failing a run
grep -v 'exit code: 0\|qualimap exit code: 255' $WORKSPACE/"$SAMPLEID".exit.log | head -n 1 > $WORKSPACE/"$SAMPLEID".error.log

aws s3 cp $WORKSPACE/ $RESULTS/ --recursive --include "*" --exclude "*$INPUT_SUFFIX"
