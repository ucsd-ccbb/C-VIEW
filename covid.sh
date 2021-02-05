#!/bin/bash

# Step 1: Map Reads and Sort
# Input: FASTQ (or FASTQ.GZ) file(s) (X.fastq or X.fastq.gz)
# Output: Sorted Untrimmed BAM (X.sorted.bam)
minimap2 -t THREADS -a -x map-ont ../ref/NC_045512.2.fas.mmi READ1.FASTQ.GZ READ2.FASTQ.GZ | samtools sort --threads THREADS -o SORTED.BAM

# Step 2: Trim Sorted BAM (resulting in unsorted trimmed BAM)
# Input: Sorted Untrimmed BAM (X.sorted.bam)
# Output: Unsorted Trimmed BAM (X.trimmed.bam)
ivar trim -e -i SORTED.BAM -b PRIMERS.bed -p TRIMMED_PREFIX

# Step 3: Sort Trimmed BAM
# Input: Unsorted Trimmed BAM (X.trimmed.bam)
# Output: Sorted Trimmed BAM (X.trimmed.sorted.bam)
samtools sort --threads THREADS -o TRIMMED_SORTED.BAM TRIMMED_PREFIX.BAM #&& rm TRIMMED_PREFIX.BAM

# Step 4: Generate Pile-Up from Trimmed Sorted BAM
# Input: Sorted Trimmed BAM (X.trimmed.sorted.bam)
# Output: Pile-up (X.pileup.txt or X.pileup.txt.gz)
samtools mpileup -A -aa -d 0 -Q 0 --reference REFERENCE.FAS TRIMMED_SORTED.BAM > PILEUP.TXT

# Step 5: Call Variants from Pile-Up
# Input: Pile-up (X.pileup.txt or X.pileup.txt.gz)
# Output: Variants (X.variants.tsv)
cat PILEUP.TXT | ivar variants -r REFERENCE.FAS -g REFERENCE.GFF -p VARIANTS.TSV -m 10

# Step 6: Call Consensus Sequence from Pile-Up
# Input: Pile-up (X.pileup.txt or X.pileup.txt.gz)
# Ouptut: Consensus Sequence (X.consensus.fas or X.consensus.fas.gz)
cat PILEUP.TXT | ivar consensus -p CONSENSUS.FAS -m 10 -n N -t 0.5

# Step 7: Call Depth (supplemental summary stats)
# Input: Sorted Trimmed BAM (X.trimmed.sorted.bam)
# Output: Depth (X.trimmed.sorted.depth.txt)
samtools depth -d 0 -Q 0 -q 0 -aa TRIMMED_SORTED.BAM > DEPTH.TXT

# Step 8: Qualimap Report (supplemental summary stats)
# Input: Sorted BAM (X.sorted.bam) and Sorted Trimmed BAM (X.trimmed.sorted.bam)
# Output: Qualimap Reports (X.sorted.stats.tar.gz and X.trimmed.sorted.stats.tar.gz)
qualimap bamqc -bam MAPPING.BAM -nt THREADS --java-mem-size=RAM -outdir STATS_DIR && tar c STATS_DIR | pigz -9 -p THREADS > STATS_DIR.tar.gz && rm -rf STATS_DIR
# Consolidating Stats
#qualimap_targz_to_TSV.py *.trimmed.sorted.stats.tar.gz > YYYY-MM-DD.trimmed.sorted.stats.tsv

# Archiving Each Sample's Files
zip -9 SAMPLE.zip SAMPLE*
