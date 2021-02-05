#!/bin/bash

export PATH=/shared/workspace/software/SD-COVID-Sequencing:/shared/workspace/software/ivar/bin:/shared/workspace/software/anaconda3/bin:$PATH

sample=359044672_S72_L002
workspace=/shared/workspace/projects/covid/data/$sample
reference_fas=/shared/workspace/software/SD-COVID-Sequencing/reference_genome/NC_045512.2.fas
reference_gff=/shared/workspace/software/SD-COVID-Sequencing/reference_genome/NC_045512.2.gff3
primers=/shared/workspace/software/SD-COVID-Sequencing/primers/swift/sarscov2_v2_primers.bed

# Step 1: Map Reads and Sort
# Input: FASTQ (or FASTQ.GZ) file(s) (X.fastq or X.fastq.gz)
# Output: Sorted Untrimmed BAM (X.sorted.bam)
# minimap2 \
# 	-t 32 -a -x map-ont \
# 	"$reference_fas".mmi \
# 	$workspace/"$sample".fastq.gz \
# 	| samtools sort \
# 	-o $workspace/"$sample".sorted.bam

parallel --jobs 32 \
	"{" time "(" minimap2 \
	-t 1 -a -x map-ont \
	"$reference_fas".mmi $workspace/"$sample"*.fastq.gz \
	"|" samtools sort \
	--threads 1 -o $workspace/{}.sorted.bam  \
	")" ";" "}" "2>" $workspace/{}.log.1.map.log ::: $sample

# Step 2: Trim Sorted BAM (resulting in unsorted trimmed BAM)
# Input: Sorted Untrimmed BAM (X.sorted.bam)
# Output: Unsorted Trimmed BAM (X.trimmed.bam)
# ivar trim \
# 	-e \
# 	-i $workspace/"$sample".sorted.bam \
# 	-b $primers \
# 	-p $workspace/"$sample".trimmed

parallel --jobs 32 \
	"{" time "(" ivar trim \
	-e -i $workspace/{}.sorted.bam \
	-b $primers \
	-p $workspace/{}.trimmed \
	")" ";" "}" ">" $workspace/{}.log.2.trim.log "2>&1" ::: $sample

# Step 3: Sort Trimmed BAM
# Input: Unsorted Trimmed BAM (X.trimmed.bam)
# Output: Sorted Trimmed BAM (X.trimmed.sorted.bam)
# samtools sort \
# 	--threads 32 \
# 	-o $workspace/"$sample".trimmed.sorted.bam \
# 	$workspace/"$sample".trimmed.bam #&& rm $workspace/"$sample".trimmed.bam

parallel --jobs 32 \
	"{" time "(" samtools sort \
	--threads 1 \
	-o $workspace/{}.trimmed.sorted.bam {}.trimmed.bam \
	"&&" rm $workspace/{}.trimmed.bam  \
	")" ";" "}" "2>" $workspace/{}.log.3.sorttrimmed.log ::: $sample

# Step 4: Generate Pile-Up from Trimmed Sorted BAM
# Input: Sorted Trimmed BAM (X.trimmed.sorted.bam)
# Output: Pile-up (X.pileup.txt or X.pileup.txt.gz)
# samtools mpileup \
# 	-A -aa -d 0 -Q 0 \
# 	--reference $reference_fas \
# 	$workspace/"$sample".trimmed.sorted.bam \
# 	> $workspace/"$sample".trimmed.sorted.pileup

parallel --jobs 32 \
	"{" time "(" samtools mpileup \
	-A -aa -d 0 -Q 0 \
	--reference $reference_fas \
	$workspace/{}.trimmed.sorted.bam \
	")" ";" "}" ">" $workspace/{}.trimmed.sorted.pileup.txt \
	"2>" $workspace/{}.log.4.pileup.log ::: $sample

# Step 5: Call Variants from Pile-Up
# Input: Pile-up (X.pileup.txt or X.pileup.txt.gz)
# Output: Variants (X.variants.tsv)
parallel --jobs 32 \
	"{" time "(" cat $workspace/{}.trimmed.sorted.pileup.txt \
	"|" ivar variants \
	-r $reference_fas \
	-g $reference_gff \
	-p $workspace/{}.trimmed.sorted.pileup.variants.tsv \
	-m 10 \
	")" ";" "}" "2>" $workspace/{}.log.5.variants.log ::: $sample

# Step 6: Call Consensus Sequence from Pile-Up
# Input: Pile-up (X.pileup.txt or X.pileup.txt.gz)
# Ouptut: Consensus Sequence (X.consensus.fas or X.consensus.fas.gz)
# cat $workspace/"$sample".trimmed.sorted.pileup \
# 	| ivar consensus \
# 	-p $workspace/"$sample".trimmed.sorted.consensus.fas \
# 	-m 10 -n N -t 0.5

parallel --jobs 32 \
	"{" time "(" cat $workspace/{}.trimmed.sorted.pileup.txt \
	"|" ivar consensus \
	-p $workspace/{}.trimmed.sorted.pileup.consensus \
	-m 10 -n N -t 0.5 ")" ";" "}" ">" \
	$workspace/{}.log.6.consensus.log "2>&1" ::: $sample

# Step 7: Call Depth (supplemental summary stats)
# Input: Sorted Trimmed BAM (X.trimmed.sorted.bam)
# Output: Depth (X.trimmed.sorted.depth.txt)
# samtools depth \
# 	-d 0 -Q 0 -q 0 -aa \
# 	$workspace/"$sample".trimmed.sorted.bam \
# 	> $workspace/"$sample".trimmed.sorted.depth.txt

parallel --jobs 32 \
	"{" time "(" samtools depth \
	-d 0 -Q 0 -q 0 -aa \
	$workspace/{}.trimmed.sorted.bam ")" ";" "}" ">" \
	$workspace/{}.trimmed.sorted.depth.txt "2>" \
	$workspace/{}.log.7.depth.log ::: $sample
# Step 8: Qualimap Report (supplemental summary stats)
# Input: Sorted BAM (X.sorted.bam) and Sorted Trimmed BAM (X.trimmed.sorted.bam)
# Output: Qualimap Reports (X.sorted.stats.tar.gz and X.trimmed.sorted.stats.tar.gz)
# qualimap bamqc \
# 	-bam $workspace/"$sample".trimmed.sorted.bam \
# 	-nt 32 \
# 	--java-mem-size=4G \
# 	-outdir $workspace/"$sample".stats && \
# 	tar c $workspace/"$sample".stats \
# 	| pigz -9 -p 32 \
# 	> $workspace/"$sample".stats.tar.gz #&& rm -rf $workspace/"$sample".stats

parallel --jobs 32 \
	"{" time "(" qualimap bamqc \
	-bam $workspace/{1}.{2}.bam \
	-nt 1 --java-mem-size=4G \
	-outdir $workspace/{1}.{2}.stats \
	"&&" tar c $workspace/{1}.{2}.stats \
	"|" pigz -9 -p 1 ">" $workspace/{1}.{2}.stats.tar.gz \
	"&&" rm -rf {1}.{2}.stats ")" ";" "}" ">" \
	$workspace/{1}.log.8.qualimap.{2}.log \
	"2>&1" ::: $sample ::: trimmed.sorted

# Consolidating Stats
python /shared/workspace/software/SD-COVID-Sequencing/qualimap_targz_to_TSV.py $workspace/*.trimmed.sorted.stats.tar.gz > $workspace/YYYY-MM-DD.trimmed.sorted.stats.tsv

# Archiving Each Sample's Files
zip -9 $workspace/"$sample".zip $workspace/"$sample"*
