#!/bin/bash

S3DOWNLOAD=$1

for SAMPLE in $(aws s3 ls $S3DOWNLOAD/  --recursive | grep fastq | awk '{print $4}' | awk -F '/|_R' '{print $4}' | sort | uniq | grep -v Undetermined); do
	qsub -v SAMPLE=$SAMPLE \
		 -v S3DOWNLOAD=$S3DOWNLOAD \
		 -N $SAMPLE \
		 -wd /shared/workspace/projects/covid/logs \
		 -pe smp 1 \
		 -S /bin/bash \
		 /shared/workspace/software/covid_sequencing_analysis_pipeline/covid.sh
done

# When all samples have finished
# qsub -v S3DOWNLOAD=$S3DOWNLOAD/results \
# 	/shared/workspace/software/covid_sequencing_analysis_pipeline/qc_summary.sh