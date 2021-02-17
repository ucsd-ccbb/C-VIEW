#!/bin/bash

S3DOWNLOAD=$1 # e.g. s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/
PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline
WORKSPACE=/scratch/$(date +"%T" | sed 's/:/-/g')

# When all samples have finished
qsub -v S3DOWNLOAD=$S3DOWNLOAD \
	 -v WORKSPACE=$WORKSPACE \
	 -wd /shared/workspace/projects/covid/logs \
	 -pe smp 1 \
	 -S /bin/bash \
     $PIPELINEDIR/qc/qc_summary.sh
