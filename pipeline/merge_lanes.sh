mkdir -p $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq

aws s3 cp $S3DOWNLOAD/ $WORKSPACE/ --recursive --exclude "*" --include "*fastq.gz"

INSPECT_DELIMITER=__
SAMPLES_LIST=()
for R1_FASTQ in $(ls $WORKSPACE | grep R1_001.fastq.gz | sort | uniq | grep -v Undetermined); do
	SEQUENCING_INFO=$(echo $R1_FASTQ | awk -F $INSPECT_DELIMITER '{print $NF}')
	SAMPLE_NUM=$(echo $SEQUENCING_INFO | awk -F _ '{print $1}')
	A_SAMPLE=$(echo $R1_FASTQ | sed "s/$SEQUENCING_INFO/$SAMPLE_NUM/g")
	SAMPLES_LIST+=($A_SAMPLE)
done

for SAMPLE in $(printf '%s\n' "${SAMPLES_LIST[@]}" | sort | uniq ); do
	LANES=$(ls $WORKSPACE/$SAMPLE* | awk -F $INSPECT_DELIMITER '{print $NF}'| awk -F '_L|_R' '{print $2}' | sort | uniq | grep 00)
	LANES_COMBINED=$(echo $LANES | sed 's/ //g')
	cat $WORKSPACE/"$SAMPLE"*R1_001.fastq.gz > $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq/"$SAMPLE"_"$LANES_COMBINED"_R1_001.fastq.gz
	cat $WORKSPACE/"$SAMPLE"*R2_001.fastq.gz > $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq/"$SAMPLE"_"$LANES_COMBINED"_R2_001.fastq.gz
done

aws s3 cp $WORKSPACE/"$SEQ_RUN"_lane_merged_fastq/ $S3DOWNLOAD/"$SEQ_RUN"_lane_merged_fastq/ --recursive --exclude "*" --include "*fastq.gz"
rm -rf $WORKSPACE