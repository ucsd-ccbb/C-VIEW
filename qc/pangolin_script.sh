#!/bin/bash

fname_prefix=$1

s3_consensus_loc=$2
s3_acceptance_loc=$3

aws s3 cp $s3_consensus_loc ./
aws s3 cp $s3_acceptance_loc ./

echo $fname_prefix

unzip consensus.zip

# filter the true samples
awk '{ if ($2 == "True") { print } }' summary.acceptance.tsv > summary.acceptance.true.tsv
awk '{print $1}' summary.acceptance.true.tsv > passQC.samples.tsv

# loop over individual .fa files, keep the ones which are in passQC.samples.tsv
touch $fname_prefix.fas # initialize the file
for f in *.fa; do 
    fshort="$(cut -d'.' -f1 <<<$f)"
    echo $fshort
    
    if (grep -qF $fshort passQC.samples.tsv); then
       #echo "Found it"
       cat $f >> $fname_prefix.fas
    fi

done

# print out file lengths for sanity check
echo $str | wc -l summary.acceptance.tsv
echo $str | wc -l passQC.samples.tsv
echo $str | wc -l $fname_prefix.fas


# note: should always update pangolin to get the most recent variant classifications... but this updates code as well as data... may be problematic
# pangolin --update
pangolin -t 32 --outfile $fname_prefix.lineage_report.csv $fname_prefix.fas

# clean up unzipped consensus files
rm *.fa

