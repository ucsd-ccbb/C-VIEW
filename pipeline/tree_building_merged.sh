aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210219_A00953_0236_BH2N77DRXY/2021-02-23_02-21-29_pe/210219_A00953_0236_BH2N77DRXY-passQC.fas ./

aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210109_A00953_0209_BHYHCVDRXX/2021-02-18_unmerged_pe/210109_A00953_0209_BHYHCVDRXX-passQC.fas ./

aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210213_A00953_0232_BHY3GLDRXX/2021-02-17_merged_pe/210213_A00953_0232_BHY3GLDRXX-passQC.fas ./

aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210116_A00953_0214_BHYGJGDRXX/2021-02-18_unmerged_pe/210116_A00953_0214_BHYGJGDRXX-passQC.fas ./

# add 3 new batches, 210223
aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210201_A00953_0224_BHY7HGDRXX/merged_lanes/2021-02-23_21-10-23_pe/210201_A00953_0224_BHY7HGDRXX-passQC.fas ./

aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210128_A00953_0221_BHYCY5DRXX/merged_lanes/2021-02-23_21-10-23_pe/210128_A00953_0221_BHYCY5DRXX-passQC.fas ./

aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/20210209/210205_A00953_0225_AHW73FDRXX/merged_lanes/2021-02-23_21-10-23_pe/210205_A00953_0225_AHW73FDRXX-passQC.fas ./

# add 2 outstanding batches, 210224
aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/2021-02-08-ARTIC/2021-02-20_00-07-14_pe/2021-02-08-ARTIC-passQC.fas ./

aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210214_A00953_0233_AHTNCMDSXY/2021-02-20_02-24-26_pe/210214_A00953_0233_AHTNCMDSXY-passQC.fas ./


# also grab the qc summary files, for help keeping track of batches

aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210219_A00953_0236_BH2N77DRXY/2021-02-23_02-21-29_pe/210219_A00953_0236_BH2N77DRXY-qc/210219_A00953_0236_BH2N77DRXY-summary.acceptance.tsv ./

aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210109_A00953_0209_BHYHCVDRXX/2021-02-18_unmerged_pe/210109_A00953_0209_BHYHCVDRXX-qc/210109_A00953_0209_BHYHCVDRXX-summary.acceptance.tsv ./

aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210213_A00953_0232_BHY3GLDRXX/2021-02-17_merged_pe/210213_A00953_0232_BHY3GLDRXX-qc/210213_A00953_0232_BHY3GLDRXX-summary.acceptance.tsv ./

aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210116_A00953_0214_BHYGJGDRXX/2021-02-18_unmerged_pe/210116_A00953_0214_BHYGJGDRXX-qc/210116_A00953_0214_BHYGJGDRXX-summary.acceptance.tsv ./

# 3 new batches, 210223
aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210201_A00953_0224_BHY7HGDRXX/merged_lanes/2021-02-23_21-10-23_pe/210201_A00953_0224_BHY7HGDRXX-qc/210201_A00953_0224_BHY7HGDRXX-summary.acceptance.tsv ./

aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210128_A00953_0221_BHYCY5DRXX/merged_lanes/2021-02-23_21-10-23_pe/210128_A00953_0221_BHYCY5DRXX-qc/210128_A00953_0221_BHYCY5DRXX-summary.acceptance.tsv ./

aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/20210209/210205_A00953_0225_AHW73FDRXX/merged_lanes/2021-02-23_21-10-23_pe/210205_A00953_0225_AHW73FDRXX-qc/210205_A00953_0225_AHW73FDRXX-summary.acceptance.tsv ./

# add 2 outstanding batches, 210224
aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/2021-02-08-ARTIC/2021-02-20_00-07-14_pe/2021-02-08-ARTIC-qc/2021-02-08-ARTIC-summary.acceptance.tsv ./

aws s3 cp s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/210214_A00953_0233_AHTNCMDSXY/2021-02-20_02-24-26_pe/210214_A00953_0233_AHTNCMDSXY-qc/210214_A00953_0233_AHTNCMDSXY-summary.acceptance.tsv ./

touch merged.fas # initialize the file
cat 210219_A00953_0236_BH2N77DRXY-passQC.fas >> merged.fas
cat 210109_A00953_0209_BHYHCVDRXX-passQC.fas >> merged.fas
cat 210213_A00953_0232_BHY3GLDRXX-passQC.fas >> merged.fas
cat 210116_A00953_0214_BHYGJGDRXX-passQC.fas >> merged.fas
# add 3 new batches, 210223
cat 210201_A00953_0224_BHY7HGDRXX-passQC.fas >> merged.fas
cat 210128_A00953_0221_BHYCY5DRXX-passQC.fas >> merged.fas
cat 210205_A00953_0225_AHW73FDRXX-passQC.fas >> merged.fas
# add 2 outstanding batches, 210224
cat 2021-02-08-ARTIC-passQC.fas >> merged.fas
cat 210214_A00953_0233_AHTNCMDSXY-passQC.fas >> merged.fas


wc -l 210219_A00953_0236_BH2N77DRXY-passQC.fas
wc -l 210109_A00953_0209_BHYHCVDRXX-passQC.fas
wc -l 210213_A00953_0232_BHY3GLDRXX-passQC.fas
wc -l 210116_A00953_0214_BHYGJGDRXX-passQC.fas
# 3 new batches, 210223
wc -l 210201_A00953_0224_BHY7HGDRXX-passQC.fas
wc -l 210128_A00953_0221_BHYCY5DRXX-passQC.fas
wc -l 210205_A00953_0225_AHW73FDRXX-passQC.fas
# add 2 outstanding batches, 210224
wc -l 2021-02-08-ARTIC-passQC.fas 
wc -l 210214_A00953_0233_AHTNCMDSXY-passQC.fas
wc -l merged.fas

# add the reference sequence
cat /shared/workspace/phylogenetics/reference_seqs/RmYN02.fas >> merged.fas
echo $str | wc -l merged.fas

# add B.1.1.7 sequence
cat /shared/workspace/phylogenetics/reference_seqs/hCoV-19_USA_CA-SEARCH-5574_2020.fasta >> merged.fas
echo $str | wc -l merged.fas

ViralMSA.py -s merged.fas -r SARS-CoV-2 -o viralmsa_out -t 32 -e brin.rosenthal@gmail.com

python /shared/workspace/phylogenetics/code/trim_msa.py -i viralmsa_out/merged.fas.aln -s 100 -e 50 -o merged.trimmed.aln

iqtree2 -T 32 -m GTR+F+G4 --polytomy -blmin 1e-9 -s merged.trimmed.aln

FastRoot.py -i merged.trimmed.aln.treefile -o merged.trimmed.aln.rooted.treefile -m OG -g "hCoV-19/bat/Yunnan/RmYN02/2019|EPI_ISL_412977|2019-06-25"

# --- metadata MAY CHANGE ----
pangolin -t 32 --outfile merged.lineage_report.csv merged.fas

awk -F "," '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' merged.lineage_report.csv > merged.metadata.txt

sed -i '1 a q2:types\tcategorical\tcategorical\tcategorical\tcategorical' merged.metadata.txt

# note currently adding batch info to metadata on local jupyter notebook... need to improve
# -------------------------

# eval "$(conda shell.bash hook)" # necessary to enter conda env from bash script
# the other way, found by amanda
#source /shared/workspace/software/anaconda3/bin activate qiime2-2020.11

conda activate qiime2-2020.11

empress tree-plot --tree merged.trimmed.aln.rooted.treefile --feature-metadata merged.metadata.batch.txt --output-dir tree-viz

conda deactivate


#aws s3 cp merged_210224 s3://ucsd-ccbb-projects/2021/20210208_COVID_sequencing/merged_210224/phylo --recursive