'''
    Makes a QC summary table from:
    1. sequencing run's multiqc_general_stats.txt file
    2. sequencing run's <seq_run>-summary.acceptance.tsv file
    # 3. Pangolin lineage_report.csv - commented out as of 3/3/21
'''

import pandas as pd
import sys

# read multiqc_general_stats.txt and <seq_run>-summary.acceptance.tsv
# into dataframes
mqc = pd.read_csv(sys.argv[1], sep="\t")
sumAcc = pd.read_csv(sys.argv[2], sep="\t")
# linRep = pd.read_csv(sys.argv[3])

# Remove QualiMap_mqc-generalstats-qualimap-general_error_rate from multiqc df
mqc = mqc.drop("QualiMap_mqc-generalstats-qualimap-general_error_rate", axis=1)

# *outer* merge all fields in acceptance df into multiqc df
combined = mqc.merge(sumAcc, left_on="Sample", right_on="fastq_id", how="outer")

# # Extract sample name from linRep[["taxon"]]
# def format_taxon(x):
#     sampleName = x.split(".")[0]
#     sampleName = sampleName.replace("Consensus_", "")
#     return sampleName

# # Extract sample names
# linRep["Sample"] = linRep["taxon"].apply(format_taxon)

# combined = combined.merge(linRep, on = "Sample", how = "outer")
combined.columns = combined.columns.str.replace("QualiMap_mqc-generalstats-qualimap-", "")
combined.columns = combined.columns.str.replace("my_genstats_mqc-generalstats-my_genstats-", "")
combined.to_csv("QCSummaryTable.csv", index=False)
