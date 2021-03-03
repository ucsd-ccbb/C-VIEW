'''
    Makes a QC summary table from:
    1. multiqc general statistics
    2. .summary.acceptance.tsv
    3. Pangolin lineage_report.csv - commented out as of 3/3/21
'''

import pandas as pd
import sys

mqc = pd.read_csv(sys.argv[1], sep = "\t")
sumAcc = pd.read_csv(sys.argv[2], sep = "\t")
# linRep = pd.read_csv(sys.argv[3])

# Remove QualiMap_mqc-generalstats-qualimap-general_error_rate
mqc = mqc.drop("QualiMap_mqc-generalstats-qualimap-general_error_rate", axis = 1)

# Take is_accepted and coverage_gte_10_reads from sumAcc
combined = mqc.merge(sumAcc[["fastq_id", "is_accepted", "coverage_gte_10_reads"]], left_on = "Sample", right_on = "fastq_id", how = "outer")

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
combined.to_csv("QCSummaryTable.csv", index = False)
