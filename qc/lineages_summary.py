# Gathers seq-run_summary.csv files and merges with pangolin lineage table
# Takes 2 arguments:
# 1. Path to directory containing per-sequence-run summary.csv files
# 2. Path of pangolin lineage file


import pandas as pd
import sys
import glob
import os

# Merge tables with a single header
combinedQCs = pd.concat(map(pd.read_csv, glob.glob(os.path.join(sys.argv[1], "*.csv"))))

# Load pangolin file
lineage_file = pd.read_csv(sys.argv[2])

def format_taxon(x):
	""" Parses sample name from taxon column in pangolin lineage report """
	sampleName = x.split(".")[0]
	sampleName = sampleName.replace("Consensus_", "")
	return sampleName

# Add Sample column to pangolin lineage report
lineage_file["Sample"] = lineage_file["taxon"].apply(format_taxon)

# Merge tables
combined = combinedQCs.merge(lineage_file, left_on = "Sample", right_on = "Sample", how = "outer")

# Write table
combined.to_csv("merged_qc_and_lineages.csv", index = False)