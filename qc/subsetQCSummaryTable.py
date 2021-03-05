'''
	Subsets QCSummaryTable.csv by run names
	Takes 2 arguments:
	1. Path to QCSummaryTable.csv
	2. Comma separated list of run names (without spaces)
	Example: python subsetQCSummaryTable.py QCSummaryWithLineage.csv test_run1,test_run3,test_run5
'''

import pandas as pd
import sys

# Read QCSummaryTable.csv
sumTable = pd.read_csv(sys.argv[1])

# Subset sumTable by runs listed in sys.argv[2]
selectRuns = sumTable[sumTable["seq_run"].isin(sys.argv[2].split(","))]

# Write table
selectRuns.to_csv("SubsetQCSummaryWithLineage.csv", index = False)