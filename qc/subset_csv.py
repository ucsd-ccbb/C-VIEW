# Filters an input csv file by stripping out any (non-header) lines that
# do not have one of the supplied allowed values in the specified column.
# Note that this only works for discrete values--
# CANNOT be used to say, e.g.,
# "keep any lines where 0 < gc_content_percent < 48.5"


import pandas as pd
import sys

# file path of the csv file to filter
csv_fp = sys.argv[1]

# file path to which to write the filtered file
output_fp = sys.argv[2]

# column on which to filter
csv_col_name = sys.argv[3]

# comma-separated list of allowed values (**without spaces**)
allowed_vals_str = sys.argv[4]

# read csv file into dataframe
orig_table = pd.read_csv(csv_fp)

# Filter by keeping only records that have an allowed value in the
# specified column
allowed_vals = allowed_vals_str.split(",")
filtered_table = orig_table[orig_table[csv_col_name].isin(allowed_vals)]

# Write table
filtered_table.to_csv(output_fp, index=False)
