#! /usr/bin/env python3

'''
Art Nasamran <cnasamran@health.ucsd.edu>
Extracts P25 insert sizes from qualimapReports
Calculates % >=Q30 from fastqc_data.txt
Takes 3 arguments:
1. qualimapReports_paths.txt
2. fastqc_data_paths.txt
3. se (single-end) or pe (paired-end)

'''

from sys import argv,stderr
from collections import Counter, defaultdict
import re

# Functions
def parseSeqQual(fastqcFile):
	""" Parses fastqc_data.txt for Per sequence quality scores and stores as a dict. """
	file = open(fastqcFile).readlines()
	data = {}
	parsing = False
	for line in file:
		line = line.strip()
		if line.startswith("#Quality"):
			parsing = True
			continue 	# Skip header
		elif line.startswith(">>END_MODULE"):
			parsing = False
		if parsing:
			key, value = line.strip().split("\t")
			value = int(float(value))
			data[key] = value
	return data

def calcQ30(table):
	""" Calculates percentage of reads >= Q30 from a dictionary in the format of quality:nSequences. """
	sumValues = sum(table.values())
	sumQ30 = 0
	for key, value in table.items():
		if int(key) >= 30:
			sumQ30 = sumQ30 + value
	GTQ30 = (sumQ30/sumValues)*100
	return(GTQ30)

def pairwise(it):
	""" From https://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list user mic_e. """
	it = iter(it)
	while True:
		try:
			yield next(it), next(it)
		except StopIteration:
			return

def getSampleName(x):
	""" Parses sample name from fastqc path. """
	name = x.split("/")
	name = name[-2]
	name = name.replace("_R1_001_fastqc", "")
	name = name.replace("_R1_fastqc", "")
	return name

## MAIN
# Load qualimapReport paths
with open(str(argv[1])) as qmapFile:
	qr_paths = [line.strip() for line in qmapFile]
qmapFile.close()

# Dict to store P25s and %>Q30
p25sQ30s = defaultdict(list)

# Extract and store insert sizes from qr_paths
for qr in qr_paths:
	with open(qr) as file:
		# Clean qr to get sample name
		x = qr.split("/")
		sample_name = x[-3]
		for line in file:
			if line.startswith("<td class=column1>P25/Median/P75</td>"):
				stats = next(file)
				clean = re.compile('<.*?>')
				stats = re.sub(clean, '', stats)
				stats = stats.split(" /")
				break
			else:
				stats = ["NA"]
		p25sQ30s[sample_name].append(stats[0])
	file.close()

# Load fastqcs
with open(str(argv[2])) as fastQCsFile:
	fastQCs_paths = [line.strip() for line in fastQCsFile]
fastQCsFile.close()

# se or pe from argv[3]
if argv[3] == "se":
	for R1 in fastQCs_paths:
		getSampleName(R1)
		S1 = parseSeqQual(R1)
		tbl = Counter(S1)
		name = getSampleName(R1)
		pctQ30 = calcQ30(tbl)
		p25sQ30s[name].append(round(pctQ30, 3)) # Add pctQ30 as the 2nd value in dict

elif argv[3] == "pe":
	# Get Q30s for both R1 and R2
	for R1, R2 in pairwise(fastQCs_paths):
		S1 = parseSeqQual(R1)
		S2 = parseSeqQual(R2)
		combined = Counter(S1) + Counter(S2)
		name = getSampleName(R1)
		pctQ30 = calcQ30(combined)
		p25sQ30s[name].append(round(pctQ30, 3)) # Add pctQ30 as the 2nd value in dict

else:
	print("Please select se or pe as argument 3.")

# Write .yaml
yamlFile = open("multiqc_custom_gen_stats.yaml", "w")
# 4 spaces per "tab" to indent a yaml
yamlFile.write('\n\ncustom_data: \n    my_genstats:\n        plot_type: "generalstats"\n        pconfig:\n            - P25 Ins. size:\n                min: 0\n                scale: "RdYlGn"\n            - Pct >=Q30:\n                max: 100\n                min: 0\n                suffix: "%"\n        data:')
for keys in p25sQ30s:
	yamlFile.write("\n            " + keys + ":\n")
	yamlFile.write("                P25 Ins. size: " + str(p25sQ30s[keys][0]) + "\n                Pct >=Q30: " + str(p25sQ30s[keys][1]))
yamlFile.close()
