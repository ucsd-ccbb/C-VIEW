#! /usr/bin/env python3

'''
Art Nasamran <cnasamran@health.ucsd.edu>
Extracts P25 insert sizes from qualimapReports
Calculates % >=Q30 from fastqc_data.txt

'''
from sys import argv,stderr
from collections import Counter, defaultdict
import re

# Load qualimapReport paths
with open(str(argv[1])) as qmapFile:
<<<<<<< HEAD
    qr_paths = [line.strip() for line in qmapFile]
=======
	qr_paths = [line.strip() for line in qmapFile]
>>>>>>> 359da712806cfe22911f6f552800ab0008ef88ca
qmapFile.close()

# Dict to store P25s and %>Q30
p25sQ30s = defaultdict(list)

# Extract and store insert sizes from qr_paths
for qr in qr_paths:
<<<<<<< HEAD
    with open(qr) as file:
        for line in file:
            if line.startswith("<td class=column1>P25/Median/P75</td>"):
                stats = next(file)
                clean = re.compile('<.*?>')
                stats = re.sub(clean, '', stats)
                stats = stats.split(" /")
                # Clean qr to get sample name
                x = qr.split("/")
                sample_name = x[-3]
                p25sQ30s[sample_name].append(stats[0]) # dict of lists
=======
	with open(qr) as file:
		for line in file:
			if line.startswith("<td class=column1>P25/Median/P75</td>"):
				stats = next(file)
				clean = re.compile('<.*?>')
				stats = re.sub(clean, '', stats)
				stats = stats.split(" /")
				# Clean qr to get sample name
				x = qr.split("/")
				sample_name = x[-3]
				p25sQ30s[sample_name].append(stats[0]) # dict of lists
>>>>>>> 359da712806cfe22911f6f552800ab0008ef88ca
file.close()

# Extract and store per sequence quality scores from fastqc_data.txt
def parseSeqQual(fastqcFile):
<<<<<<< HEAD
    """ Parses fastqc_data.txt for Per sequence quality scores and stores as a dict. """
    file = open(fastqcFile).readlines()
    data = {}
    parsing = False
    for line in file:
        line = line.strip()
        if line.startswith("#Quality"):
            parsing = True
            continue    # Skip header
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

# Load fastqcs
with open(str(argv[2])) as fastQCsFile:
    fastQCs_paths = [line.strip() for line in fastQCsFile]
=======
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

# Load fastqcs
with open(str(argv[2])) as fastQCsFile:
	fastQCs_paths = [line.strip() for line in fastQCsFile]
>>>>>>> 359da712806cfe22911f6f552800ab0008ef88ca
fastQCsFile.close()


def pairwise(it):
<<<<<<< HEAD
    """ From https://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list user mic_e """
    it = iter(it)
    while True:
        try:
            yield next(it), next(it)
        except StopIteration:
            return
=======
	""" From https://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list user mic_e """
	it = iter(it)
	while True:
		try:
			yield next(it), next(it)
		except StopIteration:
			return
>>>>>>> 359da712806cfe22911f6f552800ab0008ef88ca

# Get Q30s for both R1 and R2
# ONLY WORKS FOR PE SAMPLES
for R1, R2 in pairwise(fastQCs_paths):
<<<<<<< HEAD
    S1 = parseSeqQual(R1)
    S2 = parseSeqQual(R2)
    combined = Counter(S1) + Counter(S2)
    # Add %Q30 to list
    name = R2.split("/")
    name = name[-2]
    name = name.replace("_R2_001_fastqc", "")
    pctQ30 = calcQ30(combined)
    p25sQ30s[name].append(round(pctQ30, 3)) # Add pctQ30 as the 2nd value in dict
=======
	S1 = parseSeqQual(R1)
	S2 = parseSeqQual(R2)
	combined = Counter(S1) + Counter(S2)
	# Add %Q30 to list
	name = R2.split("/")
	name = name[-2]
	name = name.replace("_R2_001_fastqc", "")
	pctQ30 = calcQ30(combined)
	p25sQ30s[name].append(round(pctQ30, 3)) # Add pctQ30 as the 2nd value in dict
>>>>>>> 359da712806cfe22911f6f552800ab0008ef88ca

# Write .yaml
yamlFile = open("multiqc_custom_gen_stats.yaml", "w")
# 4 spaces per "tab" to indent a yaml
yamlFile.write('\n\ncustom_data: \n    my_genstats:\n        plot_type: "generalstats"\n        pconfig:\n            - P25 Ins. size:\n                min: 0\n                scale: "RdYlGn"\n            - Pct >=Q30:\n                max: 100\n                min: 0\n                suffix: "%"\n        data:')
for keys in p25sQ30s:
<<<<<<< HEAD
    yamlFile.write("\n            " + keys + ":\n")
    yamlFile.write("                P25 Ins. size: " + str(p25sQ30s[keys][0]) + "\n                Pct >=Q30: " + str(p25sQ30s[keys][1]))
yamlFile.close()
=======
	yamlFile.write("\n            " + keys + ":\n")
	yamlFile.write("                P25 Ins. size: " + str(p25sQ30s[keys][0]) + "\n                Pct >=Q30: " + str(p25sQ30s[keys][1]))
yamlFile.close()
>>>>>>> 359da712806cfe22911f6f552800ab0008ef88ca
