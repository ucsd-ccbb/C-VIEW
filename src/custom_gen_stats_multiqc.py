#! /usr/bin/env python3

# Art Nasamran <cnasamran@health.ucsd.edu>
# Extracts P25 insert sizes from qualimapReports
# Calculates % >=Q30 from q30.txt
# Calculates % alignment from subsampled_mapping_stats.tsv
# Takes 4 arguments:
# 1. qualimapReports_paths.txt
# 2. q30_reads_paths.txt
# 3. subsampled_mapping_stats_paths.txt
# 4. se (single-end) or pe (paired-end)
# 5. output file path


from sys import argv
import os
import re
import yaml

NA_VAL = "NA"
P25_KEY = "P25 Ins. size"
PCT_30_KEY = "Pct >=Q30"
UNCAPPED_READS_KEY = "Uncapped Reads"
SUBSAMPLED_MAPPED_PCT_ALIGNED_KEY = "Sub Map Pct Aligned"
SE_VALUE = "se"
PE_VALUE = "pe"


def parseSeqQual(q30File):
    """ Parses q30.txt for and stores as a list [total seqs, q30 seqs]. """
    file_obj = open(q30File)
    file = file_obj.readlines()
    data = []
    parsing = False
    for line in file:
        line = line.strip()
        if line.startswith("total sequences"):
            parsing = True
        if line.startswith("q30 reads"):
            parsing = True
        if parsing:
            value = line.strip().split(" ")[2]
            value = int(float(value))
            data.append(value)
    file_obj.close()
    return data


# Parses _subsampled_mapping_stats.tsv and stores as
# a list [Mapped seqs, Unmapped seqs]
def parseSubMapStats(subMapStatsFile):
    file_obj = open(subMapStatsFile)
    file = file_obj.readlines()
    data = []
    parsing = False
    for line in file:
        line = line.strip()
        if line.startswith("Mapped"):
            parsing = True
        if line.startswith("Unmapped"):
            parsing = True
        if parsing:
            value = line.strip().split("\t")[1]
            value = int(float(value))
            data.append(value)
    file_obj.close()
    return data


def pairwise(it):
    # From https://stackoverflow.com/questions/5389507/
    # iterating-over-every-two-elements-in-a-list user mic_e
    it = iter(it)
    while True:
        try:
            yield next(it), next(it)
        except StopIteration:
            return


def get_sequenced_pool_component_id(filepath):
    """ Parses sequenced pool component id from a file path. """

    # get the name of the folder that the file is in
    fname = os.path.basename(os.path.dirname(filepath))
    # replace qualimap name cruft if extant
    name = fname.replace(".sorted.stats", "")
    return name


def gather_p25_ins_sizes(qmap_file_list_fp):
    with open(qmap_file_list_fp) as qmapFile:
        qr_paths = [line.strip() for line in qmapFile]

    # Dict to store P25s, %>Q30, and uncapped reads
    data_dict = {}

    # Extract and store insert sizes from qr_paths
    for qr in qr_paths:
        with open(qr) as file:
            sample_name = get_sequenced_pool_component_id(qr)
            for line in file:
                if line.startswith("<td class=column1>P25/Median/P75</td>"):
                    stats = next(file)
                    clean = re.compile('<.*?>')
                    stats = re.sub(clean, '', stats)
                    stats = stats.split(" /")
                    stats[0] = int(stats[0])
                    break
                else:
                    stats = [NA_VAL]
            data_dict[sample_name] = {P25_KEY: stats[0]}
    # p25sQ30s[sample_name].append(stats[0])

    return data_dict


def insert_q30_based_values(input_dict, R1, S1, R2=None, S2=None):
    name = get_sequenced_pool_component_id(R1)

    if R2 is not None:
        second_name = get_sequenced_pool_component_id(R2)
        if name != second_name:
            raise ValueError(f"Found different sequenced pool component "
                             f"ids in R1 and R2 paths ('{R1}' and '{R2}')")

    pctQ30 = uncapped_reads = NA_VAL
    try:
        if S2 is not None:
            uncapped_reads = S1[0] + S2[0]
            pctQ30 = (S1[1] + S2[1]) / (S1[0] + S2[0]) * 100
        else:
            uncapped_reads = S1[0]
            pctQ30 = S1[1] / S1[0] * 100
    except ZeroDivisionError:
        print(f"Warning: Unable to calculate values from q30 file due "
              f"to division by zero; reporting as {NA_VAL}")

    temp_sample_dict = input_dict.get(name, dict())

    if pctQ30 != NA_VAL:
        pctQ30 = round(pctQ30, 3)
    temp_sample_dict[PCT_30_KEY] = pctQ30
    temp_sample_dict[UNCAPPED_READS_KEY] = uncapped_reads
    input_dict[name] = temp_sample_dict


def gather_pct_gte_q30(q30_file_list_fp, se_or_pe, data_dict):
    # Load q30s and total uncapped reads
    with open(q30_file_list_fp) as q30sFile:
        q30s_paths = [line.strip() for line in q30sFile]

    if se_or_pe == SE_VALUE:
        for R1 in q30s_paths:
            S1 = parseSeqQual(R1)
            insert_q30_based_values(data_dict, R1, S1)
    elif se_or_pe == PE_VALUE:
        # Get Q30s for both R1 and R2
        for R1, R2 in pairwise(q30s_paths):
            S1 = parseSeqQual(R1)
            S2 = parseSeqQual(R2)
            insert_q30_based_values(data_dict, R1, S1, R2, S2)
    else:
        print(f"Warning: Unable to run generate percent greater than q30 and "
              f"number of uncapped reads for input type '{se_or_pe}'")

    return data_dict


def _calc_pct_aligned(reads):  # reads format: [Mapped, Unmapped]
    pctAligned = NA_VAL
    try:
        pctAligned = round(reads[0] / (reads[0] + reads[1]) * 100, 3)
    except ZeroDivisionError:
        print(f"Warning: Unable to calculate values from sub map stats "
              f"file due to division by zero; reporting as {NA_VAL}")

    return pctAligned


def gather_sub_map_pct(sub_map_stats_file_list_fp, data_dict):
    with open(sub_map_stats_file_list_fp) as subMapStatsFileList:
        subMapStats_paths = [line.strip() for line in subMapStatsFileList]

    for sample in subMapStats_paths:
        name = get_sequenced_pool_component_id(sample)
        reads = parseSubMapStats(sample)  # format: [Mapped, Unmapped]
        pctAligned = _calc_pct_aligned(reads)

        temp_sample_dict = data_dict.get(name, dict())
        temp_sample_dict[SUBSAMPLED_MAPPED_PCT_ALIGNED_KEY] = pctAligned
        data_dict[name] = temp_sample_dict

    return data_dict


def generate_multiqc_dict(data_dict):
    yaml_dict = {
        "custom_data": {
            "my_genstats": {
                "plot_type": "generalstats",
                "pconfig": [
                    {P25_KEY: {
                        "min": 0,
                        "scale": "RdYlGn"
                    }},
                    {PCT_30_KEY: {
                        "max": 100,
                        "min": 0,
                        "suffix": "%"
                    }},
                    {UNCAPPED_READS_KEY: {
                        "min": 0
                    }},
                    {SUBSAMPLED_MAPPED_PCT_ALIGNED_KEY: {
                        "max": 100,
                        "min": 0,
                        "suffix": "%"
                    }}
                ],
                "data": data_dict
            }
        }
    }

    return yaml_dict


def write_custom_multiqc_yaml(arg_list):
    qmap_file_list_fp = arg_list[1]
    q30_file_list_fp = arg_list[2]
    sub_map_stats_file_list_fp = arg_list[3]
    se_or_pe = arg_list[4]
    output_fp = arg_list[5]

    data_dict = gather_p25_ins_sizes(qmap_file_list_fp)
    if se_or_pe in [SE_VALUE, PE_VALUE]:
        data_dict = gather_pct_gte_q30(q30_file_list_fp, se_or_pe, data_dict)
        data_dict = gather_sub_map_pct(sub_map_stats_file_list_fp, data_dict)
    else:
        # for non-fastq-based (single-end or paired-end) data inputs like
        # genexus bams, we won't have the inputs needed to calculate some
        # of the metrics.  Instead set those metrics to None for each sample so
        # that any code looking for those keys at least finds something
        for curr_sample in data_dict:
            data_dict[curr_sample][PCT_30_KEY] = NA_VAL
            data_dict[curr_sample][UNCAPPED_READS_KEY] = NA_VAL
            data_dict[curr_sample][SUBSAMPLED_MAPPED_PCT_ALIGNED_KEY] = NA_VAL
        # next sample
    # endif

    yaml_dict = generate_multiqc_dict(data_dict)

    with open(output_fp, 'w') as yaml_file:
        yaml.dump(yaml_dict, yaml_file, default_flow_style=False)

    return yaml_dict


if __name__ == '__main__':
    write_custom_multiqc_yaml(argv)
