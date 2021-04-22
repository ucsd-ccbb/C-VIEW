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

def parseSubMapStats(subMapStatsFile):
    """ Parses _subsampled_mapping_stats.tsv and stores as a list [Mapped seqs, Unmapped seqs]. """
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
    # p25sQ30s = defaultdict(list)
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
            data_dict[sample_name] = {"P25 Ins. size": stats[0]}
    # p25sQ30s[sample_name].append(stats[0])

    return data_dict


def insert_pct_gte_q30_in_sample_dict(data_dict, name, pctQ30):
    temp_sample_dict = data_dict.get(name, dict())
    if pctQ30 != NA_VAL:
        pctQ30 = round(pctQ30, 3)
    temp_sample_dict["Pct >=Q30"] = pctQ30
    return temp_sample_dict


def insert_uncapped_reads_in_sample_dict(data_dict, name, uncappedReads):
    temp_sample_dict = data_dict.get(name, dict())
    temp_sample_dict["Uncapped Reads"] = uncappedReads
    return temp_sample_dict


def insert_sub_map_pct_in_sample_dict(data_dict, name, pctSubMap):
    temp_sample_dict = data_dict.get(name, dict())
    temp_sample_dict["Sub Map Pct Aligned"] = pctSubMap
    return temp_sample_dict


def generate_q30_based_values(input_dict, R1, S1, R2=None, S2=None):
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
            pctQ30 = S1[1] / S1[0]
    except ZeroDivisionError:
        print(f"Warning: Unable to calculate values from q30 file due "
              f"to division by zero; reporting as {NA_VAL}")

    temp_pctQ30_dict = insert_pct_gte_q30_in_sample_dict(
        input_dict, name, pctQ30)
    temp_uncapped_reads_dict = insert_uncapped_reads_in_sample_dict(
        input_dict, name, uncapped_reads)

    return name, temp_pctQ30_dict, temp_uncapped_reads_dict


def gather_pct_gte_q30(q30_file_list_fp, se_or_pe, data_dict):
    # Load q30s and total uncapped reads
    with open(q30_file_list_fp) as q30sFile:
        q30s_paths = [line.strip() for line in q30sFile]

    if se_or_pe == "se":
        for R1 in q30s_paths:
            S1 = parseSeqQual(R1)

            name, temp_pctQ30_dict, temp_uncapped_reads_dict = \
                generate_q30_based_values(data_dict, R1, S1)

            data_dict[name] = temp_pctQ30_dict
            data_dict[name] = temp_uncapped_reads_dict
    elif se_or_pe == "pe":
        # Get Q30s for both R1 and R2
        for R1, R2 in pairwise(q30s_paths):
            S1 = parseSeqQual(R1)
            S2 = parseSeqQual(R2)

            name, temp_pctQ30_dict, temp_uncapped_reads_dict = \
                generate_q30_based_values(data_dict, R1, S1, R2, S2)

            data_dict[name] = temp_pctQ30_dict
            data_dict[name] = temp_uncapped_reads_dict
    else:
        raise ValueError(f"Unrecognized se or pe value '{se_or_pe}'")

    return data_dict

def gather_sub_map_pct(sub_map_stats_file_list_fp, data_dict):
    with open(sub_map_stats_file_list_fp) as subMapStatsFileList:
        subMapStats_paths = [line.strip() for line in subMapStatsFileList]

    for sample in subMapStats_paths:
        name = get_sequenced_pool_component_id(sample)
        reads = parseSubMapStats(sample) #format: [Mapped, Unmapped]
        pctAligned = round(reads[0]/(reads[0]+reads[1])*100, 3)

        temp_pctAligned_dict = insert_sub_map_pct_in_sample_dict(data_dict, name, pctAligned)

        data_dict[name] = temp_pctAligned_dict

    return data_dict


def generate_multiqc_dict(data_dict):
    yaml_dict = {
        "custom_data": {
            "my_genstats": {
                "plot_type": "generalstats",
                "pconfig": [
                    {"P25 Ins. size": {
                        "min": 0,
                        "scale": "RdYlGn"
                    }},
                    {"Pct >=Q30": {
                        "max": 100,
                        "min": 0,
                        "suffix": "%"
                    }},
                    {"Uncapped Reads": {
                        "min": 0
                    }},
                    {"Sub Map Pct Aligned": {
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
    data_dict = gather_pct_gte_q30(q30_file_list_fp, se_or_pe, data_dict)
    data_dict = gather_sub_map_pct(sub_map_stats_file_list_fp, data_dict)
    yaml_dict = generate_multiqc_dict(data_dict)

    with open(output_fp, 'w') as yaml_file:
        yaml.dump(yaml_dict, yaml_file, default_flow_style=False)

    return yaml_dict


if __name__ == '__main__':
    write_custom_multiqc_yaml(argv)
