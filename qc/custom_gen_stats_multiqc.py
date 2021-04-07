#! /usr/bin/env python3

# Art Nasamran <cnasamran@health.ucsd.edu>
# Extracts P25 insert sizes from qualimapReports
# Calculates % >=Q30 from q30.txt
# Takes 4 arguments:
# 1. qualimapReports_paths.txt
# 2. q30_reads_paths.txt
# 3. se (single-end) or pe (paired-end)
# 4. output file path


from sys import argv
import re
import yaml


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
            # continue  # Skip header
        if parsing:
            value = line.strip().split(" ")[2]
            value = int(float(value))
            data.append(value)
    file_obj.close()
    return data


def pairwise(it):
    """ From https://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list user mic_e. """
    it = iter(it)
    while True:
        try:
            yield next(it), next(it)
        except StopIteration:
            return


def getSampleName(x):
    """ Parses sample name from q30 path. """
    name = x.split("/")
    name = name[-2]
    name = name.replace("_R1_001", "")
    name = name.replace("_R1", "")
    return name


def gather_p25_ins_sizes(qmap_file_list_fp):
    with open(qmap_file_list_fp) as qmapFile:
        qr_paths = [line.strip() for line in qmapFile]

    # Dict to store P25s and %>Q30
    # p25sQ30s = defaultdict(list)
    data_dict = {}

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
                    stats[0] = int(stats[0])
                    break
                else:
                    stats = ["NA"]
            data_dict[sample_name] = {"P25 Ins. size": stats[0]}
    # p25sQ30s[sample_name].append(stats[0])

    return data_dict


def insert_pct_gte_q30_in_sample_dict(data_dict, name, pctQ30):
    temp_sample_dict = data_dict.get(name, dict())
    temp_sample_dict["Pct >=Q30"] = round(pctQ30, 3)
    return temp_sample_dict


def gather_pct_gte_q30(q30_file_list_fp, se_or_pe, data_dict):
    # Load q30s
    with open(q30_file_list_fp) as q30sFile:
        q30s_paths = [line.strip() for line in q30sFile]

    if se_or_pe == "se":
        for R1 in q30s_paths:
            getSampleName(R1)
            S1 = parseSeqQual(R1)
            # tbl = Counter(S1)
            name = getSampleName(R1)
            # pctQ30 = calcQ30(tbl)
            pctQ30 = S1[1]/S1[0]
            data_dict[name] = insert_pct_gte_q30_in_sample_dict(
                data_dict, name, pctQ30)
    elif se_or_pe == "pe":
        # Get Q30s for both R1 and R2
        for R1, R2 in pairwise(q30s_paths):
            S1 = parseSeqQual(R1)
            S2 = parseSeqQual(R2)
            # combined = Counter(S1) + Counter(S2)
            name = getSampleName(R1)
            # pctQ30 = calcQ30(combined)
            pctQ30 = (S1[1] + S2[1])/(S1[0] + S2[0]) * 100
            data_dict[name] = insert_pct_gte_q30_in_sample_dict(
                data_dict, name, pctQ30)
    else:
        raise ValueError(f"Unrecognized se or pe value '{se_or_pe}'")

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
    se_or_pe = arg_list[3]
    output_fp = arg_list[4]

    data_dict = gather_p25_ins_sizes(qmap_file_list_fp)
    data_dict = gather_pct_gte_q30(q30_file_list_fp, se_or_pe, data_dict)
    yaml_dict = generate_multiqc_dict(data_dict)

    with open(output_fp, 'w') as yaml_file:
        yaml.dump(yaml_dict, yaml_file, default_flow_style=False)

    return yaml_dict


if __name__ == '__main__':
    write_custom_multiqc_yaml(argv)
