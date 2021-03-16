# SARS-CoV-2 variant identification

This software implements a high-throughput data processing pipeline to identify and charaterize SARS-CoV-2 variant sequences in specimens from COVID-19 positive hosts or environments.  It is based on https://github.com/niemasd/SD-COVID-Sequencing and built for use with Amazon Web Services (AWS) EC2 machine instances and S3 data storage.

Pipeline version 1.0.1 is pre-installed on the snap-02fcbfa37e4ab6709 Amazon Web Services snapshot in region us-west-2 (Oregon).  

## Installing the pipeline
The pipeline uses the following external software programs:

* ivar 1.3.1
* minimap2 2.17-r941
* samtools 1.11
* QualiMap v.2.2.2-dev
* FastQC v0.11.9
* Pangolin (variable version)
* viralMSA 1.1.11
* IQTree 2.1.2
* FastRoot v1.5
* EMPress 1.1.0

Should one wish to set up the pipeline on a fresh instance, follow the below commands.
Create a python 3.8 conda environment and activate it, then run:

```
conda install numpy 
conda install boto3
conda install -c bioconda fastqc
conda install -c bioconda qualimap
conda install -c bioconda minimap2
conda install -c bioconda samtools

# to run unit tests, also do the following
conda install -c conda-forge wand
```

Followed by:
```
pip install multiqc
pip install nwalign3
pip install pandas
pip install seaborn
```

Then install ivar from source (see https://github.com/andersen-lab/ivar ).

It is also necessary to create a separate conda environment for Pangolin, following the instructions at 
https://github.com/cov-lineages/pangolin#install-pangolin .

Install viralMSA, IQTree, and FastRoot from github, then install EMPress via a separate qiime2 conda environment (see https://github.com/biocore/empress).

Finally, install the main (latest stable release) branch of this repository from github.

## Specifying the pipeline cluster

The pipeline is optimized to run on an AWS EC2 cluster with the following characteristics:
```
master_instance_type = t2.medium
compute_instance_type = r5d.24xlarge
cluster_type = ondemand
ebs_settings = custom
base_os = ubuntu1604
scheduler = sge
compute_root_volume_size = 500
```


## Running the pipeline

The pipeline is initiated on the head node of the cluster by calling the `run_pipeline.sh` script with an input csv file provided by the user.  This file should have a header line in the following format:

`organization,seq_run,primers,read_type,merge,variants,qc,lineage,tree_build`

It must then contain one or more data lines, each of which will trigger a run of the specified part(s) of the pipeline on a specified sequencing run dataset.

The fields are:

|Field Name|Allowed Values|Description|
|----------|--------------|-----------|
|`organization`|ucsd or helix|Specifies the organization from which all the samples in the current sequencing run are assumed to originate.  Helix sequencing runs can be combined only with data from other helix sequencing runs at the lineage and/or tree-building steps.|
|`seq_run`|artic or swift_v2|Specifies the primer set to use in trimming per-sample sorted bam files.|
|`read_type`|se or pe|Indicates whether the pipeline should be run on only R1 reads in the sequencing run or on R1 and R2 reads.|
|`merge`|true or false|Indicates whether the pipeline should attempt to merge sample read data across fastq files from multiple lanes.|
|`variants`|true or false|Indicates whether the pipeline should execute the per-sample variant-calling and consensus-sequence-building functionality on the sequencing run dataset specified in this row.|
|`qc`|true or false|Indicates whether the pipeline should execute the per-sequencing-run qc-reporting and summarization functionality on the sequencing run dataset specified in this row.|
|`lineage`|true or false|Indicates whether the pipeline should execute the lineage-calling functionality on all cumulative data available to this organization.|
|`tree_build`|true or false|Indicates whether the pipeline should execute the tree-building and tree-visualization functionality on all cumulative data available to this organization.|

An example input for `run_pipeline.sh` might look like:

```
organization,seq_run,primers,read_type,merge,variants,qc,lineage,tree_build
ucsd,210213_A00953_0232_BHY3GLDRXX,swift_v2,pe,true,true,true,true,true
```

It is also possible to run more granular elements of the pipeline directly: `run_qc_summary.sh` executes only the `qc` functionality described above on a specified sequencing run, while `run_phylogeny.sh` executes only the `lineage` and/or `tree_build` functionality on all cumulative data available to the organization. These two scripts require input files containing the field above as well as one extra final field:

|Field Name|Allowed Values|Description|
|----------|--------------|-----------|
|processing_run|string containing alphanumeric characters, hyphens, and/or underscores only|Specifies an identifier of a data processing run.  For `run_qc_summary.sh`, this specifies the per-sample data processing run to be qc'd. For `run_phylogeny.sh`, it specifies the name of this phylogeny data-processing run, which prefixed to all output artifacts.  Usually this is a timestamp, such as `2021-03-11_21-58-55`.|

An example input for either of these scripts might look like:

```
organization,seq_run,primers,read_type,merge,variants,qc,lineage,tree_build,processing_run
ucsd,210213_A00953_0232_BHY3GLDRXX,swift_v2,pe,true,true,true,true,true,2021-03-11_21-18-32
```

Finally, it is possible to run the pipeline on a single sample using the `run_sample.sh` script.  This script requires all of the fields above in the specified order, as well as an additional first field:


|Field Name|Allowed Values|Description|
|----------|--------------|-----------|
|sample|string containing alphanumeric characters, hyphens, and/or underscores only|Specifies the name of the sample to be data-processed, as provided by the sequencing center.|

An example input for `run_sample.sh` might look like:

```
sample,organization,seq_run,primers,read_type,merge,variants,qc,lineage,tree_build,processing_run
002idSEARCH-5329-SAN_L001_L002_L003_L004,ucsd,210213_A00953_0232_BHY3GLDRXX,swift_v2,pe,true,true,true,true,true,2021-03-11_21-18-32
```
