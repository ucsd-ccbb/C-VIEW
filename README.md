# SARS-CoV-2 variant identification

This software implements a high-throughput data processing pipeline to identify and charaterize SARS-CoV-2 variant sequences in specimens from COVID-19 positive hosts or environments.  It is based on https://github.com/niemasd/SD-COVID-Sequencing and built for use with Amazon Web Services (AWS) EC2 machine instances and S3 data storage.

Pipeline version 1.9.5 is pre-installed on the snap-XXXX Amazon Web Services snapshot in region us-west-2 (Oregon).  

## Installing the pipeline
The pipeline uses the following external software programs:

* ivar 1.3.1
* minimap2 2.17-r941
* samtools 1.11
* QualiMap v.2.2.2-dev
* Pangolin (variable version)
* viralMSA 1.1.11
* IQTree 2.1.2
* FastRoot v1.5
* EMPress 1.1.0
* q30 (no version info)
* samhead (no version info)
* git 2.7.4 or higher

Should one wish to set up the pipeline on a fresh AWS ubuntu instance, download the `install.sh` script from this repository, set the 
necessary variables at the top of the script, and run it. Sudo permissions are required.

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

`organization,seq_run,primers,read_type,merge,variants,qc,lineage,tree_build,read_cap,is_test`

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
|`read_cap`|a positive integer or all|Specifies the maximum number of mapped reads per sample that should be used in the per-sample variant-calling and consensus-sequence-building functionality.|
|`is_test`|true or false|Indicates whether the pipeline should execute the tree-building and tree-visualization functionality on all cumulative data available to this organization.|


An example input for `run_pipeline.sh` might look like:

```
organization,seq_run,primers,read_type,merge,variants,qc,lineage,tree_build,read_cap,is_test
ucsd,210213_A00953_0232_BHY3GLDRXX,swift_v2,pe,true,true,true,true,true,all,false
```

It is also possible to run more granular elements of the pipeline directly.  
`run_phylogeny.sh` executes only the `lineage` and/or `tree_build` functionality on all cumulative data available to the organization. It takes the inputs

`organization,lineage,tree_build,is_test`

Some of the more granular scripts require an additional `processing_run` argument:

|Field Name|Allowed Values|Description|
|----------|--------------|-----------|
|processing_run|string containing alphanumeric characters, hyphens, and/or underscores only|Specifies an identifier of a data processing run.  For `run_qc_summary.sh`, this specifies the per-sample data processing run to be qc'd. Usually this is a timestamp, such as `2021-03-11_21-58-55`.|

`run_qc_summary.sh` executes only the `qc` functionality described above on a specified sequencing run; it takes the inputs

`organization,seq_run,primers,read_type,merge,variants,qc,lineage,tree_build,processing_run,is_test`

(Note the substitution of `processing_run` for `read_cap`).

Finally, it is possible to run the pipeline on a single sample using the `run_sample.sh` script.  This script requires an additional first field:


|Field Name|Allowed Values|Description|
|----------|--------------|-----------|
|sample|string containing alphanumeric characters, hyphens, and/or underscores only|Specifies the name of the sample to be data-processed, as provided by the sequencing center.|

An example input for `run_sample.sh` might look like:

```
sample,organization,seq_run,primers,read_type,merge,variants,qc,lineage,tree_build,read_cap,processing_run
002idSEARCH-5329-SAN_L001_L002_L003_L004,ucsd,210213_A00953_0232_BHY3GLDRXX,swift_v2,pe,true,true,true,true,true,all,2021-03-11_21-18-32
```
Note that this script does not take an `is_test` argument.
