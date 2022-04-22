# C-VIEW: COVID-19 VIral Epidemiology Workflow

This software implements a high-throughput data processing pipeline to identify and charaterize SARS-CoV-2 variant sequences in specimens from COVID-19 positive hosts or environments.  It is based on https://github.com/niemasd/SD-COVID-Sequencing and built for use with Amazon Web Services (AWS) EC2 machine instances and S3 data storage.

Pipeline version 3.0.0 is pre-installed on the so-labeled Amazon Web Services snapshot in region us-west-2 (Oregon).  

## Installing the pipeline
The pipeline uses the following external software programs:

* [Minimap2 2.17-r941](https://github.com/lh3/minimap2/releases/tag/v2.17)
* [samtools 1.11](https://github.com/samtools/samtools/releases/tag/1.11)
* [Qualimap 2.2.2-dev](https://bitbucket.org/kokonech/qualimap/src/master/)
* [ivar 1.3.1](https://github.com/andersen-lab/ivar/releases/tag/v1.3.1)
* [Pangolin (variable version)](https://github.com/cov-lineages/pangolin)
* [ViralMSA 1.1.11](https://github.com/niemasd/ViralMSA/releases/tag/1.1.11)
* [q30 dev](https://github.com/artnasamran/q30)
* [samhead 1.0.0](https://github.com/niemasd/samhead/releases/tag/1.0.0)
* [pi_from_pileup 1.0.3](https://github.com/Niema-Docker/pi_from_pileup/releases/tag/1.0.3)
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

Note that, after creating a new cluster, 
the `aws cli` software must be configured on the head node with credentials for accessing the 
necessary AWS S3 resources.


## Running the pipeline

The pipeline is initiated on the head node of the cluster by calling the `run_pipeline.sh` script with an input csv file provided by the user.  This file should have a header line in the following format:

`function,organization,seq_run,merge_lanes,primer_set,fq,read_cap,sample,timestamp,istest`

It must then contain one or more data lines, each of which will trigger a run of the specified part(s) of the pipeline on a specified sequencing run dataset.

The fields are:

|Field Name|Allowed Values|Description|
|----------|--------------|-----------|
|`function`|cumulative_pipeline, pipeline, variants, sample, qc, lineages, phylogeny, cummulative_lineages, or cumulative_phylogeny | Specifies the type of functionality that should be run. See details below |
|`organization`|ucsd or helix|Specifies the organization from which all the samples in the current sequencing run are assumed to originate.  Helix sequencing runs can be combined only with data from other helix sequencing runs at the lineage and/or alignment-building steps.|
|`seq_run`|a string such as "210409_A00953_0272_AH57WJDRXY"|Specifies the sequencing center's identifier of the sequencing run to be processed, if relevant to the function provided.|
|`merge_lanes`|true or false|Indicates whether the pipeline should attempt to merge sample read data across fastq files from multiple lanes.|
|`primer_set`|artic or swift_v2|Specifies the primer set to use in trimming per-sample sorted bam files.|
|`fq`|se or pe|Indicates whether the pipeline should be run on only R1 reads in the sequencing run or on R1 and R2 reads.|
|`read_cap`|a positive integer or all|Specifies the maximum number of mapped reads per sample that should be used in the per-sample variant-calling and consensus-sequence-building functionality.|
|`sample`|a string such as "SEARCH-10003__D101802__I22__210608_A00953_0321_BH7L5LDSX2__S470_L002"|Specifies, for the sample to be processed, the part of the read one file name coming before `_R1_001.fastq.gz`.|
|`timestamp`|a string such as "2021-07-09_22-44-27"|Specifies the timestamp associated with the particular processing run that should be used.|
|`is_test`|true or false|Indicates whether data should be pulled from and written to the test S3 bucket (if true) or the production S3 bucket (if false)|

The functions supported by the pipeline are:

|Function|Description|
|--------|-----------|
|`cumulative_pipeline`|This is the primary usage. Runs variant calling, consensus sequence generation, and QC for a specified sequencing run, followed by lineage calling and alignment building on the cumulative set of all QC-passing consensus sequences ever processed by the pipeline|
|`pipeline`|Runs all pipeline functionality (including lineage calling and alignment building) for a specified sequencing run|
|`variants`|Runs variant calling and consensus sequence generation on all samples in the specified sequencing run|
|`sample`|Runs variant calling and consensus sequence generation on the specified sample in the specified sequencing run for the specified timestamp|
|`qc`|Runs QC on all outputs from the specified sequencing run processed under the specified timestamp|
|`lineages`|Runs lineage calling on all QC-passing consensus sequences in the specified sequencing run|
|`phylogeny`|Runs both lineage calling and alignment building on all QC-passing consensus sequences in the specified sequencing run|
|`cumulative_lineages`|Runs lineage calling on the cumulative set of all QC-passing consensus sequences ever processed by the pipeline|
|`cumulative_phylogeny`|Runs both lineage calling and alignment building on the cumulative set of all QC-passing consensus sequences ever processed by the pipeline|

For all functions except `sample`, some of the input fields are ignored, as shown in the table below:

|function|organization|seq_run|merge_lanes|primer_set|fq|read_cap|sample|timestamp|istest|
|--------|------------|------|----------|---------|---|-------|-----|------|------|
|cumulative_pipeline|ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|true or false|artic or swift_v2|se or pe|all or positive integer|ignored|ignored|true or false|
|pipeline|ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|true or false|artic or swift_v2|se or pe|all or positive integer|ignored|ignored|true or false|
|variants|ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|true or false|artic or swift_v2|se or pe|all or positive integer|ignored|ignored|true or false|
|sample|ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|true or false|artic or swift_v2|se or pe|all or positive integer|e.g. SEARCH-17043__D101859__L01__210409_A00953_0272_AH57WJDRXY__S82_L001_R1_001.fastq.gz|e.g. 2021-04-15_16-13-59|true or false|
|qc|ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|ignored|ignored|se or pe|ignored|ignored|e.g. 2021-04-15_16-13-59|true or false|
|lineages|ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|ignored|ignored|ignored|ignored|ignored|ignored|true or false|
|phylogeny|ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|ignored|ignored|ignored|ignored|ignored|ignored|true or false|
|cumulative_lineages|ucsd or helix|ignored|ignored|ignored|ignored|ignored|ignored|ignored|true or false|
|cumulative_phylogeny|ucsd or helix|ignored|ignored|ignored|ignored|ignored|ignored|ignored|true or false|

An example input for `run_pipeline.sh` might look like:

```
function,organization,seq_run,merge_lanes,primer_set,fq,read_cap,sample,timestamp,istest
cumulative_pipeline,ucsd,210608_A00953_0321_BH7L5LDSX2,false,swift_v2,pe,2000000,NA,NA,false
```
