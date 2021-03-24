# SARS-CoV-2 Consensus Sequence Generation and Variant Identification

This software implements a high-throughput data processing pipeline to generate consensus sequences and to identify and characterize variants in SARS-CoV-2 specimens from COVID-19 positive hosts or environments. It is based on https://github.com/niemasd/SD-COVID-Sequencing and is built for use with Amazon Web Services (AWS) [EC2 machine instances](https://aws.amazon.com/ec2/) and [S3 data storage](https://aws.amazon.com/s3/).

Pipeline version 1.1.0 is pre-installed on the snap-0444a3ef8170b3bbe Amazon Web Services snapshot in region us-west-2 (Oregon).  

## Installing the Pipeline
The pipeline uses the following external software programs:

* [iVar 1.3.1](https://github.com/andersen-lab/ivar/releases/tag/v1.3.1)
* [Minimap2 2.17-r941](https://github.com/lh3/minimap2/releases/tag/v2.17)
* [samtools 1.11](https://github.com/samtools/samtools/releases/tag/1.11)
* [Qualimap 2.2.2-dev](https://bitbucket.org/kokonech/qualimap/src/master/)
* [FastQC 0.11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip)
* [Pangolin (variable version)](https://github.com/cov-lineages/pangolin)
* [ViralMSA 1.1.11](https://github.com/niemasd/ViralMSA/releases/tag/1.1.11)
* [IQ-TREE 2.1.2](https://github.com/iqtree/iqtree2/releases/download/v2.1.2/iqtree-2.1.2-Linux.tar.gz)
* [FastRoot 1.5](https://github.com/uym2/MinVar-Rooting/releases/tag/v1.5)
* [EMPress 1.1.0](https://github.com/biocore/empress/releases/tag/v1.1.0)

Should one wish to set up the pipeline on a fresh instance, follow the below commands.
Create a Python 3.8 conda environment and activate it, then run the following:

```bash
conda install numpy 
conda install boto3
conda install -c bioconda fastqc
conda install -c bioconda qualimap
conda install -c bioconda minimap2
conda install -c bioconda samtools

# to run unit tests, also do the following
conda install -c conda-forge wand
```

Followed by the following:

```bash
pip install multiqc
pip install nwalign3
pip install pandas
pip install seaborn
```

Then install iVar from source (see https://github.com/andersen-lab/ivar).

It is also necessary to create a separate conda environment for Pangolin, following the instructions at 
https://github.com/cov-lineages/pangolin#install-pangolin.

Install ViralMSA, IQ-TREE 2, and FastRoot from GitHub, then install EMPress via a separate QIIME 2 conda environment (see https://github.com/biocore/empress).

Finally, install the main (latest stable release) branch of this repository from GitHub.

## Specifying the Pipeline Cluster

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


## Running the Pipeline

The pipeline is initiated on the head node of the cluster by calling the [`run_pipeline.sh`](pipeline/run_pipeline.sh) script with an input CSV file provided by the user. This file should have a header line in the following format:

```
organization,seq_run,primers,read_type,merge,variants,qc,lineage,tree_build
```

It must then contain one or more data lines, each of which will trigger a run of the specified part(s) of the pipeline on a specified sequencing run dataset.

The fields are as follows:

|Field Name|Allowed Values|Description|
|----------|--------------|-----------|
|`organization`|`ucsd` or `helix`|Specifies the organization from which all the samples in the current sequencing run are assumed to originate. Helix sequencing runs can be combined only with data from other helix sequencing runs at the lineage and/or tree-building steps.|
|`seq_run`|`artic` or `swift_v2`|Specifies the primer set to use in trimming per-sample sorted bam files.|
|`read_type`|`se` or `pe`|Indicates whether the pipeline should be run on only R1 reads in the sequencing run or on R1 and R2 reads.|
|`merge`|`true` or `false`|Indicates whether the pipeline should attempt to merge sample read data across fastq files from multiple lanes.|
|`variants`|`true` or `false`|Indicates whether the pipeline should execute the per-sample variant-calling and consensus-sequence-building functionality on the sequencing run dataset specified in this row.|
|`qc`|`true` or `false`|Indicates whether the pipeline should execute the per-sequencing-run qc-reporting and summarization functionality on the sequencing run dataset specified in this row.|
|`lineage`|`true` or `false`|Indicates whether the pipeline should execute the lineage-calling functionality on all cumulative data available to this organization.|
|`tree_build`|`true` or `false`|Indicates whether the pipeline should execute the tree-building and tree-visualization functionality on all cumulative data available to this organization.|

An example input for [`run_pipeline.sh`](pipeline/run_pipeline.sh) might look like the following:

```
organization,seq_run,primers,read_type,merge,variants,qc,lineage,tree_build
ucsd,210213_A00953_0232_BHY3GLDRXX,swift_v2,pe,true,true,true,true,true
```

It is also possible to run more granular elements of the pipeline directly: [`run_qc_summary.sh`](qc/run_qc_summary.sh) executes only the `qc` functionality described above on a specified sequencing run, while [`run_phylogeny.sh`](pipeline/run_phylogeny.sh) executes only the `lineage` and/or `tree_build` functionality on all cumulative data available to the organization. These two scripts require input files containing the field above as well as one extra final field:

|Field Name|Allowed Values|Description|
|----------|--------------|-----------|
|processing_run|string containing alphanumeric characters, hyphens, and/or underscores only|Specifies an identifier of a data processing run.  For [`run_qc_summary.sh`](qc/run_qc_summary.sh), this specifies the per-sample data processing run to be QC'd. For [`run_phylogeny.sh`](pipeline/run_phylogeny.sh), it specifies the name of this phylogeny data-processing run, which prefixed to all output artifacts.  Usually this is a timestamp, such as `2021-03-11_21-58-55`.|

An example input for either of these scripts might look like the following:

```
organization,seq_run,primers,read_type,merge,variants,qc,lineage,tree_build,processing_run
ucsd,210213_A00953_0232_BHY3GLDRXX,swift_v2,pe,true,true,true,true,true,2021-03-11_21-18-32
```

Finally, it is possible to run the pipeline on a single sample using the [`run_sample.sh`](pipeline/run_sample.sh) script.  This script requires all of the fields above in the specified order, as well as an additional first field:


|Field Name|Allowed Values|Description|
|----------|--------------|-----------|
|sample|string containing alphanumeric characters, hyphens, and/or underscores only|Specifies the name of the sample to be data-processed, as provided by the sequencing center.|

An example input for [`run_sample.sh`](pipeline/run_sample.sh) might look like the following:

```
sample,organization,seq_run,primers,read_type,merge,variants,qc,lineage,tree_build,processing_run
002idSEARCH-5329-SAN_L001_L002_L003_L004,ucsd,210213_A00953_0232_BHY3GLDRXX,swift_v2,pe,true,true,true,true,true,2021-03-11_21-18-32
```
