# C-VIEW: COVID-19 VIral Epidemiology Workflow

This software implements a high-throughput data processing pipeline to identify and charaterize SARS-CoV-2 variant sequences in specimens from COVID-19 positive hosts or environments.  It is based on https://github.com/niemasd/SD-COVID-Sequencing and built for use with Amazon Web Services (AWS) EC2 machine instances and S3 data storage.

# Table of Contents
1. [Installing the Pipeline](#Installing-the-Pipeline)
2. [Creating a Cluster](#Creating-a-Cluster)
3. [Running the Pipeline](#Running-the-Pipeline)


## Installing the Pipeline

**Note: Usually it will NOT be necessary to install the pipeline from scratch.**  The most
current version of the pipeline is pre-installed on the so-labeled
Amazon Web Services snapshot in region us-west-2 (Oregon),
and this snapshot can be used directly to [create a cluster](#Creating-a-Cluster).

If a fresh installation *is* required, take the following steps:

1. On AWS, launch a new ubuntu 20.04 instance
   1. Note that it MUST be version 20.04, not the latest available ubuntu version (e.g. 22.04) because 20.04 is the latest version supported by AWS ParallelCluster.
   2. Select type t2.medium
   3. Add a 35 GB root drive and a 300 GB EBS drive
   4. Set the security group to allow SSH via TCP on port 22 and all traffic via all protocols on all ports
2. `ssh` onto new instance to set up file system and mount
   1. Run `lsblk` to find the name of the 300 GB EBS drive. For the remainder, of this section, assume `lsblk` shows that the name of the 300 GB volume is `xvdb`. 
   2. Run `sudo mkfs -t xfs /dev/xvdb` to make a filesystem on the new drive 
   3. Run `sudo mkdir /shared` to create a location for the installation 
   4. Run `sudo mount /dev/xvdb /shared` to mount the 300 GB volume to the new location
   5. Run ``sudo chown `whoami` /shared`` to grant the current user permissions to the new location
3. Install anaconda and python
   1. Run `cd /shared`
   2. Run `wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh`
   3. Run `bash Anaconda3-2020.11-Linux-x86_64.sh`
   4. Answer `yes` when asked to accept the license agreement
   5. Enter `/shared/workspace/software/anaconda3` when asked for the install location
   6. Answer `yes` when asked whether to have the installer run conda init
   7. Log out of the `ssh` session and then back in to allow the conda install to take effect
4. Install C-VIEW
   1. Run `cd /shared`
   2. Download `install.sh`
   3. Run `bash install.sh`
   4. Answer yes whenever asked for permission to proceed
5. On AWS, make a snapshot of the newly installed 300 GB volume

Note that the pipeline uses the following external software programs, which are 
installed via the `install.sh` script:

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


## Creating a Cluster

The pipeline is designed to run on a version 3 or later AWS ParallelCluster. 
Begin by ensuring that ParallelCluster is installed on your local machine; if it
is not, take these steps:

1. Set up a `conda` environment and and install ParallelCluster 
   1. Run `conda create --name parallelcluster3 python=3`
   2. Run `conda activate parallelcluster3`
   3. Run `python3 -m pip install --upgrade aws-parallelcluster`
2. In the `parallelcluster3` environment, install Node Version Manager and Node.js, which are (apparently) required by AWS Cloud Development Kit (CDK)
   1. Run `curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.38.0/install.sh | bash`
   2. Run `chmod ug+x ~/.nvm/nvm.sh`
   3. Run `source ~/.nvm/nvm.sh`
   4. Run `nvm install --lts`
   5. Check your install by running `node --version` and `pcluster version`

Next, ensure you have a pem file registered with AWS and that you have run `aws configure`
locally to set up AWS command line access from your local machine.  Then 
prepare a cluster configuration yaml file using the below template:

```
Region: us-west-2
Image:
  Os: ubuntu2004
SharedStorage:
  - MountDir: /shared
    Name: custom
    StorageType: Ebs
    EbsSettings:
      Size: 300
      SnapshotId: <snapshot of current cview release, e.g. snap-09264bf02660b54ad >
HeadNode:
  InstanceType: t2.medium
  Networking:
    SubnetId: subnet-06ff527fa2d5827a3
# subnet-06ff527fa2d5827a3 is parallelcluster:public-subnet
  Ssh:
    KeyName: <name of your pem file without extension, e.g. my_key for a file named my_key.pem >
Scheduling:
  Scheduler: slurm
  SlurmQueues:
    - Name: default-queue
      ComputeSettings:
        LocalStorage:
          RootVolume:
            Size: 500
      Networking:
        SubnetIds:
          - subnet-06ff527fa2d5827a3
# subnet-06ff527fa2d5827a3 is parallelcluster:public-subnet
      ComputeResources:
        - Name: default-resource
          MaxCount: 15
          InstanceType: r5d.24xlarge
```

To create a new cluster from the command line, run

```
pcluster create-cluster \
    --cluster-name <your-cluster-name> \
    --cluster-configuration <your-config-file-name>.yaml
```

(If you experience an error referencing Node.js, you may need to once again
run `source ~/.nvm/nvm.sh` to ensure it is accessible from your shell.)  The 
cluster creation progress can be monitored from the `CloudFormation`->`Stacks` section of the AWS Console.

Once the cluster is successfully created, log in to the head node.  To avoid 
having to use its public IPv4 DNS, one can run

`pcluster ssh --cluster-name <your_cluster_name> -i /path/to/keyfile.pem`

which fills in the cluster IP address and username automatically.

From the head node, run `aws configure` to set up the head node with credentials for accessing the 
necessary AWS S3 resources.


## Running the Pipeline

The pipeline is initiated on the head node of the cluster by calling the 
`run_cview.sh` script with an input csv file provided by the user, e.g.:

`bash /shared/workspace/software/cview/pipeline/run_cview.sh /shared/runfiles/cview_test_run.csv`

This file should have a header line in the following format:

`function,organization,seq_run,merge_lanes,primer_set,fq,read_cap,sample,timestamp,istest`

It must then contain one or more data lines, each of which will trigger a run of the specified part(s) of the pipeline on a specified sequencing run dataset.

The fields are:

|Field Name| Allowed Values                                                                                                                       |Description|
|----------|--------------------------------------------------------------------------------------------------------------------------------------|-----------|
|`function`| cumulative_pipeline, pipeline, variants_qc, variants, sample, qc, lineages, phylogeny, cummulative_lineages, or cumulative_phylogeny | Specifies the type of functionality that should be run. See details below |
|`organization`| ucsd or helix                                                                                                                        |Specifies the organization from which all the samples in the current sequencing run are assumed to originate.  Helix sequencing runs can be combined only with data from other helix sequencing runs at the lineage and/or alignment-building steps.|
|`seq_run`| a string such as "210409_A00953_0272_AH57WJDRXY"                                                                                     |Specifies the sequencing center's identifier of the sequencing run to be processed, if relevant to the function provided.|
|`merge_lanes`| true or false                                                                                                                        |Indicates whether the pipeline should attempt to merge sample read data across fastq files from multiple lanes.|
|`primer_set`| artic or swift_v2                                                                                                                    |Specifies the primer set to use in trimming per-sample sorted bam files.|
|`fq`| se or pe                                                                                                                             |Indicates whether the pipeline should be run on only R1 reads in the sequencing run or on R1 and R2 reads.|
|`read_cap`| a positive integer or all                                                                                                            |Specifies the maximum number of mapped reads per sample that should be used in the per-sample variant-calling and consensus-sequence-building functionality.|
|`sample`| a string such as "SEARCH-10003__D101802__I22__210608_A00953_0321_BH7L5LDSX2__S470_L002"                                              |Specifies, for the sample to be processed, the part of the read one file name coming before `_R1_001.fastq.gz`.|
|`timestamp`| a string such as "2021-07-09_22-44-27"                                                                                               |Specifies the timestamp associated with the particular processing run that should be used.|
|`is_test`| true or false                                                                                                                        |Indicates whether data should be pulled from and written to the test S3 bucket (if true) or the production S3 bucket (if false)|

The functions supported by the pipeline are:

| Function               | Description                                                                                                                                                                                                                                                          |
|------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `cumulative_pipeline`  | This is the primary usage. Runs variant calling, consensus sequence generation, and QC for a specified sequencing run, followed by lineage calling and alignment building on the cumulative set of all QC-passing consensus sequences ever processed by the pipeline |
| `pipeline`             | Runs all pipeline functionality (including lineage calling and alignment building) for a specified sequencing run                                                                                                                                                    |
| `variants_qc`          | Runs variant calling, consensus sequence generation, and QC for a specified sequencing run                                                                                                                                                                |
| `variants`             | Runs variant calling and consensus sequence generation on all samples in the specified sequencing run                                                                                                                                                                |
| `sample`               | Runs variant calling and consensus sequence generation on the specified sample in the specified sequencing run for the specified timestamp                                                                                                                           |
| `qc`                   | Runs QC on all outputs from the specified sequencing run processed under the specified timestamp                                                                                                                                                                     |
| `lineages`             | Runs lineage calling on all QC-passing consensus sequences in the specified sequencing run                                                                                                                                                                           |
| `phylogeny`            | Runs both lineage calling and alignment building on all QC-passing consensus sequences in the specified sequencing run                                                                                                                                               |
| `cumulative_lineages`  | Runs lineage calling on the cumulative set of all QC-passing consensus sequences ever processed by the pipeline                                                                                                                                                      |
| `cumulative_phylogeny` | Runs both lineage calling and alignment building on the cumulative set of all QC-passing consensus sequences ever processed by the pipeline                                                                                                                          |

For all functions except `sample`, some of the input fields are ignored, as shown in the table below:

| function             |organization|seq_run|merge_lanes|primer_set|fq|read_cap|sample|timestamp|istest|
|----------------------|------------|------|----------|---------|---|-------|-----|------|------|
| cumulative_pipeline  |ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|true or false|artic or swift_v2|se or pe|all or positive integer|ignored|ignored|true or false|
| pipeline             |ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|true or false|artic or swift_v2|se or pe|all or positive integer|ignored|ignored|true or false|
| variants_qc          |ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|true or false|artic or swift_v2|se or pe|all or positive integer|ignored|ignored|true or false|
| variants             |ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|true or false|artic or swift_v2|se or pe|all or positive integer|ignored|ignored|true or false|
| sample               |ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|true or false|artic or swift_v2|se or pe|all or positive integer|e.g. SEARCH-17043__D101859__L01__210409_A00953_0272_AH57WJDRXY__S82_L001|e.g. 2021-04-15_16-13-59|true or false|
| qc                   |ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|ignored|ignored|se or pe|ignored|ignored|e.g. 2021-04-15_16-13-59|true or false|
| lineages             |ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|ignored|ignored|ignored|ignored|ignored|ignored|true or false|
| phylogeny            |ucsd or helix|e.g 210409_A00953_0272_AH57WJDRXY|ignored|ignored|ignored|ignored|ignored|ignored|true or false|
| cumulative_lineages  |ucsd or helix|ignored|ignored|ignored|ignored|ignored|ignored|ignored|true or false|
| cumulative_phylogeny |ucsd or helix|ignored|ignored|ignored|ignored|ignored|ignored|ignored|true or false|

An example input file might look like:

```
function,organization,seq_run,merge_lanes,primer_set,fq,read_cap,sample,timestamp,istest
cumulative_pipeline,ucsd,210608_A00953_0321_BH7L5LDSX2,false,swift_v2,pe,2000000,NA,NA,false
```
