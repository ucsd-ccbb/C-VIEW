# covid_sequencing_analysis_pipeline
AWS optimized pipeline based on https://github.com/niemasd/SD-COVID-Sequencing 

Pipeline version 0.2.0 is pre-installed on the snap-0635d84491ed5f318 Amazon Web Services snapshot in region us-east-2 (Ohio).  It uses the following external software programs:

* ivar 1.3.1
* minimap2 2.17-r941
* samtools 1.11
* QualiMap v.2.2.2-dev
* FastQC v0.11.9

Should one wish to set up the pipeline on a fresh instance, follow the below commands.
Create a conda environment and activate it, then run:

```
conda install numpy 
conda install boto3
conda install -c bioconda fastqc
conda install -c bioconda qualimap
conda install -c bioconda minimap2
conda install -c bioconda samtools
```

Followed by:
```
pip install multiqc
pip install nwalign3
```

Finally, install ivar from source (see https://github.com/andersen-lab/ivar ).

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
