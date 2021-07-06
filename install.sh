#!/bin/bash

# tested on clean ubuntu instance: ami-08962a4068733a2b6

# Provide 1. directory to install software, 2. Directory of anaconda installation ending in anaconda3
# allow user to specify where to install some software
SOFTWAREDIR=$1
ANACONDADIR=$2

# ------ this should happen pre- this script ------
# install anaconda and python 3.8 using commands below:
# wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
# bash Anaconda3-2020.11-Linux-x86_64.sh
# conda init bash # (and restart shell)
# ------------------------------------

# ------ Install basic tools ------
sudo apt-get update
sudo apt-get install python-dev-is-python3 python3-dev \
     build-essential libssl-dev libffi-dev \
     libxml2-dev libxslt1-dev zlib1g-dev
sudo apt-get install unzip
sudo apt-get install autoconf
sudo apt-get install autotools-dev
sudo apt-get install git

# ------- Make software directory ------
mkdir -p $SOFTWAREDIR

# ------- covid_sequencing_analysis_pipeline -------
git clone https://github.com/ucsd-ccbb/covid_sequencing_analysis_pipeline.git $SOFTWAREDIR/covid_sequencing_analysis_pipeline

# ------- q30 -------
git clone https://github.com/artnasamran/q30.git $SOFTWAREDIR/q30
chmod u+x SOFTWAREDIR/q30/q30.py

# ------- samhead -------
git clone https://github.com/niemasd/SD-COVID-Sequencing.git $SOFTWAREDIR/SD-COVID-Sequencing
cd $SOFTWAREDIR/SD-COVID-Sequencing/samhead
make
cd ~

# ------- covid1.1 env from pre-release README on github -------
conda create --name covid1.1
source $ANACONDADIR/bin/activate covid1.1
conda install numpy
conda install boto3
#conda install -c bioconda fastqc
conda install -c bioconda qualimap
conda install -c bioconda minimap2
conda install -c bioconda samtools

pip install multiqc
pip install nwalign3
pip install pandas
pip install seaborn
pip install fastaparser

# ------- iVar v.1.3.1 -------
source $ANACONDADIR/bin/activate covid1.1
cd $SOFTWAREDIR
wget https://github.com/andersen-lab/ivar/archive/master.zip -O ivar.zip
unzip ivar.zip
mv ivar-master ivar
rm ivar.zip

# Path to HTSlib from conda (installed by samtools)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ANACONDADIR/envs/covid1.1/lib

# Install iVar
cd ivar
./autogen.sh
# configure needs to point to location of htslib, LDFLAGS needs to point to location of environment libs
./configure --with-hts=$ANACONDADIR/envs/covid1.1 LDFLAGS=-Wl,-R$ANACONDADIR/envs/covid1.1/lib --prefix=$SOFTWAREDIR/ivar
make
make install

# ------- Pangolin -------
cd $SOFTWAREDIR
git clone https://github.com/cov-lineages/pangolin.git
cd $SOFTWAREDIR/pangolin
conda env create -f environment.yml
conda activate pangolin
pip install .
cd ~
source $ANACONDADIR/bin/activate base

# ------- ViralMSA (instructions from Niema) -------
wget "https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py"
chmod a+x ViralMSA.py
sudo mv ViralMSA.py /usr/local/bin/ViralMSA.py # optional step to install globally

# ------- biopython -------
conda install -y -c anaconda biopython

# ------- iqtree version 2.1.2 -------
cd $SOFTWAREDIR
mkdir IQTree
cd IQTree
wget https://github.com/iqtree/iqtree2/releases/download/v2.1.2/iqtree-2.1.2-Linux.tar.gz
tar -xvf iqtree-2.1.2-Linux.tar.gz
rm iqtree-2.1.2-Linux.tar.gz

# add to path
echo "export PATH=$PATH:$SOFTWAREDIR/IQTree/iqtree-2.1.2-Linux/bin"

cd ~

# ------- FastRoot: (from Niema's script: https://github.com/niemasd/ViReport/blob/master/Dockerfile#L47) -------
wget -q "https://github.com/uym2/MinVar-Rooting/archive/master.zip"
unzip -q master.zip
sudo mv MinVar-Rooting-master /usr/local/bin/MinVar-Rooting
sudo ln -s /usr/local/bin/MinVar-Rooting/FastRoot.py /usr/local/bin/FastRoot.py
rm -rf MinVar-Rooting-master master.zip
cd /usr/local/bin/MinVar-Rooting/
sudo python3 setup.py install

cd ~

# ------- empress v1.1.0 -------
wget https://data.qiime2.org/distro/core/qiime2-2020.11-py36-linux-conda.yml
conda env create -n qiime2-2020.11 --file qiime2-2020.11-py36-linux-conda.yml
source $ANACONDADIR/bin/activate qiime2-2020.11

pip install empress

source $ANACONDADIR/bin/activate base
