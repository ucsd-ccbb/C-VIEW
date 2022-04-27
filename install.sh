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

# ------- cview -------
git clone https://github.com/ucsd-ccbb/C-VIEW.git $SOFTWAREDIR/cview

# ------- q30 -------
git clone https://github.com/artnasamran/q30.git $SOFTWAREDIR/q30
chmod u+x SOFTWAREDIR/q30/q30.py

# ------- samhead -------
git clone https://github.com/niemasd/SD-COVID-Sequencing.git $SOFTWAREDIR/SD-COVID-Sequencing
cd $SOFTWAREDIR/SD-COVID-Sequencing/samhead
make
cd ~

# ------- pi metrics -----------
wget "https://raw.githubusercontent.com/Niema-Docker/pi_from_pileup/main/pi_from_pileup.cpp"
g++ -O3 --std=c++11 -o /shared/workspace/software/pi_from_pileup pi_from_pileup.cpp
rm pi_from_pileup.cpp

# ------- cview env from pre-release README on github -------
conda create --name cview
source $ANACONDADIR/bin/activate cview
conda install numpy
conda install boto3
conda install -c bioconda qualimap
conda install -c bioconda minimap2
conda install -c bioconda samtools

pip install multiqc
pip install nwalign3
pip install pandas
pip install seaborn
pip install fastaparser

# ------- iVar v.1.3.1 -------
source $ANACONDADIR/bin/activate cview
cd $SOFTWAREDIR
wget https://github.com/andersen-lab/ivar/archive/master.zip -O ivar.zip
unzip ivar.zip
mv ivar-master ivar
rm ivar.zip

# Path to HTSlib from conda (installed by samtools)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ANACONDADIR/envs/cview/lib

# Install iVar
cd ivar
./autogen.sh
# configure needs to point to location of htslib, LDFLAGS needs to point to location of environment libs
./configure --with-hts=$ANACONDADIR/envs/cview LDFLAGS=-Wl,-R$ANACONDADIR/envs/cview/lib --prefix=$SOFTWAREDIR/ivar
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

cd ~

source $ANACONDADIR/bin/activate base
