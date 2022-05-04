#!/bin/bash

# ------ BEFORE this script ------
# create a clean t2.medium ubuntu 20.04 instance
# [NB: it MUST be version 20.04, NOT the latest version (e.g. 22.04)
# because 20.04 is the latest version supported by AWS ParallelCluster]
# with a 35 GB root drive and a 300 GB data drive,
# and security group allowing:
# SSH	via TCP	on port 22
# all traffic via all protocols on all ports
#
# ssh onto new instance to set up file system and mount:
# lsblk
# for the remainder, assume lsblk shows that the name of the
# 300 gb volume is xvdb:
# sudo file -s /dev/xvdb
# sudo mkfs -t xfs /dev/xvdb
# sudo mkdir /shared
# sudo mount /dev/xvdb /shared
# sudo chown `whoami` /shared
#
# install anaconda and python:
# cd /shared
# wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
# bash Anaconda3-2020.11-Linux-x86_64.sh
# (NB: install anaconda into /shared/workspace/software/anaconda3
# and say "yes" to having the installer initialize anaconda
# by running conda init)

# THEN EXIT THE TERMINAL
# and re-enter it
# ------------------------------------

# cd /shared
# then run this script:
# bash install.sh

# ------ AFTER this script ------
# make a snapshot of the 300 gb volume.
# use cview parallelcluster template to make
# a new cluster based on this snapshot.
# ssh into head node of cluster and run
# aws configure
# to set up aws access
# ------------------------------------

SOFTWAREDIR=/shared/workspace/software
ANACONDADIR=$SOFTWAREDIR/anaconda3

# ------ Install basic tools ------
sudo apt-get update
sudo apt-get install python-dev-is-python3 python3-dev \
     build-essential libssl-dev libffi-dev \
     libxml2-dev libxslt1-dev zlib1g-dev
sudo apt-get install unzip
sudo apt-get install autoconf
sudo apt-get install autotools-dev
sudo apt-get install git
sudo apt-get install libncurses5
sudo apt-get install lib32z1
sudo apt install awscli

# ------- Make log and software directories ------
mkdir -p /shared/logs
mkdir -p /shared/runfiles
mkdir -p $SOFTWAREDIR

# ------- cview -------
git clone https://github.com/ucsd-ccbb/C-VIEW.git $SOFTWAREDIR/cview

# ------- q30 -------
git clone https://github.com/artnasamran/q30.git $SOFTWAREDIR/q30
chmod u+x $SOFTWAREDIR/q30/q30.py

# ------- samhead -------
git clone https://github.com/niemasd/samhead.git $SOFTWAREDIR/samhead
cd $SOFTWAREDIR/samhead
make
cd ~

# ------- pi metrics -----------
wget "https://raw.githubusercontent.com/Niema-Docker/pi_from_pileup/main/pi_from_pileup.cpp"
g++ -O3 --std=c++11 -o $SOFTWAREDIR/pi_from_pileup pi_from_pileup.cpp
rm pi_from_pileup.cpp

source $ANACONDADIR/bin/deactivate

# ------- iVar -------
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create --name ivar
source $ANACONDADIR/bin/activate ivar
conda install ivar
source $ANACONDADIR/bin/deactivate

# ------- cview env (not: not C-VIEW itself, which was cloned above) -------
conda env create -f $SOFTWAREDIR/cview/cview_environment.yml
source $ANACONDADIR/bin/activate cview
pip install nwalign3
source $ANACONDADIR/bin/deactivate

# ------- Pangolin -------
git clone https://github.com/cov-lineages/pangolin.git $SOFTWAREDIR/pangolin

cd $SOFTWAREDIR/pangolin
conda env create -f environment.yml

source $ANACONDADIR/bin/activate pangolin
pip install .
source $ANACONDADIR/bin/deactivate

# ------- ViralMSA (instructions from Niema) -------
mkdir $SOFTWAREDIR/viralMSA
cd $SOFTWAREDIR/viralMSA

conda create --name biopy
source $ANACONDADIR/bin/activate biopy
conda install -y -c anaconda biopython

wget "https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py"
chmod a+x ViralMSA.py

source $ANACONDADIR/bin/deactivate
