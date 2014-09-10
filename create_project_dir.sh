#! /bin/tcsh

# Create directory for project
echo "Usage: create_project_dir.sh project user group"

set proj=$1 # project name is first argument
set user=$2 # user is second argument
set grp=$3 # group is third argument

mkdir $proj
chgrp $grp $proj

cd $proj
mkdir gds
mkdir gds/samples
mkdir gds/subjects

mkdir plink
mkdir plink/samples
mkdir plink/subjects

mkdir sample_snp_annot

mkdir accessory_files
mkdir accessory_files/CC
mkdir accessory_files/SI

mkdir dbGaP
mkdir dbGaP/SI
mkdir dbGaP/To_dbGaP

mkdir imputation

mkdir chrom_anomalies
mkdir chrom_anomalies/data
mkdir chrom_anomalies/plots
mkdir chrom_anomalies/plots/QC

mkdir $user
mkdir $user/R
mkdir $user/results
