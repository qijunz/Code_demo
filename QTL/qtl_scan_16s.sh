#!/bin/bash

# untar your R installation
tar -xzf R-3.6.3.tar.gz
tar -xzf R-3.6.3-packages.tar.gz

# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH
export RHOME=$(pwd)/R
export R_LIBS=$(pwd)/R-packages

# copy data from gluster
cp /staging/qzhang333/DOdiet_input/final_genoprobs_1176.rds .
cp /staging/qzhang333/DOdiet_input/gigamuga_map_test_v2.RDS .
cp /staging/qzhang333/DOdiet_input/K_1176_DietDO.rds .

cp /staging/qzhang333/DOdiet_input/Warren_MF_DO_Study_SCFA_Rel_Abun_Attie_Phenotype_Metadata_2_16_2024_v2.csv .

# run R, with the name of your  R script
Rscript $1

# move output file to staging
mv qtl_scan_16s_add_20240216.rds /staging/qzhang333/
mv qtl_scan_16s_int_20240216.rds /staging/qzhang333/

# clear data in case return them back to home directory
rm $1
rm *.rds
rm *.RDS
rm *.csv
rm *.tar.gz