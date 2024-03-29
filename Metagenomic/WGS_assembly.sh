#!/bin/bash

# unpack pipeline 
tar -xzf SPAdes-3.15.5.tar.gz
tar -xzf Prodigal-2.6.3.tar.gz
tar -xzf python-3.8.tar.gz
tar -xzf python-3.8-packages.tar.gz

# set path
export PATH=$(pwd)/SPAdes-3.15.5/bin:$PATH
export PATH=$(pwd)/Prodigal-2.6.3:$PATH
export PATH=$(pwd)/python-3.8/bin:$PATH
export PYTHONPATH=$(pwd)/python-3.8-packages

# transfer and unpack reads
cp /staging/qzhang333/HMDP_WGS/$1_R1_trimmed_paired_mouseDNAremoved.fastq.gz .
cp /staging/qzhang333/HMDP_WGS/$1_R2_trimmed_paired_mouseDNAremoved.fastq.gz .

gzip -d *fastq.gz

# create output dir 
mkdir assembly_out_$1

# run pipeline
python assembly_pipeline.py -F $1_R1_trimmed_paired_mouseDNAremoved.fastq -R $1_R2_trimmed_paired_mouseDNAremoved.fastq -o assembly_out_$1/ -s $1

# tar output 
tar -czf assembly_out_$1.tar.gz assembly_out_$1

# move output to staging/qzhang333
mv assembly_out_$1.tar.gz /staging/qzhang333

# clear other data
rm *.fastq
rm *.tar.gz