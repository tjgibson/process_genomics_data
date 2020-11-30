#! /bin/bash

# get parameters from input --------------------------------------------------------------
file=$1
SAF_file=$2
in_dir=$3


## uncompress software
tar -xzf subread-1.6.4-Linux-x86_64.tar.gz

## extract the basename of the file by taking everything up to file extension
basename=${file%%_filtered.tar.gz}
echo ${file}
echo ${basename}

## transfer input file from staging
cp ${in_dir}${file} .
tar -xzf ${basename}_filtered.tar.gz
rm ${basename}_filtered.tar.gz

mv *_filtered/* .

## create read count table
./subread-1.6.4-Linux-x86_64/bin/featureCounts -O -f -F 'SAF' -a *.saf -o ${basename}_count_table.txt *.bam


## remove unnecessary files
rm *.bam
rm *.bam.bai
rm -r ./subread-1.6.4-Linux-x86_64
rm subread-1.6.4-Linux-x86_64.tar.gz
rm *.saf
