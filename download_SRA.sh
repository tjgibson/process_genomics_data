#! /bin/bash

# get parameters from input --------------------------------------------------------------
sra_id=$1
out_dir=$2

# uncompress software and data -----------------------------------------------------------
echo "started input data decompression"
date


tar -xzf sratoolkit.2.9.6-1-centos_linux64.tar.gz
rm sratoolkit.2.9.6-1-centos_linux64.tar.gz

echo "finished input data decompression"
date

# trim reads ------------------------------------------------------------------------------
echo "started SRA download"
date

#./sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump ${sra_id}
./sratoolkit.2.9.6-1-centos_linux64/bin/prefetch -O ./ ${sra_id}
./sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump ${sra_id}.sra

# get exit status
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
  echo "SRA download succeeded"
else
  echo "SRA download failed"
  exit 1
fi		


echo "finished SRA download"
date

# transfer output files ------------------------------------------------------------------
echo "started transferring output files"
date

gzip *.fastq

mv ./*.fastq.gz ${out_dir}
# get exit status
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
  echo "output file transfer succeeded"
  exit 0
else
  echo "output file transfer failed" >&2
  exit 1
fi		


echo "finished transferring output files"
date
