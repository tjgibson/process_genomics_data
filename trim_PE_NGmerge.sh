#! /bin/bash

# get parameters from input --------------------------------------------------------------
r1=$1
r2=$2
in_dir=$3
out_dir=$4
bn=$5

# transfer files -------------------------------------------------------------------------
echo "started file transfer"
date

cp ${in_dir}${r1} .
cp ${in_dir}${r2} .

echo "finished file transfer"
date

# uncompress software and data -----------------------------------------------------------
echo "started input data decompression"
date


tar -xzf NGmerge.tar.gz
rm NGmerge.tar.gz

echo "finished input data decompression"
date

# trim reads ------------------------------------------------------------------------------
echo "started trimming reads"
date

./NGmerge-master/NGmerge -a -e 20 -u 41 -n 8 -v -1 ./$r1 -2 ./$r2 -o ${bn}_trimmed
# get exit status
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
  echo "NGmerge succeeded"
else
  echo "NGmerge failed"
	rm ./${r1}
	rm ./${r2}
  exit 1
fi		


echo "finished trimming reads"
date
# remove input files ---------------------------------------------------------------------
rm ./${r1}
rm ./${r2}

# transfer output files ------------------------------------------------------------------
echo "started transferring output files"
date


mv ./${bn}_trimmed* ${out_dir}
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

