#! /bin/bash

# get parameters from input --------------------------------------------------------------
r1=$1
r2=$2
in_dir=$3
out_dir=$4
bn=$5

# transfer files -------------------------------------------------------------------------
cp ${in_dir}${r1} .
cp ${in_dir}${r2} .

# uncompress software and data -----------------------------------------------------------
tar -xzf NGmerge.tar.gz
rm NGmerge.tar.gz

# trim reads ------------------------------------------------------------------------------
./NGmerge-master/NGmerge -a -e 20 -u 41 -n 8 -v -1 ./$r1 -2 ./$r2 -o ${bn}_trimmed
exit_status=$?

# remove input files ---------------------------------------------------------------------
rm ./${r1}
rm ./${r2}

# transfer output files ------------------------------------------------------------------
mv ./${bn}_trimmed* ${out_dir}

# set exit status ------------------------------------------------------------------------
if [ $exit_status -eq 0 ] ; then
  echo "NGmerge succeeded"
  exit 0
else
  echo "NGmerge failed"
  exit 1
fi		
