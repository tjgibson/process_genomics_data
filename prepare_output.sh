#! /bin/bash

# create directory for output data
mkdir dag_output
mkdir dag_output/bigwigs/
mkdir dag_output/peaks/
mkdir dag_output/summary_files/
mkdir dag_output/scripts/
mkdir dag_output/input_files/
mkdir dag_output/condor_dag_files/


# move relevant files to output directory
mv *bw.tar.gz dag_output/bigwigs/

mv *macs2.tar.gz dag_output/peaks/

mv *.stats dag_output/summary_files/

mv *[0-9].err dag_output/summary_files/
mv *[0-9].log dag_output/summary_files/
mv *[0-9].out dag_output/summary_files/

mv *.sh dag_output/scripts/
mv *.sub dag_output/scripts/

mv *.txt dag_output/input_files/

mv *.dag  dag_output/condor_dag_files/
mv dag.dot  dag_output/condor_dag_files/

# compress output directory
tar -czf dag_output.tar.gz dag_output/
rm -r dag_output/