#! /bin/bash

# get parameters from input --------------------------------------------------------------
input=$1
IP=$2
in_dir=$3
out_dir=$4
ref_dir=$5
input_bn=$6
IP_bn=$7

# transfer files -------------------------------------------------------------------------
echo "started file transfer"
date

cp ${ref_dir}python38.tar.gz .
cp ${in_dir}${input} .

# get exit status for file transfer
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "input file transfer succeeded"
else
	echo "input file transfer failed" >&2
	exit 1
fi

cp ${in_dir}${IP} .

# get exit status for file transfer
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "input file transfer succeeded"
else
	echo "input file transfer failed" >&2
	exit 1
fi

echo "finished file transfer"
date

# decompress software and data -----------------------------------------------------------
echo "started file decompression"
date

tar -xzf ${input}
rm ${input}

tar -xzf ${IP}
rm ${IP}

mv *_filtered/* . 


tar -xzf python38.tar.gz
rm python38.tar.gz
export PATH=$(pwd)/python/bin:$PATH
mkdir home
export HOME=$(pwd)/home


echo "finished file decompression"
date

# call peaks using macs2 -----------------------------------------------------------------
echo "started peak calling"
date


echo "running macs2 for sample: $bn"
mkdir ./${IP_bn}_macs2

python/bin/python python/bin/macs2 callpeak -t ${IP_bn}.bam -c ${input_bn}.bam -n ${IP_bn}_macs2 --outdir ./${IP_bn}_macs2 -f BAM -g 1.2e8 --call-summits 

# get exit status for peak calling
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "peak calling succeeded"
else
	echo "peak calling failed" >&2
	exit 1
fi


echo "finished peak calling"
date

# compress output
tar -czf ${bn}_macs2.tar.gz ./${bn}_macs2


# clean up unneeded files ----------------------------------------------------------------
rm -r ./python/
rm -r ./home/
rm *.bam*
