#! /bin/bash

# get parameters from input --------------------------------------------------------------
bam=$1
in_dir=$2
out_dir=$3
bn=$4
ref_dir=$5

# transfer files -------------------------------------------------------------------------
echo "started file transfer"
date

cp ${ref_dir}python38.tar.gz .
cp ${in_dir}${bam} .

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

tar -xzf ${bam}
rm ${bam}

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


mv ./${bn}_split/*.bam* .
rmdir ./${bn}_split/

echo "running macs2 based on fragment size pileup for sample: $bn"
mkdir ./${bn}_macs2

python/bin/python python/bin/macs2 callpeak -t ${bn}_accessible.bam -n ${bn}_macs2_fragments --outdir ./${bn}_macs2 -f BAMPE --keep-dup all -g 1.2e8 --call-summits 

# get exit status for peak calling
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo " fragment peak calling succeeded"
else
	echo "fragment peak calling failed" >&2
	exit 1
fi

echo "running macs2 based on cut sites for sample: $bn"
python/bin/python python/bin/macs2 callpeak -t ${bn}_accessible.bam -n ${bn}_macs2_cutsites --outdir ./${bn}_macs2 -f BAM --keep-dup all --nomodel --shift -100 --extsize 200 -g 1.2e8 --call-summits 

# get exit status for peak calling
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "cut site peak calling succeeded"
else
	echo "cut site peak calling failed" >&2
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

