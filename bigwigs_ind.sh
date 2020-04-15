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

# use deeptools to make bigwig file for each bam file ------------------------------------
echo "started making bigwigs"
date


mv ./${bn}_filtered/*.bam* .
rmdir ./${bn}_filtered/

# make bigwig file for accessible fragments
python ./python/bin/bamCoverage --bam ${bn}.bam -o ${bn}.bw --binSize 10 -p 8

# get exit status for making bigwigs
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "making  bigwig succeeded"
else
	echo "making  bigwig failed" >&2
	exit 1
fi


# compress output
tar -czf ${bn}.bw.tar.gz ${bn}.bw



echo "finished making bigwigs"
date

# clean up unneeded files ----------------------------------------------------------------
rm -r ./python/
rm -r ./home/
rm *.bw
rm *.bam
rm *.bai
