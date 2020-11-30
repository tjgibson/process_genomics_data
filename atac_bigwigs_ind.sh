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


mv ./${bn}_split/*.bam* .
rmdir ./${bn}_split/

# make bigwig file for accessible fragments
python ./python/bin/bamCoverage --bam ${bn}_accessible.bam -o ${bn}_accessble.bw --binSize 10 -p 8
# python ./python/bin/bamCoverage --bam ${bn}_accessible.bam -o ${bn}_accessble_rpkm.bw --binSize 10 -p 8  --normalizeUsing RPKM
python ./python/bin/bamCoverage --bam ${bn}_accessible.bam -o ${bn}_accessble_cpm.bw --binSize 10 -p 8  --normalizeUsing CPM


# get exit status for making bigwigs
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "making accessible bigwig succeeded"
else
	echo "making accessible bigwig failed" >&2
	exit 1
fi

# make bigwig file for nucleosomal fragments
python ./python/bin/bamCoverage --bam ${bn}_nucleosomal.bam -o ${bn}_nucleosomal.bw --binSize 10 -p 8
# python ./python/bin/bamCoverage --bam ${bn}_nucleosomal.bam -o ${bn}_nucleosomal_rpkm.bw --binSize 10 -p 8  --normalizeUsing RPKM
python ./python/bin/bamCoverage --bam ${bn}_nucleosomal.bam -o ${bn}_nucleosomal_cpm.bw --binSize 10 -p 8  --normalizeUsing CPM

# get exit status for making bigwigs
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "making nucleosomal bigwig succeeded"
else
	echo "making nucleosomal bigwig failed" >&2
	exit 1
fi

# make bigwig file for all fragments
python ./python/bin/bamCoverage --bam ${bn}_total.bam -o ${bn}_total.bw --binSize 10 -p 8
# python ./python/bin/bamCoverage --bam ${bn}_total.bam -o ${bn}_total_rpkm.bw --binSize 10 -p 8  --normalizeUsing RPKM
python ./python/bin/bamCoverage --bam ${bn}_total.bam -o ${bn}_total_cpm.bw --binSize 10 -p 8  --normalizeUsing CPM

# get exit status for making bigwigs
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "making total bigwig succeeded"
else
	echo "making total bigwig failed" >&2
	exit 1
fi

# compress output
tar -czf ${bn}.bw.tar.gz *.bw



echo "finished making bigwigs"
date

# clean up unneeded files ----------------------------------------------------------------
rm -r ./python/
rm -r ./home/
rm *.bw
rm *.bam
rm *.bai
