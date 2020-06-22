#! /bin/bash

# get parameters from input --------------------------------------------------------------
bam=$1
in_dir=$2
out_dir=$3
bn=$4

# transfer files -------------------------------------------------------------------------
echo "started file transfer"
date

# copy input bam file
cp ${in_dir}${bam} .

# get exit status
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

# uncompress software and data -----------------------------------------------------------
echo "started file decompression"
date

tar -xzf ${bam}
rm ${bam}

gunzip fb_dmel-all-r6.26.gtf.gz
tar - xzf subread-1.6.4-Linux-x86_64.tar.gz

echo "finished file decompression"
date

# generate read count table ------------------------------------------------------------------------------
echo "started read filtering"
date

mv ${bn}_aligned/*.bam* .
rm -r ${bn}_aligned/

subread-1.6.4-Linux-x86_64/bin/featureCounts -a fb_dmel-all-r6.26.gtf  -o ${bn}_count_table.txt ${bn}.bam



# clean up unneeded files ----------------------------------------------------------------
rm *.bam
rm *.bam.bai
rm -r ./subread-1.6.4-Linux-x86_64
rm subread-1.6.4-Linux-x86_64.tar.gz
rm fb_dmel-all-r6.26.gtf*

