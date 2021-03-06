#! /bin/bash

# get parameters from input --------------------------------------------------------------
in_dir=$1
out_dir=$2
ref_dir=$3
group_name=$4
replicates=$5

# transfer files -------------------------------------------------------------------------
echo "started file transfer"
date

# transfer python
cp ${ref_dir}python38.tar.gz .

# transfer all bam files for biological replicates
for r in ${replicates[@]}; do
echo "started copying $r to wd"
date

cp ${in_dir}$r .
# get exit status for file transfer
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "input file transfer succeeded for file $r"
else
	echo "input file transfer failed for file $r" >&2
	ls ${in_dir}$r
	exit 1
fi

echo "finished copying $r to wd"
date


# decompress file and remove tar.gz file
echo "started decompressing $r"
date

tar -xzf $r
rm $r

echo "finished decompressing $r"
date


# get directory name and transfer bam files to wd
r_dir=${r%%.tar.gz}
mv ${r_dir}/* .
rmdir ${r_dir}/

done

echo "finished file transfer"
date

# decompress software and data -----------------------------------------------------------
echo "started file decompression"
date

tar -xzf python38.tar.gz
rm python38.tar.gz
export PATH=$(pwd)/python/bin:$PATH
mkdir home
export HOME=$(pwd)/home


echo "finished file decompression"
date

# merge replicate bam files --------------------------------------------------------------
echo "started merging replicate bam files"
date

ls *.bam
./samtools merge -@ 4  ${group_name}_merged.bam *.bam
./samtools index ${group_name}_merged.bam

echo "finished merging replicate bam files"
date


# create bigwigs using deeptools ---------------------------------------------------------
echo "started making bigwigs"
date

# make bigwigs from accessible fragments
python/bin/python python/bin/bamCoverage --bam ${group_name}_merged.bam -o ${group_name}.bw --binSize 10 -p 8
# python/bin/python python/bin/bamCoverage --bam ${group_name}_merged.bam -o ${group_name}_rpkm.bw --binSize 10 -p 8  --normalizeUsing RPKM
# python/bin/python python/bin/bamCoverage --bam ${group_name}_merged.bam -o ${group_name}_cpm.bw --binSize 10 -p 8  --normalizeUsing CPM

# get exit status for peak calling
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "making merged bigwigs succeeded"
else
	echo "making merged bigwigs failed" >&2
	exit 1
fi


## compress bigwigs output file
tar -czf ${group_name}_merged.bw.tar.gz *.bw



# clean up unneeded files ----------------------------------------------------------------
rm ./samtools
rm -r ./python/
rm -r ./home/
rm *.bw
rm *.bam*

