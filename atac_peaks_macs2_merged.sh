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
rm *_nucleosomal.bam*
rm *_total.bam*


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

# merge replicate bam files -----------------------------------------------------------------
echo "started merging replicate bam files"
date

./samtools merge -@ 4  ${group_name}_merged.bam *_accessible.bam
./samtools index ${group_name}_merged.bam
rm ./*_accessible.bam*

echo "finished merging replicate bam files"
date


# call peaks using macs2 -----------------------------------------------------------------
echo "started peak calling"
date


echo "running macs2 based on fragment size pileup of merged replicates for group: $group_name"
mkdir ./${group_name}_macs2

python/bin/python python/bin/macs2 callpeak -t ${group_name}_merged.bam -n ${group_name}_macs2_fragments --outdir ./${group_name}_macs2 -f BAMPE --keep-dup all -g 1.2e8 --call-summits 

# get exit status for peak calling
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo " fragment peak calling succeeded"
else
	echo "fragment peak calling failed" >&2
	exit 1
fi

echo "running macs2 based on cut sites of merged replicates for group: $group_name"
python/bin/python python/bin/macs2 callpeak -t ${group_name}_merged.bam -n ${group_name}_macs2_cutsites --outdir ./${group_name}_macs2 -f BAM --keep-dup all --nomodel --shift -100 --extsize 200 -g 1.2e8 --call-summits 

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
tar -czf ${group_name}_macs2.tar.gz ./${group_name}_macs2


# clean up unneeded files ----------------------------------------------------------------
rm ./samtools
rm -r ./python/
rm -r ./home/
rm *.bam*

