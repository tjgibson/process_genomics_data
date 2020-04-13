#! /bin/bash

# get parameters from input --------------------------------------------------------------
bam=$1
in_dir=$2
out_dir=$3
bn=$4

# transfer files -------------------------------------------------------------------------
echo "started file transfer"
date

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


gunzip sambamba-0.7.0-linux-static.gz
chmod u+x sambamba-0.7.0-linux-static
mv sambamba-0.7.0-linux-static sambamba


echo "finished file decompression"
date

# split reads based on fragment size -----------------------------------------------------
echo "started fragment separation"
date


mv ${bn}_filtered/*.bam* .
rmdir ${bn}_filtered/

 # split accessible fragments
./sambamba view -F "template_length < 100 and template_length > -100" -f bam -t 8 -o ${bn}_accessible.bam ${bn}.bam

# get exit status for file transfer
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "splitting accessible fragments succeeded"
else
	echo "splitting accessible fragments failed" >&2
	exit 1
fi

# split nucleosomal fragments
./sambamba view -F "not (template_length < 180 and template_length > -180)" -f bam -t 8 -o ${bn}_nucleosomal.bam ${bn}.bam

# get exit status for file transfer
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "splitting nucleosomal fragments succeeded"
else
	echo "splitting nucleosomal fragments failed" >&2
	exit 1
fi


# sort bam files with split fragments

# sort accessible fragments
./sambamba sort -m 12GB -t 8 -o ${bn}_sorted_accessible.bam ${bn}_accessible.bam
# get exit status for file transfer
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "sorting accessible fragments succeeded"
else
	echo "sorting accessible fragments failed" >&2
	exit 1
fi


# sort nucleosomal fragments
./sambamba sort -m 12GB -t 8 -o ${bn}_sorted_nucleosomal.bam ${bn}_nucleosomal.bam
# get exit status for file transfer
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "sorting nucleosomal fragments succeeded"
else
	echo "sorting nucleosomal fragments failed" >&2
	exit 1
fi


# rename output files
rm *.bai
rm ${bn}_accessible.bam ${bn}_nucleosomal.bam

mv ${bn}.bam ${bn}_total.bam
mv ${bn}_sorted_accessible.bam ${bn}_accessible.bam
mv ${bn}_sorted_nucleosomal.bam ${bn}_nucleosomal.bam

# index BAM files for split fragments
./sambamba index ${bn}_total.bam
./sambamba index ${bn}_accessible.bam
./sambamba index ${bn}_nucleosomal.bam



echo "finished fragment separation"
date


# transfer output files ------------------------------------------------------------------
echo "started transferring output files"
date

mkdir ${bn}_split/
mv *.bam* ./${bn}_split/

tar -czf ./${bn}_split.tar.gz ./${bn}_split/


#rm ${in_dir}${bam}

mv ./${bn}_split.tar.gz ${out_dir}
# get exit status
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "output file transfer succeeded"
else
	echo "output file transfer failed" >&2
	exit 1
fi

echo "finished transferring output files"
date

# clean up unneeded files ----------------------------------------------------------------
rm ./sambamba

# remove intermediate file from previous step --------------------------------------------
if [ $exit_status -eq 0 ] ; then
	echo "removing filtered bam files"
	rm ${in_dir}${bam}
fi




