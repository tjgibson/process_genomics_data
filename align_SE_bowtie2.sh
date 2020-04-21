#! /bin/bash

# get parameters from input --------------------------------------------------------------
r1=$1
in_dir=$2
out_dir=$3
bn=$4
ref_dir=$5
ref_genome=$6

# transfer files -------------------------------------------------------------------------
echo "started file transfer"
date

cp ${in_dir}${r1} .
cp ${ref_dir}${ref_genome} .

echo "finished file transfer"
date

# uncompress software and data -----------------------------------------------------------
echo "started file decompression"
date

tar -xzf bowtie2-2.3.5-linux-x86_64.tar.gz
rm bowtie2-2.3.5-linux-x86_64.tar.gz

tar -xzf ./${ref_genome}
rm ./${ref_genome}

echo "finished file decompression"
date
# align reads ------------------------------------------------------------------------------
echo "started alignment"
date

# run alignment
ref_name=${ref_genome%%.tar.gz}
./bowtie2-2.3.5-linux-x86_64/bowtie2 -p 8 -k 2 --very-sensitive -x ${ref_name}/${ref_name}  -U ./${r1} -S ./${bn}.sam --un-gz ./${bn}_un.fastq.gz

# get alignment exit status
exit_status=$?

# check exit status and exit if alignment failed
if [ $exit_status -eq 0 ] ; then
  echo "bowtie2 succeeded"
else
  echo "bowtie2 failed"
  rm ./samtools
  rm ./*.sam
  rm ./*.fastq.gz
  exit 1
fi		


echo "finished alignment"
date

# remove input files ---------------------------------------------------------------------
rm ./${r1}
rm ./${r2}

	
# compress, sort and index aligned reads  ------------------------------------------------


echo "started BAM compression of aligned reads"
date

# compress aligned reads
./samtools view -bh -o ./${bn}.bam ./${bn}.sam 

# get compression exit status
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
  echo "aligned read compression succeeded"
else
  echo "aligned read compression failed"
  rm ./*.sam
  rm ./*.bam
  rm ./*fastq.gz
  exit 1
fi

echo "finished BAM compression of aligned reads"
date

# remove uncompresed reads
rm ./${bn}.sam


echo "started sorting and indexing bam files"
date

# sort aligned reads
./samtools sort -o ${bn}_sorted.bam -@ 8 ${bn}.bam

# get exit status
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "read sorting succeeded"
else
	echo "read sorting  failed"
	rm ./*.sam
	rm ./*.bam
	rm ./*fastq.gz
	exit 1
fi

# index sorted reads
./samtools index ${bn}_sorted.bam

rm ${bn}.bam

echo "finished sorting and indexing bam files"
date


# transfer output files ------------------------------------------------------------------
echo "started output file transfer"
date

mkdir ${bn}_aligned

mv ./${bn}_sorted.bam* ./${bn}_aligned/
mv ./*_un.fastq.gz ./${bn}_aligned/

tar -czf ${bn}_aligned.tar.gz ./${bn}_aligned/
mv ${bn}_aligned.tar.gz ${out_dir}

# get exit status for file transfer
exit_status=$?


echo "finished output file transfer"
date

# clean up ------------------------------------------------------------------------------
rm ./samtools



