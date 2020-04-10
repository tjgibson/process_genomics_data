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

echo "finished file decompression"
date

# filter reads ------------------------------------------------------------------------------
echo "started read filtering"
date


mv ${bn}_aligned/*.bam* .
rm -r ${bn}_aligned/

 # get alignment stats for aligned reads without filtering
./samtools idxstats ${bn}_sorted.bam > ${bn}_all.stats


# filter out multi-mapping reads and get alignment stats for filtered reads
echo "filtering multi-mapping reads"
date

./samtools view -h  ${bn}_sorted.bam | grep -v "XS:i" | ./samtools view -bh -q 30 -o ./${bn}_tmp.bam
./samtools index ./${bn}_tmp.bam
./samtools idxstats ./${bn}_tmp.bam > ${bn}_filtered.stats

# remove input file
rm ./${bn}_sorted.bam*

# remove PCR and optical duplicates and get stats for dedplicated samples
echo "filtering duplicates"
date

java -jar ./picard.jar MarkDuplicates \
      I=./${bn}_tmp.bam \
      O=./${bn}_rmdup.bam \
      M=./${bn}_dup_metrics.txt \
      REMOVE_DUPLICATES=true

# get exit status of duplicate removal
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "duplicate removal succeeded"
else
	echo "duplicate removal failed" >&2
	rm ./${bn}_tmp.bam*
	exit 1
fi

# index reads with duplicates removed and get stats
./samtools index ./${bn}_rmdup.bam
./samtools idxstats ./${bn}_rmdup.bam > ${bn}_rmdup.stats

# remove input file
rm ./${bn}_tmp.bam*

# filter reads to include only reads mapping to major chromosomes (excluding contigs and mitochondrial genome)
echo "filtering reads mapping to contigs and chrM"
date
./samtools view -bh ./${bn}_rmdup.bam -L good_chroms.bed -o ${bn}_filtered.bam 

# remove input file
rm ./${bn}_rmdup.bam*

## sort and index bam file
echo "sorting filtered reads"
date
./samtools sort -o ./${bn}.bam -@ 8 ${bn}_filtered.bam 

# get exit status of sorting
exit_status=$?

# check exit status and exit if failed
if [ $exit_status -eq 0 ] ; then
	echo "read sorting succeeded"
else
	echo "read sorting failed" >&2
	rm ${bn}_filtered.bam*
	exit 1
fi


# index sorted reads
./samtools index ./${bn}.bam

rm ${bn}_filtered.bam*

echo "finished read filtering"
date


# transfer output files ------------------------------------------------------------------
echo "started transferring output files"
date

mkdir ${bn}_filtered/
mv ${bn}.bam* ./${bn}_filtered/

tar -czf ./${bn}_filtered.tar.gz ./${bn}_filtered/


rm ${in_dir}${bam}

mv ./${bn}_filtered.tar.gz ${out_dir}

echo "finished transferring output files"
date

# clean up unneeded files ----------------------------------------------------------------
rm ./samtools
rm ./picard.jar
rm ./good_chroms.bed

