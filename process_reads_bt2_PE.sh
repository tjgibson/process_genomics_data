#! /bin/bash

## detect max number of threads and set the number of threads
np=$(nproc --all)
((rp = $np / 2))

## extract the basename of the file by taking everything up to file extension
file=$1
basepath=${file%%_1.fq.gz}
tmp=${file##*/}
basename=${tmp%%_1.fq.gz}
echo ${file}
echo ${basename}
sum_file="/mnt/gluster/tjgibson2/read_processing_summaries/${basename}.summary"

## create summary file to write ouput and updates
echo "starting read processing:" > ${sum_file}
echo "sample = ${basename}" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

## transfer raw read file and reference genome file from gluster
echo "starting input file transfer from gluster:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

cp ${basepath}_1.fq.gz ./${basename}_1.fastq.gz
cp ${basepath}_2.fq.gz ./${basename}_2.fastq.gz
cp /mnt/gluster/tjgibson2/genomes/ucsc_dm6_indexed.tar.gz .
tar -xzvf ucsc_dm6_indexed.tar.gz

echo "file transfer from gluster done:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}


## uncompress trimmomatic
tar -xzvf NGmerge.tar.gz


## run trimmomatic
echo "running NGmerge:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

(NGmerge-master/NGmerge -a -e 20 -u 41 -n ${rp} -v -1 ${basename}_1.fastq.gz -2 ${basename}_2.fastq.gz -o ${basename}_trimmed) 2>> ${sum_file}

echo " DONE running NGmerge:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}


## remove NGmerge, unpaired trimmed files, and original untrimmed files
rm ./NGmerge.tar.gz
rm -r ./NGmerge-master/
rm ./${basename}_1.fastq.gz 
rm ./${basename}_2.fastq.gz


## uncompress bowtie2
tar -xzvf bowtie2-2.3.5-linux-x86_64.tar.gz


## run bowtie2
echo "running bowtie2:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

(./bowtie2-2.3.5-linux-x86_64/bowtie2 -p ${rp} -k 2 --very-sensitive --no-mixed --no-discordant -X 5000 -x ucsc_dm6/ucsc_dm6  -1 ./${basename}_trimmed_1.fastq.gz -2 ./${basename}_trimmed_2.fastq.gz -S ./${basename}_al.sam) 2>> ${sum_file}

echo "DONE running bowtie2:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

## remove bowtie2 and raw read file
rm -r ./bowtie2-2.3.5-linux-x86_64
rm ./bowtie2-2.3.5-linux-x86_64.tar.gz
rm ./${basename}_trimmed_1.fastq.gz
rm ./${basename}_trimmed_2.fastq.gz

## compress aligned reads
echo "compressing aligned reads:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

./samtools view -bh -o ./${basename}_al.bam ./${basename}_al.sam 
rm ./${basename}_al.sam

echo "DONE compressing aligned reads:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

## filter reads
echo "filtering aligned reads:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

echo "${basename} input reads:" >> ${sum_file}
./samtools view ${basename}_al.bam | wc -l >> ${sum_file}
./samtools view -h  ${basename}_al.bam | grep -v "XS:i" | ./samtools view -bh -q 30 -o ./${basename}_filtered.bam

echo "DONE filtering aligned reads:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

##Report number of filtered reads
echo "${basename} output reads:" >> ${sum_file}
./samtools view ./${basename}_filtered.bam | wc -l >> ${sum_file}


## sort and index bam file
echo "sorting aligned reads:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

#./samtools sort -o ${basename}_sorted.bam -@ ${rp} ${basename}_filtered.bam
./samtools sort -o ${basename}_sorted.bam -@ 4 ${basename}_filtered.bam

echo "DONE sorting aligned reads:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

echo "indexing aligned reads:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

./samtools index ${basename}_sorted.bam

## move filtered read file to gluster
cp ${basename}_sorted.bam /mnt/gluster/tjgibson2/aligned_reads/
cp ${basename}_sorted.bam.bai /mnt/gluster/tjgibson2/aligned_reads/

echo "DONE indexing aligned reads:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

## transfer python installation from gluster and add to path
cp /mnt/gluster/tjgibson2/software/python3.tar.gz .
tar -xzf python3.tar.gz
export PATH=$(pwd)/python/bin:$PATH
mkdir home
export HOME=$(pwd)/home

## use deeptools to make bigwig file for each bam file
echo "starting bigwig conversion:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

# python ./python/bin/bamCoverage --bam ${basename}_sorted.bam -o ${basename}.bw --binSize 10 -p ${rp}
python ./python/bin/bamCoverage --bam ${basename}_sorted.bam -o ${basename}.bw --binSize 10 -p 2

echo "DONE with bigwig conversion:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

## remove remaining files
rm ./${basename}_sorted.bam
rm ./${basename}_sorted.bam.bai
rm ./python3.tar.gz
rm -r ./python/
rm -r ./home/

## remove unfiltered reads
echo "cleaning up files and transferring output" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}


rm ./${basename}_al.bam 

## move output files to gluster
mv ${basename}.bw /mnt/gluster/tjgibson2/bigwigs/
#mv ./${basename}_un.fastq.gz /mnt/gluster/tjgibson2/aligned_reads/
#mv ./${basename}_filter.summary /mnt/gluster/tjgibson2/read_processing_summaries/
#mv ./${basename}.summary /mnt/gluster/tjgibson2/read_processing_summaries/


## remove unwanted files before exiting
rm ./${basename}_filtered.bam 
rm -r ./ucsc_dm6/
rm ucsc_dm6_indexed.tar.gz
rm samtools

echo "read processing completed!" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}
