#! /bin/bash

# Set options ---------------------------------------------------------------------------
## what is the source of the input data
input="sra" # one of 'local' or 'sra'

# what type of experiment generated your data
read_format="pe" # one of 'se', 'pe'

# should read trimming be performed
trim_reads=FALSE # logical

# what reference datasets should be used for analysis
# ref_genome_file="/mnt/gluster/tjgibson2/genomes/ucsc_dm6_indexed.tar.gz"
# genome_annotation_file="/mnt/gluster/tjgibson2/genomes/fb_dmel-all-r6.26.gtf.gz"

# what data do you want back
return_bams=FALSE # logical
return_unaligned=FALSE # logical
return_bigwigs=FALSE # logical
count_table=TRUE # logical

# should output data be transferred to gluster or output to your home directory on the submit node?
# WARNING: only transfer files back to the submit node if you are only transferring count tables, which are small. For large files such as BAMs, output files should be transferred back to gluster
transfer_to_gluster=TRUE


# setup ---------------------------------------------------------------------------------

# detect max number of threads and set number of threads to use
np=$(nproc --all)
((rp = $np / 2))

# transfer and uncompress necessary software
cp /mnt/gluster/tjgibson2/software/process_reads_software.tar.gz .
cp /mnt/gluster/tjgibson2/software/python3.tar.gz .

tar -xzf  process_reads_software.tar.gz
tar -xzf python3.tar.gz

export PATH=$(pwd)/python/bin:$PATH
mkdir home
export HOME=$(pwd)/home


# transfer genome files from gluster and decompress
cp /mnt/gluster/tjgibson2/genomes/ucsc_dm6_indexed.tar.gz .
cp /mnt/gluster/tjgibson2/genomes/fb_dmel-all-r6.26.gtf.gz .
cp /mnt/gluster/tjgibson2/genomes/good_chroms.bed .

tar -xzf ucsc_dm6_indexed.tar.gz 
gunzip fb_dmel-all-r6.26.gtf.gz

# remove tar files
rm ucsc_dm6_indexed.tar.gz 
rm fb_dmel-all-r6.26.gtf.gz

# file import for SRA files
if [ $input = "sra" ] ; then
	
	# set file extension
	ext=".fastq"
	echo "file extension = ${ext}"
	
	# retrieve files using fasterq-dump from SRA tools
	files=($@)
	basepaths=($@)
	basenames=($@)
	
	
	
	for i in ${!files[@]}; do
		f=${files[$i]}
		
		./process_reads_software/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump ${f} -e ${rp}
	done
fi


# file import for local files
if [ $input = "local" ] ; then

	# check whether file extension is .fq.gz or .fastq.gz
	if [[ $1 == *.fastq.gz ]] ; then
		ext=".fastq.gz"

	elif [[ $1 == *.fq.gz ]] ; then
		ext=".fq.gz"
	else
		echo ".fastq or .fq file extension not detected"
	fi

	echo "file extension = ${ext}"


	# extract the filename and path to the file
	files=($@)
	basepaths=()
	basenames=()
	if [ $read_format = "pe" ] ; then
		for i in ${!files[@]}; do
			f=${files[$i]}
			basepaths[$i]=${f%%_1$ext}
			tmp=${f##*/}
			basenames[$i]=${tmp%%_1$ext}
		done

	elif [ $read_format = "se" ] ; then
		for i in ${!files[@]}; do
			f=${files[$i]}
			basepaths[$i]=${f%%$ext}
			tmp=${f##*/}
			basenames[$i]=${tmp%%$ext}
		done

	else
		echo "invalid read_format: must be either 'se' or 'pe'"
	fi

	echo ${files[*]}
	echo ${basepaths[*]}
	echo ${basenames[*]}
	
	# transfer input files to working directory
	for i in ${!files[@]}; do
		f=${files[$i]}
		bp=${basepaths[$i]}
		bn=${basenames[$i]}

		if [ $read_format = "pe" ] ; then
			cp ${bp}_1${ext} . 
			cp ${bp}_2${ext} . 
		fi
	
		if [ $read_format = "se" ] ; then
			cp ${bp}${ext} . 
		fi
  
	done

else

	echo "invalid input parameter. Input should be set to either 'local' or 'sra' "

fi

# create summary file on gluster
bn=${basenames[0]}
sum_file="/mnt/gluster/tjgibson2/read_processing_summaries/${bn}.summary"


# QC ------------------------------------------------------------------------------------
echo "========================================================================================================================================" > ${sum_file}
echo "starting read processing of ${read_format} RNA-seq data:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

echo "running fastqc" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}


# use fastQC to generate HTML formatted reports for each file
(./process_reads_software/FastQC/fastqc *$ext)  2>> ${sum_file}

echo "fastqc done" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}


# trimming ------------------------------------------------------------------------------
echo "checking if trimming should be performed" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

for i in ${!files[@]}; do
	f=${files[$i]}
	bp=${basepaths[$i]}
	bn=${basenames[$i]}

	# trim reads if specified in options
	if [[ $trim_reads = TRUE ]] ; then
		
		echo "starting trimming" >> ${sum_file}
		date >> ${sum_file}
		echo "" >> ${sum_file}
		
		# use NGmerge to trim paired-end data
		if [ $read_format = "pe" ] ; then
			
			echo "trimming paired-end data with NGmerge" >> ${sum_file}
			date >> ${sum_file}
			echo "" >> ${sum_file}
			
			(./process_reads_software/NGmerge-master/NGmerge -a -e 20 -u 41 -n ${rp} -v -1 ${bn}_1$ext -2 ${bn}_2$ext -o ${bn}_trimmed) 2>> ${sum_file}
		
		fi
		
		# use trimmomatic to trim single-end data
		if [ $read_format = "se" ] ; then
			
			echo "trimming single-end data with Trimmomatic" >> ${sum_file}
			date >> ${sum_file}
			echo "" >> ${sum_file}
			
			(java -jar ./process_reads_software/Trimmomatic-0.39/trimmomatic-0.39.jar SE ${bn}${ext} ${bn}_trimmed$ext ILLUMINACLIP:./process_reads_software/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:1) 2>> ${sum_file}
		
		fi
		
		rm ${bn}${ext}
		
		echo "Trimming done" >> ${sum_file}
		date >> ${sum_file}
		echo "" >> ${sum_file}
		
	# if no trimming specified, rename file to file_trimmed format so file will be recognized in next step
	else
		echo "No trimming performed" >> ${sum_file}
		date >> ${sum_file}
		echo "" >> ${sum_file}
		
		if [ $read_format = "pe" ] ; then
		
			mv ${bn}_1$ext ${bn}_trimmed_1$ext
			mv ${bn}_2$ext ${bn}_trimmed_2$ext
				
		fi
		
		if [ $read_format = "se" ] ; then
		
			mv ${bn}${ext} ${bn}_trimmed${ext}
			
		fi
	fi
done

echo "done" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

# alignment -----------------------------------------------------------------------------	
bn=${basenames[0]}

# use hisat2 for aligning reads from RNA/cDNA libraries

echo "starting alignment with hisat2" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

if [ $read_format = "pe" ] ; then
	echo "aligning paired-end data" >> ${sum_file}
	(./process_reads_software/hisat2-2.1.0/hisat2 -k 2 -p ${rp} --no-mixed --no-discordant -x ./ucsc_dm6/ucsc_dm6 -1 ./*_trimmed_1${ext} -2 ./*_trimmed_2${ext} -S ./${bn}.sam --un-gz ./${bn}_un.fastq.gz) 2>> ${sum_file}

fi

if [ $read_format = "se" ] ; then
	echo "aligning single-end data" >> ${sum_file}
	(./process_reads_software/hisat2-2.1.0/hisat2 -p ${rp} -k 2 -x ./ucsc_dm6/ucsc_dm6 -U ./*_trimmed$ext -S ./${bn}.sam --un-gz ./${bn}_un.fastq.gz) 2>> ${sum_file}
	
fi

echo "alignment done" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

# compress aligned reads
./process_reads_software/samtools view -bh -o ./${bn}.bam ./${bn}.sam 
rm ./${bn}.sam



# filtering -----------------------------------------------------------------------------

## filter reads
echo "filtering aligned reads:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

echo "${bn} input reads:" >> ${sum_file}
./process_reads_software/samtools view ${bn}.bam | wc -l >> ${sum_file}
./process_reads_software/samtools view -h  ${bn}.bam | grep -v "XS:i" | ./process_reads_software/samtools view -bh -q 30 -o ./${bn}_filtered_tmp.bam
./process_reads_software/samtools view -bh ${bn}_filtered_tmp.bam -L good_chroms.bed -o ${bn}_filtered.bam 

##Report number of filtered reads
echo "${bn} output reads:" >> ${sum_file}
./process_reads_software/samtools view ./${bn}_filtered.bam | wc -l >> ${sum_file}


## sort and index bam file

#./samtools sort -o ${basename}_sorted.bam -@ ${rp} ${basename}_filtered.bam
./process_reads_software/samtools sort -o ${bn}_sorted.bam -@ $rp ${bn}_filtered.bam
./process_reads_software/samtools index ${bn}_sorted.bam

# cleanup intermediate files
rm ${bn}.bam 
rm ./${bn}_filtered_tmp.bam
rm ${bn}_filtered.bam

echo "DONE filtering aligned reads:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

# bigwigs -------------------------------------------------------------------------------
## use deeptools to make bigwig file for each bam file

if [[ $return_bigwigs = TRUE ]] ; then
	
	echo "starting bigwig conversion:" >> ${sum_file}
	date >> ${sum_file}
	echo "" >> ${sum_file}

	(python ./python/bin/bamCoverage --bam ${bn}_sorted.bam -o ${bn}.bw --binSize 10 -p $rp)  2>> ${sum_file}

	echo "DONE with bigwig conversion:" >> ${sum_file}
	date >> ${sum_file}
	echo "" >> ${sum_file}
	
fi


# count table ---------------------------------------------------------------------------
if [[ $count_table = TRUE ]] ; then
	(./process_reads_software/subread-1.6.4-Linux-x86_64/bin/featureCounts -a fb_dmel-all-r6.26.gtf  -o ${bn}_count_table.txt ${bn}_sorted.bam)  2>> ${sum_file}
fi


# create output to return ---------------------------------------------------------------
mkdir ${bn}_out

cp ${bn}_fastqc.html ${bn}_out

if [[ $return_bams = TRUE ]] ; then
	cp ./${bn}_sorted.bam ./${bn}_out
fi

if [[ $return_unaligned = TRUE ]] ; then
	cp ./${bn}_un.fastq.gz ./${bn}_out
fi

if [[ $return_bigwigs = TRUE ]] ; then
	cp ./${bn}.bw ./${bn}_out
fi

if [[ $count_table = TRUE ]] ; then
	cp ./${bn}_count_table.txt ./${bn}_out
fi


# clean up files- this will only remove files and not directories
rm *

# compress output directory
tar -czf ${bn}_out.tar.gz ${bn}_out/

# transfer output files
if [[ $transfer_to_gluster = TRUE ]] ; then
	mv ${bn}_out.tar.gz /mnt/gluster/tjgibson2/aligned_reads/
fi

# if transfer_to_gluster is set to FALSE, then the only single file remaining in the current directory will be the compressed output file, which will be transferred
