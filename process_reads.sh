#! /bin/bash

# Set options ---------------------------------------------------------------------------
## what is the source of the input data
input="local" # one of local or sra

# what type of experiment generated your data
assay_type="chip_w_input" # one of rna, atac, chip_w_input or chip_no_input
molecule="dna" # one of dna, rna
read_format="pe" # one of se, pe

# should read trimming be performed
trim_reads=TRUE # logical

# what reference datasets should be used for analysis
ref_genome="/mnt/gluster/tjgibson2/genomes/ucsc_dm6_indexed.tar.gz"
genome_annotation="/mnt/gluster/tjgibson2/genomes/fb_dmel-all-r6.26.gtf.gz"


# what data do you want back
return_bams=FALSE # logical
return_unaligned=FALSE # logical
return_bigwigs=TRUE # logical
call_peaks=TRUE # logical
count_table=TRUE # logical


# setup ---------------------------------------------------------------------------------

# detect max number of threads and set number of threads to use
# np=$(nproc --all)
# ((rp = $np / 2))

# file import for SRA files
if [ input = "sra" ] ; then
	
	# set file extension
	ext=".fastq"
	echo "file extension = ${ext}"
	
	# retrieve files using fasterq-dump from SRA tools
	files=($@)
	basepaths=($@)
	basenames=($@)
	
	for i in ${!files[@]}; do
		f=${files[$i]}
		
		fasterq-dump ${f} -e ${rp}
	done
	
else

	echo "invalid input parameter. Input should be set to either 'local' or 'sra' "

fi


# file import for local files
if [ input = "local" ] ; then

	# check whether file extension is .fq.gz or .fastq.gz
	if [[ $1 == *.fastq.gz ]] ; then
		ext=".fastq.gz"

	elif [[ $1 == *.fq.gz ]] ; then
		ext=".fq.gz"
	else
		echo "fastq file extension not detected"
	fi

	echo "file extension = ${ext}"


	# extract the filename and path to the file
	if [ $read_format = "se" ] ; then
		files=($1)
		basepaths=(${files[0]%%$ext})
		tmp=${files[0]##*/}
		basenames=(${tmp%%$ext})
		echo ${files[*]}
		echo ${basepaths[*]}
		echo ${basenames[*]}
  
	elif [ $read_format = "pe" ] ; then  
		files=($1)
		basepaths=(${files[0]%%_1$ext})
		tmp=${files[0]##*/}
		basenames=(${tmp%%_1$ext})
		echo ${files[*]}
		echo ${basepaths[*]}
		echo ${basenames[*]}

	else
		echo "invalid read_format: must be either 'se' or 'pe'"

	fi

	if [ $assay_type = "chip_w_input" ] ; then
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
	fi

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
# QC ------------------------------------------------------------------------------------

# use fastQC to generate HTML formatted reports for each file
./FastQC/fastqc *$ext


# trimming ------------------------------------------------------------------------------
for i in ${!files[@]}; do
	f=${files[$i]}
	bp=${basepaths[$i]}
	bn=${basenames[$i]}

	# trim reads if specified in options
	if [[ $trim_reads = TRUE ]] ; then
		
		# use NGmerge to trim paired-end data
		if [ $read_format = "pe" ] ; then

			(NGmerge-master/NGmerge -a -e 20 -u 41 -n ${rp} -v -1 ${bn}_1$ext -2 ${bn}_2$ext -o ${bn}_trimmed) 2>> ${sum_file}

	
		fi
		
		# use trimmomatic to trim single-end data
		if [ $read_format = "se" ] ; then
			(java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE ${bn}${ext} ${bn}_trimmed$ext ILLUMINACLIP:./Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:1) 2>> ${sum_file}

		fi
	# if no trimming specified, rename file to file_trimmed format so file will be recognized in next step
	else
		if [ $read_format = "pe" ] ; then
			mv ${bn}_1$ext ${bn}_trimmed_1$ext
			mv ${bn}_2$ext ${bn}_trimmed_2$ext	
		fi
		
		if [ $read_format = "se" ] ; then
			mv ${bn}${ext} ${bn}_trimmed${ext}
		fi
	
done


# alignment -----------------------------------------------------------------------------
for i in ${!files[@]}; do
	f=${files[$i]}
	bp=${basepaths[$i]}
	bn=${basenames[$i]}

	# use bowtie for aligning reads from DNA libraries
	if [ $molecule = "dna" ] ; then
		if [ $read_format = "pe" ] ; then
			(./bowtie2-2.3.5-linux-x86_64/bowtie2 -p ${rp} -k 2 --very-sensitive --no-mixed --no-discordant -X 5000 -x ucsc_dm6/ucsc_dm6  -1 ./${bn}_trimmed_1${ext} -2 ./${bn}_trimmed_2${ext} -S ./${bn}.sam --un-gz ./${bn}_un.fastq.gz) 2>> ${sum_file}

		fi
		
		# 
		if [ $read_format = "se" ] ; then
		
			(./bowtie2-2.3.5-linux-x86_64/bowtie2 -p ${rp} -k 2 --very-sensitive -x ucsc_dm6/ucsc_dm6  -U ./${bn}_trimmed$ext -S ./${bn}.sam --un-gz ./${bn}_un.fastq.gz) 2>> ${sum_file}


			fi
	fi
	
	# use hisat2 for aligning reads from RNA/cDNA libraries
	if [ $molecule = "rna" ] ; then
		if [ $read_format = "pe" ] ; then
			(./hisat2-2.1.0/hisat2 -k 2 -p ${rp} --no-mixed --no-discordant -x ./ucsc_dm6/ucsc_dm6 -1 ./${bn}_trimmed_1${ext} -2 ./${bn}_trimmed_2${ext} -S ./${bn}.sam --un-gz ./${bn}_un.fastq.gz) 2>> ${sum_file}

		fi
	
		# 
		if [ $read_format = "se" ] ; then
			(./hisat2-2.1.0/hisat2 -p ${rp} -k 2 -x ./ucsc_dm6/ucsc_dm6 -U ./${bn}_trimmed$ext -S ./${bn}.sam --un-gz ./${bn}_un.fastq.gz) 2>> ${sum_file}
			
			fi
	fi
	# compress aligned reads
	./samtools view -bh -o ./${bn}.bam ./${bn}.sam 
	rm ./${bn}.sam
done

# filtering -----------------------------------------------------------------------------

## filter reads
echo "filtering aligned reads:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

for i in ${!files[@]}; do
	f=${files[$i]}
	bp=${basepaths[$i]}
	bn=${basenames[$i]}

	echo "${bn} input reads:" >> ${sum_file}
	./samtools view ${bn}.bam | wc -l >> ${sum_file}
	./samtools view -h  ${bn}.bam | grep -v "XS:i" | ./samtools view -bh -q 30 -o ./${bn}_filtered_tmp.bam
	./samtools view -bh ${bn}_filtered_tmp.bam -L good_chroms.bed -o ${bn}_filtered.bam 

	##Report number of filtered reads
	echo "${bn} output reads:" >> ${sum_file}
	./samtools view ./${bn}_filtered.bam | wc -l >> ${sum_file}


	## sort and index bam file

	#./samtools sort -o ${basename}_sorted.bam -@ ${rp} ${basename}_filtered.bam
	./samtools sort -o ${bn}_sorted.bam -@ $rp ${bn}_filtered.bam
	./samtools index ${bn}_sorted.bam
done

echo "DONE filtering aligned reads:" >> ${sum_file}
date >> ${sum_file}
echo "" >> ${sum_file}

# bigwigs -------------------------------------------------------------------------------
## use deeptools to make bigwig file for each bam file

if [[ $return_bigwigs = TRUE ]] ; then
	for i in ${!files[@]}; do
		f=${files[$i]}
		bp=${basepaths[$i]}
		bn=${basenames[$i]}

		echo "starting bigwig conversion:" >> ${sum_file}
		date >> ${sum_file}
		echo "" >> ${sum_file}

	
		python ./python/bin/bamCoverage --bam ${bn}_sorted.bam -o ${bn}.bw --binSize 10 -p $rp

		echo "DONE with bigwig conversion:" >> ${sum_file}
		date >> ${sum_file}
		echo "" >> ${sum_file}
	done
fi

# peak calling --------------------------------------------------------------------------
if [[ $call_peaks = TRUE ]] ; then
	if [ $assay_type = "chip_w_input" ] ; then
		./python2.7/bin/python ./python2.7/bin/macs2 callpeak -t ${sample_group}*IP*_filtered_sorted_good_chroms.bam -c ${sample_group}*Input*_filtered_sorted_good_chroms.bam -n ${sample_group}_filtered --outdir ./${sample_group}_MACS2_output -f BAM -g 1.2e8 --call-summits 

fi

# count table ---------------------------------------------------------------------------
if [[ $count_table = TRUE ]] ; then
	
	# for experiments involving peak calling (ChIP or ATAC), create count table based on peaks
	if [[ $call_peaks = TRUE ]] ; then
		./subread-1.6.4-Linux-x86_64/bin/featureCounts -a peaks -o TG_count_table.txt *.bam
	
	else
		./subread-1.6.4-Linux-x86_64/bin/featureCounts -a ${gtf} -o TG_count_table.txt *.bam
	fi
fi

