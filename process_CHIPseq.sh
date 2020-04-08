#! /bin/bash

# Set options ---------------------------------------------------------------------------
## what is the source of the input data
input="local" # one of local or sra

# what type of experiment generated your data
assay_type="chip_w_input" # one of rna, atac, chip_w_input or chip_no_input
molecule="dna" # one of dna, rna
read_format="se" # one of se, pe

# should read trimming be performed
trim_reads=FALSE # logical

# what reference datasets should be used for analysis
# ref_genome="/mnt/gluster/tjgibson2/genomes/ucsc_dm6_indexed.tar.gz"
# genome_annotation="/mnt/gluster/tjgibson2/genomes/fb_dmel-all-r6.26.gtf.gz"


# what data do you want back
return_bams=TRUE # logical
return_unaligned=FALSE # logical
return_bigwigs=TRUE # logical
call_peaks=TRUE # logical

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
		
		sample_name=${basenames[0]}
		
	else
		sample_name=${basenames[0]}
	fi
	
	echo "sample_name: ${sample_name}"
	
	

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
sum_file="/mnt/gluster/tjgibson2/read_processing_summaries/${sample_name}.summary"


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
for i in ${!files[@]}; do
	f=${files[$i]}
	bp=${basepaths[$i]}
	bn=${basenames[$i]}

	# trim reads if specified in options
	if [[ $trim_reads = TRUE ]] ; then
		
		# use NGmerge to trim paired-end data
		if [ $read_format = "pe" ] ; then

			(./process_reads_software/NGmerge-master/NGmerge -a -e 20 -u 41 -n ${rp} -v -1 ${bn}_1$ext -2 ${bn}_2$ext -o ${bn}_trimmed) 2>> ${sum_file}
		
	
		fi
		
		# use trimmomatic to trim single-end data
		if [ $read_format = "se" ] ; then
			(java -jar ./process_reads_software/Trimmomatic-0.39/trimmomatic-0.39.jar SE ${bn}${ext} ${bn}_trimmed$ext ILLUMINACLIP:./process_reads_software/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:1) 2>> ${sum_file}
		
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
			(./process_reads_software/bowtie2-2.3.5-linux-x86_64/bowtie2 -p ${rp} -k 2 --very-sensitive --no-mixed --no-discordant -X 5000 -x ucsc_dm6/ucsc_dm6  -1 ./${bn}_trimmed_1${ext} -2 ./${bn}_trimmed_2${ext} -S ./${bn}.sam --un-gz ./${bn}_un.fastq.gz) 2>> ${sum_file}

		fi
		
		# 
		if [ $read_format = "se" ] ; then
		
			(./process_reads_software/bowtie2-2.3.5-linux-x86_64/bowtie2 -p ${rp} -k 2 --very-sensitive -x ucsc_dm6/ucsc_dm6  -U ./${bn}_trimmed$ext -S ./${bn}.sam --un-gz ./${bn}_un.fastq.gz) 2>> ${sum_file}


			fi
	fi
	
	# compress aligned reads
	./process_reads_software/samtools view -bh -o ./${bn}.bam ./${bn}.sam 
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

		(python ./python/bin/bamCoverage --bam ${bn}_sorted.bam -o ${bn}.bw --binSize 10 -p $rp)  2>> ${sum_file}

		echo "DONE with bigwig conversion:" >> ${sum_file}
		date >> ${sum_file}
		echo "" >> ${sum_file}

	done
fi

# remove python3 installation
rm python3.tar.gz
rm -rf python/
rm -rf home/


# peak calling --------------------------------------------------------------------------
## copy python2.7 (including MACS2) from gluster
cp /mnt/gluster/tjgibson2/software/python2.tar.gz .
tar -xzf python2.tar.gz
export PATH=$(pwd)/python/bin:$PATH
mkdir home
export HOME=$(pwd)/home


if [[ $call_peaks = TRUE ]] ; then
	echo "starting peak calling:" >> ${sum_file}
	
	sample_name=(${basenames[0]%%_input*})
	
	if [ $assay_type = "chip_w_input" ] ; then
		(./python2.7/bin/python ./python2.7/bin/macs2 callpeak -t ${sample_name}*IP*_sorted.bam -c ${sample_group}*input*_sorted.bam -n ${sample_name} --outdir ./${sample_name}_MACS2_output -f BAM -g 1.2e8 --call-summits) 2>> ${sum_file}
	fi
	
	echo "peak calling done:" >> ${sum_file}
fi

# create output to return ---------------------------------------------------------------
sample_name=(${basenames[0]})

mkdir ${sample_name}_out

mv *_fastqc.html ${sample_name}_out

mv ${sum_file} ${sample_name}_out

if [[ $return_bams = TRUE ]] ; then
	mv *_sorted.bam ./${sample_name}_out
fi

if [[ $return_unaligned = TRUE ]] ; then
	mv *_un.fastq.gz ./${sample_name}_out
fi

if [[ $return_bigwigs = TRUE ]] ; then
	mv *.bw ./${sample_name}_out
fi

if [[ $call_peaks = TRUE ]] ; then
	mv ./${sample_name}_MACS2_output/  ./${sample_name}_out
fi


# clean up files- this will only remove files and not directories
rm *

# compress output directory
tar -czf ${sample_name}_out.tar.gz ${sample_name}_out/

# transfer output files
if [[ $transfer_to_gluster = TRUE ]] ; then
	mv ${sample_name}_out.tar.gz /mnt/gluster/tjgibson2/aligned_reads/
fi

# if transfer_to_gluster is set to FALSE, then the only single file remaining in the current directory will be the compressed output file, which will be transferred

