#! /bin/bash

# get parameters -------------------------------------------------------------------------
in_files=$1
replicates_file=$2
chip_w_input=$3
chip_wo_input=$4
in_dir=/staging/tjgibson2/input_data/
out_dir=/staging/tjgibson2/output_data/
ref_dir=/staging/tjgibson2/reference_data/
ref_genome=ucsc_dm6.tar.gz

# set names of read processing scripts ---------------------------------------------------
align_name="align_SE_bowtie2"
filter_name="filter_reads"
chip_ind_peaks="peaks_macs2_ind"
ind_bigwigs="bigwigs_ind"
merged_bigwigs="bigwigs_merged"
peaks_w_input="peaks_w_input_macs2_ind"
peaks_wo_input="peaks_wo_input_macs2_ind"



# iterate over input files and create dag nodes for parallel steps --------------------
# create dag file
rm chip.dag
touch chip.dag

# initiate array for parent nodes of non-parallel steps
filter_nodes=()

# initiate line counter
i=0

# iterate over input files
while read line; do
 	((i=i+1))
 	echo $i
 	echo $line
 	
 	# check whether file extension is .fq.gz or .fastq.gz
	if [[ $line == *.fastq.gz ]] ; then
		ext=".fastq.gz"

	elif [[ $line == *.fq.gz ]] ; then
		ext=".fq.gz"
	else
		echo ".fastq or .fq file extension not detected"
	fi

	echo "file extension = ${ext}"
	
	# get basename for file
	bn=${line%%$ext}

	 
	# add align lines to dag file
	job_name=align_${i}
	echo "# align reads" >> chip.dag
	echo "JOB ${job_name} ${align_name}.sub" >> chip.dag
	echo "VARS ${job_name} r1=\"${bn}${ext}\"" >> chip.dag
	echo "VARS ${job_name} in_dir=\"${in_dir}\"" >> chip.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> chip.dag
 	echo "VARS ${job_name} bn=\"${bn}\"" >> chip.dag
	echo "VARS ${job_name} ref_dir=\"${ref_dir}\"" >> chip.dag
	echo "VARS ${job_name} ref_genome=\"${ref_genome}\"" >> chip.dag

	# add filter lines to dag file
	job_name=filter_${i}
	filter_nodes+=($job_name)
	echo "# filter reads" >> chip.dag
	echo "JOB ${job_name} ${filter_name}.sub" >> chip.dag
	echo "VARS ${job_name} bam=\"${bn}_aligned.tar.gz\"" >> chip.dag
	echo "VARS ${job_name} in_dir=\"${out_dir}\"" >> chip.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> chip.dag
 	echo "VARS ${job_name} bn=\"${bn}\"" >> chip.dag
 	
 
	# add ind_bigwig lines to dag file
	job_name=ind_bigwigs_${i}
	echo "# make bigwigs for individual replicates" >> chip.dag
	echo "JOB ${job_name} ${ind_bigwigs}.sub" >> chip.dag
	echo "VARS ${job_name} bam=\"${bn}_filtered.tar.gz\"" >> chip.dag
	echo "VARS ${job_name} in_dir=\"${out_dir}\"" >> chip.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> chip.dag
 	echo "VARS ${job_name} bn=\"${bn}\"" >> chip.dag
	echo "VARS ${job_name} ref_dir=\"${ref_dir}\"" >> chip.dag
	
	 
	# define dependencies among nodes
	echo "PARENT align_${i} CHILD filter_${i}" >> chip.dag
	echo "PARENT filter_${i} CHILD ind_bigwigs_${i}" >> chip.dag
 
done <$in_files

# iterate over replicate file and create dag lines for dependent steps -------------------
# create array to hold non-parallel node names
dep_nodes=()

# initiate line counter
i=0

# iterate over input files
while read line; do
 	((i=i+1))
 	echo $line
 	line=$(echo ${line} | sed 's/.fastq.gz/_filtered.tar.gz/g')
 	line=($line)
 	group_name=${line[0]}
 	rep_fqs=${line[@]:1}
 	reps=$rep_fqs
 	
 	echo $line
 	echo $group_name
 	echo $rep_fqs
 	echo $reps
 	
 
	# add merged_bigwig lines to dag file
	job_name=merged_bigwigs_${i}
	dep_nodes+=($job_name)
	echo "# make bigwigs for merged replicates" >> chip.dag
	echo "JOB ${job_name} ${merged_bigwigs}.sub" >> chip.dag
	echo "VARS ${job_name} in_dir=\"${out_dir}\"" >> chip.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> chip.dag
	echo "VARS ${job_name} ref_dir=\"${ref_dir}\"" >> chip.dag
	echo "VARS ${job_name} group_name=\"${group_name}\"" >> chip.dag
	echo "VARS ${job_name} replicates=\"${reps}\"" >> chip.dag
	
	 
	# define dependencies among nodes
	
done <$replicates_file

# initiate line counter
i=0

# iterate over input files
while read line; do
 	((i=i+1))
 	echo $line
 	line=$(echo ${line} | sed 's/.fastq.gz/_filtered.tar.gz/g')
 	line=($line)
 	
 	echo ${line[*]}
 	
 	chip_input_file=${line[0]}
 	chip_input_bn=${chip_input_file%%_filtered.tar.gz}
 	
 	chip_IP_file=${line[1]}
 	chip_IP_bn=${chip_IP_file%%_filtered.tar.gz}
 	
 
	# add merged_bigwig lines to dag file
	job_name=peaks_w_input_${i}
	dep_nodes+=($job_name)
	echo "# call peaks for samples with input" >> chip.dag
	echo "JOB ${job_name} ${peaks_w_input}.sub" >> chip.dag
	echo "VARS ${job_name} chip_input=\"${chip_input_file}\"" >> chip.dag
	echo "VARS ${job_name} chip_IP=\"${chip_IP_file}\"" >> chip.dag
	echo "VARS ${job_name} in_dir=\"${out_dir}\"" >> chip.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> chip.dag
	echo "VARS ${job_name} ref_dir=\"${ref_dir}\"" >> chip.dag
	echo "VARS ${job_name} chip_input_bn=\"${chip_input_bn}\"" >> chip.dag
	echo "VARS ${job_name} chip_IP_bn=\"${chip_IP_bn}\"" >> chip.dag
	
	 
	# define dependencies among nodes
	
done <$chip_w_input

# initiate line counter
i=0

# iterate over input files
while read line; do
 	((i=i+1))
 	echo $i
 	echo $line
 	
 	# check whether file extension is .fq.gz or .fastq.gz
	if [[ $line == *.fastq.gz ]] ; then
		ext=".fastq.gz"

	elif [[ $line == *.fq.gz ]] ; then
		ext=".fq.gz"
	else
		echo ".fastq or .fq file extension not detected"
	fi

	echo "file extension = ${ext}"
	
	# get basename for file
	bn=${line%%$ext}
 	
 
	# add merged_bigwig lines to dag file
	job_name=peaks_wo_input_${i}
	dep_nodes+=($job_name)
	echo "# call peaks for samples without input" >> chip.dag
	echo "JOB ${job_name} ${peaks_wo_input}.sub" >> chip.dag
	echo "VARS ${job_name} IP=\"${bn}_filtered.tar.gz\"" >> chip.dag
	echo "VARS ${job_name} in_dir=\"${out_dir}\"" >> chip.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> chip.dag
	echo "VARS ${job_name} ref_dir=\"${ref_dir}\"" >> chip.dag
	echo "VARS ${job_name} IP_bn=\"${bn}\"" >> chip.dag
	
	 
	# define dependencies among nodes
	
done <$chip_wo_input
	echo "PARENT ${filter_nodes[*]} CHILD ${dep_nodes[*]}" >> chip.dag

# create dot file describing DAG structure
echo "DOT dag.dot UPDATE" >> chip.dag