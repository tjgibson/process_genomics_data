#! /bin/bash

# get parameters -------------------------------------------------------------------------
in_files=$1
replicates_file=$2
in_dir=$3
out_dir=$4
ref_dir=$5
ref_genome=$6

# set names of read processing scripts ---------------------------------------------------
trim_name="trim_PE_NGmerge"
align_name="align_PE_bowtie2"
filter_name="filter_reads"
split_name="atac_split_fragments"
atac_ind_peaks="atac_peaks_macs2_ind"
atac_ind_bigwigs="atac_bigwigs_ind"
atac_merged_peaks="atac_peaks_macs2_merged"
atac_merged_bigwigs="atac_bigwigs_merged"



# iterate over input files and create dag nodes for parallel steps --------------------
# create dag file
touch atac.dag

# initiate array for parent nodes of non-parallel steps
split_nodes=()

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
	bn=${line%%_S[0-9]_L???_R?_???*$ext}
	pair_base=${line%%_R?_???*$ext}
	r1=${pair_base}_R1_001
	r2=${pair_base}_R2_001
	
	# add trim lines to dag file
	job_name=trim_${i}
	echo "# trim reads" >> atac.dag
	echo "JOB ${job_name} ${trim_name}.sub" >> atac.dag
	echo "VARS ${job_name} r1=\"${r1}$ext\"" >> atac.dag
	echo "VARS ${job_name} r2=\"${r2}$ext\"" >> atac.dag
	echo "VARS ${job_name} in_dir=\"${in_dir}\"" >> atac.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> atac.dag
 	echo "VARS ${job_name} bn=\"${bn}\"" >> atac.dag
 
	# add align lines to dag file
	job_name=align_${i}
	echo "# align reads" >> atac.dag
	echo "JOB ${job_name} ${align_name}.sub" >> atac.dag
	echo "VARS ${job_name} r1=\"${bn}_trimmed_1.fastq.gz\"" >> atac.dag
	echo "VARS ${job_name} r2=\"${bn}_trimmed_2.fastq.gz\"" >> atac.dag
	echo "VARS ${job_name} in_dir=\"${out_dir}\"" >> atac.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> atac.dag
 	echo "VARS ${job_name} bn=\"${bn}\"" >> atac.dag
	echo "VARS ${job_name} ref_dir=\"${ref_dir}\"" >> atac.dag
	echo "VARS ${job_name} ref_genome=\"${ref_genome}\"" >> atac.dag

	# add filter lines to dag file
	job_name=filter_${i}
	echo "# filter reads" >> atac.dag
	echo "JOB ${job_name} ${filter_name}.sub" >> atac.dag
	echo "VARS ${job_name} bam=\"${bn}_aligned.tar.gz\"" >> atac.dag
	echo "VARS ${job_name} in_dir=\"${out_dir}\"" >> atac.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> atac.dag
 	echo "VARS ${job_name} bn=\"${bn}\"" >> atac.dag
 	
	# add filter lines to dag file
	job_name=split_${i}
	split_nodes+=($job_name)
	echo "# split reads" >> atac.dag
	echo "JOB ${job_name} ${split_name}.sub" >> atac.dag
	echo "VARS ${job_name} bam=\"${bn}_filtered.tar.gz\"" >> atac.dag
	echo "VARS ${job_name} in_dir=\"${out_dir}\"" >> atac.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> atac.dag
 	echo "VARS ${job_name} bn=\"${bn}\"" >> atac.dag
 
	# add ind_bigwig lines to dag file
	job_name=ind_bigwigs_${i}
	echo "# make bigwigs for individual replicates" >> atac.dag
	echo "JOB ${job_name} ${atac_ind_bigwigs}.sub" >> atac.dag
	echo "VARS ${job_name} bam=\"${bn}_split.tar.gz\"" >> atac.dag
	echo "VARS ${job_name} in_dir=\"${out_dir}\"" >> atac.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> atac.dag
 	echo "VARS ${job_name} bn=\"${bn}\"" >> atac.dag
	echo "VARS ${job_name} ref_dir=\"${ref_dir}\"" >> atac.dag
	
	# add ind_peaks lines to dag file
	job_name=ind_peaks_${i}
	echo "# call peaks for individual replicates" >> atac.dag
	echo "JOB ${job_name} ${atac_ind_peaks}.sub" >> atac.dag
	echo "VARS ${job_name} bam=\"${bn}_split.tar.gz\"" >> atac.dag
	echo "VARS ${job_name} in_dir=\"${out_dir}\"" >> atac.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> atac.dag
 	echo "VARS ${job_name} bn=\"${bn}\"" >> atac.dag
	echo "VARS ${job_name} ref_dir=\"${ref_dir}\"" >> atac.dag
	 
	# define dependencies among nodes
	echo "PARENT trim_${i} CHILD align_${i}" >> atac.dag
	echo "PARENT align_${i} CHILD filter_${i}" >> atac.dag
	echo "PARENT filter_${i} CHILD split_${i}" >> atac.dag
	echo "PARENT split_${i} CHILD ind_bigwigs_${i} ind_peaks_${i}" >> atac.dag
 
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
 	line=$(echo ${line} | sed 's/_S[0-9]_L???_R?_???*$ext/_split.tar.gz/g')
 	line=($line)
 	group_name=${line[0]}
 	rep_fqs=${line[@]:1}
 	
 	echo $line
 	echo $group_name
 	echo $rep_fqs
 	echo $reps
 	
 
	# add merged_bigwig lines to dag file
	job_name=merged_bigwigs_${i}
	dep_nodes+=($job_name)
	echo "# make bigwigs for merged replicates" >> atac.dag
	echo "JOB ${job_name} ${atac_merged_bigwigs}.sub" >> atac.dag
	echo "VARS ${job_name} in_dir=\"${out_dir}\"" >> atac.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> atac.dag
	echo "VARS ${job_name} ref_dir=\"${ref_dir}\"" >> atac.dag
	echo "VARS ${job_name} group_name=\"${group_name}\"" >> atac.dag
	echo "VARS ${job_name} replicates=\"${reps}\"" >> atac.dag
	
	# add ind_peaks lines to dag file
	job_name=merged_peaks_${i}
	dep_nodes+=($job_name)
	echo "# call peaks for merged replicates" >> atac.dag
	echo "JOB ${job_name} ${atac_merged_peaks}.sub" >> atac.dag
	echo "VARS ${job_name} in_dir=\"${out_dir}\"" >> atac.dag
	echo "VARS ${job_name} out_dir=\"${out_dir}\"" >> atac.dag
	echo "VARS ${job_name} ref_dir=\"${ref_dir}\"" >> atac.dag
	echo "VARS ${job_name} group_name=\"${group_name}\"" >> atac.dag
	echo "VARS ${job_name} replicates=\"${reps}\"" >> atac.dag

	 
	# define dependencies among nodes
	
done <$replicates_file
	echo "PARENT ${split_nodes[*]} CHILD ${dep_nodes[*]}" >> atac.dag


# create dot file describing DAG structure
echo "DOT dag.dot UPDATE DONT-OVERWRITE" >> atac.dag