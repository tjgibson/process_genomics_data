# process_genomics_data

## Overview
This repository contains a collection of scripts used for processing data from various genomics experiments on a computing cluster using the HTCondor job manager. There are currently pipelines for processing RNA-seq, ChIP-seq and ATAC-seq data. These pipelines are intented to be used with the HTCondor DAGMan tool to submit and manage multistep workflows. There are also scripts for automated generation of the dag files required b DAGMan. 

## Current pipelines:
### RNA-seq
- Align reads to reference genome using HISAT2
- Filter out multi-mapping reads
- Generate a count table using featureCounts
- Generate bigWig files usind DeepTools

### ChIP-seq
- Align reads to reference genome using bowtie2
- Filter out multi-mapping reads
- Call peaks using MACS2
- Generate bigWig files usind DeepTools

### ATAC-seq
- Trim Tn5 adapters from reads using NGMerge
- Align reads to reference genome using bowtie2
- Filter out multi-mapping reads
- Separate aligned reads into <100bp (accessible) and >100bp (nucleosomal)
- call peaks on accessible fragments using MACS2
- Generate bigWig files usind DeepTools
