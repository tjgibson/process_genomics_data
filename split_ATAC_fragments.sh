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

 # get alignment stats for aligned reads without filtering
./sambamba view -F "template_length < 100 and template_length > -100" -f bam -t 8 -o ${bn}_accessible.bam ${bn}.bam
./sambamba view -F "not (template_length < 180 and template_length > -180)" -f bam -t 8 -o ${bn}_nucleosomal.bam ${bn}.bam

./sambamba sort -m 12GB -t 8 -o ${bn}_sorted_accessible.bam ${bn}_accessible.bam
./sambamba sort -m 12GB -t 8 -o ${bn}_sorted_nucleosomal.bam ${bn}_nucleosomal.bam



# rename output files
rm *.bai
rm ${bn}_accessible.bam ${bn}_nucleosomal.bam

mv ${bn}.bam ${bn}_total.bam
mv ${bn}_sorted_accessible.bam ${bn}_accessible.bam
mv ${bn}_sorted_nucleosomal.bam ${bn}_nucleosomal.bam

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

echo "finished transferring output files"
date

# clean up unneeded files ----------------------------------------------------------------
rm ./sambamba

