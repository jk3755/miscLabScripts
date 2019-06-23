#!/bin/bash
#
# IMPORTANT!!! REMEMBER TO REDIRECT nohup OUTPUT TO REPORT FILE WHEN LAUNCHING FOR ALIGNMENT STATS
#
# specify paths
hg38ref=/home/ubuntu/ref_genomes/hg38/hg38
fastq_dir=/home/ubuntu/atac/h508/4_good_fastq/
out_dir=/home/ubuntu/atac/h508/5_aligned_hg38_sam/
#
#
# align sample 1
bowtie2 -q -p 20 -X2000 -x $hg38ref -1 ${fastq_dir}s1_l1_r1.good.fastq -2 ${fastq_dir}s1_l1_r2.good.fastq -S ${out_dir}s1_l1_hg38.sam
bowtie2 -q -p 20 -X2000 -x $hg38ref -1 ${fastq_dir}s1_l2_r1.good.fastq -2 ${fastq_dir}s1_l2_r2.good.fastq -S ${out_dir}s1_l2_hg38.sam
bowtie2 -q -p 20 -X2000 -x $hg38ref -1 ${fastq_dir}s1_l3_r1.good.fastq -2 ${fastq_dir}s1_l3_r2.good.fastq -S ${out_dir}s1_l3_hg38.sam
bowtie2 -q -p 20 -X2000 -x $hg38ref -1 ${fastq_dir}s1_l4_r1.good.fastq -2 ${fastq_dir}s1_l4_r2.good.fastq -S ${out_dir}s1_l4_hg38.sam
#
#
# align sample 2
bowtie2 -q -p 20 -X2000 -x $hg38ref -1 ${fastq_dir}s2_l1_r1.good.fastq -2 ${fastq_dir}s2_l1_r2.good.fastq -S ${out_dir}s2_l1_hg38.sam
bowtie2 -q -p 20 -X2000 -x $hg38ref -1 ${fastq_dir}s2_l2_r1.good.fastq -2 ${fastq_dir}s2_l2_r2.good.fastq -S ${out_dir}s2_l2_hg38.sam
bowtie2 -q -p 20 -X2000 -x $hg38ref -1 ${fastq_dir}s2_l3_r1.good.fastq -2 ${fastq_dir}s2_l3_r2.good.fastq -S ${out_dir}s2_l3_hg38.sam
bowtie2 -q -p 20 -X2000 -x $hg38ref -1 ${fastq_dir}s2_l4_r1.good.fastq -2 ${fastq_dir}s2_l4_r2.good.fastq -S ${out_dir}s2_l4_hg38.sam
#
#
# align sample 3
bowtie2 -q -p 20 -X2000 -x $hg38ref -1 ${fastq_dir}s3_l1_r1.good.fastq -2 ${fastq_dir}s3_l1_r2.good.fastq -S ${out_dir}s3_l1_hg38.sam
bowtie2 -q -p 20 -X2000 -x $hg38ref -1 ${fastq_dir}s3_l2_r1.good.fastq -2 ${fastq_dir}s3_l2_r2.good.fastq -S ${out_dir}s3_l2_hg38.sam
bowtie2 -q -p 20 -X2000 -x $hg38ref -1 ${fastq_dir}s3_l3_r1.good.fastq -2 ${fastq_dir}s3_l3_r2.good.fastq -S ${out_dir}s3_l3_hg38.sam
bowtie2 -q -p 20 -X2000 -x $hg38ref -1 ${fastq_dir}s3_l4_r1.good.fastq -2 ${fastq_dir}s3_l4_r2.good.fastq -S ${out_dir}s3_l4_hg38.sam