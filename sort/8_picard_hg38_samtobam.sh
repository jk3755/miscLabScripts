#!/bin/bash
#
#
# specify paths
picard_dir=/home/ubuntu/picard/
sam_dir=/home/ubuntu/atac/h508/6_aligned_hg38_sam_rg_sorted/
out_dir=/home/ubuntu/atac/h508/7_aligned_hg38_bam/
#
#
# sample 1
java -Xmx7g -jar ${picard_dir}picard.jar SamFormatConverter \
I=${sam_dir}s1_l1_hg38_rg_sorted.sam \
O=${out_dir}s1_l1_hg38.bam &
#
java -Xmx7g -jar ${picard_dir}picard.jar SamFormatConverter \
I=${sam_dir}s1_l2_hg38_rg_sorted.sam \
O=${out_dir}s1_l2_hg38.bam &
#
java -Xmx7g -jar ${picard_dir}picard.jar SamFormatConverter \
I=${sam_dir}s1_l3_hg38_rg_sorted.sam \
O=${out_dir}s1_l3_hg38.bam &
#
java -Xmx7g -jar ${picard_dir}picard.jar SamFormatConverter \
I=${sam_dir}s1_l4_hg38_rg_sorted.sam \
O=${out_dir}s1_l4_hg38.bam &
#
#
#
# sample 2
java -Xmx7g -jar ${picard_dir}picard.jar SamFormatConverter \
I=${sam_dir}s2_l1_hg38_rg_sorted.sam \
O=${out_dir}s2_l1_hg38.bam &
#
java -Xmx7g -jar ${picard_dir}picard.jar SamFormatConverter \
I=${sam_dir}s2_l2_hg38_rg_sorted.sam \
O=${out_dir}s2_l2_hg38.bam &
#
java -Xmx7g -jar ${picard_dir}picard.jar SamFormatConverter \
I=${sam_dir}s2_l3_hg38_rg_sorted.sam \
O=${out_dir}s2_l3_hg38.bam &
#
java -Xmx7g -jar ${picard_dir}picard.jar SamFormatConverter \
I=${sam_dir}s2_l4_hg38_rg_sorted.sam \
O=${out_dir}s2_l4_hg38.bam &
#
#
# sample 3
java -Xmx7g -jar ${picard_dir}picard.jar SamFormatConverter \
I=${sam_dir}s3_l1_hg38_rg_sorted.sam \
O=${out_dir}s3_l1_hg38.bam &
#
java -Xmx7g -jar ${picard_dir}picard.jar SamFormatConverter \
I=${sam_dir}s3_l2_hg38_rg_sorted.sam \
O=${out_dir}s3_l2_hg38.bam &
#
java -Xmx7g -jar ${picard_dir}picard.jar SamFormatConverter \
I=${sam_dir}s3_l3_hg38_rg_sorted.sam \
O=${out_dir}s3_l3_hg38.bam &
#
java -Xmx7g -jar ${picard_dir}picard.jar SamFormatConverter \
I=${sam_dir}s3_l4_hg38_rg_sorted.sam \
O=${out_dir}s3_l4_hg38.bam &
#