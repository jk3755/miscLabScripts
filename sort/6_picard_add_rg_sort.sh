#!/bin/bash
#
# adds RG tags and coordinate sort the input bam/sam files
# RGID is a unique id - <experiment>.<sample>.<lane>
# RGLB is library
# RGPL is the platform - illumina
# RGPU is the instrument number - NS500277 
# RGSM is sample for marking replicates (similar to RGLB)
#
# specify paths
picard_dir=/home/ubuntu/picard/
sam_dir=/home/ubuntu/atac/h508/5_aligned_hg38_sam/
out_dir=/home/ubuntu/atac/h508/6_aligned_hg38_sam_rg_sorted/
#
#
# sample 1
java -Xmx5g -jar ${picard_dir}picard.jar AddOrReplaceReadGroups \
I=${sam_dir}s1_l1_hg38.sam \
O=${out_dir}s1_l1_hg38_rg_sorted.sam \
SORT_ORDER=coordinate \
RGID=h5308.sample1.lane1 \
RGLB=lib1 \
RGPL=illumina \
RGPU=NS500277 \
RGSM=sample1 &
#
java -Xmx5g -jar ${picard_dir}picard.jar AddOrReplaceReadGroups \
I=${sam_dir}s1_l2_hg38.sam \
O=${out_dir}s1_l2_hg38_rg_sorted.sam \
SORT_ORDER=coordinate \
RGID=h5308.sample1.lane2 \
RGLB=lib1 \
RGPL=illumina \
RGPU=NS500277 \
RGSM=sample1 &
#
java -Xmx5g -jar ${picard_dir}picard.jar AddOrReplaceReadGroups \
I=${sam_dir}s1_l3_hg38.sam \
O=${out_dir}s1_l3_hg38_rg_sorted.sam \
SORT_ORDER=coordinate \
RGID=h5308.sample1.lane3 \
RGLB=lib1 \
RGPL=illumina \
RGPU=NS500277 \
RGSM=sample1 &
#
java -Xmx5g -jar ${picard_dir}picard.jar AddOrReplaceReadGroups \
I=${sam_dir}s1_l4_hg38.sam \
O=${out_dir}s1_l4_hg38_rg_sorted.sam \
SORT_ORDER=coordinate \
RGID=h5308.sample1.lane4 \
RGLB=lib1 \
RGPL=illumina \
RGPU=NS500277 \
RGSM=sample1 &
#
#
#
# sample 2
java -Xmx5g -jar ${picard_dir}picard.jar AddOrReplaceReadGroups \
I=${sam_dir}s2_l1_hg38.sam \
O=${out_dir}s2_l1_hg38_rg_sorted.sam \
SORT_ORDER=coordinate \
RGID=h5308.sample2.lane1 \
RGLB=lib2 \
RGPL=illumina \
RGPU=NS500277 \
RGSM=sample2 &
#
java -Xmx5g -jar ${picard_dir}picard.jar AddOrReplaceReadGroups \
I=${sam_dir}s2_l2_hg38.sam \
O=${out_dir}s2_l2_hg38_rg_sorted.sam \
SORT_ORDER=coordinate \
RGID=h5308.sample2.lane2 \
RGLB=lib2 \
RGPL=illumina \
RGPU=NS500277 \
RGSM=sample2 &
#
java -Xmx5g -jar ${picard_dir}picard.jar AddOrReplaceReadGroups \
I=${sam_dir}s2_l3_hg38.sam \
O=${out_dir}s2_l3_hg38_rg_sorted.sam \
SORT_ORDER=coordinate \
RGID=h5308.sample2.lane3 \
RGLB=lib2 \
RGPL=illumina \
RGPU=NS500277 \
RGSM=sample2 &
#
java -Xmx5g -jar ${picard_dir}picard.jar AddOrReplaceReadGroups \
I=${sam_dir}s2_l4_hg38.sam \
O=${out_dir}s2_l4_hg38_rg_sorted.sam \
SORT_ORDER=coordinate \
RGID=h5308.sample2.lane4 \
RGLB=lib2 \
RGPL=illumina \
RGPU=NS500277 \
RGSM=sample2 &
#
#
#
# sample 3
java -Xmx5g -jar ${picard_dir}picard.jar AddOrReplaceReadGroups \
I=${sam_dir}s3_l1_hg38.sam \
O=${out_dir}s3_l1_hg38_rg_sorted.sam \
SORT_ORDER=coordinate \
RGID=h5308.sample3.lane1 \
RGLB=lib3 \
RGPL=illumina \
RGPU=NS500277 \
RGSM=sample3 &
#
java -Xmx5g -jar ${picard_dir}picard.jar AddOrReplaceReadGroups \
I=${sam_dir}s3_l2_hg38.sam \
O=${out_dir}s3_l2_hg38_rg_sorted.sam \
SORT_ORDER=coordinate \
RGID=h5308.sample3.lane2 \
RGLB=lib3 \
RGPL=illumina \
RGPU=NS500277 \
RGSM=sample3 &
#
java -Xmx5g -jar ${picard_dir}picard.jar AddOrReplaceReadGroups \
I=${sam_dir}s3_l3_hg38.sam \
O=${out_dir}s3_l3_hg38_rg_sorted.sam \
SORT_ORDER=coordinate \
RGID=h5308.sample3.lane3 \
RGLB=lib3 \
RGPL=illumina \
RGPU=NS500277 \
RGSM=sample3 &
#
java -Xmx5g -jar ${picard_dir}picard.jar AddOrReplaceReadGroups \
I=${sam_dir}s3_l4_hg38.sam \
O=${out_dir}s3_l4_hg38_rg_sorted.sam \
SORT_ORDER=coordinate \
RGID=h5308.sample3.lane4 \
RGLB=lib3 \
RGPL=illumina \
RGPU=NS500277 \
RGSM=sample3 &