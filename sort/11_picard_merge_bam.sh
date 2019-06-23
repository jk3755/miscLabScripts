#!/bin/bash
#
# specify paths
picard_dir=/home/ubuntu/picard/
bam_dir=/home/ubuntu/atac/h508/7_aligned_hg38_bam/
out_dir=/home/ubuntu/atac/h508/8_hg38_merged_bam/
#
#
# sample 1
java -Xmx30g -jar ${picard_dir}picard.jar MergeSamFiles \
I=${bam_dir}s1_l1_hg38.bam \
I=${bam_dir}s1_l2_hg38.bam \
I=${bam_dir}s1_l3_hg38.bam \
I=${bam_dir}s1_l4_hg38.bam \
O=${out_dir}sample1_hg38_merged.bam \
SORT_ORDER=coordinate \
ASSUME_SORTED=true \
MERGE_SEQUENCE_DICTIONARIES=true \
USE_THREADING=true &
#
#
# sample 2
java -Xmx30g -jar ${picard_dir}picard.jar MergeSamFiles \
I=${bam_dir}s2_l1_hg38.bam \
I=${bam_dir}s2_l2_hg38.bam \
I=${bam_dir}s2_l3_hg38.bam \
I=${bam_dir}s2_l4_hg38.bam \
O=${out_dir}sample2_hg38_merged.bam \
SORT_ORDER=coordinate \
ASSUME_SORTED=true \
MERGE_SEQUENCE_DICTIONARIES=true \
USE_THREADING=true &
#
#
# sample 3
java -Xmx30g -jar ${picard_dir}picard.jar MergeSamFiles \
I=${bam_dir}s3_l1_hg38.bam \
I=${bam_dir}s3_l2_hg38.bam \
I=${bam_dir}s3_l3_hg38.bam \
I=${bam_dir}s3_l4_hg38.bam \
O=${out_dir}sample3_hg38_merged.bam \
SORT_ORDER=coordinate \
ASSUME_SORTED=true \
MERGE_SEQUENCE_DICTIONARIES=true \
USE_THREADING=true &
#