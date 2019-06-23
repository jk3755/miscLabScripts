#!/bin/bash
#
# specify paths
picard_dir=/home/ubuntu/picard/
bam_dir=/home/ubuntu/atac/h508/9_deduplicated_bam/
#
#
java -Xmx50g -jar ${picard_dir}picard.jar MergeSamFiles \
I=${bam_dir}sample1_hg38_merged_deduplicated.bam \
I=${bam_dir}sample2_hg38_merged_deduplicated.bam \
I=${bam_dir}sample3_hg38_merged_deduplicated.bam \
O=${bam_dir}h508_deduplicated_merged.bam \
SORT_ORDER=coordinate \
ASSUME_SORTED=true \
MERGE_SEQUENCE_DICTIONARIES=true \
USE_THREADING=true
#
