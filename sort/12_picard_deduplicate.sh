#!/bin/bash
#
# specify paths
picard_dir=/home/ubuntu/picard/
bam_dir=/home/ubuntu/atac/h508/8_merged_bam/
out_dir=/home/ubuntu/atac/h508/9_deduplicated_bam/
#
#
# sample 1
java -Xmx30g -jar ${picard_dir}picard.jar MarkDuplicates \
I=${bam_dir}sample1_hg38_merged.bam \
O=${out_dir}sample1_hg38_merged_deduplicated.bam \
M=${out_dir}/metrics/sample1_duplication_metrics.txt \
REMOVE_DUPLICATES=true \
ASSUME_SORTED=true &
#
#
# sample 2
#java -Xmx30g -jar ${picard_dir}picard.jar MarkDuplicates \
#I=${bam_dir}sample2_hg38_merged.bam \
#O=${out_dir}sample2_hg38_merged_deduplicated.bam \
#M=${out_dir}/metrics/sample2_duplication_metrics.txt \
#REMOVE_DUPLICATES=true \
#ASSUME_SORTED=true &
#
#
# sample 3
#java -Xmx30g -jar ${picard_dir}picard.jar MarkDuplicates \
#I=${bam_dir}sample3_hg38_merged.bam \
#O=${out_dir}sample3_hg38_merged_deduplicated.bam \
#M=${out_dir}/metrics/sample3_duplication_metrics.txt \
#REMOVE_DUPLICATES=true \
#ASSUME_SORTED=true &
#