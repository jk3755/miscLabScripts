#!/bin/bash
#
# specify paths
picard_dir=/home/ubuntu/picard/
bam_dir=/home/ubuntu/atac/h508/9_hg38_dedepulicated_bam/
#
#
# sample 1
java -Xmx30g -jar ${picard_dir}picard.jar CleanSam \
I=/home/ubuntu/atac/h508/9_hg38_deduplicated_bam/sample1_hg38_merged_deduplicated.bam \
O=/home/ubuntu/atac/h508/9_hg38_deduplicated_bam/sample1_hg38_merged_deduplicated_clean.bam &
#
#
# sample 2
java -Xmx30g -jar ${picard_dir}picard.jar CleanSam \
I=/home/ubuntu/atac/h508/9_hg38_deduplicated_bam/sample2_hg38_merged_deduplicated.bam \
O=/home/ubuntu/atac/h508/9_hg38_deduplicated_bam/sample2_hg38_merged_deduplicated_clean.bam &
#
#
# sample 3
java -Xmx30g -jar ${picard_dir}picard.jar CleanSam \
I=/home/ubuntu/atac/h508/9_hg38_deduplicated_bam/sample3_hg38_merged_deduplicated.bam \
O=/home/ubuntu/atac/h508/9_hg38_deduplicated_bam/sample3_hg38_merged_deduplicated_clean.bam &
#