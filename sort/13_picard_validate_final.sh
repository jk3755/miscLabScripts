#!/bin/bash
#
# 
# use in verbose mode to get full report
# specify paths
picard_dir=/home/ubuntu/picard/
bam_dir=/home/ubuntu/atac/h508/9_hg38_deduplicated_bam/
#
#
# sample 1
java -Xmx10g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}sample1_hg38_merged_deduplicated.bam \
O=${bam_dir}/reports/sample1_dedepulicated_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
# sample 2
java -Xmx10g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}sample2_hg38_merged_deduplicated.bam \
O=${bam_dir}/reports/sample2_dedepulicated_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
# sample 3
java -Xmx10g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}sample3_hg38_merged_deduplicated.bam \
O=${bam_dir}/reports/sample3_dedepulicated_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &