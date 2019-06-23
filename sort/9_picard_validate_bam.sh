#!/bin/bash
#
# 
# use in verbose mode to get full report
# specify paths
picard_dir=/home/ubuntu/picard/
bam_dir=/home/ubuntu/atac/h508/7_aligned_hg38_bam/
out_dir=/home/ubuntu/atac/h508/7_aligned_hg38_bam/reports/
#
#
# sample 1
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}s1_l1_hg38.bam \
O=${out_dir}s1_l1_bam_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}s1_l2_hg38.bam \
O=${out_dir}s1_l2_bam_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}s1_l3_hg38.bam \
O=${out_dir}s1_l3_bam_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}s1_l4_hg38.bam \
O=${out_dir}s1_l4_bam_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
#
# sample 2
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}s2_l1_hg38.bam \
O=${out_dir}s2_l1_bam_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}s2_l2_hg38.bam \
O=${out_dir}s2_l2_bam_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}s2_l3_hg38.bam \
O=${out_dir}s2_l3_bam_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}s2_l4_hg38.bam \
O=${out_dir}s2_l4_bam_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
#
# sample 3
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}s3_l1_hg38.bam \
O=${out_dir}s3_l1_bam_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}s3_l2_hg38.bam \
O=${out_dir}s3_l2_bam_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}s3_l3_hg38.bam \
O=${out_dir}s3_l3_bam_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${bam_dir}s3_l4_hg38.bam \
O=${out_dir}s3_l4_bam_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &