#!/bin/bash
#
# 
# use in verbose mode to get full report
# specify paths
picard_dir=/home/ubuntu/picard/
sam_dir=/home/ubuntu/atac/h508/6_aligned_hg38_sam_rg_sorted/
out_dir=/home/ubuntu/atac/h508/6_aligned_hg38_sam_rg_sorted/reports/
#
#
# sample 1
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${sam_dir}s1_l1_hg38_rg_sorted.sam \
O=${out_dir}s1_l1_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &

java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${sam_dir}s1_l2_hg38_rg_sorted.sam \
O=${out_dir}s1_l2_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${sam_dir}s1_l3_hg38_rg_sorted.sam \
O=${out_dir}s1_l3_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${sam_dir}s1_l4_hg38_rg_sorted.sam \
O=${out_dir}s1_l4_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
#
# sample 2
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${sam_dir}s2_l1_hg38_rg_sorted.sam \
O=${out_dir}s2_l1_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${sam_dir}s2_l2_hg38_rg_sorted.sam \
O=${out_dir}s2_l2_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${sam_dir}s2_l3_hg38_rg_sorted.sam \
O=${out_dir}s2_l3_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${sam_dir}s2_l4_hg38_rg_sorted.sam \
O=${out_dir}s2_l4_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
#
# sample 3
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${sam_dir}s3_l1_hg38_rg_sorted.sam \
O=${out_dir}s3_l1_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${sam_dir}s3_l2_hg38_rg_sorted.sam \
O=${out_dir}s3_l2_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${sam_dir}s3_l3_hg38_rg_sorted.sam \
O=${out_dir}s3_l3_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &
#
java -Xmx5g -jar ${picard_dir}picard.jar ValidateSamFile \
I=${sam_dir}s3_l4_hg38_rg_sorted.sam \
O=${out_dir}s3_l4_report.txt \
MODE=VERBOSE \
IGNORE=MISMATCH_FLAG_MATE_UNMAPPED \
IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND &