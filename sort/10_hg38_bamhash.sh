#!/bin/bash
#
#
# Cleans the provided SAM/BAM, soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads
#
# specify paths
#bamhash_dir=/home/ubuntu/bamhash/
#bam_dir=/home/ubuntu/atac/h508/7_aligned_hg38_bam/
#fastq_dir=/home/ubuntu/atac/h508/4_good_fastq/
#
#
# sample 1
/home/ubuntu/bamhash/bamhash_checksum_bam /home/ubuntu/atac/h508/7_aligned_hg38_bam/s1_l1_hg38.bam > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s1_l1_bamhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_bam /home/ubuntu/atac/h508/7_aligned_hg38_bam/s1_l2_hg38.bam > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s1_l2_bamhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_bam /home/ubuntu/atac/h508/7_aligned_hg38_bam/s1_l3_hg38.bam > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s1_l3_bamhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_bam /home/ubuntu/atac/h508/7_aligned_hg38_bam/s1_l4_hg38.bam > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s1_l4_bamhash.txt &
#
#
# sample 2
/home/ubuntu/bamhash/bamhash_checksum_bam /home/ubuntu/atac/h508/7_aligned_hg38_bam/s2_l1_hg38.bam > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s2_l1_bamhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_bam /home/ubuntu/atac/h508/7_aligned_hg38_bam/s2_l2_hg38.bam > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s2_l2_bamhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_bam /home/ubuntu/atac/h508/7_aligned_hg38_bam/s2_l3_hg38.bam > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s2_l3_bamhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_bam /home/ubuntu/atac/h508/7_aligned_hg38_bam/s2_l4_hg38.bam > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s2_l4_bamhash.txt &
#
#
# sample 3
/home/ubuntu/bamhash/bamhash_checksum_bam /home/ubuntu/atac/h508/7_aligned_hg38_bam/s3_l1_hg38.bam > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s3_l1_bamhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_bam /home/ubuntu/atac/h508/7_aligned_hg38_bam/s3_l2_hg38.bam > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s3_l2_bamhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_bam /home/ubuntu/atac/h508/7_aligned_hg38_bam/s3_l3_hg38.bam > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s3_l3_bamhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_bam /home/ubuntu/atac/h508/7_aligned_hg38_bam/s3_l4_hg38.bam > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s3_l4_bamhash.txt &
#
#
# sample 1
/home/ubuntu/bamhash/bamhash_checksum_fastq /home/ubuntu/atac/h508/4_good_fastq/s1_l1_r1.good.fastq /home/ubuntu/atac/h508/4_good_fastq/s1_l1_r2.good.fastq > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s1_l1_fastqhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_fastq /home/ubuntu/atac/h508/4_good_fastq/s1_l2_r1.good.fastq /home/ubuntu/atac/h508/4_good_fastq/s1_l2_r2.good.fastq > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s1_l2_fastqhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_fastq /home/ubuntu/atac/h508/4_good_fastq/s1_l3_r1.good.fastq /home/ubuntu/atac/h508/4_good_fastq/s1_l3_r2.good.fastq > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s1_l3_fastqhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_fastq /home/ubuntu/atac/h508/4_good_fastq/s1_l4_r1.good.fastq /home/ubuntu/atac/h508/4_good_fastq/s1_l4_r2.good.fastq > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s1_l4_fastqhash.txt &
#
#
# sample 2
/home/ubuntu/bamhash/bamhash_checksum_fastq /home/ubuntu/atac/h508/4_good_fastq/s2_l1_r1.good.fastq /home/ubuntu/atac/h508/4_good_fastq/s2_l1_r2.good.fastq > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s2_l1_fastqhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_fastq /home/ubuntu/atac/h508/4_good_fastq/s2_l2_r1.good.fastq /home/ubuntu/atac/h508/4_good_fastq/s2_l2_r2.good.fastq > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s2_l2_fastqhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_fastq /home/ubuntu/atac/h508/4_good_fastq/s2_l3_r1.good.fastq /home/ubuntu/atac/h508/4_good_fastq/s2_l3_r2.good.fastq > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s2_l3_fastqhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_fastq /home/ubuntu/atac/h508/4_good_fastq/s2_l4_r1.good.fastq /home/ubuntu/atac/h508/4_good_fastq/s2_l4_r2.good.fastq > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s2_l4_fastqhash.txt &
#
#
# sample 3
/home/ubuntu/bamhash/bamhash_checksum_fastq /home/ubuntu/atac/h508/4_good_fastq/s3_l1_r1.good.fastq /home/ubuntu/atac/h508/4_good_fastq/s3_l1_r2.good.fastq > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s3_l1_fastqhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_fastq /home/ubuntu/atac/h508/4_good_fastq/s3_l2_r1.good.fastq /home/ubuntu/atac/h508/4_good_fastq/s3_l2_r2.good.fastq > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s3_l2_fastqhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_fastq /home/ubuntu/atac/h508/4_good_fastq/s3_l3_r1.good.fastq /home/ubuntu/atac/h508/4_good_fastq/s3_l3_r2.good.fastq > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s3_l3_fastqhash.txt &
/home/ubuntu/bamhash/bamhash_checksum_fastq /home/ubuntu/atac/h508/4_good_fastq/s3_l4_r1.good.fastq /home/ubuntu/atac/h508/4_good_fastq/s3_l4_r2.good.fastq > /home/ubuntu/atac/h508/7_aligned_hg38_bam/hashes/s3_l4_fastqhash.txt &
#