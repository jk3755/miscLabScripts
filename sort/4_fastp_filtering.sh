#!/bin/bash
#
#
# specify paths
fastp_dir=/home/ubuntu/fastp/
fastq_dir=/home/ubuntu/atac/h508/2_fastq/
out_dir=/home/ubuntu/atac/h508/4_good_fastq/
#
# 
# sample 1 processing
${fastp_dir}fastp -i ${fastq_dir}H508-1_S3_L001_R1_001.fastq -I ${fastq_dir}H508-1_S3_L001_R2_001.fastq -o ${out_dir}s1_l1_r1.good.fastq -O ${out_dir}s1_l1_r2.good.fastq &
${fastp_dir}fastp -i ${fastq_dir}H508-1_S3_L002_R1_001.fastq -I ${fastq_dir}H508-1_S3_L002_R2_001.fastq -o ${out_dir}s1_l2_r1.good.fastq -O ${out_dir}s1_l2_r2.good.fastq &
${fastp_dir}fastp -i ${fastq_dir}H508-1_S3_L003_R1_001.fastq -I ${fastq_dir}H508-1_S3_L003_R2_001.fastq -o ${out_dir}s1_l3_r1.good.fastq -O ${out_dir}s1_l3_r2.good.fastq &
${fastp_dir}fastp -i ${fastq_dir}H508-1_S3_L004_R1_001.fastq -I ${fastq_dir}H508-1_S3_L004_R2_001.fastq -o ${out_dir}s1_l4_r1.good.fastq -O ${out_dir}s1_l4_r2.good.fastq &
#
#
# sample 2 processing
${fastp_dir}fastp -i ${fastq_dir}H508-2_S2_L001_R1_001.fastq -I ${fastq_dir}H508-2_S2_L001_R2_001.fastq -o ${out_dir}s2_l1_r1.good.fastq -O ${out_dir}s2_l1_r2.good.fastq &
${fastp_dir}fastp -i ${fastq_dir}H508-2_S2_L002_R1_001.fastq -I ${fastq_dir}H508-2_S2_L002_R2_001.fastq -o ${out_dir}s2_l2_r1.good.fastq -O ${out_dir}s2_l2_r2.good.fastq &
${fastp_dir}fastp -i ${fastq_dir}H508-2_S2_L003_R1_001.fastq -I ${fastq_dir}H508-2_S2_L003_R2_001.fastq -o ${out_dir}s2_l3_r1.good.fastq -O ${out_dir}s2_l3_r2.good.fastq &
${fastp_dir}fastp -i ${fastq_dir}H508-2_S2_L004_R1_001.fastq -I ${fastq_dir}H508-2_S2_L004_R2_001.fastq -o ${out_dir}s2_l4_r1.good.fastq -O ${out_dir}s2_l4_r2.good.fastq &
#
#
# sample 3 processing
${fastp_dir}fastp -i ${fastq_dir}H508-3_S1_L001_R1_001.fastq -I ${fastq_dir}H508-3_S1_L001_R2_001.fastq -o ${out_dir}s3_l1_r1.good.fastq -O ${out_dir}s3_l1_r2.good.fastq &
${fastp_dir}fastp -i ${fastq_dir}H508-3_S1_L002_R1_001.fastq -I ${fastq_dir}H508-3_S1_L002_R2_001.fastq -o ${out_dir}s3_l2_r1.good.fastq -O ${out_dir}s3_l2_r2.good.fastq &
${fastp_dir}fastp -i ${fastq_dir}H508-3_S1_L003_R1_001.fastq -I ${fastq_dir}H508-3_S1_L003_R2_001.fastq -o ${out_dir}s3_l3_r1.good.fastq -O ${out_dir}s3_l3_r2.good.fastq &
${fastp_dir}fastp -i ${fastq_dir}H508-3_S1_L004_R1_001.fastq -I ${fastq_dir}H508-3_S1_L004_R2_001.fastq -o ${out_dir}s3_l4_r1.good.fastq -O ${out_dir}s3_l4_r2.good.fastq &