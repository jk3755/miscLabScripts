#!/bin/bash
#
# the bedpe option outputs in bedpe (paired end) format
#
# define paths
bed_dir=/home/ubuntu/bedtools/bedtools2/bin/
bam_dir=/home/ubuntu/atac/h508/9_hg38_deduplicated_bam/
out_dir=/home/ubuntu/atac/h508/10_hg38_bed/
#
${bed_dir}bedtools bamtobed -bedpe -i ${bam_dir}sample1_hg38_merged_deduplicated_clean.bam > ${out_dir}sample1.bed &
${bed_dir}bedtools bamtobed -bedpe -i ${bam_dir}sample2_hg38_merged_deduplicated_clean.bam > ${out_dir}sample2.bed &
${bed_dir}bedtools bamtobed -bedpe -i ${bam_dir}sample3_hg38_merged_deduplicated_clean.bam > ${out_dir}sample3.bed &
#