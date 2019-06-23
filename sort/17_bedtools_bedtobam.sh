#!/bin/bash
#
#
# define paths
bedtool_dir=/home/ubuntu/bedtools/bedtools2/bin/
bed_dir=/home/ubuntu/atac/h508/11_hg38_shifted_bed/
out_dir=/home/ubuntu/atac/h508/12_shifted_bam/
#
${bedtool_dir}bedtools bedtobam -i ${bed_dir}sample1_shifted.bed -g /home/ubuntu/ref_genomes/hg38/hg38.genome > ${out_dir}sample1_shifted.bam &
${bedtool_dir}bedtools bedtobam -i ${bed_dir}sample2_shifted.bed -g /home/ubuntu/ref_genomes/hg38/hg38.genome > ${out_dir}sample2_shifted.bam &
${bedtool_dir}bedtools bedtobam -i ${bed_dir}sample3_shifted.bed -g /home/ubuntu/ref_genomes/hg38/hg38.genome > ${out_dir}sample3_shifted.bam &
#