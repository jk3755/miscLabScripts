#!/bin/bash
#
# the bedpe option outputs in bedpe (paired end) format
#
# define paths
bed_dir=/home/ubuntu/bedtools/bedtools2/bin/
bam_dir=/home/ubuntu/atac/h508/9_deduplicated_bam/
ref_dir=/home/ubuntu/ref_genomes/hg38/
#
# paired end format
#${bed_dir}bedtools bamtobed -bedpe -i ${bam_dir}h508_deduplicated_namesort.bam > ${bam_dir}h508_deduplicated.bed
#
# normal bed
#${bed_dir}bedtools sort -i ${bam_dir}h508_deduplicated_nonpe.bed > ${bam_dir}h508_sorted.bed
#
#### UNIX sort can sort a bed file by chromosome first and then by start position more quickly and with less memory
sort -k 1,1 -k2,2n ${bam_dir}/h508_deduplicated_nonpe.bed > ./h508_deduplicated_nonpe_sorted.bed