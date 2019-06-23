#!/bin/bash
#
#
# convert bam file to bed
# bedtools bamtobed [OPTIONS] -i <BAM>
# -bedpe specifies to write in paired end bed format
# you may see a warning printed to console that says "marked as pe, but mate doesnt occur next to it", you can prob ignore this as it doesnt seem to skip many reads when input is namesorted
bedtools bamtobed -bedpe -i /home/ubuntu/atac/h508/9_dedup_bam/h508_all_dedup_ns.bam > /home/ubuntu/atac/h508/9_dedup_bam/shuffle/h508_all_dedup_ns.bed
#
#
# shuffle reads
# bedtools shuffle [OPTIONS] -i <BED/GFF/VCF> -g <GENOME>
# -chrom keeps the reads on the same chromosome, shuffle only the position (this keeps distribution the same across genome for null model)
# -bedpe specifies paired end bed file
# -g genome file is the chrom.sizes file
bedtools shuffle -chrom -bedpe -i /home/ubuntu/atac/h508/9_dedup_bam/shuffle/h508_all_dedup_ns.bed -g /home/ubuntu/ref_genomes/hg38/hg38.genome > /home/ubuntu/atac/h508/9_dedup_bam/shuffle/h508_all_dedup_ns_shuffled.bed
#
#
# convert bed back to bam
# bedToBam [OPTIONS] -i <BED/GFF/VCF> -g <GENOME> > <BAM>
bedtools bedtobam -i /home/ubuntu/atac/h508/9_dedup_bam/shuffle/h508_all_dedup_ns_shuffled.bed -g /home/ubuntu/ref_genomes/hg38/hg38.genome > /home/ubuntu/atac/h508/9_dedup_bam/shuffle/h508_all_dedup_ns_shuffled.bam
#