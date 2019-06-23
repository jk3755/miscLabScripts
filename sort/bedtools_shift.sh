#!/bin/bash
#
/home/jordan/bedtools2/bin/bedtools bamtobed -i lane1_sorted.bam -bedpe | awk -v OFS="\t" '{if($9=="+"){print $1,$2+4,$6+4}else if($9=="-"){print $1,$2-5,$6-5}}' > lib1_shift.bed
