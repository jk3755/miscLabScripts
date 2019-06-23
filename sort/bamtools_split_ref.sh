#!/bin/bash
#
# -mapped : split by mapped/unmapped
# -reference : split by chromosome
#
in1="/home/master/atac1/merged/ref/lib123_dp.MAPPED.bam"
in2="/home/master/atac1/merged/ref/lib456_dp.MAPPED.bam"
#
/home/master/bamtools/bin/bamtools split -in $in1 -reference &
/home/master/bamtools/bin/bamtools split -in $in2 -reference
#
