#!/bin/bash
#
# -mapped : split by mapped/unmapped
# -reference : split by chromosome
#
in1="/home/master/atac1/lib1/bam_split/lib1_dp.bam"
#
/home/master/bamtools/bin/bamtools split -in $in1 -reference
#
