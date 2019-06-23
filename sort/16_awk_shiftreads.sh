#!/bin/bash
#
# awk (gawk for gnu awk) is a programming language used for pattern matching and processing of text files
# the basic format of any awk command is as follows:
# 	awk '/search_pattern/ { action_to_take_on_matches; another_action; }' file_to_parse
#
# the BEGIN block is an awk rule that is executed only once, before first input is received
# OFS = output field separator (space by default), here it is specified as a tab
# while $0 can be used to specify the whole line, $1, $2, $3... etc are used to specify fields from 1 to n
# == is a comparison operator that means equal to
# what comes after the print command here is saying if its the positive strand reprint with adding 4 to the start and end chromosomal location,
# and otherwise do the same but with subtracting 5 from the start and end positions
#
# define paths
bed_dir=/home/ubuntu/atac/h508/10_hg38_bed/
out_dir=/home/ubuntu/atac/h508/11_hg38_shifted_bed/
#
# 
# sample 1
awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $3 + 4, $4, $5, $6; else print $1, $2 - 5, $3 - 5, $4, $5, $6}' ${bed_dir}sample1.bed > ${out_dir}sample1_shifted.bed &
awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $3 + 4, $4, $5, $6; else print $1, $2 - 5, $3 - 5, $4, $5, $6}' ${bed_dir}sample2.bed > ${out_dir}sample2_shifted.bed &
awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $3 + 4, $4, $5, $6; else print $1, $2 - 5, $3 - 5, $4, $5, $6}' ${bed_dir}sample3.bed > ${out_dir}sample3_shifted.bed &
#