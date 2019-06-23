#!/bash/bin
#
java -jar /home/master/picard/picard.jar CollectInsertSizeMetrics \
I=lib1_merge_clean.bam \
O=insertsize.txt \
H=histo.pdf \
M=0.5
