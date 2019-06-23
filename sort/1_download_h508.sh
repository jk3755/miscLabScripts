#!/bin/bash
#
#
# {access_token}=9e1acf007f014868b78778d705ce1e64
# files located at : https://basespace.illumina.com/run/67045980/JSK_COAD_ATACSEQ_020918/samples
#
# files ID's
# 
# sample 1
# https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L001_R1_001.fastq.gz?id=10172213598
# https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L001_R2_001.fastq.gz?id=10172213599
# https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L002_R1_001.fastq.gz?id=10172213596
# https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L002_R2_001.fastq.gz?id=10172213597
# https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L003_R1_001.fastq.gz?id=10172213601
# https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L003_R2_001.fastq.gz?id=10172213602
# https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L004_R1_001.fastq.gz?id=10172213600
# https://basespace.illumina.com/sample/103349250/files/tree//H508-1_S3_L004_R2_001.fastq.gz?id=10172213603
#
# sample 2
# https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L001_R1_001.fastq.gz?id=10172226239
# https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L001_R2_001.fastq.gz?id=10172226240
# https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L002_R1_001.fastq.gz?id=10172226242
# https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L002_R2_001.fastq.gz?id=10172226241
# https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L003_R1_001.fastq.gz?id=10172226243
# https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L003_R2_001.fastq.gz?id=10172226244
# https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L004_R1_001.fastq.gz?id=10172226245
# https://basespace.illumina.com/sample/103348254/files/tree//H508-2_S2_L004_R2_001.fastq.gz?id=10172226246
#
# sample 3
# https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L001_R1_001.fastq.gz?id=10172224329
# https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L001_R2_001.fastq.gz?id=10172224328
# https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L002_R1_001.fastq.gz?id=10172224331
# https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L002_R2_001.fastq.gz?id=10172224330
# https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L003_R1_001.fastq.gz?id=10172224332
# https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L003_R2_001.fastq.gz?id=10172224334
# https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L004_R1_001.fastq.gz?id=10172224333
# https://basespace.illumina.com/sample/103350250/files/tree//H508-3_S1_L004_R2_001.fastq.gz?id=10172224335
#
#
# to download files, use : 
# wget -O filename 'https://api.basespace.illumina.com/v1pre3/files/{id}/content?access_token={token}'
# replacing {id} with the specific file id from link addresses above, and {token} with the access token
# NOTE: be sure to remove the {}
#
# sample 1
wget -O H508-1_S3_L001_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213598/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L001_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213599/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L002_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213596/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L002_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213597/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L003_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213601/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L003_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213602/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L004_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213600/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-1_S3_L004_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172213603/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
#
# sample 2
wget -O H508-2_S2_L001_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226239/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L001_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226240/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L002_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226242/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L002_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226241/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L003_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226243/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L003_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226244/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L004_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226245/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-2_S2_L004_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172226246/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
#
# sample 3
wget -O H508-3_S1_L001_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224329/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L001_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224328/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L002_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224331/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L002_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224330/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L003_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224332/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L003_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224334/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L004_R1_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224333/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!
wget -O H508-3_S1_L004_R2_001.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/10172224335/content?access_token=9e1acf007f014868b78778d705ce1e64' &
wait $!