#!/bin/bash
##

# getting WTA reads with linkers

zcat 315641_4-RoCK_ROI_S3_R1_001.fastq.gz |awk 'NR%4==2' |grep -P '^[ACTG]{9,12}GTGA[ACTG]{9}GACA' > output.txt

echo "got WTA reads with linkers"

# getting reads matching CL1

while read line
do
grep -P "^[ACTG]{0,3}${line}" output.txt >> output_CLS1.txt
done < ~/whitelist_96x3/BD_CLS1.txt

echo "got reads matching CL1"

# getting reads matching CL2

while read line
do
grep -P "^[ACTG]{9,12}GTGA${line}" output_CLS1.txt >> output_CLS2.txt
done < ~/whitelist_96x3/BD_CLS2.txt

echo "got reads matching CL2"

# getting reads matching CL3

while read line
do
grep -P "^[ACTG]{9,12}GTGA[ACTG]{9}GACA${line}" output_CLS2.txt >> output_CLS3.txt
done < ~/whitelist_96x3/BD_CLS3.txt

echo "got reads matching CL3"

# outputting numbers

echo "number of reads with WTA linkers" >> wta_numbers_rockroi_second_mixing.txt
wc -l output.txt >> wta_numbers_rockroi_second_mixing.txt

echo "number of reads with CL1" >> wta_numbers_rockroi_second_mixing.txt
wc -l output_CLS1.txt >> wta_numbers_rockroi_second_mixing.txt

echo "number of reads with CL2" >> wta_numbers_rockroi_second_mixing.txt
wc -l output_CLS2.txt >> wta_numbers_rockroi_second_mixing.txt

echo "number of reads with CL3" >> wta_numbers_rockroi_second_mixing.txt
wc -l output_CLS3.txt >> wta_numbers_rockroi_second_mixing.txt

rm output.txt
rm output_CLS1.txt
rm output_CLS2.txt
rm output_CLS3.txt
