#!/bin/bash
##

# getting TSO reads with linkers

zcat 315641_1-Unmod_S4_R1_001.fastq.gz |awk 'NR%4==2' |grep -P '^[ACTG]{9}AATG[ACTG]{9}CCAC' > output.txt

echo "got TSO reads with linkers"

# getting reads matching CL1

while read line
do
grep -P "^${line}" output.txt >> output_CLS1.txt
done < ~/whitelist_96x3/BD_CLS1.txt

echo "got reads matching CL1"

# getting reads matching CL2

while read line
do
grep -P "^[ACTG]{9}AATG${line}" output_CLS1.txt >> output_CLS2.txt
done < ~/whitelist_96x3/BD_CLS2.txt

echo "got reads matching CL2"

# getting reads matching CL3

while read line
do
grep -P "^[ACTG]{9}AATG[ACTG]{9}CCAC${line}" output_CLS2.txt >> output_CLS3.txt
done < ~/whitelist_96x3/BD_CLS3.txt

echo "got reads matching CL3"

# outputting numbers

echo "number of reads with TSO linkers" >> numbers_unmod_second_mixing.txt
wc -l output.txt >> numbers_unmod_second_mixing.txt

echo "number of reads with CL1" >> numbers_unmod_second_mixing.txt
wc -l output_CLS1.txt >> numbers_unmod_second_mixing.txt

echo "number of reads with CL2" >> numbers_unmod_second_mixing.txt
wc -l output_CLS2.txt >> numbers_unmod_second_mixing.txt

echo "number of reads with CL3" >> numbers_unmod_second_mixing.txt
wc -l output_CLS3.txt >> numbers_unmod_second_mixing.txt

rm output.txt
rm output_CLS1.txt
rm output_CLS2.txt
rm output_CLS3.txt
