#!/bin/bash
##

grep -P '^[ACTG]{9}AATG[ACTG]{9}CCAC' simulation_whitelists_no_prepends.txt >> output_simulation_whitelists_no_prepends.txt
grep -P '^[ACTG]{9,12}GTGA[ACTG]{9}GACA' simulation_whitelists_no_prepends.txt >> output_simulation_whitelists_no_prepends.txt

wc -l output_simulation_whitelists_no_prepends.txt

### for WTA
# getting reads matching CL1

while read line
do
grep -P "^[ACTG]{0,3}${line}" output_simulation_whitelists_no_prepends.txt >> output_CLS1.txt
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

### for TSO
# getting reads matching CL1

while read line
do
grep -P "^${line}" output_simulation_whitelists_no_prepends.txt >> output_CLS1.txt
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

sort output_CLS3.txt| uniq| wc -l
