#!/bin/bash
##

grep -P '^[ACTG]{9}AATG[ACTG]{9}CCAC' simulation_prepends.txt >> output_simulation_prepends.txt
grep -P '^[ACTG]{9,12}GTGA[ACTG]{9}GACA' simulation_prepends.txt >> output_simulation_prepends.txt

wc -l output_simulation_prepends.txt
