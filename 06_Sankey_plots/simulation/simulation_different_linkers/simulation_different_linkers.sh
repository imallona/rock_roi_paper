#!/bin/bash
##

grep -P '^[ACTG]{9}AATG[ACTG]{9}CCAC' simulation_different_linkers.txt >> output_simulation_different_linkers.txt
grep -P '^[ACTG]{9,12}GTGA[ACTG]{9}GACA' simulation_different_linkers.txt >> output_simulation_different_linkers.txt

wc -l output_simulation_different_linkers.txt
