#!/bin/bash

for COUNT in K21 K22 K23 K31 K32 K33
do awk '{print $1 "\t" $1 "\t" $4 "\t" $2 "\t" $5 }' ../08_Counts/${COUNT}/abundance.tsv > ${COUNT}.read_cnt
done
