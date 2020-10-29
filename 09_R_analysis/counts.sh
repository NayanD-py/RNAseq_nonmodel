#!/bin/bash

awk '{print $1 "\t" $4}' ../Counts/K23/abundance.tsv > K23.counts
awk '{print $1 "\t" $4}' ../Counts/K32/abundance.tsv > K32.counts
awk '{print $1 "\t" $4}' ../Counts/U32/abundance.tsv > U32.counts
awk '{print $1 "\t" $4}' ../Counts/U13/abundance.tsv > U13.counts
awk '{print $1 "\t" $2}' ../Counts/K23/abundance.tsv > length.txt
