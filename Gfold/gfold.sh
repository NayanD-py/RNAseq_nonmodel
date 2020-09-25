#!/bin/bash
#SBATCH --job-name=gfold
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date
##########################################
## Gfold				## 
##########################################
module load gfold/1.1.4 

gfold diff -s1 K21 -s2 K31 -suf .read_cnt -o K21_vs_K31.diff

date 


