#!/bin/bash
#SBATCH --job-name=kallisto_counts
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
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
## kallisto quantification algorithm	##	
##########################################
module load kallisto/0.44.0

kallisto quant -i ../Index/Eastern_larch_index \
	-o K23 \
	-t 8 \
	../Quality_Control/trim_K23_R1.fastq ../Quality_Control/trim_K23_R2.fastq

kallisto quant -i ../Index/Eastern_larch_index \
	-o K32 \
	-t 8 \
	../Quality_Control/trim_K32_R1.fastq ../Quality_Control/trim_K32_R2.fastq

kallisto quant -i ../Index/Eastern_larch_index \
        -o U32 \
        -t 8 \
        ../Quality_Control/trim_U32_R1.fastq ../Quality_Control/trim_U32_R2.fastq

kallisto quant -i ../Index/Eastern_larch_index \
        -o U13 \
        -t 8 \
        ../Quality_Control/trim_U13_R1.fastq ../Quality_Control/trim_U13_R2.fastq


date 


