#!/bin/bash
#SBATCH --job-name=kallisto_counts
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-5]
##SBATCH --mail-type=ALL
##SBATCH --mail-user=
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

hostname
date

##########################################
## kallisto quantification algorithm	##	
##########################################

module load kallisto/0.44.0

# this is an array job. SLURM will run this script 6 times in parallel (#SBATCH --array=[0-5]) contingent on resource availability
	# each time SLURM will change the value of the variable SLURM_ARRAY_TASK_ID
	# we'll use a bash array and that variable to retrieve a diferent sample

# a bash array containing the sample IDs
LIST=($(echo K21 K22 K23 K31 K32 K33))

# get one sample ID using the slurm array task ID
SAM=${LIST[$SLURM_ARRAY_TASK_ID]}

kallisto quant \
	-i ../05_Clustering/centroids.fasta.index \
	-o ${SAM} \
	-t 8 \
	../Quality_Control/trim_${SAM}_R1.fastq.gz ../Quality_Control/trim_${SAM}_R2.fastq.gz

date 


