#!/bin/bash
#SBATCH --job-name=trinity
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 36
#SBATCH --mem=120G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-5]
##SBATCH --mail-type=ALL
##SBATCH --mail-user=
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

hostname
date

module load trinity/2.15.0
module load samtools

# this is an array job. SLURM will run this script 6 times in parallel (#SBATCH --array=[0-5]) contingent on resource availability
	# each time SLURM will change the value of the variable SLURM_ARRAY_TASK_ID
	# we'll use a bash array and that variable to retrieve a diferent sample

# a bash array containing the sample IDs
LIST=($(echo K21 K22 K23 K31 K32 K33))

# get one sample ID using the slurm array task ID
SAM=${LIST[$SLURM_ARRAY_TASK_ID]}

Trinity --seqType fq \
	--left ../02_Quality_Control/trim_${SAM}_R1.fastq.gz \
	--right ../02_Quality_Control/trim_${SAM}_R2.fastq.gz \
	--min_contig_length 300 \
	--CPU 36 \
	--max_memory 100G \
	--output trinity_${SAM} \
	--full_cleanup 
date 


