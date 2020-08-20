#!/bin/bash
#SBATCH --job-name=qc_trimmed
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# FASTQC of raw reads 
#################################################################
# create an output directory to hold fastqc output
DIR="trimmed"
mkdir -p ${DIR}_fastqc

module load fastqc/0.11.7
SAM=K21
fastqc --outdir ./${DIR}_fastqc/ trim_${SAM}_R1.fastq.gz trim_${SAM}_R2.fastq.gz
SAM=K22
fastqc --outdir ./${DIR}_fastqc/ trim_${SAM}_R1.fastq.gz trim_${SAM}_R2.fastq.gz
SAM=K23
fastqc --outdir ./${DIR}_fastqc/ trim_${SAM}_R1.fastq.gz trim_${SAM}_R2.fastq.gz
SAM=K31
fastqc --outdir ./${DIR}_fastqc/ trim_${SAM}_R1.fastq.gz trim_${SAM}_R2.fastq.gz
SAM=K32
fastqc --outdir ./${DIR}_fastqc/ trim_${SAM}_R1.fastq.gz trim_${SAM}_R2.fastq.gz
SAM=K33
fastqc --outdir ./${DIR}_fastqc/ trim_${SAM}_R1.fastq.gz trim_${SAM}_R2.fastq.gz


#################################################################
# MULTIQC of raw reads 
#################################################################
module load MultiQC/1.8

mkdir -p ${DIR}_multiqc
multiqc --outdir ${DIR}_multiqc ./${DIR}_fastqc/