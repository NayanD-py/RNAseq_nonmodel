#!/bin/bash
#SBATCH --job-name=Trimommatic
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load Trimmomatic/0.39

java -jar $Trimmomatic PE -threads 4 \
        ../Raw_Reads/U13/U13_R1.fastq \
        ../Raw_Reads/U13/U13_R2.fastq \
        trim_U13_R1.fastq singles_trim_U13_R1.fastq \
        trim_U13_R2.fastq singles_trim_U13_R2.fastq \
        ILLUMINACLIP:/isg/shared/apps/Trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10 \
        SLIDINGWINDOW:4:25 MINLEN:45

java -jar $Trimmomatic PE -threads 4 \
        ../Raw_Reads/U32/U32_R1.fastq \
        ../Raw_Reads/U32/U32_R2.fastq \
        trim_U32_R1.fastq singles_trim_U32_R1.fastq \
        trim_U32_R2.fastq singles_trim_U32_R2.fastq \
        ILLUMINACLIP:/isg/shared/apps/Trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10 \
        SLIDINGWINDOW:4:25 MINLEN:45

java -jar $Trimmomatic PE -threads 4 \
        ../Raw_Reads/K32/K32_R1.fastq \
        ../Raw_Reads/K32/K32_R2.fastq \
        trim_K32_R1.fastq singles_trim_K32_R1.fastq \
        trim_K32_R2.fastq singles_trim_K32_R2.fastq \
        ILLUMINACLIP:/isg/shared/apps/Trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10 \
        SLIDINGWINDOW:4:25 MINLEN:45


java -jar $Trimmomatic PE -threads 4 \
        ../Raw_Reads/K23/K23_R1.fastq \
        ../Raw_Reads/K23/K23_R2.fastq \
        trim_K23_R1.fastq singles_trim_K23_R1.fastq \
        trim_K23_R2.fastq singles_trim_K23_R2.fastq \
        ILLUMINACLIP:/isg/shared/apps/Trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10 \
        SLIDINGWINDOW:4:25 MINLEN:45
date
