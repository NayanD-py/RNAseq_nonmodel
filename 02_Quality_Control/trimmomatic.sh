#!/bin/bash
#SBATCH --job-name=Trimmomatic
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

SAM=K21
java -jar $Trimmomatic PE -threads 4 \
        ../01_Raw_Reads/${SAM}_R1.fastq.gz \
        ../01_Raw_Reads/${SAM}_R2.fastq.gz \
        trim_${SAM}_R1.fastq.gz singles_trim_${SAM}_R1.fastq.gz \
        trim_${SAM}_R2.fastq.gz singles_trim_${SAM}_R2.fastq.gz \
        ILLUMINACLIP:/isg/shared/apps/Trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10 \
        SLIDINGWINDOW:4:25 MINLEN:45

SAM=K22
java -jar $Trimmomatic PE -threads 4 \
        ../01_Raw_Reads/${SAM}_R1.fastq.gz \
        ../01_Raw_Reads/${SAM}_R2.fastq.gz \
        trim_${SAM}_R1.fastq.gz singles_trim_${SAM}_R1.fastq.gz \
        trim_${SAM}_R2.fastq.gz singles_trim_${SAM}_R2.fastq.gz \
        ILLUMINACLIP:/isg/shared/apps/Trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10 \
        SLIDINGWINDOW:4:25 MINLEN:45

SAM=K23
java -jar $Trimmomatic PE -threads 4 \
        ../01_Raw_Reads/${SAM}_R1.fastq.gz \
        ../01_Raw_Reads/${SAM}_R2.fastq.gz \
        trim_${SAM}_R1.fastq.gz singles_trim_${SAM}_R1.fastq.gz \
        trim_${SAM}_R2.fastq.gz singles_trim_${SAM}_R2.fastq.gz \
        ILLUMINACLIP:/isg/shared/apps/Trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10 \
        SLIDINGWINDOW:4:25 MINLEN:45

SAM=K31
java -jar $Trimmomatic PE -threads 4 \
        ../01_Raw_Reads/${SAM}_R1.fastq.gz \
        ../01_Raw_Reads/${SAM}_R2.fastq.gz \
        trim_${SAM}_R1.fastq.gz singles_trim_${SAM}_R1.fastq.gz \
        trim_${SAM}_R2.fastq.gz singles_trim_${SAM}_R2.fastq.gz \
        ILLUMINACLIP:/isg/shared/apps/Trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10 \
        SLIDINGWINDOW:4:25 MINLEN:45

SAM=K32
java -jar $Trimmomatic PE -threads 4 \
        ../01_Raw_Reads/${SAM}_R1.fastq.gz \
        ../01_Raw_Reads/${SAM}_R2.fastq.gz \
        trim_${SAM}_R1.fastq.gz singles_trim_${SAM}_R1.fastq.gz \
        trim_${SAM}_R2.fastq.gz singles_trim_${SAM}_R2.fastq.gz \
        ILLUMINACLIP:/isg/shared/apps/Trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10 \
        SLIDINGWINDOW:4:25 MINLEN:45

SAM=K33
java -jar $Trimmomatic PE -threads 4 \
        ../01_Raw_Reads/${SAM}_R1.fastq.gz \
        ../01_Raw_Reads/${SAM}_R2.fastq.gz \
        trim_${SAM}_R1.fastq.gz singles_trim_${SAM}_R1.fastq.gz \
        trim_${SAM}_R2.fastq.gz singles_trim_${SAM}_R2.fastq.gz \
        ILLUMINACLIP:/isg/shared/apps/Trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10 \
        SLIDINGWINDOW:4:25 MINLEN:45

date
