#!/bin/bash
#SBATCH --job-name=fastqer_dump_xanadu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Download fastq files from SRA 
#################################################################
module load parallel/20180122
module load sratoolkit/3.0.1


# The BioProject accession for these data is here:
	# https://www.ncbi.nlm.nih.gov/bioproject/PRJNA612058

ACCLIST=accessionlist.txt

cat $ACCLIST | xargs fasterq-dump

######################
# compress the files #
######################

ls *fastq | parallel -j 12 gzip

#############################################
# re-name the files according to sample name
#############################################


for i in {1..6}
	do 
	ACC=$(sed -n ${i}p accession_to_sampleID.txt | cut -f 1)
	SAM=$(sed -n ${i}p accession_to_sampleID.txt | cut -f 2)
	mv ${ACC}_1.fastq.gz ${SAM}_R1.fastq.gz
	mv ${ACC}_2.fastq.gz ${SAM}_R2.fastq.gz
	done
