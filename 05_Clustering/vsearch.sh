#!/bin/bash
#SBATCH --job-name=vsearch
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

##################################################
##  	Clustering using vsearch		##
##################################################
module load vsearch/2.4.3

vsearch --threads 8 --log LOGFile \
	--cluster_fast ../Coding_Regions/trinity_combine.fasta.transdecoder.cds \
	--id 0.90 \
	--centroids centroids.fasta \
	--uc clusters.uc


date 


