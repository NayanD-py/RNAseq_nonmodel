#!/bin/bash
#SBATCH --job-name=transdecoder
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

##################################################
## Combining the trinity assemblies		##
##################################################

# add a sample name prefix to each sequence ID in each assembly
SAM=K21
sed "s/>/>${SAM}_/g" ../03_Assembly/trinity_${SAM}.Trinity.fasta > ../03_Assembly/trinity_prefix_${SAM}.Trinity.fasta
SAM=K22
sed "s/>/>${SAM}_/g" ../03_Assembly/trinity_${SAM}.Trinity.fasta > ../03_Assembly/trinity_prefix_${SAM}.Trinity.fasta
SAM=K23
sed "s/>/>${SAM}_/g" ../03_Assembly/trinity_${SAM}.Trinity.fasta > ../03_Assembly/trinity_prefix_${SAM}.Trinity.fasta
SAM=K31
sed "s/>/>${SAM}_/g" ../03_Assembly/trinity_${SAM}.Trinity.fasta > ../03_Assembly/trinity_prefix_${SAM}.Trinity.fasta
SAM=K32
sed "s/>/>${SAM}_/g" ../03_Assembly/trinity_${SAM}.Trinity.fasta > ../03_Assembly/trinity_prefix_${SAM}.Trinity.fasta
SAM=K33
sed "s/>/>${SAM}_/g" ../03_Assembly/trinity_${SAM}.Trinity.fasta > ../03_Assembly/trinity_prefix_${SAM}.Trinity.fasta

# concatenate the assemblies
cat ../03_Assembly/trinity_prefix_* > ../03_Assembly/trinity_combine.fasta


##################################################
## Determine ORF using Transdecoder		##
##################################################

module load hmmer/3.2.1
module load TransDecoder/5.3.0

TransDecoder.LongOrfs -t ../03_Assembly/trinity_combine.fasta

hmmscan --cpu 16 \
       --domtblout pfam.domtblout \
       /isg/shared/databases/Pfam/Pfam-A.hmm \
       trinity_combine.fasta.transdecoder_dir/longest_orfs.pep 


TransDecoder.Predict -t ../03_Assembly/trinity_combine.fasta \
	--retain_pfam_hits pfam.domtblout \
	--cpu 16 



