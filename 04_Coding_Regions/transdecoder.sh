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
sed 's/>/>K23_/g' ../Assembly/trinity_K23.Trinity.fasta > ../Assembly/trinity_prefix_K23.Trinity.fasta
sed 's/>/>K32_/g' ../Assembly/trinity_K32.Trinity.fasta > ../Assembly/trinity_prefix_K32.Trinity.fasta
sed 's/>/>U13_/g' ../Assembly/trinity_U13.Trinity.fasta > ../Assembly/trinity_prefix_U13.Trinity.fasta
sed 's/>/>U32_/g' ../Assembly/trinity_U32.Trinity.fasta > ../Assembly/trinity_prefix_U32.Trinity.fasta

# concatenate the assemblies
cat ../Assembly/trinity_prefix_U13.Trinity.fasta \
	../Assembly/trinity_prefix_U32.Trinity.fasta \
	../Assembly/trinity_prefix_K32.Trinity.fasta \
	../Assembly/trinity_prefix_K23.Trinity.fasta >> ../Assembly/trinity_combine.fasta


##################################################
## Determine ORF using Transdecoder		##
##################################################
module load hmmer/3.2.1
module load TransDecoder/5.3.0

TransDecoder.LongOrfs -t ../Assembly/trinity_combine.fasta

hmmscan --cpu 16 \
       --domtblout pfam.domtblout \
       /isg/shared/databases/Pfam/Pfam-A.hmm \
       trinity_combine.fasta.transdecoder_dir/longest_orfs.pep 


TransDecoder.Predict -t ../Assembly/trinity_combine.fasta \
	--retain_pfam_hits pfam.domtblout \
	--cpu 16 



