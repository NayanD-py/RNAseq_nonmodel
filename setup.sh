#!/bin/bash

# run all tutorial scripts using slurm dependencies

cd 01_Raw_Reads
jid1=$(sbatch --parsable raw_data_symlinks.sh)

cd ../02_Quality_Control/
jid2=$(sbatch --parsable --dependency=afterok:$jid1 qc_raw.sh)
jid3=$(sbatch --parsable --dependency=afterok:$jid2 trimmomatic.sh)
jid4=$(sbatch --parsable --dependency=afterok:$jid3 qc_trimmed.sh)

cd ../03_Assembly/
jid5=$(sbatch --parsable --dependency=afterok:$jid4 trinity.sh)

cd ../04_Coding_Regions/
jid6=$(sbatch --parsable --dependency=afterok:$jid5 transdecoder.sh)

cd ../05_Clustering/
jid7=$(sbatch --parsable --dependency=afterok:$jid6 vsearch.sh)

cd ../06_RNAQuast/
jid8=$(sbatch --parsable --dependency=afterok:$jid7 rnaQuast.sh)

cd ../07_EnTAP/
jid9=$(sbatch --parsable --dependency=afterok:$jid7 entap.sh)

cd ../08_Counts/
jid10=$(sbatch --parsable --dependency=afterok:$jid7 kallisto_index.sh)
jid11=$(sbatch --parsable --dependency=afterok:$jid10 kallisto_counts.sh)
