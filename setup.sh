#!/bin/bash

# run all tutorial scripts using slurm dependencies

cd 01_Raw_Reads
jid1=$(sbatch --parsable raw_data_symlinks.sh)

cd ../02_Quality_Control/
jid2=$(sbatch --parsable --dependency=afterok:$jid1 qc_raw.sh)
jid3=$(sbatch --parsable --dependency=afterok:$jid2 trimmomatic.sh)
jid4=$(sbatch --parsable --dependency=afterok:$jid3 qc_trimmed.sh)



