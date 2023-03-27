#!/usr/bin/env bash
#SBATCH --job-name=mmseq2_db
#SBATCH --time 10-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 20                                # tasks (essentially threads)
#SBATCH -N 1
#SBATCH --mem 250G
#SBATCH -o logs/make_mmseqs2_nt_db.out
#SBATCH -e logs/make_mmseqs2_nt_db.err

source activate mmseqs2

date

outpath="/nobackup1b/users/chisholmlab/kve/2023_HIFI_Genome_Closing_Project/mmseqs2_db/uniprokb_20220201"

mmseqs databases UniProtKB ${outpath} /tmp

date