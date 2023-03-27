#!/usr/bin/env bash
#SBATCH --job-name=mmseq2_db
#SBATCH --time 10-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 20                                # tasks (essentially threads)
#SBATCH -N 1
#SBATCH --mem 250G
#SBATCH -o logs/make_mmseqs2_uniref90_db.out
#SBATCH -e logs/make_mmseqs2_uniref90_db.err

source activate mmseqs2

date

outpath="/nobackup1b/users/chisholmlab/kve/2023_HIFI_Genome_Closing_Project/mmseqs2_db/uniref90_20220202"

mmseqs databases UniRef90 ${outpath} /tmp

date