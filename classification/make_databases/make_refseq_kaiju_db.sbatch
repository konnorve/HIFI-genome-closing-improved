#!/usr/bin/env bash
#SBATCH --job-name=refseqdb
#SBATCH --time 10-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 20                                # tasks (essentially threads)
#SBATCH -N 1
#SBATCH --mem 250G
#SBATCH -o logs/make_refseq_db.out
#SBATCH -e logs/make_refseq_db.err

outdir=/nobackup1/chisholmlab/kve/2023_HIFI_Genome_Closing_Project/kaiju_db/refseq_20220202

mkdir -p ${outdir}

date
echo "starting refseq download"
cd ${outdir}
/home/kve/scripts/HIFI-genome-closing-improved/classification/kaiju/bin/kaiju-makedb -s refseq

date