#!/usr/bin/env bash
#SBATCH --job-name=class
#SBATCH --time 10-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 20                                # tasks (essentially threads)
#SBATCH -N 1
#SBATCH --mem 250G
#SBATCH --array=1-440%5
#SBATCH -o logs/classification.%a.out
#SBATCH -e logs/classification.%a.err

date

assembly=$(realpath /nfs/chisholmlab001/kve/2023_HIFI_Genome_Closing_Project/assembly/all_assemblies/* | sed -n ${SLURM_ARRAY_TASK_ID}p)
name=$(basename ${assembly})
name=${name%.*}

database_outdir=/nfs/chisholmlab001/kve/2023_HIFI_Genome_Closing_Project/classification/mmseqs2_databases/${name}
report_outdir=/nfs/chisholmlab001/kve/2023_HIFI_Genome_Closing_Project/classification/mmseqs2_reports

uniprot_db=/nobackup1b/users/chisholmlab/kve/2023_HIFI_Genome_Closing_Project/mmseqs2_db/GTDB
assembly_db=${database_outdir}/${name}.assembly_db
assignment_db=${database_outdir}/${name}.assignment_db
report_path=${report_outdir}/${name}.tsv

mkdir -p ${database_outdir}
mkdir -p ${report_outdir}

echo ${name}
echo ${assembly}
echo ${database_outdir}
echo ${report_path}

source activate mmseqs2

date
echo "mmseqs create assembly database"
mmseqs createdb ${assembly} ${assembly_db}

date
echo "mmseqs create taxonomy assignment database"
mmseqs taxonomy ${assembly_db} ${uniprot_db} ${assignment_db} /tmp

date
echo "mmseqs create tsv output"
mmseqs createtsv ${assembly_db} ${assignment_db} ${report_path}

date