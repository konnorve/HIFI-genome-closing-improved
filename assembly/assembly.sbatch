#!/usr/bin/env bash
#SBATCH --job-name=assem
#SBATCH --time 10-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 20                                # tasks (essentially threads)
#SBATCH -N 1
#SBATCH --mem 250G
#SBATCH --array=2-65%5
#SBATCH -o logs/assembly.%a.out
#SBATCH -e logs/assembly.%a.err

date

out_dir="/nfs/chisholmlab001/kve/2023_HIFI_Genome_Closing_Project/assembly"

batch=$(sed "${SLURM_ARRAY_TASK_ID}q;d" assembly.tsv | cut -f1 | tr -d '\r')
barcode=$(sed "${SLURM_ARRAY_TASK_ID}q;d" assembly.tsv | cut -f2 | tr -d '\r')
ccs_fastq_path=$(sed "${SLURM_ARRAY_TASK_ID}q;d" assembly.tsv | cut -f3 | tr -d '\r')
subreadset_fastq_path=$(sed "${SLURM_ARRAY_TASK_ID}q;d" assembly.tsv | cut -f4 | tr -d '\r')

echo ${batch}
echo ${barcode}
echo ${ccs_fastq_path}
echo ${subreadset_fastq_path}

flye_ccs_outdir="${out_dir}/${batch}/${barcode}/flye_ccs"
flye_subreadset_outdir="${out_dir}/${batch}/${barcode}/flye_subreadset"
metaflye_ccs_outdir="${out_dir}/${batch}/${barcode}/metaflye_ccs"
metaflye_subreadset_outdir="${out_dir}/${batch}/${barcode}/metaflye_subreadset"
canu_ccs_outdir="${out_dir}/${batch}/${barcode}/canu_ccs"
canu_subreadset_outdir="${out_dir}/${batch}/${barcode}/canu_subreadset"
spades_ccs_outdir="${out_dir}/${batch}/${barcode}/spades_ccs"
circlator_ccs_outdir="${out_dir}/${batch}/${barcode}/circlator_ccs"

mkdir -p ${flye_ccs_outdir}
mkdir -p ${flye_subreadset_outdir}
mkdir -p ${metaflye_ccs_outdir}
mkdir -p ${metaflye_subreadset_outdir}
mkdir -p ${canu_ccs_outdir}
mkdir -p ${canu_subreadset_outdir}
mkdir -p ${spades_ccs_outdir}

echo "flye assembly on CCS"

source activate flye
date
flye --version
flye --pacbio-hifi ${ccs_fastq_path} --out-dir ${flye_ccs_outdir} -t ${SLURM_NTASKS}

echo "flye assembly on subreadset"
date
flye --version
flye --pacbio-raw ${subreadset_fastq_path} --out-dir ${flye_subreadset_outdir} -t ${SLURM_NTASKS}

echo "metaflye assembly on CCS"
date
flye --version
flye --meta --pacbio-hifi ${ccs_fastq_path} --out-dir ${metaflye_ccs_outdir} -t ${SLURM_NTASKS}

echo "metaflye assembly on subreadset"
date
flye --version
flye --meta --pacbio-raw ${subreadset_fastq_path} --out-dir ${metaflye_subreadset_outdir} -t ${SLURM_NTASKS}

source deactivate

echo "canu assembly on CCS"
date
/home/kve/bin/canu-2.2/bin/canu --version
/home/kve/bin/canu-2.2/bin/canu \
   -p assembly \
   -d ${canu_ccs_outdir} \
   genomeSize=1.5m \
   -pacbio-hifi ${ccs_fastq_path} \
   useGrid=false \
   batThreads="${SLURM_NTASKS}"

echo "canu assembly on subreadset"
date
/home/kve/bin/canu-2.2/bin/canu --version
/home/kve/bin/canu-2.2/bin/canu \
   -p assembly \
   -d ${canu_subreadset_outdir} \
   genomeSize=3m \
   -pacbio ${subreadset_fastq_path} \
   useGrid=false \
   batThreads="${SLURM_NTASKS}"

echo "SPAdes assembly on CCS"

source activate spades
date
spades.py --version
spades.py -t ${SLURM_NTASKS} -s ${ccs_fastq_path} -o ${spades_ccs_outdir}

source deactivate

echo "circlator on spades output"
source activate circlator
date
circlator version
circlator all --threads ${SLURM_NTASKS} ${spades_ccs_outdir}/scaffolds.fasta ${ccs_fastq_path} ${circlator_ccs_outdir}

date