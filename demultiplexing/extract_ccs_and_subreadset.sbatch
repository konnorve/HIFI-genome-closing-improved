#!/usr/bin/env bash
#SBATCH --job-name=demultiplex
#SBATCH --time 5-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 20                                # tasks (essentially threads)
#SBATCH -N 1
#SBATCH --mem 250G
#SBATCH --array=2-10%5
#SBATCH -o logs/demultiplex.%a.out
#SBATCH -e logs/demultiplex.%a.err

export PATH="/home/kve/bin/smrtlink/smrtcmds/bin:$PATH"

barcodes_xml_path="/nfs/chisholmlab001/chisholmlab/experiment_repository/2021/HIFI_Genome_Closing_Project/analysis/raw_barcodes/Sequel_384_barcodes_v1.xml"
out_dir="/nfs/chisholmlab001/kve/2023_HIFI_Genome_Closing_Project/demultiplexing"

batch=$(sed "${SLURM_ARRAY_TASK_ID}q;d" demultiplex.tsv | cut -f1 | tr -d '\r')
sequencing_run=$(sed "${SLURM_ARRAY_TASK_ID}q;d" demultiplex.tsv | cut -f2 | tr -d '\r')
subreadset_xml=$(sed "${SLURM_ARRAY_TASK_ID}q;d" demultiplex.tsv | cut -f3 | tr -d '\r')

ccs_demux_outdir="${out_dir}/${batch}/${sequencing_run}/ccs_demux"
lima_demux_outdir="${out_dir}/${batch}/${sequencing_run}/lima_demux"

mkdir -p ${ccs_demux_outdir}
mkdir -p ${lima_demux_outdir}

echo ${batch}
echo ${sequencing_run}
echo ${subreadset_xml}
ls ${ccs_demux_outdir}
ls ${lima_demux_outdir}
ls ${barcodes_xml_path}

echo "run pacbio ccs demux workflow"

echo "smrtlink_11.1"
pbcromwell --version
pbcromwell run pb_ccs_demux -e ${subreadset_xml} -e ${barcodes_xml_path} --output-dir ${ccs_demux_outdir} --overwrite --nproc ${SLURM_NTASKS}

echo "run lima subreadset demux workflow"

lima --version
lima --same --peek-guess --split-named --num-threads ${SLURM_NTASKS} ${subreadset_xml} ${barcodes_xml_path} ${lima_demux_outdir}/demux_subreadset.bam

echo "convert lima bams to fastqs"
source activate samtools

samtools --version

for bam_file in ${lima_demux_outdir}/*.bam
do
    samtools fastq ${bam_file} > ${bam_file%.*}.fastq
done

