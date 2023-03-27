#!/usr/bin/env bash
#SBATCH --job-name=assem
#SBATCH --time 10-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 1                                # tasks (essentially threads)
#SBATCH --array=1 # -440%100
#SBATCH -o logs/classification.%a.out
#SBATCH -e logs/assembly.%a.err

date

assembly=$(realpath /nfs/chisholmlab001/kve/2023_HIFI_Genome_Closing_Project/assembly/all_assemblies/* | sed -n ${SLURM_ARRAY_TASK_ID}p)
name=$(basename ${assembly})
name=${name%.*}

out_dir="/nfs/chisholmlab001/kve/2023_HIFI_Genome_Closing_Project/classification/kaiju"
kaiju_bin="/home/kve/scripts/HIFI-genome-closing-improved/classification/kaiju/bin"
kaiju_db="/nobackup1/chisholmlab/kve/2023_HIFI_Genome_Closing_Project/kaiju_db/mar_20220201"

mkdir -p ${out_dir}

result_txt=${out_dir}/"${name}.output.txt"
result_table=${out_dir}/"${name}.output.tsv"
result_names=${out_dir}/"${name}.output_named.tsv"
result_name_class=${out_dir}/"${name}.tsv"

kaiju_index_file=${kaiju_db}/"refseq/kaiju_db_refseq.fmi"
nodes_dmp=${kaiju_db}/"nodes.dmp"
names_dmp=${kaiju_db}/"names.dmp"

${kaiju_bin}/kaiju

${kaiju_bin}/kaiju -v -t ${nodes_dmp} -f ${kaiju_index_file} -i ${assembly} -o ${result_txt} -z ${SLURM_NTASKS}

${kaiju_bin}/kaiju2table -t ${nodes_dmp} -n ${names_dmp} -r genus -o ${result_table} ${result_txt}

${kaiju_bin}/kaiju-addTaxonNames -t ${nodes_dmp} -n ${names_dmp} -i ${result_txt} -o ${result_names}

# removed unclassified contigs
sed '/^U/d' ${result_names} > ${result_name_class}

rm ${result_txt} ${result_table} ${result_names}

date