from pathlib import Path
import pandas as pd
from Bio import SeqIO

wd = Path("/nfs/chisholmlab001/kve/2023_HIFI_Genome_Closing_Project")

MIN_GENOME_LENGTH = 1000000

assemblies = (wd / "assembly" / "all_assemblies").glob("*.fasta")
classifications = (wd / "classification" / "mmseqs2_reports").glob("*.tsv")

assemblies = [
    dict(zip(
        ['batch','barcode','assembly method','assembly path'],
        assembly.name.split('.')[:3] + [assembly]
    ))
    for assembly in assemblies
]

classification_df = pd.concat([
    pd.read_table(classification,names=[
        "contig",
        "taxid",
        "taxonomic level",
        "taxonomic label",
        "# protein fragments",
        "# labeled protein fragments",
        "# fragments agree with assigned label",
        "neg log E val"
    ]).assign(**dict(zip(
        ['batch','barcode','assembly method'],
        classification.name.split('.')[:3]
    )))
    for classification in classifications
]).set_index(['batch','barcode','assembly method','contig'])

data = []
for i, assembly in enumerate(assemblies):
    print(assembly['batch'], assembly['barcode'], assembly['assembly method'], sep='\t')
    assembly_wd = wd / "assembly" / assembly['batch'] / assembly['barcode'] / assembly['assembly method']
    seq_dict = SeqIO.to_dict(SeqIO.parse(assembly['assembly path'], "fasta"))
    if 'canu' in assembly['assembly method']:
        info_df = pd.read_table(assembly_wd / "assembly.contigs.layout.tigInfo")
        info_df = info_df[info_df["tigClass"]=='contig']
        contigs = info_df["#tigID"].apply(lambda n: f"tig{n:08}")
        circularity = info_df['sugCirc'] == 'yes'
    elif 'flye' in assembly['assembly method']:
        info_df = pd.read_table(assembly_wd / "assembly_info.txt")
        contigs = info_df['#seq_name']
        circularity = info_df['circ.']=='Y'
    elif 'circlator' in assembly['assembly method']:
        info_df = pd.read_table(assembly_wd / "04.merge.circularise.log")
        info_df = info_df[info_df['#Contig'].isin(seq_dict.keys())]
        contigs = info_df['#Contig']
        circularity = info_df['circularised']==1
    else:
        raise RuntimeError('Unrecognized Assembly Method')

    df = pd.concat([contigs, circularity],axis=1)
    df.columns = ['contig','circularized']
    df['length'] = contigs.apply(lambda contig: len(seq_dict[contig]))
    df = df.assign(**assembly)
    data.append(df)

df = pd.concat(data).set_index(['batch','barcode','assembly method','contig'])
df = df.join(classification_df)
df.to_csv('contig_data.tsv',sep='\t')

df[df['circularized'] & (df['length'] >= MIN_GENOME_LENGTH)].to_csv('contig_data_filtered.tsv',sep='\t')