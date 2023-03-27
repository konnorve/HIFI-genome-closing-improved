# from pathlib import Path
import pandas as pd
# from Bio import SeqIO

df = pd.read_table("contig_data_filtered.tsv")
metadata = pd.read_table("metadata.tsv",index_col=['batch','barcode'])
metadata.columns = pd.MultiIndex.from_product([['metadata'], metadata.columns])

df = df[['batch', 'barcode', 'assembly method', 'contig', 'taxonomic label', 'length', 'assembly path']]
df = df.groupby(['batch', 'barcode', 'assembly method']).agg(lambda x: ' | '.join(x.to_list())).reset_index()
df = df.pivot(index=['batch','barcode'], columns='assembly method', values=['taxonomic label','contig','length','assembly path'])

df = metadata.join(df)

df[('statistics','num assemblers completed')] = df['contig'].count(axis=1)
df[('statistics','genome closed?')] = df[('statistics','num assemblers completed')] > 0
df[('statistics','>1 genome closed')] = df['taxonomic label'].apply(lambda r: any(['|' in x for x in r.to_list() if type(x)==str]),axis=1)

df = df.sort_index(level=['batch','barcode'])

df.to_csv("contig_data_pivoted.tsv", sep='\t')
