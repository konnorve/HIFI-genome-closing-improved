from pathlib import Path
import json
import gzip
import shutil

wd = Path("/nfs/chisholmlab001/kve/2023_HIFI_Genome_Closing_Project/demultiplexing")

(wd / 'all').mkdir(exist_ok=True, parents=True)
# (wd / 'concatenations').mkdir(exist_ok=True, parents=True)

hifi_read_paths = wd.glob("*/*/ccs_demux/outputs/*.bc*--bc*.hifi_reads.fastq.gz")
subread_paths = wd.glob("*/*/lima_demux/demux_subreadset.bc*--bc*.fastq")

hifi_reads = [
    {
        "batch" : hifi_read_path.parent.parent.parent.parent.name,
        "seq_run" : hifi_read_path.parent.parent.parent.name,
        "barcode" : hifi_read_path.name.split('.')[1].split('--')[0],
        "path" : str(hifi_read_path),
        "gzipped" : True,
        "type" : 'hifi_readset'
    }
    for hifi_read_path in hifi_read_paths
]

subreads = [
    {
        "batch" : subread_path.parent.parent.parent.name,
        "seq_run" : subread_path.parent.parent.name,
        "barcode" : subread_path.name.split('.')[1].split('--')[0],
        "path" : str(subread_path),
        "gzipped" : False,
        "type" : 'subreadset'
    }
    for subread_path in subread_paths
]

reads = hifi_reads + subreads

# combine batch 8 sequencing runs

to_concat = {}

# identify indicies in batch 8 with same barcode
i = 0
while i < len(reads):
    read = reads[i]
    # batch 8 was ran twice because of low coverage on the first run
    if read['batch'] == 'batch8':
        ident = f"{read['batch']}.comb.{read['barcode']}.{read['type']}"
        if ident in to_concat.keys():
            to_concat[ident].append(reads.pop(i))
        else:
            to_concat[ident] = [reads.pop(i)]
    elif read['batch'] == 'batch7' and read['barcode'] in ['bc1012', 'bc1022']:
        ident = f"{read['batch']}.{read['seq_run']}.bc1022.{read['type']}"
        if ident in to_concat.keys():
            to_concat[ident].append(reads.pop(i))
        else:
            to_concat[ident] = [reads.pop(i)]
    else:
        i += 1

# for ident, group in to_concat.items():
#     new_file = wd / 'all' / f"{ident}.fastq"
#     with open(new_file, 'w') as whole_file:
#         for part in group:
#             if part['gzipped']:
#                 with gzip.open(part['path'], 'rb') as file_part:
#                     for line in file_part:
#                         whole_file.write(line.decode())
#             else:
#                 with open(part['path']) as file_part:
#                     for line in file_part:
#                         whole_file.write(line)

for read in reads:
    ident = f"{read['batch']}.{read['seq_run']}.{read['barcode']}.{read['type']}"
    new_file = wd / 'all' / f"{ident}.fastq"
    if read['gzipped']:
        continue
        # with gzip.open(read['path'], 'rb') as f_in:
        #     with open(new_file, 'wb') as f_out:
        #         shutil.copyfileobj(f_in, f_out)
    else:
        new_file.hardlink_to(read['path'])