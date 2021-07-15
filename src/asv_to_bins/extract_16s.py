"""Extract 16s from scaffolds"""

import os
import pandas as pd
import numpy as np
from skbio import write as write_fa
from skbio import read as read_fa
from skbio import Sequence

def read_and_process_blast(fasta_path, data_path):
    pass
fasta_path = "../../results/original_approach/blast_fastafile-16S.fna"

def fasta_to_df(path):
    df = pd.DataFrame({i.metadata['id']:[i.values] for i in read_fa(path, format='fasta')})
    df = df.T
    df.reset_index(inplace=True)
    df.columns = ['header', 'seq']
    return df

def merge_duplicate_seqs(data:pd.DataFrame) -> pd.DataFrame:
    data.sort_values('start')
    joined = data.iloc[0]
    if len(data) <= 1:
        return joined
    for _, to_join in data.iloc[1:].iterrows():
        if to_join['start'] > joined['stop']:
            raise ValueError("Found non-overlapping duplicate in barrnap data")
        trim_point = joined['stop'] - to_join['start']
        joined['seq'] = np.concatenate([joined['seq'], to_join['seq'][trim_point:]])
        joined['stop'] = to_join['stop']
        joined['header'] = "%s-%i" % (joined['header'].split('-')[0], to_join['stop'])
    return joined


def process_barrnap(data:pd.DataFrame) -> pd.DataFrame:
    data[['name', 'start']] = data['header'].str.split(':', expand=True)
    data[['start', 'stop']] = data['start'].str.split('-', expand=True)
    data[['start', 'stop']] = data[['start', 'stop']].astype(int)
    data = data.groupby('name').apply(merge_duplicate_seqs)
    # extra assert statement to cover by back
    assert sum((data['stop'] - data['start']) != data['seq'].apply(len)) == 0, \
        "The length of the sequence dose not match the size from the indexes."
    return data


def test_barnap_procesing():
    input_df = pd.DataFrame({
        1: ['a:0-3',       'abc'],
        2: ['b:0-3',       'abc'],
        3: ['b:3-4',       'd'],
        4: ['b:3-7',       'defg'],
        5: ['c:1002-1007', 'abcde'],
        6: ['c:1003-1010', 'bcdefghij'],
        7: ['d:0-3',       'abcd'],
    }, index=['header', 'seq']).T
    expect_df = pd.DataFrame({
        'a': ['a:0-3',       'abc'],
        'b': ['b:0-7',       'abcdefg'],
        'c': ['c:1002-1010', 'abcdefghij'],
        'd': ['d:0-3',       'abcd'],
    }, index=['header', 'seq']).T
    input_df['seq'] = input_df['seq'].apply(list)
    expect_df['seq'] = expect_df['seq'].apply(list)
    output_df = process_barrnap(input_df)
    output_df = output_df[['header', 'seq']]
    pd.testing.assert_frame_equal(output_df, expect_df, check_index_type=False, check_names=False)


def read_and_process_barrnap(fasta_path):
    data = fasta_to_df(fasta_path)
    data = process_barrnap(data)
    return data

def read_and_process_blast(fasta_path, stats_path):
    pass
fasta = fasta_to_df(fasta_path)
# get headers here:
# https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
stats = pd.read_csv(stats_path, header=None, sep='\t', names=[
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
data = pd.merge(fasta, stats, right_on='sseqid', left_on='header', how='inner')
assert len(fasta) == len(data), "Something in the blast fasta dose not match" \
    "the BLASTn stats."
len(fasta)
len(data)
len(stats)
fasta[]
c
fasta_path = "../../results/original_approach/blast_fastafile-16S.fna"
stats_path = "../../results/original_approach/blast_fastafile-16S.txt"
os.system("wc -l ../../results/original_approach/rrna.gff")
os.system("cat ../../results/original_approach/blast_fastaheaders-16S.txt")
with open('../../results/original_approach/blast_fastaheaders-16S.txt') \
        as blaster:
    s = blaster.read()

s = s.split('\n')
len([i for i in headers if i in set(s)])

len(headers)
headers
os.system("grep All_IDBACoAss_bin.260_scaffold_629 ../../results/original_approach/*")
''
bla_path = "../../results/original_approach/blast_fastafile-16S.fna"
bar_path = "../../results/original_approach/barrnap_fastafile-16S.fna"



bla_data_path = "../../results/original_approach/blast_fastafile-16S.txt"
blast_fasta_path = "../../results/original_approach/blast_fastafile-16S.fna"



headers = set(['All_MegahitCoAss_bin.3_k141_331864',
           'All_IDBACoAss_bin.239_scaffold_202',
           'LM_meg_bin.137_k119_157622', 'LM_bin.60_scaffold_9911',
           'LO_scaffold.metabat-bins-_--verysensitive.129.fa_scaffold_7869',
           'LM_bin.141_scaffold_10913', 'All_IDBACoAss_bin.260_scaffold_629',
           'LO_scaffold.metabat-bins-_--verysensitive.104.fa_scaffold_2622',
           'LM_bin.91_scaffold_3243',
           'LP_scaffold.metabat-bins-_--verysensitive.56.fa_scaffold_275',
           'LP_scaffold.metabat-bins-_--verysensitive.121.fa_scaffold_13505',
           'LP_scaffold.metabat-bins-_--verysensitive.121.fa_scaffold_1582',
           'All_IDBACoAss_bin.143_scaffold_2701',
           'All_MegahitCoAss_bin.116_k141_377345',
           'KS_bin.91_scaffold_2441', 'KS_meg_bin.64_k141_265517',
           'KL_meg_bin.6_k141_250169',
           'KS_bin.40_scaffold_578', 'All_IDBACoAss_bin.59_scaffold_864',
           'LP_scaffold.metabat-bins-_--verysensitive.143.fa_scaffold_13218',
           'LM_bin.134_scaffold_10792', 'LO_meg_bin.103_k119_131170',
           'LM_meg_bin.96_k119_50248', 'All_MegahitCoAss_bin.206_k141_54051'])

bar_data = read_and_process_barrnap(bar_path)
d1 = bar_data[['name']]
d1.reset_index(inplace=True, drop=True)
d1['in_bar'] = True
d2 = pd.DataFrame({'name':list(headers)})
d2['in_hed'] = True
pd.merge(d1, d2, how='outer', on='name').sort_values('name')['name']
"""
We can read and process blas
"""
# read blast
blast_data = pd.read_csv(blast_data_path,
    sep='\t', header=None)
blast_fasta = [i for i in read_fa(blast_fasta_path, format='fasta')
               if i.metadata['id'] in headers]


os.system('man blastn')


os.system('head test.fasta ')
for i in fasta:
    write_fa(i, format='fasta', into='test.fasta')

barrnap_16s = open("../results/original_approach/barrnap_fastafile-16S.fna", "r")
{}

barrnap_16s_string = barrnap_16s.read()

#convert barrnap fasta to list
barrnap_16s_list = barrnap_16s_string.split('\n')
barrnap_16s.close()

#read in blast scaffold fasta
blast_fasta = open(, "r")
blast_fasta_string = blast_fasta.read()

process_blast
read_and_process_barrnap

test
