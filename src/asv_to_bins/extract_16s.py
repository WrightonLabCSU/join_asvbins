"""Extract 16s from scaffolds"""

import os
import pandas as pd
import numpy as np
from skbio import write as write_fa
from skbio import read as read_fa
from skbio import Sequence

# NOTE, what and why
HEADERS = set(['All_MegahitCoAss_bin.3_k141_331864',
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

def read_and_process_blast(fasta_path, data_path):
    pass

def fasta_to_df(path):
    df = pd.DataFrame({i.metadata['id']:[i.values]
                       for i in read_fa(path, format='fasta')})
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
        joined['seq'] = np.concatenate([joined['seq'],
                                        to_join['seq'][trim_point:]])
        joined['stop'] = to_join['stop']
        joined['header'] = "%s-%i" % (joined['header'].split('-')[0],
                                      to_join['stop'])
    return joined


def process_barrnap(data:pd.DataFrame) -> pd.DataFrame:
    data[['name', 'start']] = data['header'].str.split(':', expand=True)
    data[['start', 'stop']] = data['start'].str.split('-', expand=True)
    data[['start', 'stop']] = data[['start', 'stop']].astype(int)
    data = data.groupby('name').apply(merge_duplicate_seqs)
    # extra assert statement to cover by back
    assert sum((data['stop'] - data['start']) != \
               data['seq'].apply(len)) == 0, \
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
    pd.testing.assert_frame_equal(output_df, expect_df,
                                  check_index_type=False, check_names=False)


def read_and_process_barrnap(fasta_path):
    data = fasta_to_df(fasta_path)
    data = process_barrnap(data)
    return data


def read_blast(fasta_path:str, stats_path:str) -> pd.DataFrame:
    fasta = fasta_to_df(fasta_path)
    # get headers here:
    # https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    stats = pd.read_csv(stats_path, header=None, sep='\t', names=[
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    # There will still be duplicates
    data = pd.merge(fasta, stats, right_on='sseqid', left_on='header',
                    how='inner')
    return data


def read_and_process_blast(fasta_path, stats_path, headers=None):
    data = read_blast(fasta_path, stats_path)
    data = process_blast(data, headers=headers)
    return data


def get_blast_dups(data):
    """Check blast data for dups"""
    dups = data[data['header'].duplicated(keep=False)]
    dups = dups[[i for i in dups.columns
                 if i not in ['seq', 'sseqid', 'qseqid']]]
    return dups


def process_blast(data:pd.DataFrame, headers=None) -> pd.DataFrame:
    get_blast_dups(data)
    # Remove blast data dups keeping the longest string
    # NOTE unlike barrnap blast dubs don't overlap so i do it this way
    data = data.sort_values('length', ascending=False).\
        drop_duplicates('header')
    # Keep only elements of blast_fasta_list that are in headers object
    if headers is not None:
        data = pd.merge(data, pd.DataFrame({
            'header': list(headers)}), on='header')
    # Trim scaffold sequences to length of 16S from blast_16s
    # NOTE I don't know:
    #      Why the sequences need to be reversed sometimes
    #      why they need to be reversed
    # NOTE I failed to fined the meaning of the header
    data['seq'] = data.apply(
        lambda x: x["seq"][x["sstart"] - 1:x["send"]]
        if x["sstart"] < x["send"]
        else x["seq"][::-1][x["send"] - 1:x["sstart"]], axis=1)
    return data


def combine_blast_fasta(blast:pd.DataFrame, barrnap:pd.DataFrame,
                        headers=None) -> pd.DataFrame:
    barrnap = barrnap[['header', 'seq']]
    blast = blast[['header', 'seq']]
    blast['blast'] = True
    barrnap['barrnap'] = True
    # NOTE we have a choice for merge, we can convert the arays to tuples or
    # we can fix the problem after. I opt for option 2, it gives more control.
    data = pd.merge(barrnap, blast, on='header', how='outer')
    data[['blast', 'barrnap']] = data[['blast', 'barrnap']].fillna(False)
    data.rename(columns={'seq_x': 'seq_bar', 'seq_y': 'seq_bla'},
                inplace=True)
    data.apply(lambda x:
                    x['seq_bar'] if x['barrnap'] and not x['blast']
               else x['seq_bla'] if x['blast'] and not x['barrnap']
               else x['seq_bar'] if np.array_equal(x['seq_bar'], x['seq_bla'])
               else raise_(
                   Exception("Non equal duplicates in barrnap, and blast")),
        axis=1)
    print("After Merge: \n"
          " There are %i Sequences found by Barrnap.\n"
          " There are %i Sequences found by BLAST.\n"
          " There are %i Sequences found by both Barrnap and BLAST.\n" \
          % (sum(data['barrnap']),
             sum(data['blast']),
             sum(data['blast'] & data['barrnap']))
          )
    return data

# blast_fasta_path = "../../results/original_approach/blast_fastafile-16S.fna"
# blast_stats_path = "../../results/original_approach/blast_fastafile-16S.txt"
# barrnap_fasta_path = \
#     "../../results/original_approach/barrnap_fastafile-16S.fna"
# blast = read_and_process_blast(blast_fasta_path, blast_stats_path)
# barrnap = read_and_process_barrnap(barrnap_fasta_path)
# data = combine_blast_fasta(blast, barrnap)
# data = combine_blast_fasta(blast, barrnap, headers=HEADERS)
