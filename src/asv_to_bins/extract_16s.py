"""Extract 16s from scaffolds"""

import os
import pandas as pd
import numpy as np
from skbio import write as write_fa
from skbio import read as read_fa
from skbio import Sequence


def fasta_to_df(path):
    df = pd.DataFrame({i.metadata['id']:[i.values]
                       for i in read_fa(path, format='fasta')})
    df = df.T
    df.reset_index(inplace=True)
    df.columns = ['header', 'seq']
    return df


def df_to_fasta(df:pd.DataFrame, path:str):
    seqs = (Sequence(x['seq'],
                     metadata={"id": x["header"],
                               'description': x["note"]}
                     )
            for _, x in df.iterrows())
    write_fa(seqs, 'fasta', path)



def filter_fasta_from_headers(in_fasta_path, out_fasta_path, headers,
                      report_prath:str=None, show_report=False):
    headers = set(headers)
    write_fa((seq for seq in read_fa(in_fasta_path, format='fasta') if seq.metadata['id'] in headers),
             'fasta', out_fasta_path)

def merge_duplicate_seqs(data:pd.DataFrame) -> pd.DataFrame:
    data.sort_values('start')
    joined = data.iloc[0]
    if len(data) <= 1:
        return joined
    for _, to_join in data.iloc[1:].iterrows():
        if to_join['start'] > joined['stop']:
            # raise Warning("Found non-overlapping duplicate in barrnap data")
            continue
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
        7: ['d:0-3',       'abcd'], }, index=['header', 'seq']).T
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

def read_mbstats(stats_path:str) -> pd.DataFrame:
    stats = pd.read_csv(stats_path, header=None, sep='\t', names=[
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"])
    return stats

def read_mbstats_and_fasta(fasta_path:str, stats_path:str) -> pd.DataFrame:
    fasta = fasta_to_df(fasta_path)
    # get headers here:
    # https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    stats = read_mbstats(stats_path)
    data = pd.merge(fasta, stats, right_on='sseqid', left_on='header',
                    how='inner')
    return data


def read_and_process_mbstats(fasta_path:str, stats_path:str, headers=None):
    # Note get stats on balst evalue and lenther so that  we get cutofs
    # Barnnap
    data = read_mbstats_and_fasta(fasta_path, stats_path)
    data = process_mbstats(data, headers=headers)
    return data


def get_mbstats_dups(data):
    """Check mmseqs or blast stats data for dups"""
    dups = data[data['sseqid'].duplicated(keep=False)]
    dups = dups[[i for i in dups.columns
                 if i not in ['seq', 'sseqid', 'qseqid']]]
    return dups


def process_mbstats(data:pd.DataFrame, headers=None) -> pd.DataFrame:
    # TODO you may want to use this again
    # get_mbstats_dups(data)
    # Remove mbstats data dups keeping the longest string
    data = data.sort_values('length', ascending=False).\
        drop_duplicates('sseqid')
    # Keep only elements of blast_fasta_list that are in headers object
    if headers is not None:
        data = pd.merge(data, pd.DataFrame({
            'sseqid': list(headers)}), on='sseqid')
    data['seq'] = data.apply(
        lambda x: x["seq"][x["sstart"] - 1:x["send"]]
        if x["sstart"] < x["send"]
        else x["seq"][::-1][x["send"] - 1:x["sstart"]], axis=1)
    return data


def combine_fasta(mbstats:pd.DataFrame, barrnap:pd.DataFrame) -> pd.DataFrame:
    barrnap = barrnap[['header', 'seq']]
    mbstats = mbstats[['header', 'seq']]
    mbstats['mbstats'] = True
    barrnap['barrnap'] = True
    # NOTE we have a choice for merge, we can convert the arays to tuples or
    # we can fix the problem after. I opt for option 2, it gives more control.
    data = pd.merge(barrnap, mbstats, on='header', how='outer')
    data[['mbstats', 'barrnap']] = data[['mbstats', 'barrnap']].fillna(False)
    data.rename(columns={'seq_x': 'seq_bar', 'seq_y': 'seq_bla'},
                inplace=True)
    data['seq'] = data.apply(lambda x:
                    x['seq_bar'] if x['barrnap'] and not x['mbstats']
               else x['seq_bla'] if x['mbstats'] and not x['barrnap']
               else x['seq_bar'] if np.array_equal(x['seq_bar'], x['seq_bla'])
               else raise_(
                   Exception("Non equal duplicates in barrnap, and mbstats")),
        axis=1)
    data[['seq', 'note']] = data.apply(
        lambda x:
                    (x['seq_bar'], 'Barnnap')
               if x['barrnap'] and not x['mbstats']
               else (x['seq_bla'], 'BLAST')
               if x['mbstats'] and not x['barrnap']
               else (x['seq_bar'], 'Barnnap+BLAST')
               if (len(x['seq_bar']) >= len(x['seq_bla']))
               else (x['seq_bar'], 'BLAST+Barnnap'),
        axis=1, result_type='expand')
    print("After Merge: \n"
          " There are %i Sequences found by Barrnap.\n"
          " There are %i Sequences found by BLAST.\n"
          " There are %i Sequences found by both Barrnap and BLAST.\n" \
          % (sum(data['barrnap']),
             sum(data['mbstats']),
             sum(data['mbstats'] & data['barrnap']))
          )
    return data[['header', 'seq', 'note']]


def combine_mbstats_barrnap(mbstats_fasta_path:str, mbstats_stats_path:str,
                          barrnap_fasta_path:str, out_fasta_path:str):
    mbstats = read_and_process_mbstats(mbstats_fasta_path, mbstats_stats_path)
    barrnap = read_and_process_barrnap(barrnap_fasta_path)
    data = combine_fasta(mbstats, barrnap)
    df_to_fasta(data, out_fasta_path)




# DONE s1_min_pct_id=s1_min_pct_id"
# DONE s2_min_pct_id=s2_min_pct_id"
# DONE s1_min_length=s1_min_length"
# DONE s2_min_length=s2_min_length"
# DONE max_missmatch=max_missmatch"
# DONE min_length_pct=min_length_pct"
# DONE max_gaps=max_gaps"
# TODO formalize these tests
#     mbstats.loc['gapopen', 1]
#     mbstats.loc[['gapopen', 212129]]
#     mbstats.loc[212129, 'gapopen'] = 1
#     filter_data_set(mbstats)
#     filter_data_set(mbstats, max_gaps=1)
#     mbstats['mismatch']
#     filter_data_set(mbstats, max_missmatch=7)
#     mbstats['length']
#     filter_data_set(mbstats, min_length=49)
#     mbstats['sseqid']
#     mbstats['length'] / mbstats['slen']
#     filter_data_set(mbstats, min_length_pct=1)
#     mbstats['pident']
#     filter_data_set(mbstats, min_pct_id=84)
#     filter_data_set(mbstats, min_pct_id=84, max_gaps=0)

def filter_mdstats(data, min_pct_id:float=None, min_length:int=None,
                    min_length_pct:float=None, max_gaps:int=None,
                    max_missmatch:int=None):

    # NOTE For perfect matches the Alignment Length equals the DB allele Length so the percent should be length/slen.
    return data[
        data.apply(
        lambda x:
            ((x['gapopen'] <= max_gaps)
             if max_gaps is not None else True) and
            ((x['mismatch'] <= max_missmatch)
             if max_missmatch is not None else True) and
            ((x['length'] >= min_length)
             if min_length is not None else True) and
            ((x['pident'] >= min_pct_id)
             if min_pct_id is not None else True) and
            ((((x['length'] / x['slen']) * 100) >= min_length_pct)
             if min_length_pct is not None else True),
        axis=1)
        ]

def filter_to_lenth(in_data):
    """Filter to 100% length"""
    data = in_data.copy()
    return data[data['qlen'] <= data['slen']]

def filter_to_gaps(in_data):
    """Filter to 0 gaps"""
    data = in_data.copy()
    return data[data['gapopen'] <= 0]

def filter_to_mismatch(in_data, mismatch):
    """Filter to below and given number of mismatches"""
    data = in_data.copy()
    return data[data <= mismatch]

# blast_fasta_path = "../../results/original_approach/blast_fastafile-16S.fna"
# blast_stats_path = "../../results/original_approach/blast_fastafile-16S.txt"
# barrnap_fasta_path = \
#     "../../results/original_approach/barrnap_fastafile-16S.fna"
# out_test_path = "../../results/original_approach/new_output.fa"
# blast = read_and_process_blast(blast_fasta_path, blast_stats_path)
# blast_headers = read_and_process_blast(blast_fasta_path, blast_stats_path,
# headers=HEADERS)
# barrnap = read_and_process_barrnap(barrnap_fasta_path)
# data = combine_fasta(blast, barrnap)
# data = combine_fasta(blast_headers, barrnap)
# df_to_fasta(data, out_test_path)
