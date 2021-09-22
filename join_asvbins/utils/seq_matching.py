"""Extract 16s from scaffolds"""

import os
import pandas as pd
import numpy as np
from skbio import write as write_fa
from skbio import read as read_fa
from skbio import Sequence


def fasta_to_df(path, headers=None):
    if headers is not None:
        df = pd.DataFrame({seq.metadata['id']:[seq.values]
                           for seq in read_fa(path, format='fasta')
                           if seq.metadata['id'] in headers})
    else:
        df = pd.DataFrame({seq.metadata['id']:[seq.values]
                           for seq in read_fa(path, format='fasta')})
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



def filter_fasta_from_headers(in_fasta_path, out_fasta_path:str, headers,
                      report_prath:str=None, show_report=False):
    headers = set(headers)
    if out_fasta_path is not None:
        write_fa((seq for seq in read_fa(in_fasta_path, format='fasta')
                  if seq.metadata['id'] in headers),
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
        joined['barrnap_header'] = f"{joined['header']}:{joined['start']}-{joined['stop']}"
    return joined


def process_barrnap(data:pd.DataFrame) -> pd.DataFrame:
    data.rename(columns={'header': 'barrnap_header'}, inplace=True)
    data[['header', 'start']] = data['barrnap_header'].str.split(':', expand=True)
    data[['start', 'stop']] = data['start'].str.split('-', expand=True)
    data[['start', 'stop']] = data[['start', 'stop']].astype(int)
    data = data.groupby('header').apply(merge_duplicate_seqs).\
        reset_index(drop=True)
    # extra assert statement to cover by back
    assert sum((data['stop'] - data['start']) != \
               data['seq'].apply(len)) == 0, \
        "The length of the sequence dose not match the size from the indexes."
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


def get_mbstats_dups(data):
    """Check mmseqs or blast stats data for dups"""
    dups = data[data['sseqid'].duplicated(keep=False)]
    dups = dups[[i for i in dups.columns
                 if i not in ['seq', 'sseqid', 'qseqid']]]
    return dups


def process_mbdata(data:pd.DataFrame, headers=None) -> pd.DataFrame:
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
    barseqs= barrnap[['header', 'seq']].copy()
    mbseqs = mbstats[['header', 'seq']].copy()
    mbseqs['mbseqs'] = True
    barseqs['barseqs'] = True
    # NOTE we have a choice for merge, we can convert the arays to tuples or
    # we can fix the problem after. I opt for option 2, it gives more control.
    data = pd.merge(barseqs, mbseqs, on='header',
                    how='outer')
    data[['mbseqs', 'barseqs']] = data[['mbseqs', 'barseqs']].fillna(False)
    data.rename(columns={'seq_x': 'seq_bar', 'seq_y': 'seq_other'},
                inplace=True)
    def select_and_describe_seq(x):
        if x['barseqs'] and not x['mbseqs']:
            return x['seq_bar'], 'Barnnap'
        elif x['mbseqs'] and not x['barseqs']:
            return x['seq_other'], 'Other'
        elif len(x['seq_bar']) > len(x['seq_other']):
            return x['seq_bar'], 'Barnnap>Other'
        elif len(x['seq_bar']) < len(x['seq_other']):
            return x['seq_bar'], 'Other>Barnnap'
        elif len(x['seq_bar']) == len(x['seq_other']):
            return x['seq_bar'], 'Other=Barnnap'
        else:
            raise Exception("Non equal duplicates in barseqs, and mbseqs")

    data[['seq', 'note']] = data.apply(select_and_describe_seq, axis=1, result_type='expand')
    print("After Merge: \n"
         f" There are {sum(data['barseqs'])} Sequences found by Barrnap.\n"
         f" There are {sum(data['mbseqs'])} Sequences found by"
          " MMseqs/BLAST.\n"
         f" There are {sum(data['mbseqs'] & data['barseqs'])}"
          " Sequences found by both Barrnap and MMseqs/BLAST.\n"
          )
    return data[['header', 'seq', 'note']]


def filter_mdstats(data, min_pct_id:float=None, min_length:int=None,
                    min_length_pct:float=None, max_gaps:int=None,
                    max_missmatch:int=None):
    # NOTE For perfect matches the Alignment Length equals the DB allele Length so the percent should be length/slen.
    data_checks = []
    if max_gaps is not None:
        data_checks.append(lambda x: x['gapopen'] <= max_gaps)
    if max_missmatch is not None:
        data_checks.append(lambda x: x['mismatch'] <= max_missmatch)
    if min_length is not None:
        data_checks.append(lambda x: x['length'] >= min_length)
    if min_pct_id is not None:
        data_checks.append(lambda x: x['pident'] >= min_pct_id)
    if min_length_pct is not None:
        data_checks.append(lambda x:
            ((x['length'] / x['slen']) * 100) >= min_length_pct)
    return data[
        data.apply(lambda x: all([l(x) for l in data_checks]),
        axis=1)
        ]

# data.apply( lambda x: ((x['gapopen'] <= max_gaps) if max_gaps is not None else True) and ((x['mismatch'] <= max_missmatch) if max_missmatch is not None else True) and ((x['length'] >= min_length) if min_length is not None else True) and ((x['pident'] >= min_pct_id) if min_pct_id is not None else True) and ((((x['length'] / x['slen']) * 100) >= min_length_pct) if min_length_pct is not None else True), axis=1)

def combine_mbstats_barrnap(mbstats_fasta_path:str, mbstats_stats_path:str,
                            barrnap_fasta_path:str, out_fasta_path:str,
                            min_pct_id:float=None, min_length:int=None):
    mbstats = read_mbstats(mbstats_stats_path)
    mbstats = filter_mdstats(mbstats,
                             min_pct_id = min_pct_id,
                             min_length = min_length)
    mbseqs = fasta_to_df(mbstats_fasta_path, mbstats['sseqid'].values)
    mbdata = pd.merge(mbseqs, mbstats, right_on='sseqid', left_on='header',
                        how='inner')
    mbdata = process_mbdata(mbdata)
    barrnap = fasta_to_df(barrnap_fasta_path)
    barrnap = process_barrnap(barrnap)
    data = combine_fasta(mbdata, barrnap)
    df_to_fasta(data, out_fasta_path)
    df_to_fasta(data, 'test.fa')


def filter_from_mbstats(stats_file:str, fasta_file_in:str, fasta_file_out:str,
                        min_pct_id:float=None, min_length:int=None,
                        min_length_pct:float=None, max_gaps:int=None,
                        max_missmatch:int=None):
        mbstats = read_mbstats(stats_file)
        mbstats = filter_mdstats(mbstats,
                                 min_pct_id=min_pct_id,
                                 min_length=min_length,
                                 min_length_pct=min_length_pct,
                                 max_gaps=max_gaps,
                                 max_missmatch=max_missmatch
                        )
        filter_fasta_from_headers(fasta_file_in,
                                  fasta_file_out, mbstats['sseqid'].values)





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

# from sys import getsizeof
# def human_readable_size(size, decimal_places=2):
#     for unit in ['B', 'KiB', 'MiB', 'GiB', 'TiB', 'PiB']:
#         if size < 1024.0 or unit == 'PiB':
#             break
#         size /= 1024.0
#     return f"{size:.{decimal_places}f} {unit}"
# human_readable_size(getsizeof(mbseqs))

# combine_mbstats_barrnap(input['other_stats_path'],
#                         input['barrnap_fasta_path'], output['out_fasta_path'])
# combine_mbstats_barrnap("stage1_asvs_mmseqs.tab", "barrnap_fastafile-16S.fna")
#         combine_mbstats_barrnap(input['other_fasta_path'], input['other_stats_path'],
#                 input['barrnap_fasta_path'], output['out_fasta_path'])
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
