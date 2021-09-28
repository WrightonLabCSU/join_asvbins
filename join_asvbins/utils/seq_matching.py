"""Tools for extract 16s from scaffolds"""
import os
import pandas as pd
import numpy as np
from skbio import write as write_fa
from skbio import read as read_fa
from skbio import Sequence


def fasta_to_df(path, headers=None):
    """
    Convert a fasta to a dataframe

    :param path: A path to a fasta
    :param headers: Optional, headers to read from fasta
    :returns: A dataframe
    """
    if headers is not None:
        dafr = pd.DataFrame({seq.metadata['id']:[seq.values]
                           for seq in read_fa(path, format='fasta')
                           if seq.metadata['id'] in headers})
    else:
        dafr = pd.DataFrame({seq.metadata['id']:[seq.values]
                           for seq in read_fa(path, format='fasta')})
    dafr = dafr.T
    dafr.reset_index(inplace=True)
    dafr.columns = ['header', 'seq']
    return dafr


def df_to_fasta(dafr:pd.DataFrame, path:str):
    """
    Convert a dataframe to a fasta

    :param dafr: A dataframe containing 'seq', 'note' and 'header' fields
    :param path: A path to a fasta
    """
    seqs = (Sequence(x['seq'],
                     metadata={"id": x["header"],
                               'description': x["note"]}
                     )
            for _, x in dafr.iterrows())
    write_fa(seqs, 'fasta', path)


# TODO Look into how you type hint the headers
# TODO handel empty files
def filter_fasta_from_headers(in_fasta_path:str, out_fasta_path:str, headers):
    """

    Filter a fast to a list of headers

    :param in_fasta_path: Path to unfilterd fasta
    :param out_fasta_path: Path to filtered fast
    :param headers: Headers to filter by
    """
    headers = set(headers)
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


def process_barfasta(data:pd.DataFrame) -> pd.DataFrame:
    data.rename(columns={'header': 'barrnap_header'}, inplace=True)
    data[['header', 'start']] = \
        data['barrnap_header'].str.split(':', expand=True)
    data[['start', 'stop']] = data['start'].str.split('-', expand=True)
    data[['start', 'stop']] = data[['start', 'stop']].astype(int)
    data = data.groupby('header').apply(merge_duplicate_seqs).\
        reset_index(drop=True)
    # extra assert statement to cover by back
    assert sum(
        (data['stop'] - data['start']) != data['seq'].apply(len)
        ) == 0, \
        "The length of the sequence dose not match the size from the indexes."
    return data


def read_mbstats(stats_path:str) -> pd.DataFrame:
    """
    Read the tab delimited mmseqs or blast file with its very specific format.

    :param stats_path: The path to the formatted statistics
    :returns: A dataframe with proper format
    """
    stats = pd.read_csv(stats_path, header=None, sep='\t', names=[
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"])
    return stats


def get_mbstats_dups(data):
    """
    Check mmseqs or blast stats data for dups

    :param data:
    :returns:
    """
    dups = data[data['sseqid'].duplicated(keep=False)]
    dups = dups[[i for i in dups.columns
                 if i not in ['seq', 'sseqid', 'qseqid']]]
    return dups


def process_mbdata(data:pd.DataFrame, headers=None) -> pd.DataFrame:
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
    """
    Creates and then applies a filter for mmseqs or blast statistics

    :param data: Data to be filter must be in a computable blast style format
    :param min_pct_id: Optional filter
    :param min_length: Optional filter
    :param min_length_pct: Optional filter
    :param max_gaps: Optional filter
    :param max_missmatch: Optional filter
    :returns: Filtered data
    """
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


#TODO Add a handler for the case where blast/mmseqs or barrnap don't return matches
def combine_mbstats_barrnap(mbstats_fasta_path:str, mbstats_stats_path:str,
                            barrnap_fasta_path:str, out_fasta_path:str,
                            out_barstats_path:str, min_pct_id:float=None,
                            min_length:int=None, search_tool:str='Other'):
    """
    Combine the statistics from mmseqs or blast with the output from barrnap.

    :param mbstats_fasta_path: Path to mmseqs or blast fasta
    :param mbstats_stats_path: Path to mmseqs or blast statistics, see docs for format requirements
    :param barrnap_fasta_path: Path to barrnap fasta
    :param out_fasta_path: Path for the combined fasta file
    :param out_barstats_path: Path for the barrnap stats post de duplication
    :param min_pct_id: A filter for the pct identity, only for the non barrnap output
    :param min_length: A filter for min_length, only for the non barrnap output
    :param search_tool: The name of the search_tool that is not barrnap
    """
    mbstats = read_mbstats(mbstats_stats_path)
    mbstats = filter_mdstats(mbstats,
                             min_pct_id = min_pct_id,
                             min_length = min_length)
    mbseqs = fasta_to_df(mbstats_fasta_path, mbstats['sseqid'].values)
    mbdata = pd.merge(mbseqs, mbstats, right_on='sseqid', left_on='header',
                        how='inner')
    mbdata = process_mbdata(mbdata)
    barfasta = fasta_to_df(barrnap_fasta_path)
    barfasta = process_barfasta(barfasta)
    barfasta[['header', 'start', 'stop']].to_csv(out_barstats_path,
                                                 sep='\t', index=False)
    data = combine_fasta(mbdata, barfasta)
    df_to_fasta(data, out_fasta_path)


# TODO Decide if we want to split this.
def stage1_statistics(mbstats_path:str, barstats_path:str, output_path:str,
                      search_tool:str):
    mbstats = read_mbstats(mbstats_path)
    # 'corected_barrnap_stat.tsv'
    barstats = pd.read_csv(barstats_path, sep='\t')
    barstats = barstats_reformat(barstats, qname='16s')
    mbstats = mbstats_reformat(mbstats, search_tool, '16s')
    stats = pd.concat([mbstats, barstats])
    stats.to_csv(output_path, sep='\t', index=False)


def barstats_reformat(barstats_in:pd.DataFrame, qname:str)->pd.DataFrame:
    output_colums = {
        "header": f"{qname}_header",
        "start":  f"{qname}_start",
        "end":    f"{qname}_end",
    }
    barstats_out = barstats_in.copy()
    barstats_out.rename(columns=output_colums, inplace=True)
    barstats_out['search_tool'] = 'Barrnap'
    return barstats_out


def stage2_statistics(mbstats_path, output_path, search_tool):
    mbstats = read_mbstats(mbstats_path)
    stats = mbstats_reformat(mbstats, search_tool, 'asv')
    stats.to_csv(output_path, sep='\t', index=False)


def mbstats_reformat(mbstats_in:pd.DataFrame, search_tool:str, qname:str):
    output_colums = {
        "qseqid":   f"{qname}_header",
        "sseqid":    "bin_header",
        "pident":    "pident",
        "length":    "length",
        "mismatch":  "mismatch",
        "gapopen":   "gapopen",
        "qstart":   f"{qname}_start",
        "qend":     f"{qname}_end",
        "sstart":    "bin_start",
        "send":      "bin_end",
        "evalue":    "evalue",
        "bitscore":  "bitscore"
    }
    mbstats_out = mbstats_in[output_colums.keys()].copy()
    mbstats_out.rename(columns=output_colums, inplace=True)
    if search_tool == 'mmseqs':
        mbstats_out['search_tool'] = 'MMseqs2'
    elif search_tool == 'blast':
        mbstats_out['search_tool'] = 'BLAST'
    else:
        raise ValueError(f"The provided search tool name {search_tool}"
                          " is not recognized.")
    return mbstats_out


def load_barrnap_gtf(gtf_path):
    gtf_columns = ["seqname", "source", "feature", "start", "end", "score",
                   "strand", "frame", "attribute"]
    barstats = pd.read_csv(gtf_path, delimiter="\t", names=gtf_columns)
    return barstats


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





