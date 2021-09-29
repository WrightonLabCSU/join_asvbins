"""These functions are used directly by the snakemake pipline"""
import os
import pandas as pd
from join_asvbins.utils import df_to_fasta, fasta_to_df, read_mbstats, \
    filter_mdstats, mbstats_reformat,  process_mbdata, process_barfasta, \
    barstats_reformat, combine_fasta, filter_fasta_from_headers


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


def stage1_statistics(mbstats_path:str, barstats_path:str, output_path:str,
                      search_tool:str):
    mbstats = read_mbstats(mbstats_path)
    # 'corected_barrnap_stat.tsv'
    barstats = pd.read_csv(barstats_path, sep='\t')
    barstats = barstats_reformat(barstats, qname='16s')
    mbstats = mbstats_reformat(mbstats, search_tool, '16s')
    stats = pd.concat([mbstats, barstats])
    stats.to_csv(output_path, sep='\t', index=False)


def stage2_statistics(mbstats_path, output_path, search_tool):
    mbstats = read_mbstats(mbstats_path)
    stats = mbstats_reformat(mbstats, search_tool, 'asv')
    stats.to_csv(output_path, sep='\t', index=False)


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


def pullseqs_header_name_from_tab(in_fasta_path:str, out_fasta_path:str,
                                  tab_file_path:str,
                                  header_column:str='sseqid'):
       mbstats = read_mbstats(tab_file_path)
       headers = set(mbstats[header_column].values)
       filter_fasta_from_headers(in_fasta_path, out_fasta_path, headers)
