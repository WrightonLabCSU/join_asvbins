"""These functions are used directly by the snakemake pipline"""
import os
import pandas as pd
from join_asvbins.utils import df_to_fasta, fasta_to_df, read_mbstats, \
    filter_mdstats, mbstats_reformat,  process_mbdata, process_barfasta, \
    barstats_reformat, combine_fasta, filter_fasta_from_headers, read_gff


def get_stage1_mbstats_fasta(mbstats, mbstats_fasta_path):
    mbseqs = fasta_to_df(mbstats_fasta_path,
                         mbstats['sseqid'].values)
    mbdata = pd.merge(mbseqs, mbstats, right_on='sseqid',
                      left_on='header',
                        how='inner')
    mbdata = process_mbdata(mbdata)
    return mbdata


def get_stage1_barrnap_fasta(barfasta, out_barstats_path):
    barfasta = process_barfasta(barfasta)
    barfasta[['header', 'start', 'stop']].to_csv(out_barstats_path,
                                             sep='\t', index=False)
    return barfasta


#TODO Add a handler for the case where blast/mmseqs or barrnap don't return matches
def combine_mbstats_barrnap(mbstats_fasta_path:str, mbstats_stats_path:str,
                            barrnap_fasta_path:str, out_fasta_path:str,
                            out_barstats_path:str, min_pct_id:float=None,
                            min_length:int=None, search_tool:str='Other',
                            allow_empty:bool=False):
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
    # TODO add checks that these functions return empty dfs if given empty
    barfasta = fasta_to_df(barrnap_fasta_path)
    mbstats = read_mbstats(mbstats_stats_path)
    mbstats = filter_mdstats(mbstats,
                             min_pct_id = min_pct_id,
                             min_length = min_length)
    if barfasta.empty and mbstats.empty:
        raise ValueError(f"There are no hits from barrnap or {search_tool},"
                          " this is most likely caused by some irregularity in"
                          " the input files. The program can not continue.")
    if mbstats.empty:
        if not allow_empty:
            raise ValueError(f"There are no hits from {search_tool}, this is"
                             " may be caused by some irregularity the input"
                             " files. Consider if your bins are few and small,"
                             " or you could use a larger, more appropriate set"
                             " of generic 16s sequences. If you are confident"
                             " in your data you can continue by passing "
                             "--allow_empty to skip this search tool and use "
                             "only barrnap. Consider using --no_clean also to"
                             " save time.")
        else:
            data = get_stage1_barrnap_fasta(barfasta, out_barstats_path)
            df_to_fasta(data, out_fasta_path)
            return
    if barfasta.empty:
        if not allow_empty:
            raise ValueError("There are no hits from barrnap, this is may"
                             " be caused by some irregularity the input files."
                             " If you are confident in your data you can"
                             " continue by passing --allow_empty to skip this"
                             f" search tool and use only {search_tool}."
                             " Consider using --no_clean also to save time.")
        else:
            data = get_stage1_mbstats_fasta(mbstats, mbstats_fasta_path)
            df_to_fasta(data, out_fasta_path)
            return
    mbdata = get_stage1_mbstats_fasta(mbstats, mbstats_fasta_path)
    barfasta = get_stage1_barrnap_fasta(barfasta, out_barstats_path)
    data = combine_fasta(mbdata, barfasta)
    df_to_fasta(data, out_fasta_path)


def stage1_statistics(mbstats_path:str, barstats_corrected_path:str, barstats_raw_path:str,
                      output_path:str, search_tool:str):
    mbstats = read_mbstats(mbstats_path)
    # 'corected_barrnap_stat.tsv'
    barstats_corrected = pd.read_csv(barstats_corrected_path, sep='\t')
    barstats_raw = read_gff(barstats_raw_path)
    barstats = barstats_reformat(barstats_corrected, barstats_raw, 'bin')
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
