"""These functions are used directly by the snakemake pipline"""
import os
import pandas as pd
from join_asvbins.utils import df_to_fasta, fasta_to_df, read_mbstats, \
    filter_mdstats, mbstats_reformat,  process_mbdata, process_barfasta, \
    barstats_reformat, combine_fasta, filter_fasta_from_headers, read_gff, \
    get_stage1_mbstats_fasta


def combine_mbstats_barrnap(mbstats_fasta_path:str, mbstats_stats_path:str,
                            barrnap_fasta_path:str, out_fasta_path:str,
                            out_stats_path:str, barrnap_stats_path:str,
                            search_tool:str='Other', allow_empty:bool=False,
                            **filter_kargs) -> None:
    """
    :returns:

    :param out_stats_path:
    :param allow_empty:
    """
    """
    Combine the statistics from mmseqs or blast with the output from barrnap.

    :param mbstats_fasta_path: Path to mmseqs or blast fasta
    :param mbstats_stats_path: Path to mmseqs or blast statistics, see docs for format requirements
    :param barrnap_fasta_path: Path to barrnap fasta
    :param out_fasta_path: Path for the combined fasta file
    :param out_stats_path: Path for the combined statistics file
    :param barrnap_stats_path: Path for the barrnap stats gff
    :param min_pct_id: A filter for the pct identity, only for the non barrnap output
    :param min_length: A filter for min_length, only for the non barrnap output
    :param search_tool: The name of the search_tool that is not barrnap
    :param allow_empty: If true the program will continue if only one search_tool gives results
    :raises ValueError:
    """
    # TODO add checks that these functions return empty dfs if given empty
    barfasta = fasta_to_df(barrnap_fasta_path)
    mbstats = read_mbstats(mbstats_stats_path)
    mbstats = filter_mdstats(mbstats, **filter_kargs)
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
            data = process_barfasta(barfasta)
            # save_barnap_stats(barfasta, out_barstats_path)
            df_to_fasta(data, out_fasta_path)
            make_stage1_statistics(out_stats_path, search_tool, barfasta=data,
                                   barrnap_stats_path=barrnap_stats_path)
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
            make_stage1_statistics(out_stats_path, search_tool, mbstats=data)
            return
    mbdata = get_stage1_mbstats_fasta(mbstats, mbstats_fasta_path)
    barfasta = process_barfasta(barfasta)
    # save_barnap_stats(barfasta, out_barstats_path)
    data = combine_fasta(mbdata, barfasta)
    df_to_fasta(data, out_fasta_path)
    make_stage1_statistics(out_stats_path, search_tool, mbstats=mbstats,
                           barfasta=barfasta, barrnap_stats_path=barrnap_stats_path)


def make_stage1_statistics(output_path:str, search_tool:str,
                           mbstats:pd.DataFrame=None,
                           barfasta:pd.DataFrame=None, barrnap_stats_path:str=None):
    if barfasta is not None:
        barstats = read_gff(barrnap_stats_path)
        barstats_corrected = barfasta[['header', 'start', 'stop']]
        barstats = barstats_reformat(barstats_corrected, barstats, 'bin')
        if mbstats is None:
            barstats.to_csv(output_path, sep='\t', index=False)
            return
    if mbstats is not None:
        mbstats = mbstats_reformat(mbstats, search_tool, '16s')
        if barstats is None:
            mbstats.to_csv(output_path, sep='\t', index=False)
            return
    stats = pd.concat([mbstats, barstats])
    stats.to_csv(output_path, sep='\t', index=False)


# def stage1_statistics(mbstats_path:str, barstats_corrected_path:str, barstats_raw_path:str,
#                       output_path:str, search_tool:str):
#     mbstats = read_mbstats(mbstats_path)
#     # 'corected_barrnap_stat.tsv'
#     barstats_corrected = pd.read_csv(barstats_corrected_path, sep='\t')
#     barstats_raw = read_gff(barstats_raw_path)
#     barstats = barstats_reformat(barstats_corrected, barstats_raw, 'bin')
#     mbstats = mbstats_reformat(mbstats, search_tool, '16s')
#     stats = pd.concat([mbstats, barstats])
#     stats.to_csv(output_path, sep='\t', index=False)


def stage2_statistics(mbstats_path, output_path, search_tool):
    mbstats = read_mbstats(mbstats_path)
    stats = mbstats_reformat(mbstats, search_tool, 'asv')
    stats.to_csv(output_path, sep='\t', index=False)


def filter_from_mbstats(stats_file:str, fasta_file_in:str, fasta_file_out:str,
                        **filter_kargs):
                        # min_pct_id:float=None, min_length:int=None,
                        # min_len_pct:float=None, max_gaps:int=None,
                        # max_missmatch:int=None):
        mbstats = read_mbstats(stats_file)
        mbstats = filter_mdstats(mbstats, **filter_kargs)
        filter_fasta_from_headers(fasta_file_in,
                                  fasta_file_out, mbstats['sseqid'].values)


def pullseqs_header_name_from_tab(in_fasta_path:str, out_fasta_path:str,
                                  tab_file_path:str,
                                  header_column:str='sseqid'):
       mbstats = read_mbstats(tab_file_path)
       headers = set(mbstats[header_column].values)
       filter_fasta_from_headers(in_fasta_path, out_fasta_path, headers)
