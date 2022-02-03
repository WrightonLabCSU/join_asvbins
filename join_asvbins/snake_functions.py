"""These functions are used directly by the snakemake pipline"""
import os
import pandas as pd
from join_asvbins.utils import df_to_fasta, fasta_to_df, read_mbstats, \
    filter_mdstats, mbstats_reformat,  process_mbdata, process_barfasta, \
    barstats_reformat, combine_fasta, filter_fasta_from_headers, read_gff, \
    get_stage1_mbstats_fasta


def resolve_dup_gene_locs(mbstats:str, bs_name:str, bs_start:str,
                          bs_end:str, values:str, ascending:bool):
    gene_locs = {}
    def check_overlap(line:pd.Series):
        line_name = line[bs_name]
        line_start = min(line[bs_start], line[bs_end])
        line_end = max(line[bs_start], line[bs_end])
        if line_name in gene_locs:
            match_start = gene_locs[line_name]['start']
            match_end = gene_locs[line_name]['end']
            if match_start <= line_start < match_end:
                return False
            if match_start <= line_end < match_end:
                return False
        gene_locs[line_name] = {
            'start': line_start,
            'end': line_end,
        }
        return True
    mbstats.sort_values(values, inplace=True, ascending=False)
    mbstats = mbstats[mbstats.apply(check_overlap, axis=1)]
    return mbstats

def combine_mbstats_barrnap(mbstats_fasta_path:str, mbstats_stats_path:str,
                            barrnap_fasta_path:str, out_fasta_path:str,
                            out_stats_path:str, barrnap_stats_path:str,
                            search_tool:str, allow_empty:bool=False,
                            **filter_kargs) -> None:
    """
    Combine the statistics from mmseqs or blast with  barrnap.

    This function is a monster

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
    mbstats_raw = read_mbstats(mbstats_stats_path)
    mbstats_raw.drop_duplicates(inplace=True)
    mbstats = filter_mdstats(mbstats_raw, **filter_kargs)
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
                             " of generic 16S sequences. If you are confident"
                             " in your data you can continue by passing"
                             " --allow_empty to skip this search tool and use "
                             " only barrnap. Consider using --no_clean also to"
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
    data = combine_fasta(mbdata, barfasta, search_tool)
    df_to_fasta(data, out_fasta_path)
    make_stage1_statistics(out_stats_path, search_tool, mbstats=mbstats,
                           barfasta=barfasta,
                           barrnap_stats_path=barrnap_stats_path)


def make_stage1_statistics(output_path:str, search_tool:str,
                           mbstats:pd.DataFrame=None,
                           barfasta:pd.DataFrame=None,
                           barrnap_stats_path:str=None):
    if barfasta is not None:
        barstats = read_gff(barrnap_stats_path)
        barstats_corrected = barfasta[['header', 'start', 'stop']]
        barstats = barstats_reformat(barstats_corrected, barstats, 'bin')
        barstats = resolve_dup_gene_locs(barstats, bs_name='name',
                                         bs_start='bin_scaffold_start',
                                         bs_end='bin_scaffold_end',
                                         values='barrnap_e-value',
                                         ascending=True)
        if mbstats is None:
            barstats.to_csv(output_path, sep='\t', index=False, na_rep='NA')
            return
    if mbstats is not None:
        mbstats = resolve_dup_gene_locs(mbstats,
                                        bs_name='sseqid',
                                        bs_start='send', bs_end='sstart',
                                        values='bitscore', ascending=False)
        mbstats = mbstats_reformat(mbstats, search_tool, '16S')
        if barstats is None:
            mbstats.to_csv(output_path, sep='\t', index=False, na_rep='NA')
            return
    stats = pd.concat([mbstats, barstats])
    stats.to_csv(output_path, sep='\t', index=False, na_rep='NA')


def filter_from_mbstats(stats_file_in:str, fasta_file_in:str,
                        fasta_file_out:str,
                        stats_file_out:str, search_tool:str,  **filter_kargs):
        mbstats = read_mbstats(stats_file_in)
        mbstats = filter_mdstats(mbstats, **filter_kargs)
        mbstats = mbstats_reformat(mbstats, search_tool, 'ASV')
        mbstats.to_csv(stats_file_out, sep='\t', index=False, na_rep='NA')
        filter_fasta_from_headers(fasta_file_in,
                                  fasta_file_out,
                                  mbstats['bin_scaffold_header'].values)


def pullseqs_header_name_from_tab(in_fasta_path:str, out_fasta_path:str,
                                  tab_file_path:str,
                                  header_column:str='sseqid'):
       mbstats = read_mbstats(tab_file_path)
       headers = set(mbstats[header_column].values)
       filter_fasta_from_headers(in_fasta_path, out_fasta_path, headers)
