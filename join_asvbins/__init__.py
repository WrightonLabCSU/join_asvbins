"""Contains main entry point to the program and helper functions"""
import os
import io
import subprocess
import argparse
from contextlib import redirect_stdout
from snakemake import snakemake

CONFIG_VALUES = {
    "bins": None,
    "output_dir": "./",
    "asv_seqs": None,
    "blast": False,
    "allow_empty": False,
    "fasta_extention": 'fa',
    "generic_16s":
       get_package_path("data/silva_clusterd_95pct_rep_seq.fasta")
}

FILTER_VALUES = {
   "s1_min_pct_id": None,
   "s2_min_pct_id": 0.97,
   "s1_min_length": 100,
   "s2_min_length": 100,
   "min_length_pct": 1,
   "max_gaps": 3,
   "max_missmatch": 2,
   "s1_mmseqs_sensitivity": 4,
   "s2_mmseqs_sensitivity": 4,
}


def get_package_path(local_path):
    abs_snake_path = os.path.join(os.path.dirname(
        os.path.abspath(__file__)),
        local_path)
    assert os.path.exists(abs_snake_path), \
        f"Unable to locate key file path; tried {abs_snake_path}"
    return abs_snake_path


def join_asvbins(bins:str=None,
                 asv_seqs:str=CONFIG_VALUES['asv_seqs'],
                 output_dir:str=CONFIG_VALUES['output_dir'],
                 blast:bool=CONFIG_VALUES['blast'],
                 generic_16s=CONFIG_VALUES['generic_16s'],
                 s1_min_pct_id=FILTER_VALUES['s1_min_pct_id'],
                 s2_min_pct_id=FILTER_VALUES['s2_min_pct_id'],
                 s1_min_length=FILTER_VALUES['s1_min_length'],
                 s2_min_length=FILTER_VALUES['s2_min_length'],
                 s1_mmseqs_sensitivity=FILTER_VALUES["s1_mmseqs_sensitivity"],
                 s2_mmseqs_sensitivity=FILTER_VALUES["s2_mmseqs_sensitivity"],
                 min_length_pct=FILTER_VALUES['min_length_pct'],
                 max_gaps=FILTER_VALUES['max_gaps'], snake_rule:str=None,
                 fasta_extention:str='fa', bin_16s_seqs:str=None,
                 max_missmatch:int=FILTER_VALUES['max_missmatch'],
                 no_clean:bool=False, no_filter:bool=False, stats:str=None,
                 snake_args:dict={}, print_dag:bool=False,
                 keep_temp=False,
                 print_rulegraph:bool=False,
                 threads=1):
    """
    This is the main entry point of the package
    """
    output_dir = os.path.abspath(output_dir)
    if bins is not None:
        bins = os.path.abspath(bins)
    if asv_seqs is not None:
        asv_seqs = os.path.abspath(asv_seqs)
    if snake_rule is None:
        if asv_seqs is not None and bin_16s_seqs is None:
            snake_rule = 'all'
        elif asv_seqs is not None:
            snake_rule = "search2_asv_bin_matches"
        elif bin_16s_seqs is None:
            snake_rule = "search1_16s_bin_finds"
    assert len(snake_rule) > 0, "There are no tasks for join_asvbins to do."
    " Check that all arguments are logical. For example, if you provided"
    " 16s from bins but not asvs then the progam has nothing to do."
    all_locals = locals()
    config = {i:all_locals.get(i)
              for i in CONFIG_VALUES
              if all_locals.get(i) is not None}
    if not no_filter:
        # TODO when python 3.9 is more popular replace this sintax
        config = dict(config, **{i:all_locals.get(i)
                                 for i in FILTER_VALUES
                                 if all_locals.get(i) is not None})
    if print_dag or print_rulegraph:
        # NOTE you need to pass this to dot -Tpdf > name.pdf sadly.
        #      Or is? we may be able to remove the graphiz dependency
        snakemake(get_package_path('Snakefile'),
                  targets=[snake_rule], workdir=output_dir,
                  config=config, forceall=True, printdag=print_dag,
                  printrulegraph=print_rulegraph, **snake_args)
        return
    if os.path.exists(output_dir) and not no_clean:
        snakemake(get_package_path('Snakefile'), workdir=output_dir,
                  config=config, delete_all_output=True, **snake_args)
    # Note that stats='stats should work'
    snakemake(get_package_path('Snakefile'), workdir=output_dir,
              config=config, cores=threads, notemp=keep_temp, **snake_args)


