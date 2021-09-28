"""Contains main entry point to the program and helper functions"""
import os
import io
import subprocess
import argparse
from contextlib import redirect_stdout
from snakemake import snakemake


def get_package_path(local_path):
    abs_snake_path = os.path.join(os.path.dirname(
        os.path.abspath(__file__)),
        local_path)
    assert os.path.exists(abs_snake_path), \
        f"Unable to locate key file path; tried {abs_snake_path}"
    return abs_snake_path


CONFIG_VALUES = {
    "bins": None,
    "output": None,
    "asv_seqs": None,
    "blast": False,
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


# DONE alwo input of step 1 serch results externaly done
# TODO clize is a beter tool for this
def join_asvbins(output:str, bins:str=None,
                 asv_seqs:str=CONFIG_VALUES['asv_seqs'],
                 working_dir:str='./tmp', blast:bool=CONFIG_VALUES['blast'],
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
                 print_rulegraph:bool=False,
                 threads=1):
    """This is the main entry point"""
    output = os.path.abspath(output)
    assert os.path.exists(os.path.dirname(output)), \
         "The output is in a directory that dose not exists," \
        f" no such location {os.path.dirname(output)}"
    working_dir = os.path.abspath(working_dir)
    output = os.path.abspath(output)
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
                  targets=[snake_rule], workdir=working_dir,
                  config=config, forceall=True, printdag=print_dag,
                  printrulegraph=print_rulegraph, **snake_args)
        return
    if os.path.exists(working_dir) and not no_clean:
        snakemake(get_package_path('Snakefile'), workdir=working_dir,
                  config=config, delete_all_output=True, **snake_args)
    # Note that stats='stats should work'
    snakemake(get_package_path('Snakefile'), workdir=working_dir,
              config=config, cores=threads, **snake_args)


class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split('=')
            getattr(namespace, self.dest)[key] = value


def main():
    parser = argparse.ArgumentParser(description="Extract 16s from bins using "
                                    "BLAST and Barrnap.")
    parser.add_argument( "--snake_rule",  type=str, default="all",
                        help="This script is snakemake under the hood. You"
                        "can run select Snakemake rules with this argument.")
    parser.add_argument("-b", "--bins",  type=str, default=CONFIG_VALUES['bins'],
                        help="The bin that you would like to match asvs to."
                        " This can be an fna file that has all the bins"
                        " combided or a directory of bins in seperate fa"
                        " files, but you must run the rename script before"
                        " you use this tool")
    # TODO Alow for QIIME format
    parser.add_argument("-w", "--working_dir", type=str, default='./tmp',
                        help="The folder where you would like the resulting"
                        " files and output to be stored")
    parser.add_argument( "--print_dag", action='store_true',
                        help="Print the directed aciclic graph used to run the"
                        " pipline. This will be in graphiz dot format. Pipe"
                        " the output of the comand to 'dot -Tpdf > name.pdf'"
                        " to visulize")
    parser.add_argument( "--print_rulegraph", action='store_true',
                        help="Print the directed aciclic graph used to run the"
                        " pipline. This will be in graphiz dot format. Pipe"
                        " the output of the comand to 'dot -Tpdf > name.pdf'"
                        " to visulize")
    parser.add_argument( "-a", "--asv_seqs",  type=str, default=None,
                        help="The asvs you would like to atach to your bins.")
    parser.add_argument("--stats",  type=str, default=None,
                        help="Passed directly to snakemake to serve as the"
                        " output location for the stats file.")
    parser.add_argument("--generic_16s",  type=str,
                        default=CONFIG_VALUES['generic_16s'],
                        help="A set of generic_16s files that may be part of"
                        " your bins.")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="The number of threads that will be used by the "
                        "program and subprocess.")
    parser.add_argument("--blast", action='store_true',
                        help="Specifies that blast should be used instead of"
                        " mmseqs. Good if your have limited memory or don't"
                        " trust MMseqs2.")

    parser.add_argument("--bin_16s_seqs", type=str, default=None,
                        help="Provide a fasta file of 16s sequences to surve"
                        " as input to the sectond search in the sequence,"
                        " the search matching bins against asv's. If this"
                        " argument is provided then the bins argument will"
                        " be ignored and the stage on fast and stats.tab"
                        " not be made. Note that your sequences must be"
                        " trimed, before you run this program.")
    parser.add_argument("--fasta_extention", type=str, default='fa',
                        help="The extention of fasta files when providing a"
                        " directory of bins. as long as your files are in "
                        " fasta format, you can give them any extention. Also,"
                        " gunziped files are suported as long as they have the"
                        " '.gz' extention. No other formats are suported and"
                        " files that don't have the target extention are"
                        " ignored.")
    parser.add_argument("--no_clean", action='store_true',
                        help="Specifies that the directory should NOT be"
                        " cleaned of results of pass runs. If your run is"
                        " interrupted this will allow you to to pickup."
                        " where you left off. Use at your own risk.")
    parser.add_argument("--no_filter", action='store_true',
                        help="Nuclier option to remove all filters. from"
                        " the analisis.")
    parser.add_argument("--snake_args", nargs='*', action=ParseKwargs,
                        default={},
                        help="Additional args for snake make, the format is"
                        " 'snakemake_arg=value'. Remove leading dashes and"
                        " seperate entries with spaces.")
    parser.add_argument("-o", "--output", type=str, default=1,
                        help="Instead of dumping the output into the output "
                        "folder, the program should make a new folder with "
                        "this name inside the output location. This is "
                        "especially useful for testing.")
    parser.add_argument("--s1_mmseqs_sensitivity", type=int,
                        default=FILTER_VALUES['s1_mmseqs_sensitivity'],
                        help="This value is pased to mmseqs during the stage 1"
                        " search as the sensitivity setting. Increasing it"
                        " could result in more hits at the cost of time")
    parser.add_argument("--s2_mmseqs_sensitivity", type=int,
                        default=FILTER_VALUES['s2_mmseqs_sensitivity'],
                        help="This value is pased to mmseqs during the stage 2"
                        " search as the sensitivity setting. Increasing it"
                        " could result in more hits at the cost of time")
    parser.add_argument("--s1_min_pct_id", type=float,
                        default=FILTER_VALUES['s1_min_pct_id'],
                        help="This limits the percent identity that will be"
                        " tolerated in matching generic 16s scaffold to bins")
    parser.add_argument("--s2_min_pct_id", type=float,
                        default=FILTER_VALUES['s2_min_pct_id'],
                        help="This limits the percent identity that can be"
                        " used for matching ASVs to extracted 16s matches")
    #NOTE to self: We use BLAST to identify the regions in the input genome
    #that match the alleles in the ResFinder database. BLAST finds the best
    #local alignment (overlap) between a sequence in the input genome and an
    #allele in the ResFinder database. The Alignment Length is the length of
    #the alignment measured in basepairs. For perfect matches the Alignment
    #Length equals the DB allele Length.
    parser.add_argument("--s1_min_length", type=int,
                        default=FILTER_VALUES['s1_min_length'],
                        help="This limits how short the length of the match"
                        " from the ASVs to bin 16s can be for consideration.")
    parser.add_argument("--s2_min_length", type=int,
                        default=FILTER_VALUES['s2_min_length'],
                        help="This limits how short the length of the match"
                        " from the ASVs to bin 16s can be for consideration.")
    parser.add_argument("--min_length_pct", type=float,
                        default=FILTER_VALUES['min_length_pct'],
                        help="This limits how short the length of the match"
                        " from the ASVs to bin 16s can be for consideration.")
    parser.add_argument("--max_gaps", type=int,
                        default=FILTER_VALUES['max_gaps'],
                        help="This limits the number of gaps that will be"
                        " tolerated in matches from the ASVs to bin 16S.")
    parser.add_argument("--max_missmatch", type=int,
                        default=FILTER_VALUES['max_missmatch'],
                        help="This limits the number of mismatches that will"
                        " be tolerated in matches from the ASVs to bin 16S.")
    parser.set_defaults(func=join_asvbins)
    args = parser.parse_args()
    args_dict = {i: j for i, j in vars(args).items() if i != 'func'}
    args.func(**args_dict)

