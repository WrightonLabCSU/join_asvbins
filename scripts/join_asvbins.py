"""This is the main entry point of the program"""
import os
import subprocess
import argparse
from asv_to_bins.extract_16s import combine_mbstats_barrnap

def get_package_path(local_path):
    abs_snake_path = os.path.join(os.path.dirname(
        os.path.abspath(__file__)),
        local_path)
    assert os.path.exists(abs_snake_path), \
        f"Unable to locate key file path; tried {abs_snake_path}"
    return abs_snake_path

MAIN_VALUES = {
    "bins": None,
    "asv_seqs": None,
    "blast": False,
    "generic_16s":
       get_package_path("data/silva_clusterd_95pct_all_seqs.fasta")
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


# NOTE that the working_dir is acting as the output which is undified for now
# TODO Three outputs fasta with 16 from bins, fata with asv matches to bins, and 2 tab files describing them
# NOTE rename this
# TODO alwo input of step 1 serch results externaly done
def snakemake_run(bins:str, asv_seqs:str=MAIN_VALUES['asv_seqs'],
                  working_dir:str='./tmp', blast:bool=MAIN_VALUES['blast'],
                  output:str="", generic_16s=MAIN_VALUES['generic_16s'],
                  s1_min_pct_id=FILTER_VALUES['s1_min_pct_id'],
                  s2_min_pct_id=FILTER_VALUES['s2_min_pct_id'],
                  s1_min_length=FILTER_VALUES['s1_min_length'],
                  s2_min_length=FILTER_VALUES['s2_min_length'],
                  s1_mmseqs_sensitivity=FILTER_VALUES["s1_mmseqs_sensitivity"],
                  s2_mmseqs_sensitivity=FILTER_VALUES["s2_mmseqs_sensitivity"],
                  min_length_pct=FILTER_VALUES['min_length_pct'],
                  max_gaps=FILTER_VALUES['max_gaps'],
                  snake_rule:str = None,
                  max_missmatch:int=FILTER_VALUES['max_missmatch'],
                  no_clean:bool=False, no_filter:bool=False, snake_args:str="",
                  print_dag:str=None,
                  threads=1):

    working_dir = os.path.abspath(working_dir)
    bins = os.path.abspath(bins)
    if asv_seqs is not None:
        asv_seqs = os.path.abspath(asv_seqs)
    if snake_rule is None:
        if asv_seqs is not None:
            snake_rule = "all"
        else:
            snake_rule = "no_asvs"
    all_locals = locals()
    config_values = MAIN_VALUES if no_filter \
                    else dict(MAIN_VALUES, **FILTER_VALUES)
    snake_config = " ".join([f"k={all_locals[k]}"
                             for k in config_values
                             if all_locals[k] is not None])
    key_dag_args = ( # These arguments may affect the DAG so they are separate
        f" --snakefile {get_package_path('Snakefile')}"
        f" --directory {working_dir}"
        f" --config {snake_config}"
        f" {snake_args}"
    )
    if print_dag is not None:
        subprocess.run(f"snakemake {key_dag_args} --forceall --dag |"
                       " dot -Tpdf > {print_dag}.pdf",
                       check=True, shell=True)
        return
    #NOTE The folders are cleaned by default
    if os.path.exists(working_dir) and not no_clean:
        subprocess.run(f"snakemake --delete-all-output --cores {threads}"
                       f" {key_dag_args}",
                       check=True, shell=True)

    subprocess.run(f"snakemake {snake_rule} --cores {threads} {key_dag_args}",
                   check=True, shell=True)


def parse_args():
    parser = argparse.ArgumentParser(description="Extract 16s from bins using "
                                    "BLAST and Barrnap.")
    parser.add_argument( "--snake_rule",  type=str, default="all",
                        help="This script is snakemake under the hood. You"
                        "can run select Snakemake rules with this argument.")
    parser.add_argument("-b", "--bins",  type=str, default=MAIN_VALUES['bins'],
                        help="The bin that you would like to match asvs to."
                        " This can be an fna file that has all the bins"
                        " combided or a directory of bins in seperate fa"
                        " files, but you must run the rename script before"
                        " you use this tool")
    # TODO Alow for QIIME format
    parser.add_argument("-w", "--working_dir", type=str, default='./tmp',
                        help="The folder where you would like the resulting"
                        " files and output to be stored")
    parser.add_argument( "--print_dag",  type=str, default=None,
                        help="Print the DAG that is being used to run this"
                        " pipline Note that your choice of args will result"
                        " in a different DAG.")
    parser.add_argument( "-a", "--asv_seqs",  type=str, default=None,
                        help="The asvs you would like to atach to your bins.")
    parser.add_argument("--generic_16s",  type=str,
                        default=MAIN_VALUES['generic_16s'],
                        help="A set of generic_16s files that may be part of"
                        " your bins.")
    # TODO If necessary make this optional so it can connect to other command
    #      line tools or pipe designated output to other programs.
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="The number of threads that will be used by the "
                        "program and subprocess.")
    parser.add_argument("--blast", action='store_true',
                        help="Specifies that blast should be used instead of"
                        " mmseqs. Good if your have limited memory or don't"
                        " trust MMseqs2.")
    parser.add_argument("--no_clean", action='store_true',
                        help="Specifies that the directory should NOT be"
                        " cleaned of results of pass runs. If your run is"
                        " interrupted this will allow you to to pickup."
                        " where you left off. Use at your own risk.")
    parser.add_argument("--no_filter", action='store_true',
                        help="Nuclier option to remove all filters. from"
                        " the analisis.")
    parser.add_argument("--snake_args", type=str, default="",
                        help="Additional args for snake make, put them in"
                        " quots"
                        " like you are providing them to snake make")
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
    # args = parser.parse_args()
    parser.set_defaults(func=snakemake_run)
    args = parser.parse_args()
    args_dict = {i: j for i, j in vars(args).items() if i != 'func'}
    args.func(**args_dict)


if __name__ == '__main__':
    parse_args()
