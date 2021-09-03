import os
import subprocess
import argparse
from asv_to_bins.extract_16s import combine_blast_barrnap

def get_package_path(local_path):
    abs_snake_path = os.path.join(os.path.dirname(
        os.path.abspath(__file__)),
        local_path)
    assert os.path.exists(abs_snake_path), \
        f"Unable to locate the Snakemake workflow file; tried {abs_snake_path}"
    return abs_snake_path

NAME_FOR_COMBINED_BINS = "combined_bins.fna"
FILTER_VALUES = {
   "asv_seqs": get_package_path("data/silva_clusterd_95pct_all_seqs.fasta"),
   "s1_min_pct_id": None,
   "s2_min_pct_id": 0.97,
   "s1_min_length": 100,
   "s2_min_length": 100,
   "min_length_pct": 1,
   "max_gaps": 1,
   "max_missmatch": 2,
   "s1_mmseqs_sensitivity": 4,
   "s2_mmseqs_sensitivity": 4,
}

def parse_args():
    parser = argparse.ArgumentParser(description="Extract 16s from bins using "
                                    "BLAST and Barrnap.")
    parser.add_argument( "--snake_rule",  type=str, default="all",
                        help="This script is snakemake under the hood. You"
                        "can run select Snakemake rules with this argument.")
    parser.add_argument("-b", "--bins",  type=str, default=None,
                        help="The bin that you would like to match asvs to."
                        " This can be an fna file that has all the bins"
                        " combided or a directory of bins in seperate fa files,"
                        " but you must run the rename script before you use"
                        " this tool")
    # TODO Alow for QIIME format
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="The folder where you would like the resulting"
                        " files and output to be stored")
    # TODO add these features:
    #
    # Exact matches
    # No gaps
    # 0, 1, 2 mismatches
    # 100 % length
    parser.add_argument( "--print_dag",  type=str, default=None,
                        help="Print the DAG that is being used to run this"
                        " pipline Note that your choice of args will result"
                        " in a different DAG.")
    parser.add_argument( "-a", "--asv_seqs",  type=str, default=None,
                        help="The asvs you would like to atach to your bins.")
    parser.add_argument("--generic_16s",  type=str,
                        default=FILTER_VALUES['asv_seqs'],
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
    parser.add_argument("--snake_args", type=str, default="",
                        help="Additional args for snake make, put them in quots"
                        " like you are providing them to snake make")
    parser.add_argument("-n", "--name", type=str, default=1,
                        help="Instead of dumping the output into the output "
                        "folder, the program should make a new folder with "
                        "this name inside the output location. This is "
                        "especially useful for testing.")
    #TODO
    parser.add_argument("--s1_mmseqs_sensitivity", type=int,
                        default=FILTER_VALUES['s1_mmseqs_sensitivity'],
                        help="This value is pased to mmseqs during the stage 1"
                        " search as the sensitivity setting. Increasing it"
                        " could result in more hits at the cost of time")
    #TODO
    parser.add_argument("--s2_mmseqs_sensitivity", type=int,
                        default=FILTER_VALUES['s2_mmseqs_sensitivity'],
                        help="This value is pased to mmseqs during the stage 2"
                        " search as the sensitivity setting. Increasing it"
                        " could result in more hits at the cost of time")
    #TODO
    parser.add_argument("--s1_min_pct_id", type=float,
                        default=FILTER_VALUES['s1_min_pct_id'],
                        help="This limits the percent identity that will be"
                        " tolerated in matching generic 16s scaffold to bins")
    #TODO
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
    #TODO
    parser.add_argument("--s1_min_length", type=int,
                        default=FILTER_VALUES['s1_min_length'],
                        help="This limits how short the length of the match"
                        " from the ASVs to bin 16s can be for consideration.")
    #TODO
    parser.add_argument("--s2_min_length", type=int,
                        default=FILTER_VALUES['s2_min_length'],
                        help="This limits how short the length of the match"
                        " from the ASVs to bin 16s can be for consideration.")
    #TODO
    parser.add_argument("--min_length_pct", type=float,
                        default=FILTER_VALUES['min_length_pct'],
                        help="This limits how short the length of the match"
                        " from the ASVs to bin 16s can be for consideration.")
    #TODO
    parser.add_argument("--max_gaps", type=int,
                        default=FILTER_VALUES['max_gaps'],
                        help="This limits the number of gaps that will be"
                        " tolerated in matches from the ASVs to bin 16S.")
    #TODO
    parser.add_argument("--max_missmatch", type=int,
                        default=FILTER_VALUES['max_missmatch'],
                        help="This limits the number of mismatches that will"
                        " be tolerated in matches from the ASVs to bin 16S.")
    args = parser.parse_args()
    return args



def snakemake_run(bins_folder:str, all_bin_seqs_path:str,
                  asv_seqs_path:str,
                  working_dir, tool='mmseqs',
                  generic_16s_data=FILTER_VALUES['asv_seqs'],
                  s1_min_pct_id=FILTER_VALUES['s1_min_pct_id'],
                  s2_min_pct_id=FILTER_VALUES['s2_min_pct_id'],
                  s1_min_length=FILTER_VALUES['s1_min_length'],
                  s2_min_length=FILTER_VALUES['s2_min_length'],
                  s1_mmseqs_sensitivity=FILTER_VALUES["s1_mmseqs_sensitivity"],
                  s2_mmseqs_sensitivity=FILTER_VALUES["s2_mmseqs_sensitivity"],
                  min_length_pct=FILTER_VALUES['min_length_pct'],
                  max_gaps=FILTER_VALUES['max_gaps'],
                  run_rule:str = "all",
                  max_missmatch=FILTER_VALUES['max_missmatch'],
                  clean=True, snake_args="", print_dag:str=None
                  threads=1):
    key_dag_args = ( # These arguments may affect the DAG so they are separate
        f" --snakefile {get_package_path('Snakefile')}"
        f" --directory {working_dir}"
        " --config"
            f" bins_folder=\'{bins_folder}\'"
            f" all_bin_seqs_path=\'{all_bin_seqs_path}\'"
            f" asv_seqs_path=\'{asv_seqs_path}\'"
            f" generic_16s_data=\'{generic_16s_data}\'"
            f" tool=\'{tool}\'"
            f" s1_min_pct_id=s1_min_pct_id"
            f" s2_min_pct_id=s2_min_pct_id"
            f" s1_min_length=s1_min_length"
            f" s2_min_length=s2_min_length"
            f" s1_mmseqs_sensitivity=s1_mmseqs_sensitivity"
            f" s2_mmseqs_sensitivity=s2_mmseqs_sensitivity"
            f" min_length_pct=min_length_pct"
            f" max_gaps=max_gaps"
            f" {snake_args}"
    )
    if print_dag is not None:
        subprocess.run(f"snakemake {key_dag_args} --forceall --dag | dot -Tpdf > {print_dag}.pdf",
                       check=True, shell=True)
        return
    if os.path.exists(working_dir) and clean:
        subprocess.run(f"snakemake --delete-all-output {key_dag_args}",
                       check=True, shell=True)

    subprocess.run(f"snakemake {run_rule} --cores {threads} {key_dag_args}",
                   check=True, shell=True)

def main():
    args = parse_args()
    working_dir = os.path.abspath(args.output)
    print(working_dir)
    if args.name is not None:
        working_dir = os.path.join(working_dir, args.name)
    if os.path.isdir(args.bins): # this is a dir of fa files
        bins_folder = os.path.abspath(args.bins)
        all_bin_seqs_path = NAME_FOR_COMBINED_BINS
    else:
        bins_folder = None
        all_bin_seqs_path = os.path.abspath(args.bins)

    if args.asv_seqs is not None:
        asv_seqs_path = os.path.abspath(args.asv_seqs)
    else:
        asv_seqs_path = None
    if args.snake_rule is not None:
        run_rule = args.snake_rule
    else:
        if asv_seqs_path is not None:
            run_rule = "all"
        else:
            run_rule = "no_asvs"
    tool = 'blast' if args.blast else 'mmseqs'
    clean = not args.no_clean
    snakemake_run(bins_folder, all_bin_seqs_path, asv_seqs_path, working_dir, tool,
                  generic_16s_data=args.generic_16s_data,
                  s1_min_pct_id=args.s1_min_pct_id,
                  s2_min_pct_id=args.s2_min_pct_id,
                  s1_min_length=args.s1_min_length,
                  s2_min_length=args.s2_min_length,
                  s1_mmseqs_sensitivity=args.s1_mmseqs_sensitivity,
                  s2_mmseqs_sensitivity=args.s2_mmseqs_sensitivity,
                  min_length_pct=args.min_length_pct,
                  snake_args=args.snake_args, run_rule=run_rule, clean=clean,
                  threads=args.threads, print_dag=args.print_dag)


if __name__ == '__main__':
    main()
