import os
import subprocess
import argparse
from asv_to_bins.extract_16s import combine_blast_barrnap

DEFAULT_NAME_FOR_COMBINED_BINS = "combined_bins"

def parse_args():
    parser = argparse.ArgumentParser(description="Extract 16s from bins using "
                                    "BLAST and Barrnap.")
    parser.add_argument( "-r", "--rule",  type=str, default=None,
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
                        help="The folder where you would like the resulting "
                        "files and output to be stored")
    # TODO add these features:
    #
    # Exact matches
    # No gaps
    # 0, 1, 2 mismatches
    # 100 % length

    parser.add_argument( "-a", "--asvs",  type=str, default=None,
                        help="The asvs you would like to atach to your bins.")
    # TODO If necessary make this optional so it can connect to other command
    #      line tools or pipe designated output to other programs.
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="The number of threads that will be used by the "
                        "program and subprocess.")
    parser.add_argument("-n", "--name", type=str, default=1,
                        help="Instead of dumping the output into the output "
                        "folder, the program should make a new folder with "
                        "this name inside the output location. This is "
                        "especially useful for testing.")
    args = parser.parse_args()
    return args


def get_snake_path():
    abs_snake_path = os.path.join(os.path.dirname(
        os.path.abspath(__file__)),
        "Snakefile")
    assert os.path.exists(abs_snake_path), "Unable to locate the Snakemake workflow file; tried %s" % (abs_snake_path)
    return abs_snake_path

def snakemake_run(run_rule, bins_folder, bins_all_seqs_path, asv_seqs_path, working_dir, threads=1):
    cmd = ("snakemake {run_rule}"
           " --snakefile {snakefile}"
           " --directory {working_dir}"
           " --config"
               " bins_folder=\'{bins_folder}\'"
               " bins_all_seqs_path=\'{bins_all_seqs_path}\'"
               " asv_seqs_path=\'{asv_seqs_path}\'"
           " --cores {threads}").format(
               run_rule=run_rule,
               snakefile=get_snake_path(),
               working_dir=working_dir,
               bins_folder=bins_folder,
               bins_all_seqs_path=bins_all_seqs_path,
               asv_seqs_path=asv_seqs_path,
               text="My text",
               threads=threads,
           )
    subprocess.run(cmd, check=True, shell=True)


def main():
    args = parse_args()
    working_dir = os.path.abspath(args.output)
    print(working_dir)
    if args.name is not None:
        working_dir = os.path.join(working_dir, args.name)
    if os.path.isdir(args.bins): # this is a dir of fa files
        bins_folder = os.path.abspath(args.bins)
        bins_all_seqs_path = DEFAULT_NAME_FOR_COMBINED_BINS + ".fna"
    else:
        bins_folder = None
        bins_all_seqs_path = os.path.abspath(args.bins)

    if args.asvs is not None:
       asv_seqs_path = os.path.abspath(args.asvs)
    else:
        asv_seqs_path = None
    if args.rule is not None:
        run_rule = args.rule
    else:
        if asv_seqs_path is not None:
            run_rule = "all"
        else:
            run_rule = "no_asvs"

    snakemake_run(run_rule, bins_folder, bins_all_seqs_path, asv_seqs_path, working_dir, threads=args.threads)

if __name__ == '__main__':
    main()
