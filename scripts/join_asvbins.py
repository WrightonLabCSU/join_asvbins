"""This is the main entry point of the program"""
import argparse

from asv_to_bins import join_asvbins, CONFIG_VALUES, FILTER_VALUES


class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split('=')
            getattr(namespace, self.dest)[key] = value


if __name__ == '__main__':
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
    # args = parser.parse_args()
    parser.set_defaults(func=join_asvbins)
    args = parser.parse_args()
    args_dict = {i: j for i, j in vars(args).items() if i != 'func'}
    args.func(**args_dict)

