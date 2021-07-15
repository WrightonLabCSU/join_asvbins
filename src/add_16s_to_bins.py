from subprocess import Popen, PIPE
import argparse
import os
from glob import glob
import pytest

# TODO
def parse_args():
    parser = argparse.ArgumentParser(description="Extract 16s from bins using "
                                    "BLAST and Barrnap.")
    parser.add_argument("input", type=str, default=None,
                        help="A name for the run")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="The number of threads that will be used by the "
                        "program and subprocess.")
    args = parser.parse_args()
    return args


def barrnap_extact_16s(threads:int, input_path:str, input_type='fna'):
    if input_type == 'fna':
        Popen(["barrnap", "--threads", str(threads), input_path],
              stdout=PIPE)
    elif input_type == 'fa':
        cat = Popen(["cat"] + glob(os.path.join(input_path, "*.fa")),
              stdout=PIPE)
        # you would hope that if barrnap can take standerd inpute it could output
        barrnap =Popen(["barrnap", "--threads", str(threads),
               ], stdout=PIPE, stdin=cat.stdout)
        barrnap = Popen(["barrnap", "--threads", str(threads),
               ], stdout=PIPE, stdin=cat.stdout)
        cat.wait()
        barrnap.wait()
        return barrnap.communicate()
    else:
        raise ValueError("The file type can only be 'fna' or 'fa'.")

barrnap_extact_16s(10, '../data/salmonella/', 'fa')
os.system('cat ../data/salmonella/*.fa')
input_path = "../data/salmonella/"
def validate_input(input_path):
    if os.path.isdir(input_path): # this is a dir of fa files
        # TODO check for fa files
        return 'fa'
    # else this is a combined fna file
    # TODO check not empty
    return 'fna'


def main():
    args = parse_args()
    input_type = validate_input(args.input_path)
    barrnap_extact_16s(args.threads, args.input, input_type)

if __name__ == '__main__':
    main()

def test_vaid_input_type():
    assert 'fa' == validate_input('../data/salmonella/'), "Problem with file or folder detection"


#  barrnap --threads 10 All_bins.fna > rrna.gff
#  #gives you rrna.gff
#  grep '16S' rrna.gff > 16S-gff.gff
#  bedtools getfasta -fi All_bins.fna -bed 16S-gff.gff -fo 16S-fasta.fna
#  grep ">" 16S-fasta.fna|sed 's/>//g' > 16S-id.txt
#  xargs samtools faidx 16S-fasta.fna < 16S-id.txt > barrnap_fastafile-16S.fna
#
#
#  Ran blast on all bins
#
#  #pull a full length 16S fasta from silva-  e. coli.
#  makeblastdb -dbtype nucl -in ../All_bins.fna -out All_bins.fna.fa_BLASTDB
#  blastn -db All_bins.fna.fa_BLASTDB -out blast_fastafile-16S.txt -query ecoli_16S.txt -outfmt 6
#
#  #create list of fasta headers
#  awk '{print $2}' blast_fastafile-16S.txt > blast_fastaheaders-16S.txt
#
#  #use pullseq to pull the scaffolds with 16S hits
#  pullseq_header_name.py -i ../All_bins.fna -o blast_fastafile-16S.fna -e F -n blast_fastaheaders-16S.txt
