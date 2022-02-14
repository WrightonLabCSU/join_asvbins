"""Here we can test the steps of the pipeline in a methodically"""
import pytest

import subprocess
from join_asvbins import join_asvbins
import os


MINI_BINS = os.path.join('tests', 'data', 'mini_bins')
MINI_ASV = os.path.join('tests', 'data', 'mini_salmonella_asv.fa')
MINI_16S = os.path.join('tests', 'data', 'mini_silva.fa')
EXPECTED_CANDIDATE_FNA = os.path.join('tests', 'data', 'expected_output', 'candidate_sequences.fna')
MINI_BINS_BARRNAP = os.path.join(MINI_BINS, 'mini_salmonella_barrnap.fa')
MINI_BINS_MMSEQS  = os.path.join(MINI_BINS, 'mini_salmonella_mmseqsblast.fa')



# allow_empty
# keep_temp
# no_clean
# no_filter
# keep_temp
def test_run_blast(tmp_path):
    output_path = os.path.join(tmp_path, 'test_run_blast')
    join_asvbins(
        bins = MINI_BINS,
        asv_seqs = MINI_ASV,
        output_dir = output_path,
        generic_16S = MINI_16S,
        threads=4,
        verbosity=6,
        blast=True
    )
    assert os.path.exists(os.path.join(output_path, "candidate_sequences.fna"))
    assert os.path.exists(os.path.join(output_path, "candidate_statistics.tsv"))
    assert os.path.exists(os.path.join(output_path, "match_sequences.fna"))
    assert os.path.exists(os.path.join(output_path, "match_statistics.tsv"))
    assert len(os.listdir(output_path)) == 5 # acounts for .snakemake folder


def test_run_mmseqs(tmp_path):
    output_path = os.path.join(tmp_path, 'test_run_mmseqs')
    join_asvbins(
        bins = MINI_BINS,
        asv_seqs = MINI_ASV,
        output_dir = output_path,
        generic_16S = MINI_16S,
        threads=4,
        verbosity=6
    )
    assert os.path.exists(os.path.join(output_path, "candidate_sequences.fna"))
    assert os.path.exists(os.path.join(output_path, "candidate_statistics.tsv"))
    assert os.path.exists(os.path.join(output_path, "match_sequences.fna"))
    assert os.path.exists(os.path.join(output_path, "match_statistics.tsv"))
    assert len(os.listdir(output_path)) == 5 # acounts for .snakemake folder


def test_get_matches(tmp_path):
    output_path = os.path.join(tmp_path, 'test_get_matches')
    join_asvbins(
        candidate_16S_seqs=EXPECTED_CANDIDATE_FNA,
        output_dir = output_path,
        asv_seqs = MINI_ASV,
        generic_16S = MINI_16S,
        threads=4,
        verbosity=6
    )
    assert os.path.exists(os.path.join(output_path, "match_sequences.fna"))
    assert os.path.exists(os.path.join(output_path, "match_statistics.tsv"))
    assert len(os.listdir(output_path)) == 3 # acounts for .snakemake folder


def test_get_candidates(tmp_path):
    output_path = os.path.join(tmp_path, 'test_get_candidates')
    join_asvbins(
        bins = MINI_BINS,
        output_dir = output_path,
        generic_16S = MINI_16S,
        threads=4,
        verbosity=6
    )
    assert os.path.exists(os.path.join(output_path, "candidate_sequences.fna"))
    assert os.path.exists(os.path.join(output_path, "candidate_statistics.tsv"))
    assert len(os.listdir(output_path)) == 3 # acounts for .snakemake folder

