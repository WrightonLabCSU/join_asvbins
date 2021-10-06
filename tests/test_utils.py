import random
from itertools import combinations
import pytest
import pathlib
import pandas as pd
from join_asvbins.utils import process_barfasta, filter_mdstats, \
    fasta_to_df, df_to_fasta, filter_fasta_from_headers, check_overlap, \
    get_stage1_barrnap_fasta


def test_filter_mdstats():
    """Test test_filter_mdstats"""
    test_values = {
    "max_gaps": 2,
    "max_missmatch": 2,
    "min_length": 300,
    "min_pct_id": 50,
    "min_len_pct": 50
    }
    input_df = pd.DataFrame({
       "gapopen"  : [  3,   2,   1,   0,   2,   2],
       "mismatch" : [  2,   3,   1,   1,   2,   2],
       "length"   : [400, 400, 148, 300, 300, 300],
       "qlen"     : [500, 800, 200, 600, 601, 600],
       "pident"   : [100,  50,  60,  40,  80,  60],
    }, index= ["max_gaps", "max_missmatch", "min_length", "min_pct_id",
               "min_len_pct", "None"])
    for row in range(len(test_values)):
        for i in combinations(test_values, row):
            out = filter_mdstats(input_df, **{j:test_values[j] for j in i})
            out = set(out.index)
            expect = {'None'}.union({k for k in test_values if k not in i})
            assert len(out - expect) == 0, \
               f"The test case/cases {out - expect} wrongfully NOT filtered"
            assert len(expect - out) == 0, \
                f"The test case/cases {expect - out} wrongfully filtered"


def test_barnap_procesing():
    """Test test_barnap_procesing"""
    input_df = pd.DataFrame({
        1: ['a:0-3',       'abc'],
        2: ['b:0-3',       'abc'],
        3: ['b:3-4',       'd'],
        4: ['b:3-7',       'defg'],
        5: ['c:1002-1007', 'abcde'],
        6: ['c:1003-1012', 'bcdefghij'],
        7: ['d:0-4',       'abcd'], }, index=['header', 'seq']).T
    expect_df = pd.DataFrame({
        0: ['a', 'a:0-3',       'abc'],
        1: ['b', 'b:0-7',       'abcdefg'],
        2: ['c', 'c:1002-1012', 'abcdefghij'],
        3: ['d', 'd:0-4',       'abcd'],
    }, index=['header', 'barrnap_header', 'seq']).T
    input_df['seq'] = input_df['seq'].apply(list)
    expect_df['seq'] = expect_df['seq'].apply(list)
    output_df = process_barfasta(input_df)
    output_df = output_df[['header','barrnap_header', 'seq']]
    pd.testing.assert_frame_equal(output_df, expect_df,
                                  check_index_type=False, check_names=False)


@pytest.fixture(scope="session")
def temp_fasta_protien_100(tmpdir_factory):
    """Creates a fake protein file"""
    codes = "APBQCRDSETFUGVHWIYKZLXMN"
    path = str(tmpdir_factory.mktemp('fasta').join('temp.fna'))
    headers = {f">sequence-{i}:{i+10}" for i in range(100)}
    with open(path, 'w') as out :
        for i in headers:
            out.write(i + "\n")
            seqlen = int(i.split(':')[1])
            out.write(
                "".join(random.choice(codes) for j in range(seqlen)) + "\n")
    return path


def test_fasta_to_df_no_filter(temp_fasta_protien_100):
    """Test fasta_to_df_no_filter"""
    data = fasta_to_df(temp_fasta_protien_100)
    headers = {f"sequence-{i}:{i+10}" for i in range(100)}
    assert set(data['header'].values) == headers, \
        "The header was not correctly read"
    expect_len = data['header'].str.split(':', expand=True)[1].astype(int)
    output_len = data['seq'].apply(len)
    assert expect_len.equals(output_len), \
        "The header was not correctly matched"


def test_empty_fasta_to_df(tmp_path):
    empty_file = tmp_path / "emptyfile"
    empty_file.touch()
    empty_df = fasta_to_df(empty_file)
    assert empty_df.empty, "The program will fail if barrnap fails"

def test_empty_process_barfasta():
empty_df = process_barfasta(pd.DataFrame())
    assert empty_df.empty, "The program will fail if barrnap fails"

def get_stage1_mbstats_fasta(mbstats, mbstats_fasta_path):
def save_barnap_stats(barfasta, out_barstats_path):
def (data:pd.DataFrame) -> pd.DataFrame:


def test_filter_fasta_from_headers(temp_fasta_protien_100, tmp_path):
    """Test filter_fasta_from_headers"""
    headers = {f"sequence-{i}:{i+10}" for i in range(100)}
    sample_headers = set(random.sample(headers, 20))
    filtered_path = str(tmp_path / 'filter.fa')
    filter_fasta_from_headers(temp_fasta_protien_100,
                              filtered_path, sample_headers)
    data = fasta_to_df(filtered_path)
    assert set(data['header'].values) == sample_headers, \
        "The header was not correctly read"
    expect_len = data['header'].str.split(':', expand=True)[1].astype(int)
    output_len = data['seq'].apply(len)
    assert expect_len.equals(output_len), \
        "The header was not correctly matched"


@pytest.fixture()
def overlapp_df() -> pd.DataFrame:
    input_df = pd.DataFrame({
    # TODO ask about:                               this   this
        "pass":     [  True,   True, False, False, False, False, False],
        "pass_nol": [  True,   True,  True,  True,  True, False, False],
        "qstart":   [  1070,   1057,     1,  1323,  1490,   175,  1042],
        "qend":     [     1,      1,    76,  1498,  1457,  1498,    40],
        "qlen":     [  1410,   1500,  1498,  1498,  1490,  1498,  1408],
        "sstart":   [     2,      3, 75940,     2, 33132, 56104, 81459],
        "send":     [  1105,   1069, 76015,   179, 33165, 57427, 82460],
        "slen":     [131290, 131290, 76015, 38213, 33168, 57427, 90450],
        "length":   [  1104,   1067,    76,   178,    34,  1323,  1003]
    })
    return input_df


def test_check_overlap(overlapp_df):
    assert overlapp_df['pass_nol'].equals(overlapp_df.apply(check_overlap, axis=1)), \
        "Can't properly filter for overlapping 16s"


def test_overlap_filter(overlapp_df):
    min_len_with_overlap = 1000
    min_len_pct_no_overlap = 95
    ouput_df = filter_mdstats(overlapp_df,
                   min_len_with_overlap=min_len_with_overlap,
                   min_len_pct_no_overlap=min_len_pct_no_overlap)
    assert overlapp_df[overlapp_df['pass']].equals(ouput_df)
