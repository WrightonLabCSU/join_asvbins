import random
from itertools import combinations
import pytest
import pandas as pd
from pathlib import Path
from join_asvbins.snake_functions import combine_mbstats_barrnap
from join_asvbins.utils import MBSTATS_NAMES

# TODO Enable stats for howmayn bins had finds and how many 16s where founds STAGE 1
# TODO Enable stats for how many matches where founds STAGE 1


def test_combine_mbstats_barrnap_all_empty(tmp_path):
    """
    Test combine_mbstats_barrnap_all_empty

    :param tmp_path:
    """
    arguments = {
        "mbstats_fasta_path": None,
        "mbstats_stats_path": None,
        "barrnap_fasta_path": None,
        "barrnap_stats_path": None,
        "out_fasta_path": None,
        "out_stats_path": None,
        "search_tool": "test"
    }
    for arg in arguments:
        arguments[arg] = tmp_path / arg
        if not arg.startswith('out'):
            arguments[arg].touch()
        arguments[arg] = str(arguments[arg])
    with pytest.raises(ValueError, match=r'.*no hits*.') as e_info:
        combine_mbstats_barrnap(**arguments)


@pytest.fixture(scope="session")
def temp_fasta_protien_10(tmpdir_factory):
    """Creates a fake protein file"""
    codes = "APBQCRDSETFUGVHWIYKZLXMN"
    path = str(tmpdir_factory.mktemp('fasta').join('temp.fna'))
    headers = {f">sequence-{i}:{i+10}" for i in range(10)}
    with open(path, 'w') as out :
        for i in headers:
            out.write(i + "\n")
            seqlen = int(i.split(':')[1])
            out.write(
                "".join(random.choice(codes) for j in range(seqlen)) + "\n")
    return path


def test_combine_mbstats_barrnap_all_empty_filter(tmp_path):
    """
    Test combine_mbstats_barrnap_all_empty

    :param tmp_path:
    """
    arguments = {
        "mbstats_fasta_path": None,
        "mbstats_stats_path": None,
        "barrnap_fasta_path": None,
        "barrnap_stats_path": None,
        "out_fasta_path": None,
        "out_stats_path": None,
        "search_tool": "test"
    }
    for arg in arguments:
        arguments[arg] = tmp_path / arg
        if not arg.startswith('out'):
            arguments[arg].touch()
        arguments[arg] = str(arguments[arg])

    tmp_mbstats = pd.DataFrame({i:[None, None, None] for i in MBSTATS_NAMES},
                               index=[0, 1, 2])
    tmp_mbstats['length'] = [2,2,2]
    tmp_mbstats.to_csv(arguments['mbstats_stats_path'], index=False, sep='\t',
                       header=False)
    with pytest.raises(ValueError, match=r'.*no hits*.') as e_info:
        combine_mbstats_barrnap(**arguments, min_length=100)


def test_barrnap_empty_combine_mbstats_barrnap_all(tmp_path):
    """
    Test combine_mbstats_barrnap_all_empty

    :param tmp_path:
    """
    arguments = {
        "mbstats_fasta_path": None,
        "mbstats_stats_path": None,
        "barrnap_fasta_path": None,
        "barrnap_stats_path": None,
        "out_fasta_path": None,
        "out_stats_path": None,
        "search_tool": "test"
    }
    for arg in arguments:
        arguments[arg] = tmp_path / arg
        if not arg.startswith('out'):
            arguments[arg].touch()
        arguments[arg] = str(arguments[arg])

    tmp_mbstats = pd.DataFrame({i:[None, None, None] for i in MBSTATS_NAMES},
                               index=[0, 1, 2])
    tmp_mbstats['length'] = [200, 200, 200]
    tmp_mbstats.to_csv(arguments['mbstats_stats_path'], index=False, sep='\t',
                       header=False)
    with pytest.raises(ValueError, match=r'.*no hits*.') as e_info:
        combine_mbstats_barrnap(**arguments, min_length=100)
