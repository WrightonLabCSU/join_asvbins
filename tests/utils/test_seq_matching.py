import pytest
import pandas as pd

from join_asvbins.utils.seq_matching import process_barrnap

def test_barnap_procesing():
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
    output_df = process_barrnap(input_df)
    output_df = output_df[['header','barrnap_header', 'seq']]
    pd.testing.assert_frame_equal(output_df, expect_df,
                                  check_index_type=False, check_names=False)


