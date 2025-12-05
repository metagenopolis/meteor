# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Test phylogeny"""

from ..session import Component
from ..phylogeny import Phylogeny
import pytest
from pathlib import Path
import pickle
import pandas as pd

@pytest.fixture
def phylogeny_builder(datadir: Path, tmp_path: Path) -> Phylogeny:
    meteor = Component
    meteor.tmp_dir = tmp_path
    meteor.tree_dir = tmp_path / "tree"
    meteor.tree_dir.mkdir()
    meteor.threads = 1
    meteor.DEFAULT_GAP_CHAR = "?"
    return Phylogeny(meteor, [Path(datadir / "msp_0864.fasta")], 0.5, 4, False)


def test_compute_site_info(phylogeny_builder: Phylogeny):
    """Test compute_site_info"""
    assert phylogeny_builder.compute_site_info(["?GCT", "A?CT", "A?CT", "A?CT"]) == [
        0.25,
        0.75,
        0.0,
        0.0,
    ]
    assert phylogeny_builder.compute_site_info(["?G", "A?", "A?", "A?"]) == [0.25, 0.75]
    assert phylogeny_builder.compute_site_info(["?CG", "A?T", "A?G", "AC?"]) == [
        0.25,
        0.5,
        0.25,
    ]
    assert phylogeny_builder.compute_site_info(["?G", "?A", "?T", "GC"]) == [0.75, 0]


def test_clean_sites(phylogeny_builder: Phylogeny, datadir: Path, tmpdir: Path):
    msp = tmpdir / "msp_0864_clean.fasta"
    msp_expected_file = datadir / "msp_0864_dict.pck"
    with open(msp_expected_file, "rb") as msp_file:
        msp_expected = pickle.load(msp_file)
        with msp.open("w") as f:
            result_dict, _ = phylogeny_builder.clean_sites(
                phylogeny_builder.msp_file_list[0], f
            )
            for key in msp_expected:
                assert key in result_dict
                assert len(result_dict[key]) == len(msp_expected[key])
                assert result_dict[key] == msp_expected[key]


def test_remove_edge_labels(phylogeny_builder: Phylogeny):
    assert phylogeny_builder.remove_edge_labels("edge.0:1") == ":1"
    assert (
        phylogeny_builder.remove_edge_labels("(A:0.1, B:0.2)edge.1:0.3")
        == "(A:0.1, B:0.2):0.3"
    )

# def test_write_distance_matrix_to_tsv(phylogeny_builder: Phylogeny, tmp_path: Path):
#     """Test _write_distance_matrix_to_tsv method"""
    
#     # Create a mock distance matrix object similar to cogent3's DistanceMatrix
#     class MockDistanceMatrix:
#         def __init__(self):
#             self.names = ['seq1', 'seq2', 'seq3']
#             self._distances = {
#                 ('seq1', 'seq1'): 0.0,
#                 ('seq1', 'seq2'): 0.1,
#                 ('seq1', 'seq3'): 0.2,
#                 ('seq2', 'seq1'): 0.1,
#                 ('seq2', 'seq2'): 0.0,
#                 ('seq2', 'seq3'): 0.15,
#                 ('seq3', 'seq1'): 0.2,
#                 ('seq3', 'seq2'): 0.15,
#                 ('seq3', 'seq3'): 0.0,
#             }
        
#         def get_distance(self, name1, name2):
#             return self._distances.get((name1, name2), 0.0)
    
#     # Create mock distance matrix
#     mock_dists = MockDistanceMatrix()
    
#     # Define output file
#     dist_file = tmp_path / "msp_0864.tsv"
    
#     # Call the method
#     phylogeny_builder._write_distance_matrix_to_tsv(mock_dists, dist_file)
    
#     # Verify file exists
#     assert dist_file.exists()
    
#     # Read and verify content
#     with dist_file.open('r') as f:
#         lines = f.readlines()
    
#     # Check header line
#     expected_header = "\tseq1\tseq2\tseq3\n"
#     assert lines[0] == expected_header
    
#     # Check data lines
#     expected_lines = [
#         "seq1\t0.0\t0.1\t0.2\n",
#         "seq2\t0.1\t0.0\t0.15\n", 
#         "seq3\t0.2\t0.15\t0.0\n"
#     ]
    
#     for i, expected_line in enumerate(expected_lines, 1):
#         assert lines[i] == expected_line
    
#     # Verify total number of lines (header + 3 data lines)
#     assert len(lines) == 4

def test_generate_pairwise_comparison_table(phylogeny_builder: Phylogeny, tmp_path: Path):
    """Test _generate_pairwise_comparison_table method"""
    # Create mock sequences
    sequences = {
        "sample1": "ATCGATCGATCG",
        "sample2": "ATCGATCGATCG",
        "sample3": "ATCGNTCGATCG"
    }
    
    # Create mock distance matrix
    class MockDistanceMatrix:
        def __init__(self):
            self.names = ['sample1', 'sample2', 'sample3']
            
        def get_distance(self, name1, name2):
            # Simple distance calculation for testing
            if name1 == name2:
                return 0.0
            return 0.1
    
    dists = MockDistanceMatrix()
    
    # Define output file
    output_file = tmp_path / "test_comparison.tsv"
    
    # Call the method
    phylogeny_builder._generate_pairwise_comparison_table(sequences, dists, output_file)
    
    # Verify file exists
    assert output_file.exists()
    
    # Read and verify content
    df = pd.read_csv(output_file, sep='\t')
    
    # Check that we have the expected columns
    expected_columns = ['sample1', 'sample2', 'total_length', 'coverage_overlap_nq', 
                       'coverage_overlap_all', 'compared_bases_count', 'coverage_sample1_pct',
                       'coverage_sample2_pct', 'mean_depth_sample1', 'mean_depth_sample2', 'distance']
    assert all(col in df.columns for col in expected_columns)
    
    # Check that we have the expected number of rows (3 samples, so 3 pairs: 1-2, 1-3, 2-3)
    # Actually, with 3 samples we should have 3 pairs: (0,1), (0,2), (1,2)
    assert len(df) == 3
    
    # Check that distances are properly recorded
    assert all(df['distance'] >= 0)


def test_execute(phylogeny_builder: Phylogeny):
    phylogeny_builder.execute()
    # result = phylogeny_builder.meteor.tree_dir / "msp_0864.tree"
    # assert result.exists()
    matrix = phylogeny_builder.meteor.tree_dir / "msp_0864.tsv"
    assert matrix.exists()
