"""Tests unitaires de la CGR, de la FCGR et des distances."""

import numpy as np
import pytest

from cgr_signatures.cgr import cgr_points
from cgr_signatures.fcgr import count_kmers, fcgr_matrix, kmer_cell, kmer_label_grid
from cgr_signatures.phylogeny import distance_matrix, upgma_tree


def test_cgr_first_point_is_midway_to_corner():
    """Le premier point est au milieu entre le centre et le coin."""
    points = cgr_points("A")
    assert points[0] == pytest.approx((0.25, 0.75))


def test_cgr_ignores_ambiguous_bases():
    """Les N ne produisent pas de point."""
    assert len(cgr_points("ANNT")) == 2


def test_k1_cells_match_corners():
    """Pour k=1, chaque base occupe le quadrant de son coin."""
    assert kmer_cell("A") == (0, 0)
    assert kmer_cell("T") == (0, 1)
    assert kmer_cell("C") == (1, 0)
    assert kmer_cell("G") == (1, 1)


@pytest.mark.parametrize("k", [1, 2, 3, 4, 5])
def test_each_kmer_has_its_own_cell(k):
    """La correspondance k-mer -> case est une bijection."""
    grid = kmer_label_grid(k)
    assert grid.shape == (2 ** k, 2 ** k)
    assert len(set(grid.ravel())) == 4 ** k


def test_count_kmers_simple():
    """Comptage des 2-mers de AAAC : AA deux fois, AC une fois."""
    counts = count_kmers("AAAC", 2)
    assert counts.sum() == 3
    assert counts[0] == 2       # AA = 0*4 + 0
    assert counts[1] == 1       # AC = 0*4 + 1


def test_count_kmers_skips_ambiguous():
    """Les k-mers contenant N ne sont pas comptés."""
    assert count_kmers("ACNGT", 2).sum() == 2   # AC et GT


def test_point_after_k_steps_falls_in_kmer_cell():
    """Le point CGR après un k-mer tombe dans la case FCGR de ce k-mer."""
    rng = np.random.default_rng(1)
    sequence = "".join(rng.choice(list("ACGT"), 500))
    k = 4
    points = cgr_points(sequence)
    size = 2 ** k
    for i in range(k - 1, len(sequence)):
        x, y = points[i]
        cell = (min(int((1 - y) * size), size - 1),
                min(int(x * size), size - 1))
        assert cell == kmer_cell(sequence[i - k + 1:i + 1])


def test_fcgr_is_normalized():
    """La matrice normalisée somme à 1."""
    rng = np.random.default_rng(0)
    sequence = "".join(rng.choice(list("ACGT"), 10_000))
    matrix = fcgr_matrix(sequence, 5)
    assert matrix.shape == (32, 32)
    assert matrix.sum() == pytest.approx(1.0)


def test_identical_genomes_have_zero_distance():
    """Deux signatures identiques sont à distance nulle."""
    rng = np.random.default_rng(0)
    sequence = "".join(rng.choice(list("ACGT"), 5_000))
    vectors = np.array([fcgr_matrix(sequence, 3).ravel()] * 2
                       + [fcgr_matrix("AT" * 2_500, 3).ravel()])
    dist, condensed = distance_matrix(vectors, metric="cosine")
    assert dist[0, 1] == pytest.approx(0.0, abs=1e-12)
    assert dist[0, 2] > 0.1
    assert upgma_tree(condensed).shape == (2, 4)


def test_empty_signature_raises():
    """Une signature vide provoque une erreur explicite."""
    vectors = np.array([[0.0, 0.0], [1.0, 0.0]])
    with pytest.raises(ValueError):
        distance_matrix(vectors)
