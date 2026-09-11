"""The reindexing operator search."""

from __future__ import annotations

import numpy as np
import pytest

import fixtures
from mxeq import reindex


def test_candidate_pool_holds_both_lattice_groups():
    pool = reindex.candidate_operators()
    assert len(pool) >= 48
    identity = np.eye(3, dtype=np.int64)
    assert any(np.array_equal(m, identity) for m in pool)
    # Every candidate must be invertible over the integers.
    assert all(abs(round(float(np.linalg.det(m)))) == 1 for m in pool)
    # Closed under transposition, which is what lets the search be blind to
    # the row-vector / column-vector convention.
    keys = {m.tobytes() for m in pool}
    assert all(m.T.astype(np.int64).tobytes() in keys for m in pool)


def test_planted_reindexing_operator_is_recovered():
    operator = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]], dtype=np.int64)
    rng = np.random.default_rng(0)
    hkl_a = fixtures._unique_hkl(400, rng)
    hkl_b = hkl_a.astype(np.int64) @ operator
    found = reindex.find_operator(hkl_a, hkl_b)
    assert np.array_equal(found.operator, operator)
    assert found.fraction == pytest.approx(1.0)
    assert not found.is_identity


def test_identical_indices_give_the_identity():
    rng = np.random.default_rng(1)
    hkl = fixtures._unique_hkl(200, rng)
    found = reindex.find_operator(hkl, hkl)
    assert found.is_identity
    assert found.fraction == pytest.approx(1.0)


def test_unindexed_rows_do_not_flatter_the_search():
    # (0, 0, 0) agrees under every operator. If it were counted, a table of
    # mostly-unindexed rows would score highly for whichever came first.
    hkl = np.zeros((100, 3), dtype=np.int64)
    hkl[:5] = [[1, 2, 3], [4, 5, 6], [1, 1, 1], [2, 0, 1], [0, 3, 2]]
    found = reindex.find_operator(hkl, hkl)
    assert found.n_tested == 5


def test_unrelated_indices_score_low():
    rng = np.random.default_rng(2)
    a = fixtures._unique_hkl(300, rng)
    b = fixtures._unique_hkl(300, np.random.default_rng(99))
    assert reindex.find_operator(a, b).fraction < 0.1


def test_metric_compatibility_rejects_an_operator_the_cell_cannot_admit():
    # A tetragonal cell: swapping a and b is fine, swapping a and c is not.
    real = np.diag([50.0, 50.0, 120.0])
    rng = np.random.default_rng(3)
    hkl = fixtures._unique_hkl(300, rng)
    swap_ab = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]], dtype=np.int64)
    swap_ac = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]], dtype=np.int64)
    assert reindex.find_operator(hkl, hkl @ swap_ab, real).metric_compatible
    assert not reindex.find_operator(hkl, hkl @ swap_ac, real).metric_compatible


def test_space_group_lookup_from_hall():
    assert reindex.space_group_name(fixtures.HALL) is not None
    assert reindex.point_group_operators(fixtures.HALL)
    assert reindex.space_group_name("not a hall symbol") is None
