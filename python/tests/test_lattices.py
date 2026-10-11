"""mxi_index --max-lattices on Graeme's two crystals (ins2_27), and what it
does with spots more than one lattice indexes.

Two lattices, 99.6 per cent of the spots between them -- 13617 and 7579, each
spot assigned globally to the lattice it lies nearest -- refined with one
detector. A spot as near a second lattice's point as that lattice's own spots
lie to theirs -- within --ambiguity-sigma of each lattice's RMS scatter -- is
ambiguous: by default left out, unindexed, so that refinement pulls neither
model towards it; with --ambiguous closest, given to the nearer (Graeme: leave
them out, configurably). 2, 6 and 17 at 1.5, 2 and 3. Before spots were
assigned globally, the first lattice had taken some 1700 of the second's by
chance within the indexing tolerance, its scatter was ten times its own, and
1308 were called ambiguous.

Needs MXI_INDEX and MXI_TWO_LATTICES, the strong spots mxi_find wrote for
ins2_27 with its imported models (strong.rflx).
"""

import os
import subprocess

import numpy as np
import pytest

from mxeq import refl, rflx

INDEX = os.environ.get("MXI_INDEX")
STRONG = os.environ.get("MXI_TWO_LATTICES")
pytestmark = pytest.mark.skipif(
    not (INDEX and STRONG), reason="set MXI_INDEX and MXI_TWO_LATTICES"
)


def index(tmp_path, name, *words):
    r = subprocess.run(
        [INDEX, STRONG, "--max-lattices", "3", "-o", name, *words],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert r.returncode == 0, r.stderr
    return r.stdout, refl.load(str(tmp_path / name))


def indexed(table):
    return ~np.all(
        np.asarray(table.columns["miller_index"]).reshape(-1, 3) == 0, axis=1
    )


def ambiguous_count(stdout):
    for line in stdout.splitlines():
        if "more than one lattice" in line:
            return int(line.split()[0])
    return 0


def test_two_lattices_one_detector(tmp_path):
    out, table = index(tmp_path, "closest.rflx", "--ambiguous", "closest")
    # Each lattice its own spots, the second's not taken by the first: some
    # 1700 were, before spots were assigned globally.
    assert np.sum(np.asarray(table.columns["id"])[indexed(table)] == 1) > 7000
    d = rflx.experiments_dict(str(tmp_path / "closest.rflx"))
    assert (
        len(d["experiment"]) == 2 and len(d["detector"]) == 1 and len(d["crystal"]) == 2
    )
    assert indexed(table).mean() > 0.99
    ids = np.asarray(table.columns["id"])[indexed(table)]
    assert set(np.unique(ids)) == {0, 1}
    assert table.identifiers == {
        i: e["identifier"] for i, e in enumerate(d["experiment"])
    }


def test_ambiguous_spots_left_out_by_default_or_given_to_the_nearer(tmp_path):
    out_x, left_out = index(tmp_path, "exclude.rflx")
    out_c, given = index(tmp_path, "closest.rflx", "--ambiguous", "closest")
    n = ambiguous_count(out_x)
    assert n > 0 and n == ambiguous_count(out_c)
    # The same spots either way; left out, they are the difference in what is indexed.
    assert indexed(given).sum() - indexed(left_out).sum() == n
    assert np.all(indexed(given)[indexed(left_out)])


def test_fewer_are_ambiguous_as_the_threshold_tightens(tmp_path):
    counts = [
        ambiguous_count(index(tmp_path, f"s{s}.rflx", "--ambiguity-sigma", s)[0])
        for s in ("1.5", "2", "3")
    ]
    assert counts[0] <= counts[1] <= counts[2] and counts[0] < counts[2]
