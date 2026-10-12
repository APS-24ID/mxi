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
ins2_27 with its imported models (strong.rflx); and for integration,
MXI_INTEGRATE and MXI_TWO_LATTICES_REFINED, its two lattices as mxi_refine left
them; for symmetry, MXI_SYMMETRY and MXI_TWO_LATTICES_INTEGRATED, the two
integrated; and for scaling, MXI_SCALE and MXI_TWO_LATTICES_SYMMETRIZED, the two
as mxi_symmetry left them.
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


INTEGRATE = os.environ.get("MXI_INTEGRATE")
REFINED = os.environ.get("MXI_TWO_LATTICES_REFINED")


@pytest.mark.skipif(
    not (INTEGRATE and REFINED), reason="set MXI_INTEGRATE and MXI_TWO_LATTICES_REFINED"
)
def test_two_lattices_integrate_each_with_its_own_profile_sharing_the_rest(tmp_path):
    # Each lattice integrated with a profile model of its own, estimated from
    # its own spots, and the other lattice's reflections masked out of its
    # backgrounds, as DIALS's 3d_threaded masks a reflection's neighbours.
    # Joined, the models that are the same are one again: each run wrote its
    # own copy, and the first join gave two detectors and two imagesets.
    r = subprocess.run(
        [INTEGRATE, REFINED, "--threads", "2", "-o", "integrated.rflx"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert r.returncode == 0, r.stderr[-2000:]
    d = rflx.experiments_dict(str(tmp_path / "integrated.rflx"))
    for kind in ("detector", "beam", "scan", "imageset"):
        assert len(d[kind]) == 1, kind
    assert len(d["crystal"]) == 2 and len(d["profile"]) == 2
    assert d["profile"][0]["sigma_m"] != d["profile"][1]["sigma_m"]
    table = refl.load(str(tmp_path / "integrated.rflx"))
    ids = np.asarray(table.columns["id"])
    assert all(40000 < np.sum(ids == l) < 45000 for l in (0, 1))
    # Each lattice's backgrounds masked of the other's reflections (and its
    # own): overlapped reflections, both kinds, in both lattices -- none, when
    # the neighbours were found before the list they index was sorted.
    flags = np.asarray(table.columns["flags"]).astype(np.int64)
    for l in (0, 1):
        assert np.sum((flags[ids == l] >> 11) & 1) > 1000  # background
        assert np.sum((flags[ids == l] >> 12) & 1) > 100  # foreground
    # Every row on the one image set: numbered before the joined list shared
    # it, the second lattice's named an image set no longer there.
    assert set(np.asarray(table.columns["imageset_id"]).tolist()) == {0}


SYMMETRY = os.environ.get("MXI_SYMMETRY")
INTEGRATED = os.environ.get("MXI_TWO_LATTICES_INTEGRATED")


@pytest.mark.skipif(
    not (SYMMETRY and INTEGRATED),
    reason="set MXI_SYMMETRY and MXI_TWO_LATTICES_INTEGRATED",
)
def test_cosym_brings_the_two_crystals_to_one_setting_and_finds_i23(tmp_path):
    # Two crystals of insulin, I23 in an I432 lattice: each indexed in a basis
    # of its own choosing, the second's half the time the first's twin. cosym
    # (Gildea & Winter, 2018) reindexes the second into the first's setting;
    # pooled, they are I m -3, and I 2 3 by absences. Pooled as they were,
    # without it, C 1 2/m 1.
    def symmetry(*words):
        r = subprocess.run(
            [SYMMETRY, INTEGRATED, "-o", "s.rflx", *words],
            capture_output=True,
            text=True,
            cwd=tmp_path,
        )
        assert r.returncode == 0, r.stderr
        return r.stdout

    out = symmetry()
    assert "One setting for 2 crystals (cosym)" in out
    line = next(l for l in out.splitlines() if l.strip().startswith("experiment 1:"))
    assert float(line.split("CC with the first")[1].split()[0]) > 0.9
    assert "Recommended space group: I 2 3" in out
    assert "Recommended space group: I 2 3" not in symmetry("--no-cosym")


SCALE = os.environ.get("MXI_SCALE")
SYMMETRIZED = os.environ.get("MXI_TWO_LATTICES_SYMMETRIZED")


@pytest.mark.skipif(
    not (SCALE and SYMMETRIZED), reason="set MXI_SCALE and MXI_TWO_LATTICES_SYMMETRIZED"
)
def test_the_two_crystals_scale_each_with_its_own_model_and_agree(tmp_path):
    # A scale model a crystal; scaled, each crystal merged apart agrees with
    # the other (0.987 on Graeme's two), the second's intensities near the
    # first's.
    r = subprocess.run(
        [SCALE, SYMMETRIZED, "-o", "scaled.rflx"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert r.returncode == 0, r.stderr
    assert "2 experiments, a model each" in r.stdout
    t = refl.load(str(tmp_path / "scaled.rflx"))
    ids = np.asarray(t.columns["id"])
    hkl = np.asarray(t.columns["miller_index"]).reshape(-1, 3)
    scaled = np.asarray(t.columns["intensity.scale.value"]) / np.asarray(
        t.columns["inverse_scale_factor"]
    )
    ok = np.isfinite(scaled)
    # m-3 in the conventional setting: 23's rotations, with inversion.
    rotations = []
    for p in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        for s in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)):
            m = np.zeros((3, 3), int)
            for i in range(3):
                m[i, p[i]] = s[i]
            rotations.append(m)

    def merged(lattice):
        out = {}
        for h, x in zip(hkl[(ids == lattice) & ok], scaled[(ids == lattice) & ok]):
            key = max(
                tuple(int(v) for v in sign * (m @ h))
                for m in rotations
                for sign in (1, -1)
            )
            out.setdefault(key, []).append(x)
        return {k: np.mean(v) for k, v in out.items()}

    a, b = merged(0), merged(1)
    common = [k for k in a if k in b]
    x = np.array([a[k] for k in common])
    y = np.array([b[k] for k in common])
    assert np.corrcoef(x, y)[0, 1] > 0.97
    strong = x > np.median(x)
    assert 0.9 < np.median(y[strong] / x[strong]) < 1.1


EXPORT = os.environ.get("MXI_EXPORT")
SCALED = os.environ.get("MXI_TWO_LATTICES_SCALED")


@pytest.mark.skipif(
    not (EXPORT and SCALED), reason="set MXI_EXPORT and MXI_TWO_LATTICES_SCALED"
)
def test_the_second_crystals_batches_start_a_block_of_a_thousand(tmp_path):
    # Graeme: a sweep whose images overlap another's moved up by whole blocks,
    # the next power of ten above the largest image number -- 1 to 600, then
    # 1001 to 1600 -- where dials' rule gives 701 to 1300.
    gemmi = pytest.importorskip("gemmi")
    r = subprocess.run(
        [EXPORT, SCALED, "-o", "scaled.mtz"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert r.returncode == 0, r.stderr
    numbers = sorted(
        b.number for b in gemmi.read_mtz_file(str(tmp_path / "scaled.mtz")).batches
    )
    assert numbers == list(range(1, 601)) + list(range(1001, 1601))


@pytest.mark.skipif(
    not (INTEGRATE and REFINED), reason="set MXI_INTEGRATE and MXI_TWO_LATTICES_REFINED"
)
def test_one_pass_over_the_images_is_each_lattice_alone(tmp_path):
    # The lattices of a sweep integrated in one pass over its images, each its
    # own profile model and reference profiles, each masking the others; the
    # same, row for row, as a run a lattice (--each-alone) -- but partial_id,
    # the row's number, and the fitted values to the rounding of their sums,
    # which gather the same reflections in another order.
    def integrate(name, *words):
        r = subprocess.run(
            [INTEGRATE, REFINED, "--threads", "2", "-o", name, *words],
            capture_output=True,
            text=True,
            cwd=tmp_path,
        )
        assert r.returncode == 0, r.stderr[-2000:]
        return refl.load(str(tmp_path / name))

    one, alone = integrate("one.rflx"), integrate("alone.rflx", "--each-alone")
    assert one.nrows == alone.nrows and sorted(one.columns) == sorted(alone.columns)
    assert one.identifiers == alone.identifiers

    def order(t):
        h = np.asarray(t.columns["miller_index"]).reshape(-1, 3)
        z = np.asarray(t.columns["xyzcal.px"]).reshape(-1, 3)[:, 2]
        return np.lexsort((z, h[:, 2], h[:, 1], h[:, 0], np.asarray(t.columns["id"])))

    a, b = order(one), order(alone)
    for k in one.columns:
        if k == "partial_id":
            continue
        x = np.asarray(one.columns[k]).reshape(one.nrows, -1)[a]
        y = np.asarray(alone.columns[k]).reshape(alone.nrows, -1)[b]
        if k.startswith("intensity.prf"):  # differing by 2e-11 at most
            assert np.allclose(x, y, rtol=1e-9, atol=1e-9, equal_nan=True), k
        elif k.startswith("profile."):  # the correlation, by 1e-8 at most
            assert np.allclose(x, y, rtol=1e-7, atol=1e-7, equal_nan=True), k
        else:
            assert np.array_equal(x, y, equal_nan=x.dtype.kind == "f"), k
    d = rflx.experiments_dict(str(tmp_path / "one.rflx"))
    assert len(d["detector"]) == 1 and len(d["imageset"]) == 1
    assert len(d["crystal"]) == 2 and len(d["profile"]) == 2
