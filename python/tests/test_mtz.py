"""mxeq mtz: two unmerged MTZ files of one sweep, the second indexed the
other way of I 2 3's two, ten per cent brighter, its batches numbered from
1001 and a column different -- each found."""

import numpy as np
import pytest

gemmi = pytest.importorskip("gemmi")

from mxeq import mtz as mxeq_mtz  # noqa: E402

CELL = (77.89, 77.89, 77.89, 90, 90, 90)
SWAP = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])  # (k, h, -l): I 2 3's other way


def unmerged(path, hkl, batch, i, sigi, extra=None):
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.SpaceGroup("I 2 3")
    mtz.set_cell_for_all(gemmi.UnitCell(*CELL))
    ds = mtz.add_dataset("insulin")
    ds.wavelength = 0.953738
    columns = [("M/ISYM", "Y"), ("BATCH", "B"), ("I", "J"), ("SIGI", "Q")]
    columns += [(extra, "R")] if extra else []
    for label, kind in columns:
        mtz.add_column(label, kind)
    # Each row as an exporter writes it: its index in the asymmetric unit and
    # M/ISYM saying how to get back (gemmi's ensure_asu moves the indices but
    # leaves ISYM 0 -- the measured indices lost).
    asu = gemmi.ReciprocalAsu(mtz.spacegroup)
    ops = mtz.spacegroup.operations()
    stored = [asu.to_asu([int(x) for x in h], ops) for h in hkl]
    in_asu = np.array([list(h) for h, _ in stored])
    isym = np.array([n for _, n in stored])
    rows = [in_asu[:, 0], in_asu[:, 1], in_asu[:, 2], isym, batch, i, sigi]
    if extra:
        rows.append(np.ones(len(hkl)))
    mtz.set_data(np.column_stack(rows).astype(np.float32))
    mtz.write_to_file(str(path))


@pytest.fixture
def pair(tmp_path):
    rng = np.random.default_rng(7)
    sg = gemmi.SpaceGroup("I 2 3")
    ops = sg.operations()
    asu = gemmi.ReciprocalAsu(sg)
    cell = gemmi.UnitCell(*CELL)
    unique = [
        (h, k, l)
        for h in range(0, 26)
        for k in range(0, 26)
        for l in range(-25, 26)
        if (h + k + l) % 2 == 0
        and (h, k, l) != (0, 0, 0)
        and asu.is_in([h, k, l])
        and cell.calculate_d([h, k, l]) > 3.0
    ]
    truth = rng.exponential(1000.0, len(unique))
    rotations = [np.array(op.rot) // op.DEN for op in ops.sym_ops]
    hkl, batch, which = [], [], []
    for n, h in enumerate(unique):
        for _ in range(4):
            r = rotations[rng.integers(len(rotations))]
            sign = 1 if rng.random() < 0.5 else -1
            hkl.append(sign * (np.array(h) @ r))
            batch.append(int(rng.integers(1, 101)))
            which.append(n)
    hkl, batch, which = np.array(hkl), np.array(batch), np.array(which)
    sigma = np.sqrt(truth[which] + 50.0)
    i_a = truth[which] + rng.normal(0.0, 1.0, len(which)) * sigma
    i_b = 1.1 * (truth[which] + rng.normal(0.0, 0.3, len(which)) * sigma)
    a = tmp_path / "dials.mtz"
    b = tmp_path / "mxi.mtz"
    unmerged(a, hkl, batch, i_a, sigma, extra="BG")
    unmerged(b, hkl @ SWAP, batch + 1000, i_b, 1.1 * sigma, extra="LP")
    return a, b


def test_the_other_indexing_is_found_by_the_observations(pair):
    a = mxeq_mtz.load(str(pair[0]), "DIALS")
    b = mxeq_mtz.load(str(pair[1]), "mxi")
    m, n_best, n_identity, n_sample = mxeq_mtz.choose_indexing(a, b)
    # All but the few the fixture drew twice -- one index in one batch twice,
    # matched once.
    assert n_best >= 0.999 * n_sample and n_identity < n_sample // 4
    # Taking mxi's indices through M gives DIALS's: (k, h, -l) is its own
    # inverse, so M is it, or it composed with an operation of I 2 3 that
    # moves no observation between batches -- the identity only, here.
    assert np.array_equal(m, SWAP)


def test_the_report_finds_the_scale_the_columns_and_the_batches(pair, tmp_path):
    a = mxeq_mtz.load(str(pair[0]), "DIALS")
    b = mxeq_mtz.load(str(pair[1]), "mxi")
    text = mxeq_mtz.compare(a, b, csv=str(tmp_path / "m.csv"))
    assert "only in DIALS: BG(R)" in text and "only in mxi: LP(R)" in text
    assert "under h' = h M" in text
    matched = next(l for l in text.splitlines() if " matched: " in l)
    assert "100.0% of DIALS's" in matched and "100.0% of mxi's" in matched
    ratio = next(l for l in text.splitlines() if l.strip().startswith("intensity: CC"))
    assert abs(float(ratio.split()[-1]) - 1.1) < 0.01
    sigma = next(l for l in text.splitlines() if l.strip().startswith("sigma:"))
    assert abs(float(sigma.split()[-1]) - 1.1) < 1e-3
    overall = text.splitlines()[-3]  # the shells' last row, all of them
    cc = float(overall.split()[3])
    assert cc > 0.99
    assert (tmp_path / "m.csv").read_text().count("\n") > 1000


def test_a_merged_file_is_refused_by_name(tmp_path):
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.SpaceGroup("I 2 3")
    mtz.set_cell_for_all(gemmi.UnitCell(*CELL))
    mtz.add_dataset("x")
    mtz.add_column("I", "J")
    mtz.add_column("SIGI", "Q")
    mtz.set_data(np.array([[1, 1, 0, 10, 1]], dtype=np.float32))
    path = tmp_path / "merged.mtz"
    mtz.write_to_file(str(path))
    with pytest.raises(ValueError, match="no M/ISYM column: is it unmerged"):
        mxeq_mtz.load(str(path), "merged")


def test_the_command(pair, capsys):
    from mxeq.cli import main

    code = main(["mtz", str(pair[0]), str(pair[1]), "--labels", "DIALS", "mxi"])
    assert code == 0
    assert "The merged data, merged here the same way" in capsys.readouterr().out
