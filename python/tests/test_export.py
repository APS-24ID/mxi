"""mxi_export: an unmerged MTZ as dials.export writes one (docs/export.md).

The insulin fixture, scaled, exported and read back with gemmi: every row's I
and SIGI the scaled intensity and its sigma divided by the inverse scale, its
SCALEUSED, XDET, YDET, ROT and FRACTIONCALC the table's, every batch a header
with U orthonormal and the scan axis on +Z -- the Cambridge frame. And two
sweeps whose images overlap -- both numbered from one -- the second's batches
moved to the next number ending in 01, as dials' offsets have it. Needs
MXI_EXPORT, MXI_SCALE, MXI_SCALE_EXPT and MXI_SCALE_REFL; and gemmi.
"""

import copy
import json
import os
import subprocess

import numpy as np
import pytest

from mxeq import refl

gemmi = pytest.importorskip("gemmi")
EXPORT = os.environ.get("MXI_EXPORT")
SCALE = os.environ.get("MXI_SCALE")
EXPT = os.environ.get("MXI_SCALE_EXPT")
REFL = os.environ.get("MXI_SCALE_REFL")
needs = pytest.mark.skipif(
    not (EXPORT and SCALE and EXPT and REFL),
    reason="set MXI_EXPORT, MXI_SCALE, MXI_SCALE_EXPT, MXI_SCALE_REFL",
)


def run(tmp_path, *words):
    r = subprocess.run(list(words), capture_output=True, text=True, cwd=tmp_path)
    assert r.returncode == 0, r.stderr + r.stdout[-1000:]
    return r.stdout


def scaled(tmp_path, expt, reflections):
    run(
        tmp_path,
        SCALE,
        expt,
        reflections,
        "-o",
        "scaled.refl",
        "--output-expt",
        "scaled.expt",
    )
    return tmp_path / "scaled.expt", tmp_path / "scaled.refl"


@needs
def test_the_scaled_intensities_and_headers_are_as_the_table_and_models_give(tmp_path):
    expt, table_path = scaled(tmp_path, EXPT, REFL)
    out = run(tmp_path, EXPORT, str(expt), str(table_path), "-o", "out.mtz")
    assert "scaled intensities" in out
    m = gemmi.read_mtz_file(str(tmp_path / "out.mtz"))
    labels = [c.label for c in m.columns]
    assert labels[:9] == [
        "H",
        "K",
        "L",
        "M/ISYM",
        "BATCH",
        "I",
        "SIGI",
        "SCALEUSED",
        "SIGSCALEUSED",
    ]
    assert labels[-6:] == ["FRACTIONCALC", "XDET", "YDET", "ROT", "LP", "QE"]

    t = refl.load(str(table_path))
    c = t.columns
    hkl = np.asarray(c["miller_index"]).reshape(-1, 3)
    z = np.asarray(c["xyzobs.px.value"], float).reshape(-1, 3)[:, 2]
    batch = np.floor(z).astype(int) + 1
    isf = np.asarray(c["inverse_scale_factor"], float)
    expect = {}
    for i in range(t.nrows):
        expect.setdefault((tuple(hkl[i]), int(batch[i])), []).append(i)

    m.switch_to_original_hkl()
    data = np.array(m, copy=False)
    col = {name: labels.index(name) for name in labels}
    matched = 0
    for row in data[:2000]:
        key = ((int(row[0]), int(row[1]), int(row[2])), int(row[col["BATCH"]]))
        candidates = expect.get(key, [])
        assert candidates, key
        hit = [
            i
            for i in candidates
            if abs(c["intensity.scale.value"][i] / isf[i] - row[col["I"]])
            <= 1e-4 * max(1, abs(row[col["I"]]))
        ]
        assert hit, key
        i = hit[0]
        sigma = np.sqrt(c["intensity.scale.variance"][i]) / isf[i]
        assert row[col["SIGI"]] == pytest.approx(sigma, rel=1e-5)
        assert row[col["SCALEUSED"]] == pytest.approx(isf[i], rel=1e-6)
        xyz = np.asarray(c["xyzcal.px"], float).reshape(-1, 3)[i]
        assert row[col["XDET"]] == pytest.approx(xyz[0], rel=1e-6)
        assert row[col["YDET"]] == pytest.approx(xyz[1], rel=1e-6)
        assert row[col["FRACTIONCALC"]] == pytest.approx(c["partiality"][i], rel=1e-6)
        matched += 1
    assert matched == min(2000, m.nreflections)

    numbers = {b.number for b in m.batches}
    assert set(data[:, col["BATCH"]].astype(int)) <= numbers
    for b in m.batches[:: max(1, len(m.batches) // 10)]:
        f = list(b.floats)
        U = np.array(f[6:15]).reshape(3, 3).T
        assert np.allclose(U @ U.T, np.eye(3), atol=1e-5)
        assert np.allclose(f[38:41], [0, 0, 1], atol=1e-6)
        assert f[80] == -1.0


@needs
def test_overlapping_sweeps_batches_move_to_the_next_epoch(tmp_path):
    # The sweep split in two, both halves numbered from image one: the second
    # overlaps the first, and goes to batches 201 to 350.
    e = json.load(open(EXPT))
    n = e["scan"][0]["image_range"][1]
    half = n // 2
    out = copy.deepcopy(e)
    for k in ("beam", "detector", "goniometer", "scan", "imageset"):
        out[k] = [copy.deepcopy(e[k][0]), copy.deepcopy(e[k][0])]
    c0 = e["crystal"][0]
    varying = "A_at_scan_points" in c0
    if varying:
        out["crystal"] = [copy.deepcopy(c0), copy.deepcopy(c0)]
    out["experiment"] = [
        copy.deepcopy(e["experiment"][0]),
        copy.deepcopy(e["experiment"][0]),
    ]
    for i, (lo, hi) in enumerate(((0, half), (half, n))):
        s = out["scan"][i]
        s["image_range"] = [1, hi - lo]
        for key, values in s["properties"].items():
            if isinstance(values, list) and len(values) == n:
                s["properties"][key] = values[lo:hi]
        if varying:
            out["crystal"][i]["A_at_scan_points"] = c0["A_at_scan_points"][lo : hi + 1]
        x = out["experiment"][i]
        for k in ("beam", "detector", "goniometer", "scan", "imageset"):
            x[k] = i
        x["crystal"] = i if varying else 0
        x["identifier"] = f"half-{i}"
        x.pop("profile", None)
    out.pop("profile", None)
    json.dump(out, open(tmp_path / "two.expt", "w"))
    t = refl.load(REFL)
    for name in ("xyzcal.px", "xyzobs.px.value"):
        xyz = np.asarray(t.columns[name], float).reshape(-1, 3).copy()
        second = xyz[:, 2] >= half
        xyz[second, 2] -= half
        t.columns[name] = xyz
    zc = np.asarray(refl.load(REFL).columns["xyzcal.px"], float).reshape(-1, 3)[:, 2]
    t.columns["id"] = (zc >= half).astype(np.int32)
    t.types["id"] = "int"
    t.identifiers = {0: "half-0", 1: "half-1"}
    refl.write(str(tmp_path / "two.refl"), t)
    expt, table = scaled(tmp_path, "two.expt", "two.refl")
    run(tmp_path, EXPORT, str(expt), str(table), "-o", "two.mtz")
    m = gemmi.read_mtz_file(str(tmp_path / "two.mtz"))
    numbers = [b.number for b in m.batches]
    assert numbers[:half] == list(range(1, half + 1))
    assert numbers[half:] == list(range(201, 201 + n - half))
    labels = [c.label for c in m.columns]
    batches = set(np.array(m, copy=False)[:, labels.index("BATCH")].astype(int))
    assert batches <= set(numbers) and any(b >= 201 for b in batches)


@needs
def test_integrated_intensities_are_corrected_as_dials_export_corrects_them(tmp_path):
    # Profile-fitted by LP / QE; summed by LP / QE / partiality, as a summed
    # intensity is the part recorded and a profile-fitted one the whole.
    run(tmp_path, EXPORT, EXPT, REFL, "-o", "int.mtz")
    m = gemmi.read_mtz_file(str(tmp_path / "int.mtz"))
    labels = [x.label for x in m.columns]
    assert labels[5:9] == ["IPR", "SIGIPR", "I", "SIGI"]
    t = refl.load(REFL)
    c = t.columns
    hkl = np.asarray(c["miller_index"]).reshape(-1, 3)
    z = np.asarray(c["xyzobs.px.value"], float).reshape(-1, 3)[:, 2]
    batch = np.floor(z).astype(int) + 1
    expect = {}
    for i in range(t.nrows):
        expect.setdefault((tuple(hkl[i]), int(batch[i])), []).append(i)
    m.switch_to_original_hkl()
    data = np.array(m, copy=True)
    col = {name: labels.index(name) for name in labels}
    checked = 0
    for row in data[:2000]:
        key = ((int(row[0]), int(row[1]), int(row[2])), int(row[col["BATCH"]]))
        for i in expect.get(key, []):
            k = c["lp"][i] / c["qe"][i]
            ipr = c["intensity.prf.value"][i] * k
            if abs(ipr - row[col["IPR"]]) > 1e-4 * max(1.0, abs(ipr)):
                continue
            assert row[col["SIGIPR"]] == pytest.approx(
                np.sqrt(c["intensity.prf.variance"][i]) * k, rel=1e-5
            )
            ks = k / c["partiality"][i]
            assert row[col["I"]] == pytest.approx(
                c["intensity.sum.value"][i] * ks, rel=1e-5, abs=1e-3
            )
            assert row[col["SIGI"]] == pytest.approx(
                np.sqrt(c["intensity.sum.variance"][i]) * ks, rel=1e-5
            )
            assert row[col["LP"]] == pytest.approx(c["lp"][i], rel=1e-6)
            checked += 1
            break
    assert checked == min(2000, m.nreflections)
