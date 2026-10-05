"""mxi_integrate on several sweeps: each integrated alone, the tables and the
models joined.

The 300 image insulin sweep split at image 150 into two experiments -- each its
own beam, detector, goniometer, scan, image set, and crystal, the crystal being
scan-varying and so cut at its 150th point -- and the reflections given the id
of the half their frame is in. Integrated with the profile model given, each
row's id must match its frame, each sweep must have its own profile model, and
a reflection well inside either half must sum to exactly what the sweep
integrated whole gives it: summation depends on the box and the background, not
on the reference profiles, so a difference would be a sweep's geometry -- its
angles above all -- wrong. Needs MXI_INTEGRATE, MXI_POSTREFINE_EXPT and
MXI_POSTREFINE_REFL.
"""

import copy
import json
import os
import subprocess

import numpy as np
import pytest

from mxeq import refl

BINARY = os.environ.get("MXI_INTEGRATE")
EXPT = os.environ.get("MXI_POSTREFINE_EXPT")
REFL = os.environ.get("MXI_POSTREFINE_REFL")
needs = pytest.mark.skipif(
    not (BINARY and EXPT and REFL),
    reason="set MXI_INTEGRATE, MXI_POSTREFINE_EXPT, MXI_POSTREFINE_REFL",
)
SPLIT = 150
PROFILE = ["--sigma-b", "0.0273", "--sigma-m", "0.1286", "--threads", "2"]


def split(tmp_path):
    e = json.load(open(EXPT))
    n = e["scan"][0]["image_range"][1]
    out = copy.deepcopy(e)
    for k in ("beam", "detector", "goniometer", "scan", "imageset", "crystal"):
        out[k] = [copy.deepcopy(e[k][0]), copy.deepcopy(e[k][0])]
    out["experiment"] = [
        copy.deepcopy(e["experiment"][0]),
        copy.deepcopy(e["experiment"][0]),
    ]
    for i, (lo, hi) in enumerate(((0, SPLIT), (SPLIT, n))):
        s = out["scan"][i]
        s["image_range"] = [lo + 1, hi]
        for key, values in s["properties"].items():
            if isinstance(values, list) and len(values) == n:
                s["properties"][key] = values[lo:hi]
        out["imageset"][i]["single_file_indices"] = list(range(lo, hi))
        c = out["crystal"][i]
        c["A_at_scan_points"] = c["A_at_scan_points"][lo : hi + 1]
        x = out["experiment"][i]
        for k in ("beam", "detector", "goniometer", "scan", "imageset", "crystal"):
            x[k] = i
        x["identifier"] = f"half-{i}"
    out.pop("profile", None)
    json.dump(out, open(tmp_path / "two.expt", "w"))
    t = refl.load(REFL)
    z = np.asarray(t.columns["xyzobs.px.value"], float).reshape(-1, 3)[:, 2]
    t.columns["id"] = (z >= SPLIT).astype(np.int32)
    t.types["id"] = "int"
    t.identifiers = {0: "half-0", 1: "half-1"}
    refl.write(str(tmp_path / "two.refl"), t)


def integrate(tmp_path, expt, reflections, name):
    run = subprocess.run(
        [
            BINARY,
            expt,
            reflections,
            *PROFILE,
            "-o",
            name + ".refl",
            "--output-expt",
            name + ".expt",
        ],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr + run.stdout[-2000:]
    return (
        refl.load(str(tmp_path / (name + ".refl"))),
        json.load(open(tmp_path / (name + ".expt"))),
        run.stdout,
    )


@needs
def test_two_sweeps_integrate_as_the_sweep_whole_does(tmp_path):
    split(tmp_path)
    whole, _, _ = integrate(tmp_path, EXPT, REFL, "whole")
    two, models, out = integrate(
        tmp_path, str(tmp_path / "two.expt"), str(tmp_path / "two.refl"), "two"
    )
    assert "Integrating 2 sweeps" in out
    assert len(models["experiment"]) == 2 and len(models["crystal"]) == 2
    assert [x["profile"] for x in models["experiment"]] == [0, 1]
    assert [x["scan"] for x in models["experiment"]] == [0, 1]
    assert two.identifiers == {0: "half-0", 1: "half-1"}
    ids = np.asarray(two.columns["id"]).ravel()
    zc = np.asarray(two.columns["xyzcal.px"], float).reshape(-1, 3)[:, 2]
    assert np.all(zc[ids == 0] <= SPLIT + 1e-6) and np.all(zc[ids == 1] >= SPLIT - 1e-6)
    assert not os.path.exists(
        tmp_path / "two.refl.sweep0.refl"
    )  # the parts tidied away

    # Inside either half, well away from the split and the ends: the same sum.
    def keyed(t):
        hkl = np.asarray(t.columns["miller_index"]).reshape(-1, 3)
        z = np.asarray(t.columns["xyzcal.px"], float).reshape(-1, 3)[:, 2]
        bb = np.asarray(t.columns["bbox"]).reshape(-1, 6)
        flags = np.asarray(t.columns["flags"]).astype(np.int64)
        s = np.asarray(t.columns["intensity.sum.value"], float)
        inside = (
            ((bb[:, 5] <= SPLIT - 2) | (bb[:, 4] >= SPLIT + 2))
            & (bb[:, 4] >= 2)
            & (bb[:, 5] <= 298)
        )
        good = inside & ((flags & (1 << 8)) != 0)
        return {
            (tuple(h), round(float(zz), 1)): v
            for h, zz, v, g in zip(hkl, z, s, good)
            if g
        }

    a, b = keyed(whole), keyed(two)
    common = sorted(set(a) & set(b))
    assert len(common) > 0.9 * min(len(a), len(b))
    x = np.array([a[k] for k in common])
    y = np.array([b[k] for k in common])
    assert np.allclose(x, y, rtol=1e-9, atol=1e-9), np.max(np.abs(x - y))
