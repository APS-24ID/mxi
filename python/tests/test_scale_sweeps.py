"""mxi_scale on several sweeps: a model each, sharing the merged intensities.

The symmetrized 300 image insulin sweep split at image 150 into two
experiments -- each its own scan and, the crystal being scan-varying, its own
part of it -- with every reflection given the id of the half its frame is in,
and scaled as two: the intensities are the same, so the merging statistics must
be what scaling them as one sweep gives, near enough; every scaled row of both
sweeps must have its scale. Needs MXI_SCALE, MXI_SCALE_EXPT and MXI_SCALE_REFL.
"""

import copy
import json
import os
import re
import subprocess

import numpy as np
import pytest

from mxeq import equivalents, refl

BINARY = os.environ.get("MXI_SCALE")
EXPT = os.environ.get("MXI_SCALE_EXPT")
REFL = os.environ.get("MXI_SCALE_REFL")
needs = pytest.mark.skipif(
    not (BINARY and EXPT and REFL),
    reason="set MXI_SCALE, MXI_SCALE_EXPT, MXI_SCALE_REFL",
)
SPLIT = 150


def split(tmp_path):
    e = json.load(open(EXPT))
    n = e["scan"][0]["image_range"][1]
    out = copy.deepcopy(e)
    for k in ("beam", "detector", "goniometer", "scan", "imageset"):
        out[k] = [copy.deepcopy(e[k][0]), copy.deepcopy(e[k][0])]
    c = e["crystal"][0]
    varying = "A_at_scan_points" in c
    if varying:
        out["crystal"] = [copy.deepcopy(c), copy.deepcopy(c)]
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
        if "single_file_indices" in out["imageset"][i]:
            out["imageset"][i]["single_file_indices"] = list(range(lo, hi))
        if varying:
            out["crystal"][i]["A_at_scan_points"] = c["A_at_scan_points"][lo : hi + 1]
        x = out["experiment"][i]
        for k in ("beam", "detector", "goniometer", "scan", "imageset"):
            x[k] = i
        x["crystal"] = i if varying else 0
        x["identifier"] = f"half-{i}"
        x.pop("profile", None)
    out.pop("profile", None)
    json.dump(out, open(tmp_path / "two.expt", "w"))
    t = refl.load(REFL)
    z = np.asarray(t.columns["xyzcal.px"], float).reshape(-1, 3)[:, 2]
    t.columns["id"] = (z >= SPLIT).astype(np.int32)
    t.types["id"] = "int"
    t.identifiers = {0: "half-0", 1: "half-1"}
    refl.write(str(tmp_path / "two.refl"), t)


def scale(tmp_path, expt, reflections, name):
    run = subprocess.run(
        [
            BINARY,
            expt,
            reflections,
            "-o",
            name + ".refl",
            "--output-expt",
            name + ".expt",
        ],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    overall = {}
    for label in ("I/sigma", "CC half"):
        m = re.search(
            r"^" + re.escape(label) + r"\s+([-0-9.]+)", run.stdout, flags=re.M
        )
        assert m, label
        overall[label] = float(m.group(1))
    return overall, run.stdout, refl.load(str(tmp_path / (name + ".refl")))


@needs
def test_two_sweeps_of_the_same_intensities_scale_as_one(tmp_path):
    split(tmp_path)
    whole, _, _ = scale(tmp_path, EXPT, REFL, "whole")
    two, out, table = scale(
        tmp_path, str(tmp_path / "two.expt"), str(tmp_path / "two.refl"), "two"
    )
    assert "2 sweeps, a model each" in out
    assert "sweep 0:" in out and "sweep 1:" in out
    assert abs(two["CC half"] - whole["CC half"]) <= 0.002, (two, whole)
    assert abs(two["I/sigma"] / whole["I/sigma"] - 1.0) < 0.05, (two, whole)
    ids = np.asarray(table.columns["id"]).ravel()
    flags = np.asarray(table.columns["flags"]).astype(np.int64)
    isf = np.asarray(table.columns["inverse_scale_factor"], float)
    scaled = (flags & equivalents.SCALED) != 0
    for i in (0, 1):
        assert (scaled & (ids == i)).sum() > 1000
        assert np.all(isf[scaled & (ids == i)] > 0)
