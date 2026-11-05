"""mxi triage end to end, on a scaled file made of intensities drawn as
Wilson's statistics have them -- a P 4 crystal, 60 x 60 x 80 A, to 2 A,
each reflection observed as two of its symmetry mates -- untwinned, and
twinned by k,h,-l at 0.3: the command merges them, and says each as it is."""

import json
import math
import os
import re
import subprocess

import numpy as np
import pytest

BUILD = os.path.dirname(os.environ.get("MXI_FIND", "")) or os.environ.get(
    "MXI_BUILD", ""
)
MXI = os.path.join(BUILD, "mxi")

pytestmark = pytest.mark.skipif(not os.path.exists(MXI), reason="mxi is not built")

SCALED = 1 << 26  # flag::kScaled


def p4_mates(h, k, l):
    """The point group 4's equivalents of h, k, l, Friedel's with them."""
    rot = [(h, k, l), (-k, h, l), (-h, -k, l), (k, -h, l)]
    return rot + [(-a, -b, -c) for a, b, c in rot]


def scaled_file(tmp_path, twin=0.0, seed=7):
    from mxeq import refl
    from test_verify import indexed

    expt, _ = indexed(tmp_path)  # an indexed list from the planted series
    doc = json.load(open(expt))
    doc["crystal"][0].update(
        real_space_a=[60.0, 0, 0],
        real_space_b=[0, 60.0, 0],
        real_space_c=[0, 0, 80.0],
        space_group_hall_symbol=" P 4",
    )
    json.dump(doc, open(tmp_path / "scaled.expt", "w"))
    identifier = doc["experiment"][1]["identifier"]

    rng = np.random.default_rng(seed)
    truth = {}
    for h in range(-31, 32):
        for k in range(-31, 32):
            for l in range(-41, 42):
                if (h, k, l) == (0, 0, 0):
                    continue
                d = 1 / math.sqrt((h * h + k * k) / 3600 + l * l / 6400)
                if d < 2.0:
                    continue
                key = max(p4_mates(h, k, l))
                if key in truth:
                    continue
                eps = 4 if h == 0 and k == 0 else 1
                mean = 1000 * math.exp(-20 / (2 * d * d)) * eps
                centric = l == 0
                truth[key] = (
                    mean * rng.normal() ** 2 if centric else mean * rng.exponential(),
                    mean,
                    d,
                )
    rows = []
    for key, (i, mean, d) in truth.items():
        if twin:
            h, k, l = key
            i = (1 - twin) * i + twin * truth[max(p4_mates(k, h, -l))][0]
        mates = sorted(set(p4_mates(*key)))
        for m in (mates[0], mates[len(mates) // 2]):  # two observations each
            s = 0.02 * mean
            rows.append((m, i + s * rng.normal(), s * s))
    n = len(rows)
    t = refl.ReflectionTable(nrows=n)
    t.columns["id"] = np.ones(n, dtype=np.int32)
    t.columns["miller_index"] = np.array([r[0] for r in rows], dtype=np.int32)
    t.columns["flags"] = np.full(n, SCALED, dtype=np.uint64)
    t.columns["intensity.scale.value"] = np.array([r[1] for r in rows])
    t.columns["intensity.scale.variance"] = np.array([r[2] for r in rows])
    t.types.update(
        {
            "id": "int",
            "miller_index": "cctbx::miller::index<>",
            "flags": "std::size_t",
            "intensity.scale.value": "double",
            "intensity.scale.variance": "double",
        }
    )
    t.identifiers = {1: identifier}
    refl.write(str(tmp_path / "scaled.refl"), t)
    run = subprocess.run(
        [MXI, "convert", "scaled.expt", "scaled.refl", "-o", "scaled.rflx"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    return tmp_path / "scaled.rflx", len(truth)


def triage(tmp_path, *words):
    run = subprocess.run(
        [MXI, "triage", "scaled.rflx", *words],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    return run.stdout


def h_alpha(out):
    row = re.search(r"^  (\S+)\s+M\s+[\d.]+\s+\d+\s+([\d.]+)\s+([\d.]+)$", out, re.M)
    assert row, out
    return row.group(1), float(row.group(2)), float(row.group(3))


def test_untwinned_data_said_untwinned(tmp_path):
    _, unique = scaled_file(tmp_path)
    out = triage(tmp_path)
    assert (
        f"{unique} unique reflections merged, in P 4" in out
    ), "every reflection, merged"
    assert "The twin tests: 10.00 to 2.0" in out
    law, r_obs, alpha = h_alpha(out)
    assert law in ("k,h,-l", "-k,-h,-l", "-h,k,-l", "h,-k,-l") and alpha < 0.05
    assert "no twinning is suspected" in out
    assert (tmp_path / "mxi_triage.log").exists()


def test_twinned_data_said_twinned(tmp_path):
    scaled_file(tmp_path, twin=0.3)
    out = triage(tmp_path)
    _, _, alpha = h_alpha(out)
    assert abs(alpha - 0.3) < 0.04, alpha
    assert "twinning could be why" in out
    assert "The twin tests: 10.00 to 3.00" in triage(tmp_path, "--d-min", "3")


def test_it_needs_a_scaled_file(tmp_path):
    from test_verify import indexed

    indexed(tmp_path)
    subprocess.run(
        [MXI, "convert", "indexed.expt", "indexed.refl", "-o", "indexed.rflx"],
        capture_output=True,
        cwd=tmp_path,
    )
    run = subprocess.run(
        [MXI, "triage", "indexed.rflx"], capture_output=True, text=True, cwd=tmp_path
    )
    assert run.returncode == 1 and "not a scaled file" in run.stderr
