"""mxi_scale --absorption-level: dials.scale's levels, each a degree for the
absorption surface and a restraint on its harmonics.

low (the default): degree 4, 5e5; medium: 6, 5e4; high: 6, 5e3. --l-max still
sets the degree; an unknown level, or a level with --no-absorption, refused. The
fixture's 30 degree sweep is told it is 90 -- 0.3 degrees an image -- so that it
is wide enough for an absorption surface: the geometry is then wrong, and only
the model mxi_scale builds is being checked. Needs MXI_SCALE, MXI_SCALE_EXPT and
MXI_SCALE_REFL.
"""

import json
import os
import subprocess

import pytest

BINARY = os.environ.get("MXI_SCALE")
EXPT = os.environ.get("MXI_SCALE_EXPT")
REFL = os.environ.get("MXI_SCALE_REFL")
needs = pytest.mark.skipif(
    not (BINARY and EXPT and REFL),
    reason="set MXI_SCALE, MXI_SCALE_EXPT, MXI_SCALE_REFL",
)


def wide(tmp_path):
    e = json.load(open(EXPT))
    s = e["scan"][0]
    n = s["image_range"][1] - s["image_range"][0] + 1
    start = s["properties"]["oscillation"][0]
    s["properties"]["oscillation"] = [start + 0.3 * i for i in range(n)]
    json.dump(e, open(tmp_path / "wide.expt", "w"))
    return str(tmp_path / "wide.expt")


def scale(tmp_path, expt, *extra):
    return subprocess.run(
        [BINARY, expt, REFL, *extra, "-o", "s.refl", "--output-expt", "s.expt"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )


@needs
@pytest.mark.parametrize(
    "extra, degree, weight",
    [
        ((), 4, "500000"),
        (("--absorption-level", "medium"), 6, "50000"),
        (("--absorption-level", "high"), 6, "5000"),
        (("--absorption-level", "high", "--l-max", "4"), 4, "5000"),
    ],
)
def test_each_level_is_dials_scales_degree_and_restraint(
    tmp_path, extra, degree, weight
):
    r = scale(tmp_path, wide(tmp_path), *extra)
    assert r.returncode == 0, r.stderr
    level = extra[1] if extra else "low"
    assert (
        f"absorption level {level}: degree {degree}, each harmonic restrained by {weight}"
        in r.stdout
    )


@needs
@pytest.mark.parametrize(
    "extra",
    [
        ("--absorption-level", "extreme"),
        ("--absorption-level", "low", "--no-absorption"),
    ],
)
def test_a_level_unknown_or_contradicted_is_refused(tmp_path, extra):
    r = scale(tmp_path, EXPT, *extra)
    assert r.returncode == 2
