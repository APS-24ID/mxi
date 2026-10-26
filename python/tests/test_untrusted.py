"""Untrusted regions, as mxi_mask draws them, through spot finding.

The planted series (mxeq.fixtures.nxmx: a grid of spots on a 384 x 256 frame),
its experiment list given an image set naming the master and, the second
time, a rectangle over part of the grid: mxi_find finds no spot inside it,
and the same spots outside it as without. And mxi_untrusted sets, lists and
clears a mask file's regions in a .rflx.
"""

import json
import os
import subprocess
import sys

import numpy as np
import pytest

from mxeq import refl, rflx

BUILD = os.environ.get("MXI_BUILD", "")
FIND = os.environ.get("MXI_FIND") or os.path.join(BUILD, "mxi_find")
UNTRUSTED = os.path.join(os.path.dirname(FIND), "mxi_untrusted")

pytestmark = pytest.mark.skipif(
    not (os.path.exists(FIND) and os.path.exists(UNTRUSTED)),
    reason="mxi_find and mxi_untrusted are not built",
)

RECTANGLE = [100, 50, 260, 170]  # x0, y0, x1, y1


def planted(tmp_path):
    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(tmp_path), "6"],
        check=True,
        capture_output=True,
    )
    master = str(tmp_path / "series.nxs")
    from mxeq.fixtures import nxmx

    nxmx.write_expt(str(tmp_path / "plain.expt"), 6)
    with open(tmp_path / "plain.expt") as f:
        doc = json.load(f)
    doc["experiment"][0]["imageset"] = 0
    doc["imageset"] = [{"__id__": "ImageSequence", "template": master}]
    return doc


def find(tmp_path, doc, name):
    expt = tmp_path / f"{name}.expt"
    with open(expt, "w") as f:
        json.dump(doc, f)
    out = tmp_path / f"{name}.refl"
    run = subprocess.run(
        [FIND, str(expt), "-o", str(out)], capture_output=True, text=True, cwd=tmp_path
    )
    assert run.returncode == 0, run.stderr + run.stdout[-2000:]
    xyz = np.asarray(refl.load(str(out))["xyzobs.px.value"])
    return xyz, run.stdout


def inside(xyz):
    x0, y0, x1, y1 = RECTANGLE
    return (xyz[:, 0] >= x0) & (xyz[:, 0] <= x1) & (xyz[:, 1] >= y0) & (xyz[:, 1] <= y1)


def test_spot_finding_leaves_out_an_untrusted_region(tmp_path):
    doc = planted(tmp_path)
    plain, _ = find(tmp_path, doc, "plain")
    assert inside(plain).sum() >= 5, "the rectangle covers planted spots"
    doc["imageset"][0]["untrusted"] = [{"panel": 0, "rectangle": RECTANGLE}]
    masked, out = find(tmp_path, doc, "masked")
    assert "Untrusted regions: 1 rectangle" in out
    assert inside(masked).sum() == 0, "no spot inside"
    # Outside, the same spots: a margin of a few pixels, where a spot astride
    # the edge loses part of itself.
    x0, y0, x1, y1 = RECTANGLE
    far = (
        (plain[:, 0] < x0 - 8)
        | (plain[:, 0] > x1 + 8)
        | (plain[:, 1] < y0 - 8)
        | (plain[:, 1] > y1 + 8)
    )
    assert len(masked) >= far.sum()


def test_mxi_untrusted_sets_lists_and_clears(tmp_path):
    doc = planted(tmp_path)
    expt = tmp_path / "a.expt"
    with open(expt, "w") as f:
        json.dump(doc, f)
    convert = os.path.join(os.path.dirname(FIND), "mxi_convert")
    # An experiment list as a .rflx, as mxi_import writes one.
    first = tmp_path / "a.rflx"
    subprocess.run(
        [convert, str(expt), "-o", str(first)],
        check=True,
        capture_output=True,
        cwd=tmp_path,
    )
    mask = tmp_path / "mask.json"
    mask.write_text(
        json.dumps([{"circle": [50, 60, 10]}, {"polygon": [[0, 0], [9, 0], [0, 9]]}])
    )
    run = subprocess.run(
        [UNTRUSTED, str(first), "--set", str(mask)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    assert "1 circle, 1 polygon" in run.stdout
    shapes = rflx.experiments_dict(str(first))["imageset"][0]["untrusted"]
    assert shapes[0]["circle"] == [50, 60, 10]
    run = subprocess.run(
        [UNTRUSTED, str(first), "--clear"], capture_output=True, text=True, cwd=tmp_path
    )
    assert run.returncode == 0 and "none" in run.stdout
    assert "untrusted" not in rflx.experiments_dict(str(first))["imageset"][0]
