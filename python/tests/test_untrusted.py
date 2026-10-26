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
    # Set: a new file, masked.rflx, the input as it was (Graeme, 22 October).
    run = subprocess.run(
        [UNTRUSTED, str(first), "--set", str(mask)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    assert "1 circle, 1 polygon" in run.stdout
    masked = tmp_path / "masked.rflx"
    shapes = rflx.experiments_dict(str(masked))["imageset"][0]["untrusted"]
    assert shapes[0]["circle"] == [50, 60, 10]
    assert "untrusted" not in rflx.experiments_dict(str(first))["imageset"][0]
    # Clear, into a file named: the masked one left as it was too.
    cleared = tmp_path / "cleared.rflx"
    run = subprocess.run(
        [UNTRUSTED, str(masked), "--clear", "-o", str(cleared)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0 and "none" in run.stdout
    assert "untrusted" not in rflx.experiments_dict(str(cleared))["imageset"][0]
    assert "untrusted" in rflx.experiments_dict(str(masked))["imageset"][0]


def test_shoeboxes_are_read_and_the_outliers_exported(tmp_path):
    # The planted series' strong spots stand in for an integrated file: each
    # box's shape its bbox's, its brightest pixel at the spot's centroid; and
    # mxeq truth's export writes the boxes of the observations asked for.
    from mxeq import truth

    doc = planted(tmp_path)
    expt = tmp_path / "s.expt"
    with open(expt, "w") as f:
        json.dump(doc, f)
    strong = tmp_path / "strong.rflx"
    run = subprocess.run(
        [FIND, str(expt), "-o", str(strong)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    table = rflx.load_table(str(strong))
    xyz = np.asarray(table["xyzobs.px.value"], float)
    rows = list(range(0, table.nrows, 5))
    boxes = rflx.load_shoeboxes(str(strong), rows)
    for r, b, data in zip(rows, boxes["bbox"], boxes["data"]):
        x0, x1, y0, y1, z0, z1 = b
        assert data.shape == (z1 - z0, y1 - y0, x1 - x0)
        k, j, i = np.unravel_index(np.argmax(data), data.shape)
        assert (
            abs(x0 + i + 0.5 - xyz[r, 0]) <= 2.5
            and abs(y0 + j + 0.5 - xyz[r, 1]) <= 2.5
        )

    # Observations at those spots, the first three made the worst.
    n = table.nrows
    obs = truth.Observations(
        hkl=np.zeros((n, 3), int),
        intensity=np.ones(n),
        variance=np.ones(n),
        cell=(10, 10, 10, 90, 90, 90),
        hall=" P 1",
        x=xyz[:, 0],
        y=xyz[:, 1],
        image=xyz[:, 2],
        experiment=np.asarray(table["id"], float),
    )
    z = np.zeros(n)
    z[[0, 1, 2]] = [30.0, -20.0, 10.0]
    res = truth.Result(
        matched=n,
        unmatched=0,
        k=1.0,
        b=0.0,
        r_unmerged=0.0,
        var_z=1.0,
        robust_var_z=1.0,
        tails={},
        between=0.0,
        within=0.0,
        z=z,
        zbar=z,
        model=np.full(n, 100.0),
        rows=np.arange(n),
    )
    out = tmp_path / "outliers.npz"
    got = truth.export_shoeboxes(obs, res, str(strong), str(out), worst=3, normal=5)
    assert got["worst"] == 3 and got["normal"] == 5 and got["unmatched"] == 0
    saved = np.load(out)
    assert list(saved["integrated_row"][:3]) == [0, 1, 2]
    assert list(saved["z"][:3]) == [30.0, -20.0, 10.0]
    first = rflx.load_shoeboxes(str(strong), [0])["data"][0]
    assert np.array_equal(
        saved["data"][saved["offsets"][0] : saved["offsets"][1]], first.ravel()
    )
