"""The detector's own pixel_mask, honoured as dials.import honours it.

A defective pixel recording large counts was data to mxi: on a small-molecule
sweep of Graeme's, every shoebox over one had a background a hundred thousand
times too high, and two observations came out at -25 sigma, which DIALS,
masking the pixel from the master's pixel_mask, never had. Now every frame read
has the pixels the mask lists made the bad-pixel marker, which spot finding,
integration and the maximum projection all take for no measurement.

Three 64x64 frames, a hot cluster of four pixels reading 5000 counts on every
one: with a pixel_mask listing them, mxi_find finds nothing there; with an all
zero pixel_mask, it does find them -- so the test cannot pass for want of a
hot pixel. Needs mxi_find (MXI_FIND) and hdf5plugin.
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

from mxeq import refl

h5py = pytest.importorskip("h5py")
hdf5plugin = pytest.importorskip("hdf5plugin")
BINARY = os.environ.get("MXI_FIND") or shutil.which("mxi_find")
pytestmark = pytest.mark.skipif(
    not BINARY, reason="mxi_find is not built or not on PATH"
)
# Three frames: on more, a pixel equally bright on all of them has its peak
# frame more than 2 px from its centroid, and the spot finder's peak-centroid
# filter -- DIALS's, max_separation 2 -- would remove it unmasked too.
FRAMES, SIZE = 3, 64
HOT = (slice(40, 42), slice(50, 52))  # y, x


def series(directory, mask_hot):
    data = directory / "data_000001.h5"
    master = directory / "master.nxs"
    rng = np.random.default_rng(0)
    with h5py.File(data, "w") as f:
        d = f.create_dataset(
            "data",
            shape=(FRAMES, SIZE, SIZE),
            dtype="u2",
            chunks=(1, SIZE, SIZE),
            **hdf5plugin.Bitshuffle(nelems=0, cname="lz4"),
        )
        for i in range(FRAMES):
            frame = rng.poisson(0.1, (SIZE, SIZE)).astype("u2")
            frame[HOT] = 5000
            d[i] = frame
    with h5py.File(master, "w") as f:
        layout = h5py.VirtualLayout(shape=(FRAMES, SIZE, SIZE), dtype="u2")
        layout[...] = h5py.VirtualSource(data.name, "data", shape=(FRAMES, SIZE, SIZE))
        g = f.create_group("entry")
        g.attrs["NX_class"] = "NXentry"
        dg = g.create_group("data")
        dg.attrs["NX_class"] = "NXdata"
        dg.create_virtual_dataset("data", layout)
        det = g.create_group("instrument").create_group("detector")
        det.attrs["NX_class"] = "NXdetector"
        mask = np.zeros((SIZE, SIZE), dtype="u4")
        if mask_hot:
            mask[HOT] = 8  # DECTRIS's flag for a noisy pixel
        det.create_dataset("pixel_mask", data=mask)
    return master


def find(directory, master):
    run = subprocess.run(
        [BINARY, str(master), "-o", "strong.refl"],
        capture_output=True,
        text=True,
        cwd=directory,
    )
    assert run.returncode == 0, run.stderr
    t = refl.load(str(directory / "strong.refl"))
    xyz = (
        np.asarray(t.columns["xyzobs.px.value"], float).reshape(-1, 3)
        if t.nrows
        else np.zeros((0, 3))
    )
    on_hot = (np.abs(xyz[:, 0] - 50.5) < 2) & (np.abs(xyz[:, 1] - 40.5) < 2)
    return int(on_hot.sum()), run.stdout


def test_pixels_the_detector_masks_are_not_measured(tmp_path):
    found, out = find(tmp_path, series(tmp_path, mask_hot=True))
    assert found == 0
    assert "4 pixels masked by the detector's own pixel_mask" in out


def test_without_the_mask_the_hot_pixels_are_found(tmp_path):
    found, out = find(tmp_path, series(tmp_path, mask_hot=False))
    assert found >= 1
    assert "masked by the detector" not in out
