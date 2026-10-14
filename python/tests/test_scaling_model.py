"""mxi_scale writes its scale model into the scaled experiments (outstanding
item 26): an mxi_physical entry for each lattice, the experiment pointing at it,
enough to evaluate g without mxi. Needs MXI_SCALE, MXI_SCALE_EXPT and
MXI_SCALE_REFL."""

import os
import subprocess

import h5py
import hdf5plugin  # noqa: F401
import numpy as np
import pytest

from mxeq import rflx

SCALE = os.environ.get("MXI_SCALE")
EXPT = os.environ.get("MXI_SCALE_EXPT")
REFL = os.environ.get("MXI_SCALE_REFL")
needs = pytest.mark.skipif(
    not (SCALE and EXPT and REFL),
    reason="set MXI_SCALE, MXI_SCALE_EXPT, MXI_SCALE_REFL",
)


def spline(points, t):
    # As the entry's configuration says: uniform cubic B-spline, n + 1 segments.
    p = np.asarray(points, float)
    n = len(p)
    if n < 2:
        return np.full_like(t, p[0])
    t = np.clip(t, 0.0, 1.0)
    u = t * (n + 1)
    i = np.minimum(np.floor(u), n).astype(int)
    f = u - i
    w = [
        (1 - 3 * f + 3 * f**2 - f**3) / 6,
        (4 - 6 * f**2 + 3 * f**3) / 6,
        (1 + 3 * f + 3 * f**2 - 3 * f**3) / 6,
        f**3 / 6,
    ]
    return sum(w[k] * p[np.clip(i + k - 2, 0, n - 1)] for k in range(4))


@needs
def test_the_scale_model_is_written_and_gives_back_the_scales(tmp_path):
    out = tmp_path / "scaled.rflx"
    r = subprocess.run(
        [SCALE, EXPT, REFL, "-o", str(out)], capture_output=True, text=True
    )
    assert r.returncode == 0, r.stderr
    d = rflx.experiments_dict(str(out))
    models = d["scaling_model"]
    pointed = [e for e in d["experiment"] if e.get("scaling_model") is not None]
    assert models and len(pointed) == len(models)
    for e in d["experiment"]:
        assert (e.get("scaling_model") is None) == (e.get("crystal") is None)
    for m in models:
        assert m["__id__"] == "mxi_physical"
        for part in ("scale", "decay", "absorption"):
            assert len(m[part]["parameters"]) == m[part]["n_parameters"]
            assert all(s >= 0 for s in m[part]["est_standard_devs"])
        assert {"a", "b", "a_background", "split"} <= set(m["error_model"])
    # Without absorption, g from the entry alone against what scaling wrote, for
    # every reflection it scaled.
    table = h5py.File(out, "r")["reflections/0"]
    ids = table["id"][()]
    flags = table["flags"][()].astype(np.int64)
    z = table["xyzcal.px"][()][:, 2]
    dd = table["d"][()]
    g_file = table["inverse_scale_factor"][()]
    checked = 0
    for k, e in enumerate(d["experiment"]):
        if e.get("scaling_model") is None:
            continue
        m = models[e["scaling_model"]]
        if m["absorption"]["lmax"] != 0:
            continue
        first, last = d["scan"][e["scan"]]["image_range"]
        rows = (ids == k) & (((flags >> 26) & 1) == 1)
        t = (z[rows] - (first - 1)) / (last - first + 1)
        g = spline(m["scale"]["parameters"], t) * np.exp(
            spline(m["decay"]["parameters"], t) / (2 * dd[rows] ** 2)
        )
        assert np.allclose(g, g_file[rows], rtol=1e-9)
        checked += rows.sum()
    if checked == 0:
        pytest.skip("every lattice has absorption: g not evaluated here")
