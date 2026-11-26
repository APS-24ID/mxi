"""mxi merge: a merged MTZ from scaled data, as dials.merge writes one (26
November). Run with MXI_SYMMETRY_DATA naming an integrated.rflx (insulin's,
say): scaled, merged, and the file read back."""

import os
import subprocess

import numpy as np
import pytest

gemmi = pytest.importorskip("gemmi")

BUILD = os.environ.get("MXI_BUILD")
MXI = os.path.join(BUILD, "mxi") if BUILD else None
DATA = os.environ.get("MXI_SYMMETRY_DATA")

pytestmark = pytest.mark.skipif(
    not (MXI and DATA), reason="set MXI_BUILD and MXI_SYMMETRY_DATA"
)


def run(*args, cwd):
    r = subprocess.run([MXI, *args], capture_output=True, text=True, cwd=cwd)
    assert r.returncode == 0, r.stderr[-2000:]
    return r.stdout


def test_the_merged_mtz_is_dials_merges_shape_and_self_consistent(tmp_path):
    run("scale", DATA, "-o", "scaled.rflx", cwd=tmp_path)
    run("merge", "scaled.rflx", "-o", "merged.mtz", cwd=tmp_path)
    m = gemmi.read_mtz_file(str(tmp_path / "merged.mtz"))
    labels = [c.label for c in m.columns]
    assert labels == [
        "H",
        "K",
        "L",
        "FreeR_flag",
        "IMEAN",
        "SIGIMEAN",
        "N",
        "F",
        "SIGF",
        "I(+)",
        "SIGI(+)",
        "I(-)",
        "SIGI(-)",
        "N(+)",
        "N(-)",
        "F(+)",
        "SIGF(+)",
        "F(-)",
        "SIGF(-)",
        "DANO",
        "SIGDANO",
    ]
    d = {c.label: np.array(c, copy=False) for c in m.columns}
    both = np.isfinite(d["I(+)"]) & np.isfinite(d["I(-)"])
    assert both.any()
    assert np.all(d["N"][both] == d["N(+)"][both] + d["N(-)"][both])
    w1, w2 = 1 / d["SIGI(+)"][both] ** 2, 1 / d["SIGI(-)"][both] ** 2
    mean = (w1 * d["I(+)"][both] + w2 * d["I(-)"][both]) / (w1 + w2)
    assert np.allclose(mean, d["IMEAN"][both], rtol=1e-4, atol=1e-3)
    strong = d["IMEAN"] / d["SIGIMEAN"] > 10
    ratio = np.nanmedian(d["F"][strong] ** 2 / d["IMEAN"][strong])
    assert 0.99 < ratio < 1.0, ratio
    assert np.nanmin(d["F"]) > 0
    free = np.mean(d["FreeR_flag"] == 0)
    assert 0.03 < free < 0.07, free
    dano = np.isfinite(d["DANO"])
    assert np.allclose(d["DANO"][dano], (d["F(+)"] - d["F(-)"])[dano], atol=1e-3)


def test_each_crystal_merged_alone(tmp_path):
    run("scale", DATA, "-o", "scaled.rflx", cwd=tmp_path)
    out = run("merge", "scaled.rflx", "--each", "-o", "per.mtz", cwd=tmp_path)
    assert "per_0.mtz" in out
    assert (tmp_path / "per_0.mtz").exists()


def _table(text):
    lines = text.splitlines()
    start = lines.index("Merging statistics, Friedel mates merged")
    out = []
    for line in lines[start:]:
        if not line.strip() or line.startswith("  columns:"):
            break
        out.append(line)
    return out


def test_the_merge_reproduces_scalings_merging_statistics(tmp_path):
    """Graeme, 26 November: since we are merging, it should at least reproduce
    the merging statistics -- the same table, line for line, as mxi scale's,
    from the same observations."""
    scaled = run("scale", DATA, "-o", "scaled.rflx", cwd=tmp_path)
    merged = run("merge", "scaled.rflx", "-o", "merged.mtz", cwd=tmp_path)
    table = _table(scaled)
    assert len(table) > 5
    assert _table(merged) == table
