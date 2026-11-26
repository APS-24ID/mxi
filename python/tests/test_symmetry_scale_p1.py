"""mxi symmetry --scale-p1: scaled in P1 first, for the symmetry's judgement
alone (25 November). Four of 30 A71EV2AZ sweeps, P 61, came out monoclinic or
triclinic: their three- and six-fold mates, measured 120 and 60 degrees apart,
correlated 0.62 to 0.71 where the identity's Friedel pairs did 0.92 -- a
sweep's systematic changes parting symmetry mates. On well-behaved data the
judgement is the same with it, and the table written stays unscaled.

Run with MXI_SYMMETRY_DATA naming an integrated.rflx (insulin's, say).
"""

import os
import subprocess

import pytest

from mxeq import refl

BUILD = os.environ.get("MXI_BUILD")
MXI = os.path.join(BUILD, "mxi") if BUILD else None
DATA = os.environ.get("MXI_SYMMETRY_DATA")

pytestmark = pytest.mark.skipif(
    not (MXI and DATA), reason="set MXI_BUILD and MXI_SYMMETRY_DATA"
)


def recommended(stdout):
    for line in stdout.splitlines():
        if line.startswith("Recommended space group:"):
            return line.split(":", 1)[1].strip()
    return None


def test_scaled_in_p1_first_the_same_group_and_the_table_unscaled(tmp_path):
    plain = subprocess.run(
        [MXI, "symmetry", DATA, "-o", "plain.rflx"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    p1 = subprocess.run(
        [MXI, "symmetry", DATA, "--scale-p1", "-o", "p1.rflx"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert plain.returncode == 0 and p1.returncode == 0, p1.stderr[-2000:]
    assert "Scaled in P1 first (--scale-p1; scale and decay)" in p1.stdout
    assert recommended(p1.stdout) == recommended(plain.stdout)
    table = refl.load(str(tmp_path / "p1.rflx"))
    assert "inverse_scale_factor" not in table.columns, "the file written unscaled"


def test_a_scaled_files_intensities_are_divided_by_their_scale_factors(tmp_path):
    """intensity.scale.value is the corrected intensity before scaling, as in
    DIALS; the scaled one is it over inverse_scale_factor. Until 26 November
    mxi symmetry took it as scaled already: --scale-p1 never applied its
    scales. Every other observation's value and factor multiplied by ten, its
    variance by a hundred -- its scaled intensity unchanged -- and the
    symmetry's scoring must not change at all."""
    import shutil

    import h5py

    scaled = subprocess.run(
        [MXI, "scale", DATA, "-o", "scaled.rflx"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert scaled.returncode == 0, scaled.stderr[-2000:]
    shutil.copy(tmp_path / "scaled.rflx", tmp_path / "tenfold.rflx")
    with h5py.File(tmp_path / "tenfold.rflx", "r+") as f:
        g = f["reflections/0"]
        for name, factor in (
            ("intensity.scale.value", 10.0),
            ("intensity.scale.variance", 100.0),
            ("inverse_scale_factor", 10.0),
        ):
            values = g[name][...]
            values[::2] *= factor
            g[name][...] = values

    def scoring(name):
        run = subprocess.run(
            [MXI, "symmetry", name, "-o", "sym_" + name],
            capture_output=True,
            text=True,
            cwd=tmp_path,
        )
        assert run.returncode == 0, run.stderr[-2000:]
        keep = []
        for line in run.stdout.splitlines():
            if line.startswith("E(CC; S)") or line.startswith("Recommended"):
                keep.append(line)
            elif line.startswith("  ") and "|(" in line:  # each element
                keep.append(line)
        return keep

    original = scoring("scaled.rflx")
    assert original, "the scoring printed"
    assert scoring("tenfold.rflx") == original
