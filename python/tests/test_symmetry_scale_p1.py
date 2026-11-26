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
