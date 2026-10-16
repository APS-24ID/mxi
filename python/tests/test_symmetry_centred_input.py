"""mxi_symmetry on data already in a centred setting: taken to its primitive
cell first, the analysis as on an indexer's primitive cell. Taken as it was,
I 2 3's conventional cell was analysed as a primitive cubic lattice -- every
Patterson group P, P 2 3 recommended (Graeme, a cluster of 36 insulin
crystals; docs/outstanding.md item 70)."""

import os
import subprocess

import pytest

SYMMETRY = os.environ.get("MXI_SYMMETRY")
EXPT = os.environ.get("MXI_SCALE_EXPT")
REFL = os.environ.get("MXI_SCALE_REFL")
pytestmark = pytest.mark.skipif(
    not (SYMMETRY and EXPT and REFL),
    reason="set MXI_SYMMETRY, MXI_SCALE_EXPT and MXI_SCALE_REFL",
)


def test_symmetrised_insulin_comes_back_as_i_2_3(tmp_path):
    # Insulin symmetrised in I 2 3 by dials.symmetry, through mxi_symmetry
    # again.
    out = subprocess.run(
        [SYMMETRY, EXPT, REFL, "-o", str(tmp_path / "again.rflx")],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert out.returncode == 0, out.stderr
    assert "The input is in a centred setting, I 2 3" in out.stdout, out.stdout
    assert "Recommended space group: I 2 3" in out.stdout, out.stdout
