"""mxi integrate stops on a frame it cannot read.

A frame the writer never received is an unallocated chunk. Integration no
longer clears its shoeboxes' values before filling them from the frames --
every voxel of a box lies on a frame of the scan, and is written -- so a box
spanning a missing frame would hold whatever its memory held. Graeme has never
seen an NXmx file with one, over a great deal of data, and would rather mxi
said so and stopped than integrated around it (13 November). The planted
series with --holey: one frame never written.
"""

import os
import subprocess
import sys

import pytest

BUILD = os.path.dirname(os.environ.get("MXI_FIND", "")) or os.environ.get(
    "MXI_BUILD", ""
)
MXI = os.path.join(BUILD, "mxi")

# The planted series has no shoeboxes to estimate a profile model from.
MODEL = ["--sigma-b", "0.6", "--sigma-m", "0.1", "-j", "2"]

pytestmark = pytest.mark.skipif(not os.path.exists(MXI), reason="mxi is not built")


def test_a_missing_frame_stops_integration_and_says_which(tmp_path):
    from test_verify import indexed

    indexed(tmp_path)
    whole = subprocess.run(
        [MXI, "integrate", "indexed.expt", "indexed.refl", "-o", "whole.refl", *MODEL],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert whole.returncode == 0, whole.stderr[-2000:]
    # The same series written again with a frame missing.
    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(tmp_path), "6", "--holey"],
        check=True,
        capture_output=True,
    )
    holey = subprocess.run(
        [MXI, "integrate", "indexed.expt", "indexed.refl", "-o", "holey.refl", *MODEL],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert holey.returncode == 1, holey.stdout[-1000:] + holey.stderr[-1000:]
    assert "could not be read" in holey.stderr
    assert "series.nxs" in holey.stderr
    assert not (tmp_path / "holey.refl").exists()
